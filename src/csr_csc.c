#include "../include/csr_csc.h"
#include "../include/lib_poisson1D.h"

void set_CSR_poisson1D(int N, double **values, int **col_indices,
                       int **row_ptr) {
  int nz = 3 * N - 2;
  *values = (double *)malloc(nz * sizeof(double));
  *col_indices = (int *)malloc(nz * sizeof(int));
  *row_ptr = (int *)malloc((N + 1) * sizeof(int));

  int idx = 0;
  for (int i = 0; i < N; i++) {
    (*row_ptr)[i] = idx;
    if (i > 0) { // Subdiagonal
      (*values)[idx] = -1.0;
      (*col_indices)[idx] = i - 1;
      idx++;
    }

    (*values)[idx] = 2.0;
    (*col_indices)[idx] = i;
    idx++;
    if (i < N - 1) {
      (*values)[idx] = -1.0;
      (*col_indices)[idx] = i + 1;
      idx++;
    }
  }
  (*row_ptr)[N] = idx;
}

void set_CSC_poisson1D(int N, double **values, int **row_indices,
                       int **col_ptr) {
  int nz = 3 * N - 2;
  *values = (double *)malloc(nz * sizeof(double));
  *row_indices = (int *)malloc(nz * sizeof(int));
  *col_ptr = (int *)malloc((N + 1) * sizeof(int));

  int idx = 0;
  for (int j = 0; j < N; j++) {
    (*col_ptr)[j] = idx;
    if (j > 0) {
      (*values)[idx] = -1.0;
      (*row_indices)[idx] = j - 1;
      idx++;
    }

    (*values)[idx] = 2.0;
    (*row_indices)[idx] = j;
    idx++;
    if (j < N - 1) {
      (*values)[idx] = -1.0;
      (*row_indices)[idx] = j + 1;
      idx++;
    }
  }
  (*col_ptr)[N] = idx;
}

void dcsrmv(int N, double *values, int *col_indices, int *row_ptr, double *x,
            double *y) {
  for (int i = 0; i < N; i++) {
    y[i] = 0.0;
    for (int j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
      y[i] += values[j] * x[col_indices[j]];
    }
  }
}

void dcscmv(int N, double *values, int *row_indices, int *col_ptr, double *x,
            double *y) {
  for (int i = 0; i < N; i++) {
    y[i] = 0.0;
  }
  for (int j = 0; j < N; j++) {
    for (int k = col_ptr[j]; k < col_ptr[j + 1]; k++) {
      y[row_indices[k]] += values[k] * x[j];
    }
  }
}


// /*The following functions are a work in progress*/

void richardson_alpha_csr(int N, double *values, int *col_indices, int *row_ptr,
                          double *RHS, double *X, double *alpha_rich,
                          double *tol, int *maxit, double *resvec, int *nbite) {
    double *b = (double *)calloc(N, sizeof(double));
    double norm_b = cblas_dnrm2(N, RHS, 1); // Compute ||b||_2
    printf("Computed RHS norm: %lf\n", norm_b);

    if (norm_b == 0.0) {
        printf("RHS norm is zero. Exiting Richardson method.\n");
        free(b);
        return;
    }

    memcpy(b, RHS, N * sizeof(double));
    dcsrmv(N, values, col_indices, row_ptr, X, b);
    for (int i = 0; i < N; i++) {
        b[i] = RHS[i] - b[i];
    }

    resvec[0] = cblas_dnrm2(N, b, 1) / norm_b;
    printf("Initial Residual: %lf\n", resvec[0]);

    while (resvec[*nbite] > *tol && *nbite < *maxit) {
        cblas_daxpy(N, *alpha_rich, b, 1, X, 1); // X = X + alpha_rich * b

        memcpy(b, RHS, N * sizeof(double));
        dcsrmv(N, values, col_indices, row_ptr, X, b);
        for (int i = 0; i < N; i++) {
            b[i] = RHS[i] - b[i];
        }

        (*nbite)++;
        resvec[*nbite] = cblas_dnrm2(N, b, 1) / norm_b;
    }

    free(b);
}



void richardson_alpha_csc(int N, double *values, int *row_indices, int *col_ptr,
                          double *RHS, double *X, double *alpha_rich,
                          double *tol, int *maxit, double *resvec, int *nbite) {
  double *b = (double *)calloc(N, sizeof(double));
  double norm_b = cblas_dnrm2(N, RHS, 1); // Compute ||b||_2

  memcpy(b, RHS, N * sizeof(double));
  dcscmv(N, values, row_indices, col_ptr, X, b); // b = b - Ax
  for (int i = 0; i < N; i++) {
    b[i] = RHS[i] - b[i];
  }

  resvec[0] = cblas_dnrm2(N, b, 1) / norm_b;

  while (resvec[*nbite] > *tol && *nbite < *maxit) {
    cblas_daxpy(N, *alpha_rich, b, 1, X, 1); // X = X + alpha_rich * b

    memcpy(b, RHS, N * sizeof(double));
    dcscmv(N, values, row_indices, col_ptr, X, b); // b = b - Ax
    for (int i = 0; i < N; i++) {
      b[i] = RHS[i] - b[i];
    }

    (*nbite)++;
    resvec[*nbite] = cblas_dnrm2(N, b, 1) / norm_b;
  }

  free(b);
}

void richardson_MB_csr(int N, double *AB_values, int *AB_col_indices,
                       int *AB_row_ptr, double *RHS, double *X,
                       double *MB_values, int *MB_col_indices, int *MB_row_ptr,
                       double *tol, int *maxit, double *resvec, int *nbite) {
  double *rk = (double *)calloc(N, sizeof(double));
  double norm_B = cblas_dnrm2(N, RHS, 1); // L2 norm
  double inv_norm = 1.0 / norm_B;

  for ((*nbite) = 0; (*nbite) < (*maxit); (*nbite)++) {
    memcpy(rk, RHS, N * sizeof(double));
    dcsrmv(N, AB_values, AB_col_indices, AB_row_ptr, X, rk);
    for (int i = 0; i < N; i++) {
      rk[i] = RHS[i] - rk[i];
    }

    double norm = cblas_dnrm2(N, rk, 1) * inv_norm;
    resvec[*nbite] = norm;

    if (norm <= (*tol)) {
      break;
    }

    dcsrmv(N, MB_values, MB_col_indices, MB_row_ptr, rk, rk);
    cblas_daxpy(N, 1.0, rk, 1, X, 1); // X = X + rk
  }

  free(rk);
}

void richardson_MB_csc(int N, double *AB_values, int *AB_row_indices,
                       int *AB_col_ptr, double *RHS, double *X,
                       double *MB_values, int *MB_row_indices, int *MB_col_ptr,
                       double *tol, int *maxit, double *resvec, int *nbite) {
  double *rk = (double *)calloc(N, sizeof(double));
  double norm_B = cblas_dnrm2(N, RHS, 1); // L2 norm
  double inv_norm = 1.0 / norm_B;

  for ((*nbite) = 0; (*nbite) < (*maxit); (*nbite)++) {
    memcpy(rk, RHS, N * sizeof(double));
    dcscmv(N, AB_values, AB_row_indices, AB_col_ptr, X, rk);
    for (int i = 0; i < N; i++) {
      rk[i] = RHS[i] - rk[i];
    }

    double norm = cblas_dnrm2(N, rk, 1) * inv_norm;
    resvec[*nbite] = norm;

    if (norm <= (*tol)) {
      break;
    }
    

    dcscmv(N, MB_values, MB_row_indices, MB_col_ptr, rk, rk);
    cblas_daxpy(N, 1.0, rk, 1, X, 1); // X = X + rk
  }

  free(rk);
}
