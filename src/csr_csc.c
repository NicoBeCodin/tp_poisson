#include <stdlib.h>

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
    y[i] = 0.0; // Initialize y
  }
  for (int j = 0; j < N; j++) {
    for (int k = col_ptr[j]; k < col_ptr[j + 1]; k++) {
      y[row_indices[k]] += values[k] * x[j];
    }
  }
}
