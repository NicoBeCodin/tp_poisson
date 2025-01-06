/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"
#include <cblas.h>
#include <math.h>

void eig_poisson1D(double *eigval, int *la) {
  *eigval = 2.0 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

double eigmax_poisson1D(int *la) {
  return -2.0 * cos((*la) * M_PI / (*la + 1)) + 2.0;
}

double eigmin_poisson1D(int *la) { return -2.0 * cos(M_PI / (*la + 1)) + 2.0; }

double richardson_alpha_opt(int *la) {
  return 2.0 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

// void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich,
// int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec,
// int *nbite){
//   double* V = (double*)calloc(*la, sizeof(double));
//   const double normb = (1/cblas_dnrm2(*la, RHS, 1));

//   //Copie du vecteur V
//   cblas_dcopy(*la, RHS, 1.0, V, 1.0);
//   cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB,
//   *lab, X, 1.0, 1.0,V,1.0);
//   //residu

//   double residu = cblas_dnrm2(*la, V,1)/ normb;
//   resvec[*nbite] = residu;

//   while (residu> *tol && *maxit > *nbite){
//     cblas_daxpy(*la, *alpha_rich, V, 1.0,X,1.0);
//     //residu
//     cblas_dcopy(*la, RHS, 1.0, V, 1.0);
//     cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB,
//     *lab, X, 1.0, 1.0, V, 1.0);

//     resvec[*nbite] =residu;
//     residu = cblas_dnrm2(*la, V, 1)* normb;
//     ++*nbite;
//   }
//   free(V);

// }

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich,
                      int *lab, int *la, int *ku, int *kl, double *tol,
                      int *maxit, double *resvec, int *nbite) {

  double *b = (double *)calloc(*la, sizeof(double));
  double norm_b = cblas_dnrm2(*la, RHS, 1); // Compute ||b||_2

  memcpy(b, RHS, *la * sizeof(double));
  cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab,
              X, 1, 1.0, b, 1);

  resvec[0] = cblas_dnrm2(*la, b, 1) / norm_b;

  while (resvec[*nbite] > *tol && *nbite < *maxit) {

    cblas_daxpy(*la, *alpha_rich, b, 1, X, 1);

    memcpy(b, RHS, *la * sizeof(double));
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab,
                X, 1, 1.0, b, 1);

    (*nbite)++;
    resvec[*nbite] = cblas_dnrm2(*la, b, 1) / norm_b;
  }

  free(b);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,
                               int *ku, int *kl, int *kv) {
  for (int i = 0; i < *la; i++) {
    MB[i * (*lab) + 1] = AB[i * (*lab) + 1];
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,
                                     int *ku, int *kl, int *kv) {
  for (int i = 0; i < *la; i++) {
    MB[*lab * i + 1] = AB[i * (*lab) + 1];
    MB[*lab * i + 2] = AB[i * (*lab) + 2];
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab,
                   int *la, int *ku, int *kl, double *tol, int *maxit,
                   double *resvec, int *nbite) {
  double residual;
  double residual_norm;
  for (int iter = 0; iter < *maxit; iter++) {
    residual = 0;
    residual_norm = 0;
    for (int i = 0; i < *la; ++i) {
      double tmp = 0;
      for (int j = (*kl) - 1; j < *ku; ++j) {
        if (j >= 0 && j < *lab) {
          tmp += AB[*lab * i + j] * X[j - (*kl) + 1];
        }
      }
      double tmp2 = tmp + RHS[i] - (*kl) * (*ku) * X[i + 1];
      double divisor =
          (*ku) * (*ku) * MB[i + 1] - (*kl) * (*ku) * (1 - MB[i + 1]);

      if (divisor != 0.0) {
        double oldX = X[i + 1];
        X[i + 1] = (1 - MB[i + 1]) * oldX + MB[i + 1] * tmp2 / divisor;
        double new_residual = X[i + 1] - oldX;
        residual += new_residual * new_residual;
        residual_norm += oldX * oldX;
      }
    }
    resvec[iter] = sqrt(residual / residual_norm);
    if (resvec[iter < *tol]) {
      *nbite = iter + 1;
      return;
    }
    // Debugging
    printf("Iteration %d, Residual: %e\n", iter, resvec[iter]);
  }
  *nbite = *maxit;
}
