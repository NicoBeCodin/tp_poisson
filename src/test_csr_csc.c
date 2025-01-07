#include "../include/csr_csc.h"
#include "lib_poisson1D_writers.c"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


//To try out the different csr and csc functions
int main() {
  int N = 5; // Size of the Poisson 1D problem
  double *csr_values, *csc_values;
  int *csr_col_indices, *csr_row_ptr, *csc_row_indices, *csc_col_ptr;

  set_CSR_poisson1D(N, &csr_values, &csr_col_indices, &csr_row_ptr);

  set_CSC_poisson1D(N, &csc_values, &csc_row_indices, &csc_col_ptr);

  // CSR
  printf("CSR Representation:\n");
  printf("Values: ");
  for (int i = 0; i < 3 * N - 2; i++) {
    printf("%lf ", csr_values[i]);
  }
  printf("\nCol Indices: ");
  for (int i = 0; i < 3 * N - 2; i++) {
    printf("%d ", csr_col_indices[i]);
  }
  printf("\nRow Ptr: ");
  for (int i = 0; i <= N; i++) {
    printf("%d ", csr_row_ptr[i]);
  }
  printf("\n");

  // CSC
  printf("CSC Representation:\n");
  printf("Values: ");
  for (int i = 0; i < 3 * N - 2; i++) {
    printf("%lf ", csc_values[i]);
  }
  printf("\nRow Indices: ");
  for (int i = 0; i < 3 * N - 2; i++) {
    printf("%d ", csc_row_indices[i]);
  }
  printf("\nCol Ptr: ");
  for (int i = 0; i <= N; i++) {
    printf("%d ", csc_col_ptr[i]);
  }
  printf("\n");

  // matrix-vector multiplication
  double x[] = {1, 2, 3, 4, 5}; // Input vector
  double y_csr[5] = {0};
  double y_csc[5] = {0};

  dcsrmv(N, csr_values, csr_col_indices, csr_row_ptr, x, y_csr);
  dcscmv(N, csc_values, csc_row_indices, csc_col_ptr, x, y_csc);

  printf("Matrix-Vector Multiplication (CSR):\n");
  for (int i = 0; i < N; i++) {
    printf("%lf ", y_csr[i]);
  }
  printf("\n");

  printf("Matrix-Vector Multiplication (CSC):\n");
  for (int i = 0; i < N; i++) {
    printf("%lf ", y_csc[i]);
  }
  printf("\n");

  // Solve using Richardson methods
  double RHS[5] = {1, 2, 3, 4, 5};
  double SOL[5] = {0};
  double alpha_rich = 0.5;
  double tol = 1e-6;
  int maxit = 1000;
  double resvec[1000];
  int nbite = 0;

  clock_t start = clock();
  richardson_alpha_csr(N, csr_values, csr_col_indices, csr_row_ptr, RHS, SOL,
                       &alpha_rich, &tol, &maxit, resvec, &nbite);
  clock_t end = clock();
  printf("Time for CSR Richardson method: %ld CPU cycles\n",
         (long)(end - start));

  for (int i = 0; i < N; i++)
    SOL[i] = 0;

printf("nbite : %d\n", nbite);
  write_vec(resvec, &nbite, "RESVEC_csr.dat");

  start = clock();
  richardson_alpha_csc(N, csc_values, csc_row_indices, csc_col_ptr, RHS, SOL,
                       &alpha_rich, &tol, &maxit, resvec, &nbite);
  end = clock();
  printf("Time for CSC Richardson method: %ld CPU cycles\n",
         (long)(end - start));

  write_vec(resvec, &nbite, "RESVEC_csc.dat");

  printf("nbite : %d\n", nbite);

  free(csr_values);
  free(csr_col_indices);
  free(csr_row_ptr);
  free(csc_values);
  free(csc_row_indices);
  free(csc_col_ptr);

  return 0;
}
