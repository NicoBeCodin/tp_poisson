#ifndef CSR_CSC_H
#define CSR_CSC_H

#include <stdlib.h>

// Function to set the CSR representation of the 1D Poisson matrix
void set_CSR_poisson1D(int N, double **values, int **col_indices, int **row_ptr);

// Function to set the CSC representation of the 1D Poisson matrix
void set_CSC_poisson1D(int N, double **values, int **row_indices, int **col_ptr);

// Function to perform matrix-vector multiplication in CSR format
void dcsrmv(int N, double *values, int *col_indices, int *row_ptr, double *x, double *y);

// Function to perform matrix-vector multiplication in CSC format
void dcscmv(int N, double *values, int *row_indices, int *col_ptr, double *x, double *y);

// Function to perform Richardson iteration using CSR matrix format
void richardson_alpha_csr(int N, double *values, int *col_indices, int *row_ptr, 
                          double *RHS, double *SOL, double *alpha_rich, 
                          double *tol, int *maxit, double *resvec, int *nbite);

// Function to perform Richardson iteration using CSC matrix format
void richardson_alpha_csc(int N, double *values, int *row_indices, int *col_ptr, 
                          double *RHS, double *SOL, double *alpha_rich, 
                          double *tol, int *maxit, double *resvec, int *nbite);


#endif // CSR_CSC_H
