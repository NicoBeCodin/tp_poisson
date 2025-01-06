#include "../include/csr_csc.h"
#include <stdio.h>

int main() {
    int N = 5; // Size of the Poisson 1D problem
    double *csr_values, *csc_values;
    int *csr_col_indices, *csr_row_ptr, *csc_row_indices, *csc_col_ptr;

    // Create the CSR representation
    set_CSR_poisson1D(N, &csr_values, &csr_col_indices, &csr_row_ptr);

    // Create the CSC representation
    set_CSC_poisson1D(N, &csc_values, &csc_row_indices, &csc_col_ptr);

    // Print CSR
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

    // Print CSC
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

    // Perform matrix-vector multiplication
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

    // Free allocated memory
    free(csr_values);
    free(csr_col_indices);
    free(csr_row_ptr);
    free(csc_values);
    free(csc_row_indices);
    free(csc_col_ptr);

    return 0;
}
