/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "../include/lib_poisson1D.h"
#include "atlas_headers.h"


#include "lib_poisson1D_writers.c" //This is temporary, I need to try out my functions


#include <cblas.h>

#define TRF 0
#define TRI 1
#define SV 2

int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info = 1;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  double relres;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);//done
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);//done
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  printf("Trying out the Poisson 1D function...\n");
  
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

  //On utlise la bilbiotheque cblas pour effectuer la multiplcation matrice vecteur
  cblas_dgbmv(CblasColMajor, CblasTrans, la, la, kl, ku, 1, AB + 1, lab, EX_SOL, 1, 0,  RHS, 1);
  
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  //On utilise dnrm2 pour calculer la norme du vecteur
  double norm_2 = cblas_dnrm2(la, RHS, 1);
  //Puisque t0 = -5 et t1 = 5, la norme est sensé tendre vers 0
  printf("Norm 2 of RHS is : %lf \n", norm_2);

  printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));

  /* LU Factorization */
  if (IMPLEM == TRF) {
    //utilisation de la factorisation LU avec LAPCK
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
    if (info!=0){
      printf("problem with dgbtrf_ ...\n");

    } else{
      printf("Correct execution of dgbtrf_\n");
    }
  }
  
  /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  if (IMPLEM == TRI) {
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  }

  if (IMPLEM == TRI || IMPLEM == TRF){
    /* Solution (Triangular) */
    if (info==0){
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      if (info!=0){printf("\n INFO DGBTRS = %d\n",info);} else{
        printf("dgbtrs_ solving succesful!\n");
      }
    }else{
      printf("\n INFO = %d\n",info);
    }
  }

  /* It can also be solved with dgbsv */
  if (IMPLEM == SV) {
    int ldb = 2* kl +ku +1;
    int *ipiv_sv = (int* )malloc(sizeof(int)*la);
    double *AB_sv = (double*)malloc(sizeof(double)*ldb*la);
    memcpy(AB_sv, AB,sizeof(double)*ldb*la); 

    dgbsv_(&la, &kl, &ku, &NRHS, AB_sv, &ldb, ipiv_sv, RHS, &la, &info);

    if (info==0){
      printf("dgbsv succesful!\n");
    }else{
      printf("dgbsv failed info= %d\n", info);
    }
  free(ipiv_sv);
  free(AB_sv);

  }

  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  relres = relative_forward_error(RHS, EX_SOL, &la);
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
