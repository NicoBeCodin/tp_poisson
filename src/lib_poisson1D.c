/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"
#include <cblas.h>

//Stockage GB en priorité colonne pour la matrice de poisson 1D
void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int i, j,k;
  for(j=0;j<(*la); j++){
    k = j * (*lab);
    if (*kv >=0){
      for (i=0; i<*kv; i++){
        AB[k+i] = 0.0;
      }
      
    }
    AB[k+*kv] = -1.0;
    AB[k+*kv + 1] = 2.0;
    AB[k+*kv + 2] = -1.0;
    }
  AB[0]=0.0;
  if (*kv == 1){
    AB[1] = 0;
  }
  AB[(*lab)* (*la) - 1] =0.0;
  
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int i, j, k;
  for (j=0;j<(*la);j++){
    k = j*(*lab);
    if (*kv>=0){
      for (i=0;i< *kv;i++){
	AB[k+i]=0.0;
      }
    }
    AB[k+ *kv]=0.0;
    AB[k+ *kv+1]=1.0;
    AB[k+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  RHS[0] = *BC0;
  RHS[(*la)-1] = *BC1;
  for(int i= 1; i<(*la) - 1; i++){
    RHS[i] = 0.0;
  }

}  

//pour calculer la solution exacte de manière linéaire
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  for (int i=0; i<*la; i++, EX_SOL++, X++){
    *EX_SOL = *BC0 + *X * (*BC1 - *BC0);
  }
}  

void set_grid_points_1D(double* x, int* la){
  double point_size = 1.0/(1.0 * (*la + 1.0));
  for(int i =0; i<(*la); ++i){
    x[i] = (i+1) *point_size;
  }

}

double relative_forward_error(double* x, double* y, int* la){
    double err = 0.0;
    double norm_x = 0.0;
    
    for (int i = 0; i < *la; i++) {
        err += pow((x[i] - y[i]), 2); 
        norm_x += pow(x[i], 2);
    }
    //calcul norme L2 de x
    norm_x = sqrt(norm_x);
    
    return (sqrt(err) / norm_x);
}

int indexABCol(int i, int j, int *lab){
  return i + j*(*lab);
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}

void eig_poisson1D(double* eigval, int *la){
  *eigval = 2.0/(eigmin_poisson1D(la) + eigmax_poisson1D(la));

}

double eigmax_poisson1D(int *la){
  return -2.0*cos((*la)*M_PI/(*la+1)) +2.0;
}

double eigmin_poisson1D(int *la){
  return -2.0*cos(M_PI/(*la+1)) +2.0;
}

double richardson_alpha_opt(int *la){
  return 2.0/(eigmin_poisson1D(la) + eigmax_poisson1D(la));
}


void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich,
                      int *lab, int *la, int *ku, int *kl, double *tol,
                      int *maxit, double *resvec, int *nbite) {


    double *b = (double *)calloc(*la, sizeof(double));
    double norm_b = cblas_dnrm2(*la, RHS, 1); // Compute ||b||_2

    memcpy(b, RHS, *la * sizeof(double));
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, b, 1);


    resvec[0] = cblas_dnrm2(*la, b, 1) / norm_b;


    while (resvec[*nbite] > *tol && *nbite < *maxit) {

        cblas_daxpy(*la, *alpha_rich, b, 1, X, 1);


        memcpy(b, RHS, *la * sizeof(double));
        cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, b, 1);

        (*nbite)++;
        resvec[*nbite] = cblas_dnrm2(*la, b, 1) / norm_b;
    }

    free(b);
}



void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv) {
  for (int i = 0; i < *la; i++) {
    for (int j = 0; j < *lab; j++) {
      int ind = (*lab)*i;
      MB[ind + j] = 0.0;
    }
    // Fill the diagonal
    int ind = ((*lab)*i) + (*ku);
    MB[ind] = AB[ind];
  }
}


void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  for (int i=0; i<*la; i++){
    MB[*lab * i +1] = AB[i * (*lab) + 1];
    MB[*lab * i + 2] = AB[i*(*lab)+2];
  }
}


void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite) {
  double* rk = malloc(sizeof(double) * (*la));
  double norm_B = cblas_dnrm2(*la, RHS, 1); // L2 norm
  double inv_norm = 1 / norm_B;
  double norm = 0.0;
  int* ipiv = malloc(sizeof(int) * (*la)); // Pivots
  int info = 0;// 0->success
  int NHRS = 1; // Right-hand sides
  int kuu = (*ku)-1;



  dgbtrf_(la, la, kl, &kuu, MB, lab, ipiv, &info);
  for ((*nbite) = 0; (*nbite) < (*maxit); (*nbite)++) {

    for (int i = 0; i < (*la); i++) {
      rk[i] = RHS[i];
    }

    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, rk, 1);

    norm = cblas_dnrm2(*la, rk, 1) * inv_norm;
    resvec[(*nbite)] = norm;

    dgbtrs_("N", la, kl, &kuu, &NHRS, MB, lab, ipiv, rk, la, &info, 1);

    cblas_daxpy(*la, 1, rk, 1, X, 1);

    if (norm <= (*tol))
      break;
  }

  free(rk);
  free(ipiv);

}
