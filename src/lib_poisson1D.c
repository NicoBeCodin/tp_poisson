/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "../include/lib_poisson1D.h"
#include <cblas.h>

//Stockage GB en priorité colonne pour la mtrice de poisson 1D
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
