/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

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
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
}  

void set_grid_points_1D(double* x, int* la){
}

double relative_forward_error(double* x, double* y, int* la){
}

int indexABCol(int i, int j, int *lab){
  return 0;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}
