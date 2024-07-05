#include "utilities.h"
#include "linsolvers.h"
#include <stdlib.h>

void solve_gssor(int nx, int ny, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **var, int max_iter, double tol, double alp_relax, double **result)
{
  int i, j, ip, jp, iter;
  double l2err, arrmax1, arrmax2, err_ref, rel_err;
  
  double padd[nx+2][ny+2];
  double padold[nx+2][ny+2];
  alp_relax = 1.20;

//Initialize
  for (i = 0; i < nx; i++) 
    for (j = 0; j < ny; j++) 
        result[i][j] = 0.0;

  for (i = 0; i < nx+2; i++) {
    for (j = 0; j < ny+2; j++) {
      padd[i][j] = 0.0;
      padold[i][j] = 0.0;
        } }
    
 // copy to padded array
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            ip = i + 1;
            jp = j + 1;
            padd[ip][jp] = var[i][j];
        }
    }
 //------------------------------------------------------------------
 // perform iterations
    for (iter = 0; iter < max_iter; iter++) {
    
        // Copy the updated values of upad back to upadold
        for (i = 0; i < nx + 2; i++) {
            for (j = 0; j < ny + 2; j++) {
                padold[i][j] = padd[i][j];
            } }
    
        for (i = 0; i < nx; i++) {
            for (j = 0; j < ny; j++) {
                ip = i + 1;
                jp = j + 1;
                            
                padd[ip][jp] = (1.0 - alp_relax) * padold[ip][jp] + (alp_relax * (b[i][j] + aE[i][j] * padd[ip+1][jp] + aW[i][j] * padd[ip-1][jp] + aN[i][j] * padd[ip][jp+1] + aS[i][j] * padd[ip][jp-1]) / aP[i][j]);
                /*
                if (ip ==1 && jp ==7)
                  printf("\n%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", aP[i][j], aE[i][j], aW[i][j], aN[i][j] ,aS[i][j], b[i][j], padd[ip+1][jp], padd[ip-1][jp], padd[ip][jp+1], padd[ip][jp-1], padold[ip][jp], padd[ip][jp]); */
            }
        }
/*    // check for convergence
    l2err = get_l2err_norm(nx+2, ny+2, upadd, upadold);
    arrmax1 = get_max_of_array(nx+2, ny+2, upadd );
    arrmax2 = get_max_of_array(nx+2, ny+2, upadold);
    err_ref = fmax(arrmax1, arrmax2);  err_ref = fmax(err_ref, 1.0e-6);
    rel_err = l2err/err_ref;
    
    //printf("   > %d %9.5e  %9.5e  %9.5e\n", iter, l2err, err_ref, rel_err);
   if(rel_err < tol)
     break;*/
}  
  
  // Copy the final solution from upad to u
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            ip = i + 1;
            jp = j + 1;
            result[i][j] = padd[ip][jp];
        }
    }
  //printf("   > %d %9.5e  %9.5e  %9.5e\n", iter, l2err, err_ref, rel_err);
}