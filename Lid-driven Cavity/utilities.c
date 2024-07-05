#include "utilities.h"

double get_max_of_array(int nx, int ny, double **arr)
{
  int i, j;
  double arrmax, val;

  arrmax = arr[0][0];
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     val = fabs(arr[i][j]);
     if(arrmax < val)
       arrmax = val;
   }
  return arrmax;
}

double get_l2err_norm(int nx, int ny, double **arr1, double **arr2)
{
  double l2err = 0.0, val;
  int i, j;

  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     val = arr1[i][j] - arr2[i][j];
     l2err += val*val;
   }
   //printf("l2err = %lf\n", l2err);
  l2err = l2err/((double) (nx*ny));
  l2err = sqrt(l2err);

  return l2err;
}

double get_rel_l2err_norm(int nx, int ny, double **arr1, double **arr2)
{

  double l2err, arrmax1, arrmax2, err_ref, rel_err;

  l2err = get_l2err_norm(nx, ny, arr1, arr2);
  arrmax1 = get_max_of_array(nx, ny, arr1);
  arrmax2 = get_max_of_array(nx, ny, arr2);
  err_ref = fmax(arrmax1, arrmax2);  err_ref = fmax(err_ref, 1.0e-6);
  rel_err = l2err/err_ref;
  return rel_err;
}

void output_soln_u(int nx, int ny, int iter, double *x, double *y, double **T, char* strin)
{
  int i, j;
  FILE* fp1, *fp2;
  char fname1[100], fname2[100];

  sprintf(fname1, "T_xy_%s_%03d_%03d_%04d.dat", strin, nx, ny, iter);
  sprintf(fname2, "T_xy_%s_%03d_%03d_%04d_plotting.dat", strin, nx, ny, iter);
  
  fp1 = fopen(fname1, "w");
  fp2 = fopen(fname2, "w");
  
  for(i=0; i<nx; i++) { 
   for(j=0; j<ny; j++)  
      {
     fprintf(fp1, "%lf\t", T[i][j]); //, Tex[i][j]);
     fprintf(fp2, "%lf, %lf, %lf\n", x[i], y[j], T[i][j]); 
     }
     fprintf(fp1,"\n");
     }
  fclose(fp1);
  fclose(fp2);
  
  printf("Done writing %s for stamp = %d to file %s\n", strin, iter, fname1);
  }

  void output_soln_v(int nx, int ny, int iter, double *x, double *y, double **T, char* strin)
{
  int i, j;
  FILE* fp1, *fp2;
  char fname1[100], fname2[100];

  sprintf(fname1, "T_xy_%s_%03d_%03d_%04d.dat", strin, nx, ny, iter);
  sprintf(fname2, "T_xy_%s_%03d_%03d_%04d_plotting.dat", strin, nx, ny, iter);
  
  fp1 = fopen(fname1, "w");
  fp2 = fopen(fname2, "w");
  
  for(j=0; j<ny; j++) {
  for(i=0; i<nx; i++) { 
     fprintf(fp1, "%lf\t", T[i][j]); //, Tex[i][j]);
     fprintf(fp2, "%lf, %lf, %lf\n", y[j], x[i], T[i][j]); 
     }
     fprintf(fp1,"\n");
     }
  fclose(fp1);
  fclose(fp2);
  
  printf("Done writing %s for stamp = %d to file %s\n", strin, iter, fname1);
  }

void copy_arr_2d(int nx, int ny, double **arrin, double **arrout)
{

  int i, j;
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
     arrout[i][j] = arrin[i][j];
}

void print_arr_2d(int nx, int ny, double **arr)
{

  int i, j;
  printf("------------------------------\n");
  for(i=0; i<nx; i++)
  {
   for(j=0; j<ny; j++)
     {
       printf("%0.4e\t",arr[i][j]);
     }
     printf("\n");
     }
     //printf("%3d %3d %.6e\n",i,j,arr[i][j]);
  printf("------------------------------\n");
}
