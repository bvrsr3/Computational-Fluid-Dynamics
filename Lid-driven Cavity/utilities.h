#include<stdio.h>
#include<math.h>

double get_max_of_array(int nx, int ny, double **arr);
double get_l2err_norm(int nx, int ny, double **arr1, double **arr2);
void output_soln_u(int nx, int ny, int iter, double *x, double *y, double **T, char *strin);
void output_soln_v(int nx, int ny, int iter, double *x, double *y, double **T, char* strin);
double get_rel_l2err_norm(int nx, int ny, double **arr1, double **arr2);
void print_arr_2d(int nx, int ny, double **arr);
void copy_arr_2d(int nx, int ny, double **arrin, double **arrout);
