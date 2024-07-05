#include "functions.h"
#include "face_calc.h"
#include <stdio.h>


void get_u_face(int nx, int ny, double *dxf,double *dyf,double **aPu, double **u, double **p, double **du, double **uface, int type)
{

  int i, j;

if (type == 0) {     // Linear Interpolation
   
for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
      du[i][j] = 0.0;

for (i = 1; i <= nx-1; ++i) 
   for (j = 0; j < ny; ++j) 
      uface[i][j] = (u[i-1][j] + u[i][j]) * 0.5;

i=0; 
   for (j = 0; j < ny; ++j)  
      uface[i][j] = 0.0;  
 
i=nx; 
   for (j = 0; j < ny; ++j)
      uface[i][j] = 0.0; 
}     // end of LIP

else if (type == 1) {     // RC

// calculation of d's  
for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
      du[i][j] = dyf[j]/aPu[i][j];       //first time du's are calculated based on linear interpolation of face velocities

  // start populating the coefficients 
  // ------ Step 1 :: interior points ------

for(i=2; i<=nx-2; i++)  
   for(j=0; j<ny; j++) 
      uface[i][j] = 0.5*(u[i-1][j] + u[i][j]) + 0.5*(du[i-1][j] + du[i][j])*(p[i-1][j] - p[i][j]) - 0.25*du[i-1][j]*(p[i-2][j] - p[i][j]) - 0.25*du[i][j]*(p[i-1][j] - p[i+1][j]);

i = 1;
for(j=0; j<ny; j++) {
      uface[i][j] = 0.5*(u[i-1][j] + u[i][j]) + 0.5*(du[i-1][j] + du[i][j])*(p[i-1][j] - p[i][j]) - 0.25*du[i-1][j]*(p[0][j] - p[i][j]) - 0.25*du[i][j]*(p[i-1][j] - p[i+1][j]);
      }
      
i = nx-1;
for(j=0; j<ny; j++) {
      uface[i][j] = 0.5*(u[i-1][j] + u[i][j]) + 0.5*(du[i-1][j] + du[i][j])*(p[i-1][j] - p[i][j]) - 0.25*du[i-1][j]*(p[i-2][j] - p[i][j]) - 0.25*du[i][j]*(p[i-1][j] - p[nx-1][j]);
      }

i=0;
   for (j = 0; j < ny; ++j) {
      uface[i][j] = 0.0;
} 

i=nx;
   for (j = 0; j < ny; ++j) {
      uface[i][j] = 0.0;
}
}     // end of RC

}

void get_v_face(int nx, int ny, double *dxf, double *dyf, double **dv, double **vface, int type, double **v, double **p, double **aPv)
{
  int i, j;

if (type == 0) {     // Linear Interpolation

for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
      dv[i][j] = 0.0;

for (i = 0; i < nx; ++i) 
   for (j = 1; j <= ny-1; ++j )  
         vface[i][j] = (v[i][j-1] + v[i][j]) * 0.5;
   
j=0;
   for (i = 0; i < nx; ++i) 
      vface[i][j] = 0.0;

j=ny;
   for (i = 0; i < nx; ++i) 
      vface[i][j] = 0.0;

}     // end of LIP


else if (type == 1) {     // RC

  // calculation of d's  
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
      dv[i][j] = dxf[i]/aPv[i][j];       //first time du's are calculated based on linear interpolation 

  // start populating the coefficients 
  // ------ Step 1 :: interior points ------

for(i=0; i<nx; i++)
   for(j=2; j<=ny-2; j++)  
      vface[i][j] = (v[i][j-1] + v[i][j])*0.5 + 0.5*(dv[i][j-1] + dv[i][j])*(p[i][j-1] - p[i][j]) - 0.25*dv[i][j-1]*(p[i][j-2] - p[i][j]) - 0.25*dv[i][j]*(p[i][j-1] - p[i][j+1]);
      

j = 1;
for(i=0; i<nx; i++) 
      vface[i][j] = (v[i][j-1] + v[i][j])*0.5 + 0.5*(dv[i][j-1] + dv[i][j])*(p[i][j-1] - p[i][j]) - 0.25*dv[i][j-1]*(p[i][0] - p[i][j]) - 0.25*dv[i][j]*(p[i][j-1] - p[i][j+1]);
      
      
j = ny-1;
for(i=0; i<nx; i++)
      vface[i][j] =  (v[i][j-1] + v[i][j])*0.5 + 0.5*(dv[i][j-1] + dv[i][j])*(p[i][j-1] - p[i][j]) - 0.25*dv[i][j-1]*(p[i][j-2] - p[i][j]) - 0.25*dv[i][j]*(p[i][j-1] - p[i][ny-1]);
      

j=0;
   for(i=0; i<nx; i++)
      vface[i][j] = 0.0;


j=ny;
   for(i=0; i<nx; i++) 
      vface[i][j] = 0.0;

}     // end of RC

}  // end of function

void face_corrections(int nx, int ny, double **ppr, double **du, double **dv, double **ufacecorr, double **vfacecorr)
{
  int i, j;

// uface corrections 
for (j = 0; j < ny; ++j)
   for (i = 1; i <= nx-1; ++i) 
      ufacecorr[i][j] = (du[i-1][j] + du[i][j])*0.5*(ppr[i-1][j] - ppr[i][j]);

i = 0;
   for (j = 0; j < ny; ++j) 
      ufacecorr[i][j] = 0.0;

i = nx;
   for (j = 0; j < ny; ++j)
      ufacecorr[i][j] = 0.0;

// vface corrections
for (i = 0; i < nx; ++i)   
   for (j = 1; j <= ny-1; ++j) 
         vfacecorr[i][j] = (dv[i][j-1] + dv[i][j])*0.5*(ppr[i][j-1] - ppr[i][j]);
   
j = 0;
  for (i = 0; i < nx; ++i) 
      vfacecorr[i][j] = 0.0;

j = ny;
   for (i = 0; i < nx; ++i)
      vfacecorr[i][j] = 0.0;

}     





