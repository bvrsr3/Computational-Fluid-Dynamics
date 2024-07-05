#include "functions.h"
#include <stdio.h>


void grid(int nx, double xst, double xen, double *xc, double *xf, double *dxc, double *dxf)
{
  int i;
  double dxunif;

  // uniform mesh for now; 
  // can use stretching factors to place xc nodes later
  dxunif = (xen-xst)/(double)nx;

  // xc[i] s are inside the CVs
  for(i=0; i<nx; i++)
    xc[i] = ((double)i + 0.5) * dxunif;

  // xf[i] s are at the faces; mid-way betwee xc[i]s
  for(i=0; i<=nx; i++)
    xf[i] = ((double)i) * dxunif;

  // dxc[i] s are spacing between adjacent xcs; ends are different 
  dxc[0] = xc[0]-xf[0];
  for(i=1; i<=nx-1; i++)
    dxc[i] = xc[i] - xc[i-1];
  dxc[nx] = xf[nx]- xc[nx-1];

  // dxf[i] s are spacing between adjacent xes
  for(i=0; i<nx; i++)
    dxf[i] = xf[i+1] - xf[i];

/*
  // debug -- print xc
  printf("--xc--\n");
  for(i=0; i<nx; i++)
    printf("%d %lf\n", i, xc[i]);

//  // debug -- print xf
  printf("--xf--\n");
  for(i=0; i<nx+1; i++)
    printf("%d %lf\n", i, xf[i]);

  // debug -- print dxc
  printf("--dxc--\n");
  for(i=0; i<nx+1; i++)
    printf("%d %lf\n", i, dxc[i]);

//  // debug -- print dxf
  printf("--dxf--\n");
  for(i=0; i<nx; i++)
    printf("%d %lf\n", i, dxf[i]);     */
}

void set_initial_guess(int nx, int ny, double **u, double **v, double **p)
{
  int i, j;

  // set u guess
  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      u[i][j] = 0.0;

  // set v guess
  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      v[i][j] = 0.0;

// set p guess
  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      p[i][j] = 1.0;
}

void calc_diffusivity(int nx, int ny, double **kdiff, double *kleft, double *krght, double *ktop, double *kbot)
{
  int i, j;

  // calculate diffusivity (may be dependent on T)
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
     kdiff[i][j] = k_const; // set to uniform for now

 // left boundary -- Non-homogeneous Dirichlet
 for(j=0; j<ny; j++)
	kleft[j] = k_const;
	
 // right boundary -- Homogeneous Dirichlet
 for(j=0; j<ny; j++)
	krght[j] = k_const;

 // top boundary -- Non-homogeneous Dirichlet
 for(i=0; i<nx; i++)
	ktop[i] = k_const;

 // bottom boundary -- Homogeneous Dirichlet
 for(i=0; i<nx; i++)
	kbot[i] = k_const;

}

void calc_sources_u(int nx, int ny, double *dxf, double *dyf, double **p, double **bu, double **Spu)
{
  int i, j;
   
    // calculate Sp source
  for(i=0; i<nx; i++) {
   for(j=0; j<ny; j++)
   {
     // set the source and multiply by volume of CV
     Spu[i][j] = 0.0 * dxf[i]*dyf[j];
     // set the source and multiply by volume of CV
     bu[i][j] = 0.0 * dxf[i]*dyf[j];
   }}
   
  //Pressures - interior
  for(i=1; i<nx-1; i++)
   for(j=0; j<ny; j++)
     bu[i][j] = (p[i-1][j] - p[i+1][j]) * dyf[j] *0.5;
  
  // Left boundary
  i = 0;
      for(j=0; j<ny; j++) 
         bu[i][j] = (p[0][j] - p[i+1][j]) * dyf[j] * 0.5;
         
  //right boundaries      
  i = nx-1;     
      for(j=0; j<ny; j++)    
         bu[i][j] = (p[i-1][j] - p[nx-1][j]) * dyf[j] * 0.5;   
 
}

void calc_sources_v(int nx, int ny, double *dxf, double *dyf, double **p, double **bv, double **Spv)
{
  int i, j;

    // calculate Sp source
  for(i=0; i<nx; i++) {
   for(j=0; j<ny; j++)
   {
     // set the source and multiply by volume of CV
     Spv[i][j] = 0.0 * dxf[i]*dyf[j];
     // set the source and multiply by volume of CV
     bv[i][j] = 0.0 * dxf[i]*dyf[j];
   } }
   
  //Pressures - interior
  for(i=0; i<nx; i++)
   for(j=1; j<ny-1; j++)
     bv[i][j] = (p[i][j-1] - p[i][j+1]) * dxf[i] * 0.5;
  
  // Top boundary
  j = ny-1;
      for(i=0; i<nx; i++) 
         bv[i][j] = (p[i][j-1] - p[i][ny-1]) * dxf[i] * 0.5;
         
  //Bot boundaries      
  j = 0;     
      for(i=0; i<nx; i++)     
         bv[i][j] = (p[i][0] - p[i][j+1]) * dxf[i] * 0.5; 
}


void set_boundary_conditions(int nx, int ny, double *Fuleft, double *Furght, double *Fvtop, double *Fvbot)
{
  int i, j;
  // left and right  boundary -- stationary wall

  for(j=0; j<ny; j++)  {	
		Fuleft[j] = 0.0;
  	Furght[j] = 0.0; 
	}

  // top and bot  boundary -- moving lid
  for(i=0; i<nx; i++) {	
    Fvtop[i] = 0.0;
    Fvbot[i] = 0.0;
  }
}

void get_u_coeffs(int nx, int ny, double *xc, double *yc, double *dxc, double *dxf, double *dyc, double *dyf, double **aPu, double **aEu, double **aWu, double **aNu, double **aSu, double **bu, double **Spu, double **kdiff, double **rho, double **u, double **v, double **p, double *Fuleft, double *Furght, double *Fvtop, double *Fvbot,double *kleft, double *krght, double *ktop, double *kbot, double **uface, double **vface, int convec_scheme, double u_bc_e, double u_bc_w, double u_bc_n, double u_bc_s,double v_bc_e, double v_bc_w,double v_bc_n, double v_bc_s, double **du)                                                                      
{
  int i, j;
  double bcE, bcW, bcN, bcS, kw, ke, kn, ks, ifac, Fe, Fw, Fn, Fs;
  double De, Dw, Dn, Ds, rhoe, rhow, rhon, rhos, ue, uw, vn, vs;

  calc_diffusivity(nx, ny, kdiff, kleft, krght, ktop, kbot);
  set_boundary_conditions(nx, ny, Fuleft, Furght, Fvtop, Fvbot);
  calc_sources_u(nx, ny, dxf, dyf, p, bu, Spu);

  // start populating the coefficients

  // ------ Step 1 :: interior points ------
  for(i=1; i<nx-1; i++) {
   for(j=1; j<ny-1; j++)
   {
      // east-west -- diffusion
      ke = k_const;   De = ke * dyf[j] / dxc[i+1];
      kw = k_const;   Dw = kw * dyf[j] / dxc[i];

      // north-south -- diffusion
      kn = k_const;   Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;   Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = rho_const;      Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;      Fw = rhow * uface[i][j] * dyf[j];

      // north-south -- convection
      rhon = rho_const;     Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;     Fs = rhos * vface[i][j] * dxf[i];

      //#include "convective_scheme_coeffs_u.h"
      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spu, convec_scheme, aPu, aEu, aWu, aNu, aSu);

    //printf("%d %d %lf %lf %lf %lf \n", i, j, ke, De, dyc[j], dxc[i]);
    //printf("%d %d %lf %lf %lf %lf \n", i, j, kw, Dw, dyc[j], dxc[i-1]);
    //printf("%d %d %lf %lf %lf %lf \n", i, j, kn, Dn, dyf[j], dxf[i]);
    //printf("%d %d %lf %lf %lf %lf \n", i, j, ks, Ds, dyf[j-1], dxf[i]);
    //printf("%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", i, j, aPu[i][j], aEu[i][j], aWu[i][j], aNu[i][j], aSu[i][j], bu[i][j], De, Dw, Dn, Ds, Fe, Fw, Fn, Fs);

      // source --  has already been populated
   } }

  // ------ Step 2  :: 1-boundary points ------
  //// ------ Step 2a :: left boundary ----------
  i = 0;
   for(j=1; j<ny-1; j++)
   {	
      // east-west -- diffusion
      ke = k_const;   De = ke * dyf[j] / dxc[i+1];
      kw = k_const;	  Dw = 0.0;

      // north-south -- diffusion
      kn = k_const;   Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;   Ds = ks * dxf[i] / dyc[j];
      
      // east-west -- convection
      rhoe = rho_const;    Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;    Fw = 0.0;

      // north-south -- convection
      rhon = rho_const;     Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;     Fs = rhos * vface[i][j] * dxf[i];

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spu, convec_scheme, aPu, aEu, aWu, aNu, aSu);

     if(convec_scheme==1)  // UW
       ifac = -1.0;
     else if (convec_scheme==2) // CD
       ifac = -1.0;
       
       bcW = (kleft[j]/dxc[i] - ifac * Fuleft[j]) * dyf[j];
       bu[i][j]  += bcW * u_bc_w;        
       Spu[i][j] = -1*bcW;
       aPu[i][j] -= Spu[i][j]; 
   }
  //// ------ Step 2a :: left boundary done ---

  // ------ Step 2b :: right boundary ----------
  i = nx-1;
   for(j=1; j<ny-1; j++)
   {
      // east-west
      ke = k_const;   De = 0.0;
      kw = k_const;   Dw = kw * dyf[j] / dxc[i];    // since size of dxc is nx

      // north-south
      kn = k_const;   Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;   Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = rho_const;   Fe = 0.0;
      rhow = rho_const;   Fw = rhow * uface[i][j] * dyf[j];

      // north-south -- convection
      rhon = rho_const;     Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;     Fs = rhos * vface[i][j] * dxf[i];

     convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spu, convec_scheme, aPu, aEu, aWu, aNu, aSu);

     if(convec_scheme==1)  // UW
       ifac = 0.0;
     else if (convec_scheme==2) // CD
       ifac = 1.0;
     
       bcE = (krght[j]/dxc[i+1] - ifac * Furght[j]) * dyf[j];
       bu[i][j]  += bcE * u_bc_e;
       Spu[i][j] = -1*bcE;
       aPu[i][j] -= Spu[i][j]; 
   }
  //// ------ Step 2b :: right boundary done ---

  // ------ Step 2c :: top boundary ----------
  j = ny-1;
  for(i=1; i<nx-1; i++)
   {
      // east-west
      ke = k_const;   De = ke * dyf[j] / dxc[i+1];
      kw = k_const;   Dw = kw * dyf[j] / dxc[i];

      // north-south
      kn = k_const;   Dn = 0.0;
      ks = k_const;   Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = rho_const;     Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;     Fw = rhow * uface[i][j] * dyf[j];

      // north-south -- convection
      rhon = rho_const;      Fn = 0.0;
      rhos = rho_const;      Fs = rhos * vface[i][j] * dxf[i];

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spu, convec_scheme, aPu, aEu, aWu, aNu, aSu);

      if(convec_scheme==1)  // UW
       ifac = 0.0;
      else if (convec_scheme==2) // CD
       ifac = 1.0;
       
       bcN = (ktop[i]/dyc[j+1] - ifac * Fvtop[i]) * dxf[i];
       bu[i][j]  += bcN * u_bc_n;   
       Spu[i][j] = -1*bcN;
       aPu[i][j] -= Spu[i][j]; 
   }
  // ------ Step 2c :: top boundary done ---

  // ------ Step 2d :: bottom boundary ----------
  j = 0;
  for(i=1; i<nx-1; i++)
   {
      // east-west
      ke = k_const;   De = ke * dyf[j] / dxc[i+1];
      kw = k_const;   Dw = kw * dyf[j] / dxc[i];
    
      // north-south
      kn = k_const;   Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;   Ds = 0.0;
    
      // east-west -- convection
      rhoe = rho_const;     Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;     Fw = rhow * uface[i][j] * dyf[j];
    
      // north-south -- convection
      rhon = rho_const;     Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;     Fs = 0.0;
    
      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spu, convec_scheme, aPu, aEu, aWu, aNu, aSu);

      if(convec_scheme==1)  // UW
       ifac = -1.0;
      else if (convec_scheme==2) // CD
       ifac = -1.0;
       
       bcS = (kbot[i]/dyc[j] - ifac * Fvbot[i]) * dxf[i];
       bu[i][j]  += bcS * u_bc_s;         
       Spu[i][j] = -1*bcS;
       aPu[i][j] -= Spu[i][j]; 
        
   }
  // ------ Step 2d :: bottom boundary done ---

  //// ------ Step 3  :: 2-boundary points ------
  //// ------ Step 3a :: bottom-left boundary ----------
  i = 0; j = 0;
      // east-west
      ke = k_const;   De = ke * dyf[j] / dxc[i+1];
      kw = k_const;   Dw = 0.0;

      // north-south
      kn = k_const;   Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;   Ds = 0.0;

      // east-west -- convection
      rhoe = rho_const ;    Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;     Fw = 0.0;

      // north-south -- convection
      rhon = rho_const;       Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;       Fs = 0.0;

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spu, convec_scheme, aPu, aEu, aWu, aNu, aSu);

      if(convec_scheme==1)  // UW
        ifac = -1.0;
      else if (convec_scheme==2) // CD
        ifac = -1.0;
      bcW = (kleft[j]/dxc[i] - ifac * Fuleft[j]) * dyf[j];
      bu[i][j]  += bcW * u_bc_w;   
      Spu[i][j] = -1*bcW;
      aPu[i][j] -= Spu[i][j];     
      
      if(convec_scheme==1)  // UW
       ifac = -1.0;
      else if (convec_scheme==2) // CD
       ifac = -1.0;
       bcS = (kbot[i]/dyc[j] - ifac * Fvbot[i]) * dxf[i];
       bu[i][j]  += bcS * u_bc_s;         
       Spu[i][j] = -1*bcS;
       aPu[i][j] -= Spu[i][j]; 
      
  // ------ Step 3a :: bottom-left boundary done ---

  // ------ Step 3b :: bottom-right boundary ----------
  j = 0; i = nx-1;
      // east-west
      ke = k_const;       De = 0.0;
      kw = k_const;       Dw = kw * dyf[j] / dxc[i];

      // north-south
      kn = k_const;       Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;       Ds = 0.0; 
      // east-west -- convection
      rhoe = rho_const;    Fe = 0.0;
      rhow = rho_const;    Fw = rhow * uface[i][j] * dyf[j];

      // north-south -- convection
      rhon = rho_const;    Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;    Fs = 0.0;

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spu, convec_scheme, aPu, aEu, aWu, aNu, aSu);

      if(convec_scheme==1)  // UW
        ifac = 0.0;
      else if (convec_scheme==2) // CD
        ifac = 1.0;
       bcE = (krght[j]/dxc[i+1] - ifac * Furght[j]) * dyf[j];
       bu[i][j]  += bcE * u_bc_e;        
       Spu[i][j] = -1*bcE;
       aPu[i][j] -= Spu[i][j];;

       if(convec_scheme==1)  // UW
       ifac = -1.0;
      else if (convec_scheme==2) // CD
       ifac = -1.0;
       bcS = (kbot[i]/dyc[j] - ifac * Fvbot[i]) * dxf[i];
       bu[i][j]  += bcS * u_bc_s;         
       Spu[i][j] = -1*bcS;
       aPu[i][j] -= Spu[i][j]; 

  // ------ Step 3b :: bottom-right boundary done ---

  //// ------ Step 3c :: top-left boundary ----------
  j = ny-1; i = 0;
      // east-west
      ke = k_const;       De = ke * dyf[j] / dxc[i+1];
      kw =  k_const;      Dw = 0.0;

      // north-south
      kn =  k_const;      Dn = 0.0;
      ks =  k_const;      Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = rho_const;    Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;    Fw = 0.0;

      // north-south -- convection
      rhon = rho_const;    Fn = 0.0;
      rhos = rho_const;    Fs = rhos * vface[i][j] * dxf[i];

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spu, convec_scheme, aPu, aEu, aWu, aNu, aSu);

      if(convec_scheme==1)  // UW
       ifac = 0.0;
      else if (convec_scheme==2) // CD
       ifac = 1.0;
       bcN = (ktop[i]/dyc[j+1] - ifac * Fvtop[i]) * dxf[i];
       bu[i][j]  += bcN * u_bc_n;   
       Spu[i][j] = -1*bcN;
       aPu[i][j] -= Spu[i][j]; 

      if(convec_scheme==1)  // UW
        ifac = -1.0;
      else if (convec_scheme==2) // CD
        ifac = -1.0;
       bcW = (kleft[j]/dxc[i] - ifac * Fuleft[j]) * dyf[j];
       bu[i][j]  += bcW * u_bc_w;        
       Spu[i][j] = -1*bcW;
       aPu[i][j] -= Spu[i][j];

  //     //printf("%d %d %lf %lf %lf %lf %lf\n", i, j, kleft[j], dxc[i], Fuleft[j], dyf[j], bcW);
  //     //printf("%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", i, j, aE[i][j], aW[i][j], aN[i][j], aS[i][j], Fe, Fw, Fn, Fs, bcW, bcS, aP[i][j], bu[i][j]);

  //// ------ Step 3c :: top-left boundary done ---

  //// ------ Step 3d :: top-right boundary ----------
  j = ny-1;  i = nx-1;
      // east-west
      ke = k_const;       De = 0.0;
      kw = k_const;       Dw = kw * dyf[j] / dxc[i];

      // north-south
      kn = k_const;       Dn = 0.0;
      ks = k_const;       Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = rho_const;   Fe = 0.0;
      rhow = rho_const;   Fw = rhow * uface[i][j] * dyf[j];

      // north-south -- convection
      rhon = rho_const;   Fn = 0.0;
      rhos = rho_const;   Fs = rhos * vface[i][j] * dxf[i];

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spu, convec_scheme, aPu, aEu, aWu, aNu, aSu);

      if(convec_scheme==1)  // UW
       ifac = 0.0;
      else if (convec_scheme==2) // CD
       ifac = 1.0;
  
       bcN = (ktop[i]/dyc[j+1] - ifac * Fvtop[i]) * dxf[i];
       bu[i][j]  += bcN * u_bc_n;   
       Spu[i][j] = -1*bcN;
       aPu[i][j] -= Spu[i][j]; 

      if(convec_scheme==1)  // UW
        ifac = 0.0;
      else if (convec_scheme==2) // CD
        ifac = 1.0;
       bcE = (krght[j]/dxc[i+1] - ifac * Furght[j]) * dyf[j];
       bu[i][j]  += bcE * u_bc_e;        
       Spu[i][j] = -1*bcE;
       aPu[i][j] -= Spu[i][j];

  //// ------ Step 3d :: top-right boundary done ---
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
      du[i][j] = dyf[j]/aPu[i][j];
  //// debug
//  for(i=0; i<nx; i++)
//  for(j=0; j<ny; j++)
//   {
//    printf("%d %d %lf %lf %lf %lf %lf %lf\n", i, j, aPu[i][j], aEu[i][j], aWu[i][j], aNu[i][j], aSu[i][j], bu[i][j]);
//  }
}

void get_v_coeffs(int nx, int ny, double *xc, double *yc, double *dxc, double *dxf, double *dyc, double *dyf, double **aPv, double **aEv, double **aWv, double **aNv, double **aSv, double **bv, double **Spv, double **kdiff, double **rho, double **u, double **v, double **p, double *Fuleft, double *Furght, double *Fvtop, double *Fvbot, double *kleft, double *krght, double *ktop, double *kbot, double **uface, double **vface, int convec_scheme, double u_bc_e, double u_bc_w, double u_bc_n, double u_bc_s,double v_bc_e, double v_bc_w,double v_bc_n, double v_bc_s, double **dv)
{

  int i, j;
  double bcE, bcW, bcN, bcS, kw, ke, kn, ks, ifac, Fe, Fw, Fn, Fs;
  double De, Dw, Dn, Ds, rhoe, rhow, rhon, rhos, ue, uw, vn, vs;

  // calculate diffusivity at [xc, yc]; may be dependent on T
  calc_diffusivity(nx, ny, kdiff, kleft, krght, ktop, kbot);
  set_boundary_conditions(nx, ny, Fuleft, Furght, Fvtop, Fvbot);
  calc_sources_v(nx, ny, dxf, dyf, p, bv, Spv);
           
  // start populating the coefficients

  // ------ Step 1 :: interior points ------
  for(i=1; i<nx-1; i++)
   for(j=1; j<ny-1; j++)
   {
      ke = k_const;   De = ke * dyf[j] / dxc[i+1];
      kw = k_const;   Dw = kw * dyf[j] / dxc[i];

      // north-south -- diffusion
      kn = k_const;   Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;   Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = rho_const;     Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;     Fw = rhow * uface[i][j] * dyf[j];

      // north-south -- convection
      rhon = rho_const;     Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;     Fs = rhos * vface[i][j] * dxf[i];	

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spv, convec_scheme, aPv, aEv, aWv, aNv, aSv);

    //printf("%d %d %lf %lf %lf %lf \n", i, j, ke, De, dyc[j], dxc[i]);
    //printf("%d %d %lf %lf %lf %lf \n", i, j, kw, Dw, dyc[j], dxc[i-1]);
    //printf("%d %d %lf %lf %lf %lf \n", i, j, kn, Dn, dyf[j], dxf[i]);
    //printf("%d %d %lf %lf %lf %lf \n", i, j, ks, Ds, dyf[j-1], dxf[i]);
    //printf("%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", i, j, aP[i][j], aE[i][j], aW[i][j], aN[i][j], aS[i][j], bu[i][j], De, Dw, Dn, Ds, Fe, Fw, Fn, Fs);
      // source --  has already been populated
   }

  // ------ Step 2  :: 1-boundary points ------
  // ------ Step 2a :: left boundary ----------
  i = 0;
   for(j=1; j<ny-1; j++)
   {
      // east-west
      ke = k_const;       De = ke * dyf[j] / dxc[i+1];
      kw = k_const;       Dw = 0.0;

      // north-south
      kn =k_const;    Dn = kn * dxf[i] / dyc[j+1];
      ks =k_const;    Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = rho_const;           Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;           Fw = 0.0;

      // north-south -- convection
      rhon = rho_const;     Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;     Fs = rhos * vface[i][j] * dxf[i];

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spv, convec_scheme, aPv, aEv, aWv, aNv, aSv);

      if(convec_scheme==1)  // UW
        ifac = -1.0;
      else if (convec_scheme==2) // CD
        ifac = -1.0;
       bcW = (kleft[j]/dxc[i] - ifac * Fuleft[j]) * dyf[j];
       bv[i][j]  += bcW * v_bc_w;        
       Spv[i][j] = -1*bcW;
       aPv[i][j] -= Spv[i][j]; 
   }
  // ------ Step 2a :: left boundary done ---

  // ------ Step 2b :: right boundary ----------
  i = nx-1;
   for(j=1; j<ny-1; j++)
   {
      // east-west
      ke = k_const;   De = 0.0;
      kw = k_const;   Dw = kw * dyf[j] / dxc[i];

      // north-south
      kn = k_const;   Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;   Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = rho_const;     Fe = 0.0;
      rhow = rho_const;     Fw = rhow * uface[i][j] * dyf[j];
      
      // north-south -- convection
      rhon = rho_const;     Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;     Fs = rhos * vface[i][j] * dxf[i];

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spv, convec_scheme, aPv, aEv, aWv, aNv, aSv);

       if(convec_scheme==1)  // UW
        ifac = 0.0;
       else if (convec_scheme==2) // CD
        ifac = 1.0;
       bcE = (krght[j]/dxc[i+1] - ifac * Furght[j]) * dyf[j];
       bv[i][j]  += bcE * v_bc_e;        
       Spv[i][j] = -1*bcE;
       aPv[i][j] -= Spv[i][j]; 
   }
  // ------ Step 2b :: right boundary done ---

  //// ------ Step 2c :: top boundary ----------
  j = ny-1;
  for(i=1; i<nx-1; i++)
   {
      // east-west
      ke = k_const;   De = ke * dyf[j] / dxc[i+1];
      kw = k_const;   Dw = kw * dyf[j] / dxc[i];
    
      // north-south
      kn = k_const;   Dn = 0.0;
      ks = k_const;       Ds = ks * dxf[i] / dyc[j];
    
      // east-west -- convection
      rhoe = rho_const;     Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;    Fw = rhow * uface[i][j] * dyf[j];
    
      // north-south -- convection
      rhon = rho_const;     Fn = 0.0;
      rhos = rho_const;           Fs = rhos * vface[i][j] * dxf[i];
    
      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spv, convec_scheme, aPv, aEv, aWv, aNv, aSv);

       if(convec_scheme==1)  // UW
        ifac = 0.0;
       else if (convec_scheme==2) // CD
        ifac = 1.0;
       bcN = (ktop[i]/dyc[j+1] - ifac * Fvtop[i]) * dxf[i];
       bv[i][j]  += bcN * v_bc_n;        
       Spv[i][j] = -1*bcN;
       aPv[i][j] -= Spv[i][j]; 
   }
  //// ------ Step 2c :: top boundary done ---
    
  //// ------ Step 2d :: bottom boundary ----------
  j = 0;
  for(i=1; i<nx-1; i++)
   {
      // east-west
      ke = k_const;   De = ke * dyf[j] / dxc[i+1];
      kw = k_const;   Dw = kw * dyf[j] / dxc[i];
    
      // north-south
      kn = k_const;       Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;       Ds = 0.0;
    
      // east-west -- convection
      rhoe = rho_const;     Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;     Fw = rhow * uface[i][j] * dyf[j];
    
      // north-south -- convection
      rhon = rho_const;           Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;     Fs = 0.0;
    
      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spv, convec_scheme, aPv, aEv, aWv, aNv, aSv);

       if(convec_scheme==1)  // UW
        ifac = -1.0;
       else if (convec_scheme==2) // CD
        ifac = -1.0;
       bcS = (kbot[i]/dyc[j] - ifac * Fvbot[i]) * dxf[i];
       bv[i][j]  += bcS * v_bc_s;
       Spv[i][j] = -1*bcS;
       aPv[i][j] -= Spv[i][j]; 
   }
  //// ------ Step 2d :: top boundary done ---

  //// ------ Step 3  :: 2-boundary points ------
  //// ------ Step 3a :: bottom-left boundary ----------
  i = 0; j = 0;
      // east-west
      ke = k_const;       De = ke * dyf[j] / dxc[i+1];
      kw = k_const;       Dw = 0.0;

      // north-south
      kn = k_const;       Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;       Ds = 0.0;

      // east-west -- convection
      rhoe = rho_const;     Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;     Fw = 0.0;

      // north-south -- convection
      rhon = rho_const;     Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;     Fs = 0.0;

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spv, convec_scheme, aPv, aEv, aWv, aNv, aSv);

       if(convec_scheme==1)  // UW
        ifac = -1.0;
       else if (convec_scheme==2) // CD
        ifac = -1.0;
       bcW = (kleft[j]/dxc[i] - ifac * Fuleft[j]) * dyf[j];
       bv[i][j]  += bcW * v_bc_w;        
       Spv[i][j] = -1*bcW;
       aPv[i][j] -= Spv[i][j]; 

       if(convec_scheme==1)  // UW
        ifac = -1.0;
       else if (convec_scheme==2) // CD
        ifac = -1.0;
       bcS = (kbot[i]/dyc[j] - ifac * Fvbot[i]) * dxf[i];
       bv[i][j]  += bcS * v_bc_s;
       Spv[i][j] = -1*bcS;
       aPv[i][j] -= Spv[i][j];  

  //// ------ Step 3a :: bot-left boundary done ---

  //// ------ Step 3b :: bot-right boundary ----------
  j = 0; i = nx-1;
      // east-west
      ke = k_const;         De = 0.0;
      kw = k_const;         Dw = kw * dyf[j] / dxc[i];

      // north-south
      kn = k_const;         Dn = kn * dxf[i] / dyc[j+1];
      ks = k_const;         Ds = 0.0;

      // east-west -- convection
      rhoe = rho_const;     Fe = 0.0;
      rhow = rho_const;     Fw = rhow * uface[i][j] * dyf[j];

      // north-south -- convection
      rhon = rho_const;     Fn = rhon * vface[i][j+1] * dxf[i];
      rhos = rho_const;     Fs = 0.0;

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spv, convec_scheme, aPv, aEv, aWv, aNv, aSv);

       if(convec_scheme==1)  // UW
        ifac = 0.0;
       else if (convec_scheme==2) // CD
        ifac = 1.0;
       bcE = (krght[j]/dxc[i+1] - ifac * Furght[j]) * dyf[j];
       bv[i][j]  += bcE * v_bc_e;        
       Spv[i][j] = -1*bcE;
       aPv[i][j] -= Spv[i][j]; 

       if(convec_scheme==1)  // UW
        ifac = -1.0;
       else if (convec_scheme==2) // CD
        ifac = -1.0;
       bcS = (kbot[i]/dyc[j] - ifac * Fvbot[i]) * dxf[i];
       bv[i][j]  += bcS * v_bc_s;
       Spv[i][j] = -1*bcS;
       aPv[i][j] -= Spv[i][j]; 

  //// ------ Step 3b :: bot-right boundary done ---

  //// ------ Step 3c :: top-left boundary ----------
  j = ny-1; i = 0;
      // east-west
      ke = k_const;       De = ke * dyf[j] / dxc[i+1];
      kw = k_const;       Dw = 0.0;

      // north-south
      kn = k_const;       Dn = 0.0;
      ks = k_const;       Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = rho_const;     Fe = rhoe * uface[i+1][j] * dyf[j];
      rhow = rho_const;     Fw = 0.0;

      // north-south -- convection
      rhon = rho_const;           Fn = 0.0;
      rhos = rho_const;           Fs = rhos * vface[i][j] * dxf[i];

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spv, convec_scheme, aPv, aEv, aWv, aNv, aSv);

       if(convec_scheme==1)  // UW
        ifac = 0.0;
       else if (convec_scheme==2) // CD
        ifac = 1.0;
       bcN = (ktop[i]/dyc[j+1] - ifac * Fvtop[i]) * dxf[i];
       bv[i][j]  += bcN * v_bc_n;        
       Spv[i][j] = -1*bcN;
       aPv[i][j] -= Spv[i][j]; 

       if(convec_scheme==1)  // UW
        ifac = -1.0;
       else if (convec_scheme==2) // CD
        ifac = -1.0;
       bcW = (kleft[j]/dxc[i] - ifac * Fuleft[j]) * dyf[j];
       bv[i][j]  += bcW * v_bc_w;        
       Spv[i][j] = -1*bcW;
       aPv[i][j] -= Spv[i][j]; 

  //     //printf("%d %d %lf %lf %lf %lf %lf\n", i, j, kleft[j], dxc[i], Fvleft[j], dyf[j], bcW);
  //     //printf("%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", i, j, aE[i][j], aW[i][j], aN[i][j], aS[i][j], Fe, Fw, Fn, Fs, bcW, bcS, aP[i][j], bu[i][j]);

  //// ------ Step 3c :: top-left boundary done ---

  //// ------ Step 3d :: top-right boundary ----------
  j = ny-1;  i = nx-1;
      // east-west
      ke = k_const;       De = 0.0;
      kw = k_const;       Dw = kw * dyf[j] / dxc[i];

      // north-south
      kn = k_const;       Dn = 0.0;
      ks = k_const;       Ds = ks * dxf[i] / dyc[j];

      // east-west -- convection
      rhoe = rho_const;     Fe = 0.0;
      rhow = rho_const;     Fw = rhow * uface[i][j] * dyf[j];


      // north-south -- convection
      rhon = rho_const;      Fn = 0.0;
      rhos = rho_const;      Fs = rhos * vface[i][j] * dxf[i];

      convective_scheme_coeffs(i, j, De, Dw, Dn, Ds, Fe, Fw, Fn, Fs, Spv, convec_scheme, aPv, aEv, aWv, aNv, aSv);

       if(convec_scheme==1)  // UW
        ifac = 0.0;
       else if (convec_scheme==2) // CD
        ifac = 1.0;
       bcN = (ktop[i]/dyc[j+1] - ifac * Fvtop[i]) * dxf[i];
       bv[i][j]  += bcN * v_bc_n;        
       Spv[i][j] = -1*bcN;
       aPv[i][j] -= Spv[i][j];  

       if(convec_scheme==1)  // UW
        ifac = 0.0;
       else if (convec_scheme==2) // CD
        ifac = 1.0;
       bcE = (krght[j]/dxc[i+1] - ifac * Furght[j]) * dyf[j];
       bv[i][j]  += bcE * v_bc_e;        
       Spv[i][j] = -1*bcE;
       aPv[i][j] -= Spv[i][j]; 
  //// ------ Step 3d :: bottom-right boundary done ---


  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
      dv[i][j] = dxf[i]/aPv[i][j];


  //// debug
  //printf("/n V coeff");
//  for(i=0; i<nx; i++)
//   for(j=0; j<ny; j++)
//   {
//    printf("%d %d %lf %lf %lf %lf %lf %lf\n", i, j, aPv[i][j], aEv[i][j], aWv[i][j], aNv[i][j], aSv[i][j], bv[i][j]);
//  }
}


void get_p_coeffs(int nx, int ny, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **aPu, double **aPv, double **du, double **dv, double *dxf, double *dyf, double **uface, double **vface)
{

  int i, j;
  
  //// start populating the coefficients

  // ------ Step 1 :: interior ------
  for(i=1; i<nx-1; i++)
   for(j=1; j<ny-1; j++)
   {
      aE[i][j] = (du[i][j] + du[i+1][j]) * dyf[j] * 0.5;      aW[i][j] = (du[i][j] + du[i-1][j]) * dyf[j] * 0.5;
      aN[i][j] = (dv[i][j] + dv[i][j+1]) * dxf[i] * 0.5;      aS[i][j] = (dv[i][j] + dv[i][j-1]) * dxf[i] * 0.5;
  
      b[i][j] = -1.0*((uface[i+1][j] - uface[i][j])* dyf[j] + (vface[i][j+1] - vface[i][j])* dxf[i]);      // Linear interpolation for vel* terms
    }
    
    //top left corner
    j = ny-1; i = 0;
      aE[i][j] = (du[i][j] + du[i+1][j]) * dyf[j] * 0.5;      aW[i][j] = 0.0;
      aN[i][j] = 0.0;                                         aS[i][j] = (dv[i][j] + dv[i][j-1]) * dxf[i] * 0.5;
  
      b[i][j] = -1.0*((uface[i+1][j] - uface[i][j])* dyf[j] + (vface[i][j+1] - vface[i][j])* dxf[i]); 
    
    //bot left corner (to be replaced by 1.0 -- Nuemann BC)
    i = 0; j = 0;
      aE[i][j] = (du[i][j] + du[i+1][j]) * dyf[j] * 0.5;      aW[i][j] = 0.0;
      aN[i][j] = (dv[i][j] + dv[i][j+1]) * dxf[i] * 0.5;      aS[i][j] = 0.0;
  
      b[i][j] = -1.0*((uface[i+1][j] - uface[i][j])* dyf[j] + (vface[i][j+1] - vface[i][j])* dxf[i]);   
    
    // bot right corner
    j = 0; i = nx-1;
      aE[i][j] = 0.0;                                         aW[i][j] = (du[i][j] + du[i-1][j]) * dyf[j] * 0.5;
      aN[i][j] = (dv[i][j] + dv[i][j+1]) * dxf[i] * 0.5;      aS[i][j] = 0.0;
  
      b[i][j] = -1.0*((uface[i+1][j] - uface[i][j])* dyf[j] + (vface[i][j+1] - vface[i][j])* dxf[i]);   
    
    // top right corner
    j = ny-1; i = nx-1;
      aE[i][j] = 0.0;                                         aW[i][j] = (du[i][j] + du[i-1][j]) * dyf[j] * 0.5;
      aN[i][j] = 0.0;                                         aS[i][j] = (dv[i][j] + dv[i][j-1]) * dxf[i] * 0.5;
  
      b[i][j] = -1.0*((uface[i+1][j] - uface[i][j])* dyf[j] + (vface[i][j+1] - vface[i][j])* dxf[i]);  

    //left boundaries
    i = 0;
      for(j=1; j<ny-1; j++) 
      {  
      aE[i][j] = (du[i][j] + du[i+1][j]) * dyf[j] * 0.5;      aW[i][j] = 0.0;
      aN[i][j] = (dv[i][j] + dv[i][j+1]) * dxf[i] * 0.5;      aS[i][j] = (dv[i][j] + dv[i][j-1]) * dxf[i] * 0.5;
  
      b[i][j] = -1.0*((uface[i+1][j] - uface[i][j])* dyf[j] + (vface[i][j+1] - vface[i][j])* dxf[i]); 
            
       }     

    //right boundaries      
    i = nx-1;     
      for(j=1; j<ny-1; j++)
      {
      aE[i][j] = 0.0;                                         aW[i][j] = (du[i][j] + du[i-1][j]) * dyf[j] * 0.5;
      aN[i][j] = (dv[i][j] + dv[i][j+1]) * dxf[i] * 0.5;      aS[i][j] = (dv[i][j] + dv[i][j-1]) * dxf[i] * 0.5;
  
      b[i][j] = -1.0*((uface[i+1][j] - uface[i][j])* dyf[j] + (vface[i][j+1] - vface[i][j])* dxf[i]); 
      }
            
    //top boundaries
    j = ny-1;            
      for(i=1; i<nx-1; i++)
      {
      aE[i][j] = (du[i][j] + du[i+1][j]) * dyf[j] * 0.5;      aW[i][j] = (du[i][j] + du[i-1][j]) * dyf[j] * 0.5;
      aN[i][j] = 0.0;                                         aS[i][j] = (dv[i][j] + dv[i][j-1]) * dxf[i] * 0.5;
  
      b[i][j] = -1.0*((uface[i+1][j] - uface[i][j])* dyf[j] + (vface[i][j+1] - vface[i][j])* dxf[i]); 
      }
      
    //bot boundaries 
    j = 0;           
      for(i=1; i<nx-1; i++)
      {
      aE[i][j] = (du[i][j] + du[i+1][j]) * dyf[j] * 0.5;      aW[i][j] = (du[i][j] + du[i-1][j]) * dyf[j] * 0.5;
      aN[i][j] = (dv[i][j] + dv[i][j+1]) * dxf[i] * 0.5;      aS[i][j] = 0.0;
  
      b[i][j] = -1.0*((uface[i+1][j] - uface[i][j])* dyf[j] + (vface[i][j+1] - vface[i][j])* dxf[i]); 
      }   

  // ------ Step 3 :: aP for all points ------
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
      aP[i][j] = aE[i][j] + aW[i][j] + aN[i][j] + aS[i][j];

  // ------ Step 4 :: Neumann conditions on all boundaries makes the matrix singular ------
  // ------           Set one reference value ------------
  aP[0][0] = 1.0; aE[0][0] = 0.0; aW[0][0] = 0.0; aN[0][0] = 0.0; aS[0][0] = 0.0; b[0][0] = 0.0;

  //// debug
//  for(i=0; i<nx; i++)
//   for(j=0; j<ny; j++)
//   {
//    printf("%d %d %lf %lf %lf %lf %lf %lf\n", i, j, aP[i][j], aE[i][j], aW[i][j], aN[i][j], aS[i][j], b[i][j]);
//  }
 
}

void correct_uvp(int nx, int ny, double **du, double **dv, double **ppr, double **ucorr, double **vcorr, double **pst, double relax_p, double **u, double **v, double **pnew)
{
  int i, j;

//Initialize
  for(i=0; i<nx; i++) {
   for(j=0; j<ny; j++) {
     pnew[i][j] = 0.0;
     ucorr[i][j] = 0.0;
     vcorr[i][j] = 0.0;
   }}

  //Pressure correction
  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
     pnew[i][j] = pst[i][j] + relax_p*ppr[i][j];                  // PRESSURE FIELD UPDATED
  
//  fill_ghost_values(nx, ny, pLeft, pRght, pBot, pTop, pst, uLeft, uRght, uBot, uTop, vLeft, vRght, vBot, vTop, u, v, u_bc_e, u_bc_w, u_bc_n, u_bc_s, v_bc_e, v_bc_w, v_bc_n, v_bc_s);

  //correction at interiors
  for(i=1; i<nx-1; i++) {
   for(j=1; j<ny-1; j++)  {
     ucorr[i][j] = 0.5*du[i][j]*(ppr[i-1][j] - ppr[i+1][j]); 
     vcorr[i][j] = 0.5*dv[i][j]*(ppr[i][j-1] - ppr[i][j+1]);
     }  }
  //correction at boundaries
  // left
   i = 0;
      for(j=1; j<ny-1; j++) 
      {
         ucorr[i][j] = 0.5*du[i][j]*(ppr[0][j] - ppr[i+1][j]);
         vcorr[i][j] = 0.5*dv[i][j]*(ppr[i][j-1] - ppr[i][j+1]);
       }
    
  //right boundaries      
    i = nx-1;     
      for(j=1; j<ny-1; j++)
      {     
         ucorr[i][j] = 0.5*du[i][j]*(ppr[i-1][j] - ppr[nx-1][j]);
         vcorr[i][j] = 0.5*dv[i][j]*(ppr[i][j-1] - ppr[i][j+1]);
       }
       
    //top boundaries
    j = ny-1;            
      for(i=1; i<nx-1; i++)
      {
         ucorr[i][j] = 0.5*du[i][j]*(ppr[i-1][j] - ppr[i+1][j]);
         vcorr[i][j] = 0.5*dv[i][j]*(ppr[i][j-1] - ppr[i][ny-1]);
      }
      
     //bot boundaries 
    j = 0;           
      for(i=1; i<nx-1; i++)
      {
       ucorr[i][j] = 0.5*du[i][j]*(ppr[i-1][j] - ppr[i+1][j]);
       vcorr[i][j] = 0.5*dv[i][j]*(ppr[i][0] - ppr[i][j+1]);
      } 
       
   // top left
   j = ny-1; i = 0;
         ucorr[i][j] = 0.5*du[i][j]*(ppr[0][j] - ppr[i+1][j]);
         vcorr[i][j] = 0.5*dv[i][j]*(ppr[i][j-1] - ppr[i][ny-1]);
   
   //bot left
   j = 0; i = 0;
           ucorr[i][j] = 0.5*du[i][j]*(ppr[0][j] - ppr[i+1][j]);
           vcorr[i][j] = 0.5*dv[i][j]*(ppr[i][0] - ppr[i][j+1]);
           
    //top right 
    j = ny-1; i = nx-1; 
         ucorr[i][j] = 0.5*du[i][j]*(ppr[i-1][j] - ppr[nx-1][j]);
          vcorr[i][j] = 0.5*dv[i][j]*(ppr[i][j-1] - ppr[i][ny-1]);
         
     //bot right
     j = 0; i = nx-1;
          ucorr[i][j] = 0.5*du[i][j]*(ppr[i-1][j] - ppr[nx-1][j]);
          vcorr[i][j] = 0.5*dv[i][j]*(ppr[i][0] - ppr[i][j+1]);  
     
}

void velocity_update(int nx, int ny, double **pst, double **u, double **ust, double **ucorr, double **v, double **vst, double **vcorr, double **uface, double **ustface, double **ufacecorr, double **vface, double **vstface, double **vfacecorr, double relax_u, double relax_v, double **unew, double **vnew, double **ufacenew, double **vfacenew)
{
  int i, j;
  
  //Initialize
  for(i=0; i<nx; i++) {
   for(j=0; j<ny; j++) {
    unew[i][j] = 0.0;
    vnew[i][j] = 0.0;
   }}

  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++){
      unew[i][j] = (1-relax_u)*u[i][j] + relax_u*(ust[i][j] + ucorr[i][j]); 
      vnew[i][j] = (1-relax_v)*v[i][j] + relax_v*(vst[i][j] + vcorr[i][j]);
    }
  
  for(i=0; i<nx+1; i++)
    for(j=0; j<ny; j++)
      ufacenew[i][j] = (1-relax_u)*uface[i][j] + relax_u*(ustface[i][j] + ufacecorr[i][j]);
    
  for(i=0; i<nx; i++)
    for(j=0; j<ny+1; j++)
      vfacenew[i][j] = (1-relax_v)*vface[i][j] + relax_v*(vstface[i][j] + vfacecorr[i][j]);
    

}
