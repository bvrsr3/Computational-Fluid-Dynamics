#include "utilities.h"
#include "functions.h"
#include "linsolvers.h"
#include "face_calc.h"

#include<stdlib.h>

int main()
{

  int nx, ny;
  double *xf, *xc, *dxf, *dxc, *yf, *yc, *dyc, *dyf;
  double **aPu, **aEu, **aWu, **aNu, **aSu, **bu, **Spu, **ust, **du;
  double **aPv, **aEv, **aWv, **aNv, **aSv, **bv, **Spv, **vst, **dv;
  double **aPp, **aEp, **aWp, **aNp, **aSp, **bp, **Spp, **pst, **ppr;
  double **T, **kdiff, **rho, **u, **v, **p;
  double xst, xen, yst, yen, Re;
  double *kleft, *krght, *ktop, *kbot;
  double *Fuleft, *Furght;
  double *Fvtop, *Fvbot;
  double **uface, **vface;
  int i, j, max_iter_u,max_iter_v,max_iter_p, max_simple_iters, isim;
  double relax_u, relax_v, relax_p, relax_gssor, tol, l2u, l2v, l2p, l2max;
  double **ucorr, **vcorr, **ufacecorr, **vfacecorr, **ustface, **vstface;
  double u_bc_w, u_bc_e, u_bc_n, u_bc_s, v_bc_w, v_bc_e, v_bc_n, v_bc_s;
  double **unew, **vnew, **pnew, **ufacenew, **vfacenew;

  FILE* fp;  

  // read inputs
  fp = fopen("input.in", "r");
  fscanf(fp, "%d %d\n", &nx, &ny);
  fscanf(fp, "%lf %lf\n", &xst, &xen);
  fscanf(fp, "%lf %lf\n", &yst, &yen);
  fscanf(fp, "%lf\n", &Re);
  fscanf(fp, "%lf %lf %lf %lf\n", &relax_gssor, &relax_p, &relax_u, &relax_v);
  fscanf(fp, "%d %d %d %d %lf\n", &max_simple_iters, &max_iter_u, &max_iter_v, &max_iter_p, &tol);
  fclose(fp);

  printf("Inputs are: %d %d\n %lf %lf\n %lf %lf\n %lf \n%lf %lf %lf %lf\n %d %d %d %d %lf\n", nx, ny, xst, xen, yst, yen, Re, 
                                relax_gssor, relax_p, relax_u, relax_v, max_simple_iters, max_iter_u, max_iter_v, max_iter_p, tol);
  
  k_const = 1.0/Re;  // will be used in functions.c
  rho_const = 1.0;  // will be used in functions.c

  // allocate memory
  printf("\n > Allocating Memory -- \n");
  xc  = (double *)malloc((nx    ) * sizeof(double));   // CV centers (interiors)
  xf  = (double *)malloc((nx + 1) * sizeof(double));   // CV faces
  dxc = (double *)malloc((nx + 1) * sizeof(double));   // spacing betw centers (interiors)
  dxf = (double *)malloc((nx    ) * sizeof(double));   // spacing betw faces

  yc  = (double *)malloc((ny    ) * sizeof(double));   // CV centers (interiors)
  yf  = (double *)malloc((ny + 1) * sizeof(double));   // CV faces
  dyc = (double *)malloc((ny + 1) * sizeof(double));   // spacing betw centers (interiors)
  dyf = (double *)malloc((ny    ) * sizeof(double));   // spacing betw faces

  Fuleft  = (double *)malloc( ny * sizeof(double));   // left   Dirichlet condition
  Furght  = (double *)malloc( ny * sizeof(double));   // right  Dirichlet condition
  
  Fvtop = (double *)malloc( nx * sizeof(double));   // top    Dirichlet condition
  Fvbot  = (double *)malloc( nx * sizeof(double));   // right  Dirichlet condition
  
  kleft  = (double *)malloc( ny * sizeof(double));   // left   Dirichlet condition
  krght  = (double *)malloc( ny * sizeof(double));   // right  Dirichlet condition
  ktop   = (double *)malloc( nx * sizeof(double));   // top    Dirichlet condition
  kbot   = (double *)malloc( nx * sizeof(double));   // bottom Neumann   condition
  

  printf("   >> Done allocating 1D arrays -- \n");


  // allocate 2D arrays dynamically
  // -- for rho --
  rho = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    rho[i] = (double *)malloc(ny*sizeof(double));

  // -- for u --
  u = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    u[i] = (double *)malloc(ny*sizeof(double));

  // -- for v --
  v = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    v[i] = (double *)malloc((ny)*sizeof(double));

  // -- for p --
  p = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    p[i] = (double *)malloc(ny*sizeof(double));

  // -- for ust --
  ust = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    ust[i] = (double *)malloc(ny*sizeof(double));

  // -- for vst --
  vst = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    vst[i] = (double *)malloc((ny)*sizeof(double));

  // -- for pst --
  pst = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    pst[i] = (double *)malloc(ny*sizeof(double));


//------for u-momentum---------------------------------
  // -- for aPu --
  aPu = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    aPu[i] = (double *)malloc(ny*sizeof(double));

  // -- for aE --
  aEu = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    aEu[i] = (double *)malloc(ny*sizeof(double));

  // -- for aW --
  aWu = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    aWu[i] = (double *)malloc(ny*sizeof(double));

  // -- for aN --
  aNu = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    aNu[i] = (double *)malloc(ny*sizeof(double));

  // -- for aS --
  aSu = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    aSu[i] = (double *)malloc(ny*sizeof(double));

  // -- for b --
  bu = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    bu[i] = (double *)malloc(ny*sizeof(double));

  // -- for Sp --
  Spu = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    Spu[i] = (double *)malloc(ny*sizeof(double));

  // -- for du --
  du = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    du[i] = (double *)malloc(ny*sizeof(double));

  uface = (double **)malloc((nx+1)*sizeof(double *));
  for(i=0; i<nx+1; i++)
    uface[i] = (double *)malloc((ny)*sizeof(double));
  
  vface = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    vface[i] = (double *)malloc((ny+1)*sizeof(double));
    
//------for u-momentum---------------------------------


//------for v-momentum---------------------------------
  // -- for aPv --
  
  aPv = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aPv[i] = (double *)malloc((ny)*sizeof(double));

  // -- for aE --
  aEv = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aEv[i] = (double *)malloc((ny)*sizeof(double));

  // -- for aW --
  aWv = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aWv[i] = (double *)malloc((ny)*sizeof(double));

  // -- for aN --
  aNv = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aNv[i] = (double *)malloc((ny)*sizeof(double));
  
  // -- for aS --
  aSv = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aSv[i] = (double *)malloc((ny)*sizeof(double));

  // -- for b --
  bv = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    bv[i] = (double *)malloc((ny)*sizeof(double));

  // -- for Sp --
  Spv = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    Spv[i] = (double *)malloc((ny)*sizeof(double));

  // -- for dv --
  dv = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    dv[i] = (double *)malloc((ny)*sizeof(double));

//------for v-momentum---------------------------------

//------for pressure  ---------------------------------
  // -- for aPp --
  aPp = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aPp[i] = (double *)malloc(ny*sizeof(double));

  // -- for aE --
  aEp = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aEp[i] = (double *)malloc(ny*sizeof(double));

  // -- for aW --
  aWp = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aWp[i] = (double *)malloc(ny*sizeof(double));

  // -- for aN --
  aNp = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aNp[i] = (double *)malloc(ny*sizeof(double));

  // -- for aS --
  aSp = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    aSp[i] = (double *)malloc(ny*sizeof(double));

  // -- for b --
  bp = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    bp[i] = (double *)malloc(ny*sizeof(double));


  ucorr = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    ucorr[i] = (double *)malloc(ny*sizeof(double));
  
  // -- for v --
  vcorr = (double **)malloc(nx*sizeof(double *));
  for(i=0; i<nx; i++)
    vcorr[i] = (double *)malloc((ny)*sizeof(double));
    
    
      ufacecorr = (double **)malloc((nx+1)*sizeof(double *));
  for(i=0; i<nx+1; i++)
    ufacecorr[i] = (double *)malloc((ny)*sizeof(double));
  
  // -- for v --
  vfacecorr = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    vfacecorr[i] = (double *)malloc((ny+1)*sizeof(double));
    
      ustface = (double **)malloc((nx+1)*sizeof(double *));
  for(i=0; i<nx+1; i++)
    ustface[i] = (double *)malloc((ny)*sizeof(double));
  
  // -- for v --
  vfacenew = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    vfacenew[i] = (double *)malloc((ny+1)*sizeof(double));

      ufacenew = (double **)malloc((nx+1)*sizeof(double *));
  for(i=0; i<nx+1; i++)
    ufacenew[i] = (double *)malloc((ny)*sizeof(double));
  
  // -- for v --
  vstface = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    vstface[i] = (double *)malloc((ny+1)*sizeof(double));

      vnew = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    vnew[i] = (double *)malloc(ny*sizeof(double));

      unew = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    unew[i] = (double *)malloc(ny*sizeof(double));

      pnew = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    pnew[i] = (double *)malloc(ny*sizeof(double));

//------for v-momentum---------------------------------

  // -- for kdiff --
  kdiff = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    kdiff[i] = (double *)malloc(ny*sizeof(double));

  // -- for work arrays wrk1, wrk2 --
  ppr = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    ppr[i] = (double *)malloc(ny*sizeof(double));

  printf("   >> Done allocating 2D arrays -- \n");
  printf(" > Done allocating memory -------- \n");
  
  //-------------------------------------------------------------------------
  
    // set the bc values
  u_bc_w = 0.0; u_bc_e = 0.0; u_bc_s = 0.0; u_bc_n = 1.0;
  v_bc_w = 0.0; v_bc_e = 0.0; v_bc_s = 0.0; v_bc_n = 0.0;
  
  // initialize the grid
  grid(nx, xst, xen, xc, xf, dxc, dxf);  // -- along x --
  grid(ny, yst, yen, yc, yf, dyc, dyf);  // -- along y --
  printf("\n > Done setting up grid ---------- \n");

  set_initial_guess(nx, ny, u, v, pst);  // initial condition
  printf("\n > Done setting up initial guess -- \n");
  // set ust, vst, pst equal to u, v, p
  
  copy_arr_2d(nx, ny, u, ust);	
  copy_arr_2d(nx, ny, v, vst);	
  
  //Face vel
  int type = 0;  //0 for Linear interpolation, 1 for RC

  for(i=0; i<nx; i++)
   for(j=0; j<ny; j++)
   {
     rho[i][j] = rho_const; // set to uniform for now
   }
  
  get_u_face(nx, ny, dxf, dyf, aPu, u, pst, du, uface, type);
  get_v_face(nx, ny, dxf, dyf, dv, vface, type, v, pst, aPv);
    
    copy_arr_2d(nx+1, ny, uface, ustface);
    copy_arr_2d(nx, ny+1, vface, vstface);

  // Start SIMPLE loop
  int advec_scheme = 1; //1 for upwind, 2 for CD2, 3 for power law; Only upwind is working now
  

  for(isim=0; isim<max_simple_iters; isim++)
  {
    get_u_coeffs(nx, ny, xc, yc, dxc, dxf, dyc, dyf, aPu, aEu, aWu, aNu, aSu, bu, Spu, kdiff, rho, u, v, pst, Fuleft, Furght, Fvtop, Fvbot, kleft, krght, ktop, kbot, uface, vface, advec_scheme, u_bc_e, u_bc_w,  u_bc_n, u_bc_s, v_bc_e, v_bc_w, v_bc_n, v_bc_s, du);                                                            
    //print_arr_2d(nx, ny, aPu);
    //print_arr_2d(nx, ny, bu);
  
    //print_arr_2d(nx, ny, u);
    solve_gssor(nx, ny, aPu, aEu, aWu, aNu, aSu, bu, u, max_iter_u, tol, relax_gssor, ust);
    //print_arr_2d(nx, ny, ust);
    
    get_v_coeffs(nx, ny, xc, yc, dxc, dxf, dyc, dyf, aPv, aEv, aWv, aNv, aSv, bv, Spv, kdiff, rho, u, v, pst, Fuleft, Furght, Fvtop, Fvbot, kleft, krght, ktop, kbot, uface, vface, advec_scheme, u_bc_e, u_bc_w,  u_bc_n, u_bc_s, v_bc_e, v_bc_w, v_bc_n, v_bc_s, dv);
    //print_arr_2d(nx, ny, aPv);
    //print_arr_2d(nx, ny, bv);
    
    solve_gssor(nx, ny, aPv, aEv, aWv, aNv, aSv, bv, v, max_iter_v, tol, relax_gssor, vst);
    //print_arr_2d(nx, ny, vst);

    //print_arr_2d(nx, ny, du);
    //print_arr_2d(nx, ny, dv);
    //print_arr_2d(nx, ny, pst);
    type = 1; //RC Interpolation to obtain ust_faces
    get_u_face(nx, ny, dxf, dyf, aPu, ust, pst, du, ustface, type);
    get_v_face(nx, ny, dxf, dyf, dv, vstface, type, vst, pst, aPv);
    
    //print_arr_2d(nx+1, ny, ustface);
    //print_arr_2d(nx, ny+1, vstface);
//------------------------------------------------------
    // get p-coeffs
    get_p_coeffs(nx, ny, aPp, aEp, aWp, aNp, aSp, bp, aPu, aPv, du, dv, dxf, dyf, ustface, vstface);
    //print_arr_2d(nx, ny, aPp);
    //print_arr_2d(nx, ny, bp);

    //--------------------------------------------------------------------------------------------------------
    // solve for p'
    solve_gssor(nx, ny, aPp, aEp, aWp, aNp, aSp, bp, ppr, max_iter_p, tol, relax_gssor, ppr);
    //print_arr_2d(nx, ny, ppr);
    
    // correct p, and obtain correction velocities
    correct_uvp(nx, ny, du, dv, ppr, ucorr, vcorr, pst, relax_p, ust, vst, pnew);
    
    //print_arr_2d(nx, ny, pnew);
    //print_arr_2d(nx, ny, ppr);
    
    face_corrections(nx, ny, ppr, du, dv, ufacecorr, vfacecorr);
    
    //Correct velocity and faces u,v: Corrections have been addded to ust, vst
    velocity_update(nx, ny, pst, u, ust, ucorr, v, vst, vcorr, uface, ustface, ufacecorr, vface, vstface, vfacecorr, relax_u, relax_v, unew, vnew, ufacenew, vfacenew); 
    //fill_ghost_values(nx, ny, pLeft, pRght, pBot, pTop, pnew, uLeft, uRght, uBot, uTop, vLeft, vRght, vBot, vTop, ust, vst, u_bc_e, u_bc_w, u_bc_n, u_bc_s, v_bc_e, v_bc_w, v_bc_n, v_bc_s);
    //set_boundary_conditions(nx, ny, uLeft, uRght, uBot, uTop, vLeft, vRght, vBot, vTop, pnew, pLeft, pRght, pBot, pTop, Fuleft, Furght, Fvtop, Fvbot, ust, vst, u_bc_e, u_bc_w, u_bc_n,                               u_bc_s ,v_bc_e, v_bc_w, v_bc_n, v_bc_s);


    // check for convergence
    l2u = get_rel_l2err_norm(nx , ny, u, unew);
    l2v = get_rel_l2err_norm(nx , ny, v, vnew);
    l2p = get_rel_l2err_norm(nx , ny, pst, pnew);	
    l2max = fmax(fmax(l2u, l2v), l2p);

    //print_arr_2d(nx, ny, u);

    printf("--> SIMPLE :: Iter No: %d, l2u = %.6e, l2v = %.6e, l2p = %.6e\n", isim, l2u, l2v, l2p);

    if(l2max < tol)
    {
      printf("--> SIMPLE :: Converged to %.6e tolerance after %d iterations\n", tol, isim);
      break;
    }

    // prepare for the next iteration
    copy_arr_2d(nx, ny, unew, u);
    copy_arr_2d(nx, ny, vnew, v);
    copy_arr_2d(nx, ny, pnew, pst);
    copy_arr_2d(nx+1, ny, ufacenew, uface);
    copy_arr_2d(nx, ny+1, vfacenew, vface);


  //print_arr_2d(nx, ny, u);
  //print_arr_2d(nx, ny, v);
  //print_arr_2d(nx, ny, pst);

  }


  output_soln_u(nx, ny, isim, xc, yc, u, "Uvel");
  output_soln_v(nx, ny, isim, xc, yc, v, "Vvel");
 // output_soln(nx, ny, isim, xc, yc, p, "Pressure");
 
  
  //output_1d_soln(nx, ny, simple_iter, Re, l2u, "u_norm");
  //output_1d_soln(nx, ny, simple_iter, Re, l2v, "v_norm");
  //output_1d_soln(nx, ny, simple_iter, Re, l2p, "p_norm");
  
  // free memory
   // ----1D arrays ---
   free(yf);  free(yc);    free(dyf);    free(dyc);
   free(xf);  free(xc);    free(dxf);    free(dxc); 
   free(kleft); free(krght); free(ktop);  free(kbot);
   free(Fuleft); free(Furght);
   free(Fvtop);  free(Fvbot);
  // // --- Done 1D arrays ---

  // // ----2D arrays --- 
   for(i=0; i<nx; i++)     free(rho[i]);    free(rho);
   for(i=0; i<nx; i++)     free(u[i]);      free(u);
   for(i=0; i<nx; i++)     free(v[i]);      free(v);
   for(i=0; i<nx; i++)     free(p[i]);      free(p);
   for(i=0; i<nx; i++)     free(unew[i]);     free(unew);
   for(i=0; i<nx; i++)     free(vnew[i]);     free(vnew);
   for(i=0; i<nx; i++)     free(pnew[i]);     free(pnew);
   for(i=0; i<nx; i++)     free(ust[i]);    free(ust);
   for(i=0; i<nx; i++)     free(vst[i]);    free(vst);
   for(i=0; i<nx; i++)     free(pst[i]);    free(pst);
   for(i=0; i<nx+1; i++)     free(ufacenew[i]);      free(ufacenew);
   for(i=0; i<nx; i++)     free(vfacenew[i]);      free(vfacenew);


   
   for(i=0; i<nx; i++)     free(kdiff[i]);   free(kdiff);
   for(i=0; i<nx; i++)     free(aPu[i]);     free(aPu);
   for(i=0; i<nx; i++)     free(aEu[i]);     free(aEu);
   for(i=0; i<nx; i++)     free(aWu[i]);     free(aWu);
   for(i=0; i<nx; i++)     free(aNu[i]);     free(aNu);
   for(i=0; i<nx; i++)     free(aSu[i]);     free(aSu);
   for(i=0; i<nx; i++)     free(Spu[i]);     free(Spu);
   for(i=0; i<nx; i++)     free(bu[i]);      free(bu);
   for(i=0; i<nx; i++)     free(du[i]);      free(du);
   for(i=0; i<nx; i++)     free(ppr[i]);    free(ppr);
   
   for(i=0; i<nx; i++)     free(aPv[i]);     free(aPv);
   for(i=0; i<nx; i++)     free(aEv[i]);     free(aEv);
   for(i=0; i<nx; i++)     free(aWv[i]);     free(aWv);
   for(i=0; i<nx; i++)     free(aNv[i]);     free(aNv);
   for(i=0; i<nx; i++)     free(aSv[i]);     free(aSv);
   for(i=0; i<nx; i++)     free(Spv[i]);     free(Spv);
   for(i=0; i<nx; i++)     free(bv[i]);      free(bv);
   for(i=0; i<nx; i++)     free(dv[i]);      free(dv);
   
   for(i=0; i<nx; i++)     free(aPp[i]);     free(aPp);
   for(i=0; i<nx; i++)     free(aEp[i]);     free(aEp);
   for(i=0; i<nx; i++)     free(aWp[i]);     free(aWp);
   for(i=0; i<nx; i++)     free(aNp[i]);     free(aNp);
   for(i=0; i<nx; i++)     free(aSp[i]);     free(aSp);
   for(i=0; i<nx; i++)     free(bp[i]);      free(bp);


   for(i=0; i<nx; i++)     free(ucorr[i]);      free(ucorr);
   for(i=0; i<nx; i++)     free(vcorr[i]);      free(vcorr);
   for(i=0; i<nx+1; i++)     free(ufacecorr[i]);      free(ufacecorr);
   for(i=0; i<nx; i++)     free(vfacecorr[i]);      free(vfacecorr);

   // // --- Done 2D arrays ---
  printf("\n > Done freeing up memory --------- \n");
  return 0;
}

