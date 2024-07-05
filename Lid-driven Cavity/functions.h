

void grid(int nx, double xst, double xen, double *xc, double *xf, double *dxc, double *dxf);

void set_initial_guess(int nx, int ny, double **u, double **v, double **p);

void calc_diffusivity(int nx, int ny, double **kdiff, double *kleft, double *krght, double *ktop, double *kbot);

void calc_sources_u(int nx, int ny, double *dxf, double *dyf, double **p, double **bu, double **Spu);

void calc_sources_v(int nx, int ny, double *dxf, double *dyf, double **p, double **bv, double **Spv);

void set_boundary_conditions(int nx, int ny, double *Fuleft, double *Furght, double *Fvtop, double *Fvbot);

void get_u_coeffs(int nx, int ny, double *xc, double *yc, double *dxc, double *dxf, double *dyc, double *dyf, double **aPu, double **aEu, double **aWu, double **aNu, double **aSu, double **bu, double **Spu, double **kdiff, double **rho, double **u, double **v, double **p, double *Fuleft, double *Furght, double *Fvtop, double *Fvbot, double *kleft, double *krght, double *ktop, double *kbot, double **uface, double **vface, int convec_scheme, double u_bc_e, double u_bc_w, double u_bc_n, double u_bc_s,double v_bc_e, double v_bc_w,double v_bc_n, double v_bc_s, double **du);                                                                  

void get_v_coeffs(int nx, int ny, double *xc, double *yc, double *dxc, double *dxf, double *dyc, double *dyf, double **aPv, double **aEv, double **aWv, double **aNv, double **aSv, double **bv, double **Spv, double **kdiff, double **rho, double **u, double **v, double **p, double *Fuleft, double *Furght, double *Fvtop, double *Fvbot, double *kleft, double *krght, double *ktop, double *kbot, double **uface, double **vface, int convec_scheme, double u_bc_e, double u_bc_w, double u_bc_n, double u_bc_s,double v_bc_e, double v_bc_w,double v_bc_n, double v_bc_s, double **dv);

void get_p_coeffs(int nx, int ny, double **aP, double **aE, double **aW, double **aN, double **aS, double **b, double **aPu, double **aPv, double **du, double **dv, double *dxf, double *dyf, double **uface, double **vface);

void correct_uvp(int nx, int ny, double **du, double **dv, double **ppr, double **ucorr, double **vcorr, double **pst, double relax_p, double **u, double **v, double **pnew);

void velocity_update(int nx, int ny, double **pst, double **u, double **ust, double **ucorr, double **v, double **vst, double **vcorr, double **uface, double **ustface, double **ufacecorr, double **vface, double **vstface, double **vfacecorr, double relax_u, double relax_v, double **unew, double **vnew, double **ufacenew, double **vfacenew); 

void convective_scheme_coeffs(int i, int j, double De, double Dw, double Dn, double Ds, double Fe, double Fw, double Fn, double Fs, double **Sp, int convec_scheme, double **aP, double **aE, double **aW, double **aN, double **aS);


// global variables
double k_const, rho_const;
