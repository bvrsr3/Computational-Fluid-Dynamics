
void get_u_face(int nx, int ny, double *dxf, double *dyf, double **aPu, double **u, double **p, double **du, double **uface, int type);

void get_v_face(int nx, int ny, double *dxf, double *dyf, double **dv, double **vface, int type, double **v, double **p, double **aPv);

void face_corrections(int nx, int ny, double **ppr, double **du, double **dv, double **ufacecorr, double **vfacecorr);




