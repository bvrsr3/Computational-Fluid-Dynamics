#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void grid(int nx, int nxglob, int istglob, int ienglob, double xstglob, double xenglob, double *x, double *dx)
{
  int i, iglob;

  *dx = (xenglob - xstglob)/(double)(nxglob-1);
  for(i=0; i<nx; i++) { 
     iglob = istglob + i; 
     x[i] = xstglob + (double)iglob* (*dx); 
  }
}

void enforce_bcs(int nx, int nxglob, int istglob, int ienglob, double *x, double *T)
{
 if(istglob==0)
   T[0] = -1.0;
 if(ienglob==(nxglob-1))
   T[nx-1] = 1.0;
}

void set_initial_condition(int nx, int nxglob, int istglob, int ienglob, double *x, double *T, int rank)
{
  int i;

  for(i=0; i<nx; i++)
    T[i] = tanh((x[i]-0.5)/0.05);

  enforce_bcs(nx,nxglob,istglob,ienglob,x,T); //ensure BCs are satisfied at t = 0
/*  
  FILE* fid;
  char debugfname[100];
  sprintf(debugfname, "debug_ic_%04d.dat", rank);
  fid = fopen(debugfname, "w");
  for (i = 0; i<nx; ++i)
   fprintf(fid, "%lf\n",T[i]);
  fclose(fid);
*/ 
}

void get_rhs(int nx, int nxglob, int istglob, int ienglob, double dx, double xleftghost, double xrghtghost, double *x, double *T, double *rhs, int rank, int it)
{
  int i;
  double dxsq = dx*dx;
  /*    
  FILE* fid;
  char debugfname[100];
  sprintf(debugfname, "debug_rhs_%04d.dat", rank);  
  fid = fopen(debugfname, "a");
  fprintf(fid, "Iter %d;", it);
  */
  for(i=0; i<nx; ++i) {
      
   if (i == 0)
        rhs[i] = (T[i+1]+xleftghost-2.0*T[i])/dxsq;
             
   else if (i == nx-1)
        rhs[i] = (xrghtghost+T[i-1]-2.0*T[i])/dxsq;
    
   else  rhs[i] = (T[i+1]+T[i-1]-2.0*T[i])/dxsq;
    
   //fprintf(fid, "%lf\t", rhs[i]);
    }
  /*
  fprintf(fid, "\n");
  fclose(fid);    
  */
}

double halo_exchange_1d_x(int rank, int size, int nx, double *x, double *T, double *xleftghost, double *xrghtghost, int it)
{
  MPI_Status status;
  int left_nb, rght_nb;
/*  
  FILE* fid;
  char debugfname[100];
  sprintf(debugfname, "debug_halo_%04d.dat", rank);
  fid = fopen(debugfname, "a");
*/
      left_nb = rank -1;
      rght_nb = rank +1;
  
      if (left_nb <= -1) { 
      left_nb = MPI_PROC_NULL; 
      *xleftghost = 0.0; 
        }
      else if (rght_nb >= size) { 
      rght_nb = MPI_PROC_NULL; 
      *xrghtghost = 0.0;
        }
      
/*  fprintf(fid, "--Start Communication--\n"); */

      double t_comm;
      //Communicating to right nb
      double t_start = MPI_Wtime(); 
      MPI_Send(&T[nx-1], 1, MPI_DOUBLE, rght_nb, 100, MPI_COMM_WORLD); 
      MPI_Recv(xrghtghost, 1, MPI_DOUBLE, rght_nb, 101, MPI_COMM_WORLD, &status);  
              
     //COmmunicating to left nb
      MPI_Send(&T[0], 1, MPI_DOUBLE, left_nb, 101, MPI_COMM_WORLD); 
      MPI_Recv(xleftghost, 1, MPI_DOUBLE, left_nb, 100, MPI_COMM_WORLD, &status);  
      double t_end = MPI_Wtime();
      t_comm = t_end - t_start;
      
      if (left_nb == -1) {  
      *xleftghost = 0.0; 
        }
      else if (rght_nb == size) {
      *xrghtghost = 0.0;
        }
        
        //printf("%lf\n",t_comm);
        return t_comm;
  /*
      fprintf(fid, "Iter: %d; Rank %d has sent %lf to rank %d (rnb)\n ", it, rank, T[nx-1], rght_nb);
      fprintf(fid, "Iter: %d; Rank %d has received %lf from rank %d (rnb)\n ", it, rank, *xrghtghost, rght_nb);
      fprintf(fid, "Iter: %d; Rank %d has sent %lf to rank %d (lnb)\n", it, rank, T[0], left_nb);
      fprintf(fid, "Iter: %d; Rank %d has received %lf from rank %d (lnb)\n", it, rank, *xleftghost, left_nb);
 
  fprintf(fid, "--Done Communication--\n");
  fclose(fid);  
  */
}

double timestep_Euler(int rank, int size, int nx, int nxglob, int istglob, int ienglob, double dt, double dx, double *x, double *T, double *rhs, int it, double *avg_tcomm)
{
  int i;
  double xleftghost, xrghtghost, t_comm;

  xleftghost = 0.0;  xrghtghost = 0.0;
  t_comm = halo_exchange_1d_x(rank, size, nx, x, T, &xleftghost, &xrghtghost, it);
  *avg_tcomm = *avg_tcomm + t_comm;
  
  // compute rhs
  get_rhs(nx,nxglob,istglob,ienglob,dx,xleftghost,xrghtghost,x,T,rhs, rank, it);

  // (Forward) Euler scheme
  for(i=0; i<nx; i++)
    T[i] = T[i] + dt*rhs[i];

  // set Dirichlet BCs
  enforce_bcs(nx,nxglob,istglob,ienglob,x,T);
  return *avg_tcomm;

}

void get_exct_soln(int nx, int nxglob, int istglob, int ienglob, double *x, double *T_exct, double tcurr)
{
  int i;

  for(i=0; i<nx; i++)
    T_exct[i] = erf((x[i]-0.5)/(2*tcurr));
}

void get_l2norm(int rank, int nx, double *T_exct, double *T, double l2norm, double tcurr, int it)
{
  int i;
  double *error;
  double sum_error_sq = 0.0;
  error = (double *)malloc(nx*sizeof(double));
  
  FILE* fp;
  char fname[100];
  sprintf(fname, "L2_%04d_%04d.dat", it, rank);
  fp = fopen(fname, "w");
  for (i=0; i<nx; ++i) {
    error[i] = pow(T_exct[i]-T[i],2);
    sum_error_sq = sum_error_sq + error[i];
  }  
  l2norm = (sqrt(sum_error_sq))/nx;
  fprintf(fp,"\n %15.13e",l2norm);
  fclose(fp);
}

void output_soln(int rank, int nx, int it, double tcurr, double *x, double *T)
{
  int i;
  FILE* fp;
  char fname[100];
  sprintf(fname, "T_x_Gathered_%04d_%04d.dat", it, rank);
  fp = fopen(fname, "w");
  
  for(i=0; i<nx; i++)
    fprintf(fp, "%lf %15.13e\n", x[i], T[i]);
  fclose(fp);
}

int main(int argc, char** argv)
{

  int nx;
  double *x, *T, *rhs, tst, ten, xst, xen, dx, dt, tcurr, xlen, l2norm, *T_exct, *buffer_T, *buffer_x, avg_tcomm;
  int i, it, num_time_steps, it_print;
  char debugfname[100];
  FILE* fid, *fp;

  int rank, size;
  int nxglob, istglob, ienglob;
  double xstglob, xenglob; 
 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // read inputs only on rank0 bcast to others
  if(rank==0)
  {
    fp = fopen("input.in", "r");
    fscanf(fp, "%d\n", &nxglob);
    fscanf(fp, "%lf %lf\n", &xstglob, &xenglob);
    fscanf(fp, "%lf %lf\n", &tst, &ten);
    fclose(fp);

    // calculate global/local variables
    nx = nxglob/size;                       // nx = 5 on each rank
    xlen = (xenglob-xstglob)/(double)size;  // xlen = 0.5 on each rank
  }
  
  // Use MPI_Bcast to send all the variables read or calculated above
  // nxglob, xstglob, xenglob, tst, ten, nx, xlen
    MPI_Bcast (&nxglob, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&nx, 1, MPI_INT, 0, MPI_COMM_WORLD);    
    
    MPI_Bcast (&xstglob, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&xenglob, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&tst, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&ten, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&xlen, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    
    istglob = rank * (nxglob/size);         // istglob = 0 for rank0; 5 for rank1
    ienglob = (rank+1) * (nxglob/size)-1;     // ienglob = 4 for rank0; 9 for rank1
    xst = xstglob + rank*xlen;              // xst = 0.0 for rank0; 0.5 for rank1
    xen = xst + xlen;                       // xen = 0.5 for rank0; 1.0 for rank1
  

  x = (double *)malloc(nx*sizeof(double));     // nx is 'local'
  T = (double *)malloc(nx*sizeof(double));
  T_exct = (double *)malloc(nx*sizeof(double));
  rhs = (double *)malloc(nx*sizeof(double));
  

  grid(nx,nxglob,istglob,ienglob,xstglob,xenglob,x,&dx); // initialize the grid
  set_initial_condition(nx,nxglob,istglob,ienglob,x,T,rank);  // initial condition

  /*
  sprintf(debugfname, "debug_%04d.dat", rank);
  fid = fopen(debugfname, "w");
  fprintf(fid, "\n--Debug-1- %d %d %d %d %d\n", rank, nx, nxglob, istglob, ienglob);
  fprintf(fid, "--Debug-2- %d %lf %lf %lf %lf %lf\n", rank, xst, xen, xstglob, xenglob, xlen);
  fprintf(fid, "\n--Writing grid points--\n");
  for(i=0; i<nx; i++)
    fprintf(fid, "%d %d %d %lf\n", rank, i, i+istglob, x[i]);
    
  fprintf(fid, "--Done writing grid points--\n");
  fclose(fid);  
  */

  // prepare for time loop
  dt = 1.25e-9; 
  num_time_steps = (int)((ten-tst)/dt) + 1; // why add 1 to this?
  it_print = num_time_steps/20;             // write out approximately 10 intermediate results

  // start time stepping loop
  double avg_time_tot = 0.0;
  double time_tot;
  
  avg_tcomm = 0.0;
  double t_start = MPI_Wtime();
  for(it=0; it<num_time_steps; it++)
  {
    tcurr = tst + (double)it * dt;    
    avg_tcomm = timestep_Euler(rank,size,nx,nxglob,istglob,ienglob,dt,dx,x,T,rhs,it,&avg_tcomm);    // get T^{n+1} from T^{n} using Forward-Euler method
    //get_exct_soln(nx,nxglob, istglob,ienglob,x,T_exct,tcurr);
    //get_l2norm(rank, nx, T_exct, T, l2norm, tcurr, it);
    //output soln every it_print time steps
    //if(it%it_print==0)
    //  output_soln(rank,nx,it,tcurr,x,T);
  }
    double t_end = MPI_Wtime();
    time_tot = t_end-t_start;
    
    avg_time_tot = time_tot/(it-1);    // Average total time
    avg_tcomm = (avg_tcomm)/(it-1);    // Average communication time
   
  //Gathering all data into root processor
  buffer_T = (double *)malloc(nxglob*sizeof(double));
  buffer_x = (double *)malloc(nxglob*sizeof(double));
  
  MPI_Gather(T, nx, MPI_DOUBLE, buffer_T, nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
  MPI_Gather(x, nx, MPI_DOUBLE, buffer_x, nx, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
  
  if (rank ==0) {
  output_soln(rank,nxglob,it,tcurr,buffer_x,buffer_T);
  printf("\n Avg Total time: %lf; Avg commn time: %lf; Average calc time: %lg",avg_time_tot, avg_tcomm, avg_time_tot-avg_tcomm);
  
  // Writing an output dat file to store time values
  FILE* fp;
  char fname[100];
  sprintf(fname, "Time_measured_%d_.dat",rank);
  fp = fopen(fname, "a");
  fprintf(fp, "%d %d %10.8e %10.8e %10.8e\n", nxglob, size, avg_time_tot, avg_tcomm, avg_time_tot-avg_tcomm);
  fclose(fp);
  
  printf("\n---- Execution Completed. Files have been written----\n");
  }
  
  free(rhs);
  free(T);
  free(x);
  free(T_exct);
  free(buffer_x);
  free(buffer_T);
  
  MPI_Finalize();
  return 0;
}
