/********************************************************************
 
 modified by huyao for data compression

 This benchmark test program is measuring a cpu performance
 of floating point operation by a Poisson equation solver.

 If you have any question, please ask me via email.
 written by Ryutaro HIMENO, November 26, 2001.
 Version 3.0
 ----------------------------------------------
 Ryutaro Himeno, Dr. of Eng.
 Head of Computer Information Division,
 RIKEN (The Institute of Pysical and Chemical Research)
 Email : himeno@postman.riken.go.jp
 ---------------------------------------------------------------
 You can adjust the size of this benchmark code to fit your target
 computer. In that case, please chose following sets of
 (mimax,mjmax,mkmax):
 small : 33,33,65
 small : 65,65,129
 midium: 129,129,257
 large : 257,257,513
 ext.large: 513,513,1025
 This program is to measure a computer performance in MFLOPS
 by using a kernel which appears in a linear solver of pressure
 Poisson eq. which appears in an incompressible Navier-Stokes solver.
 A point-Jacobi method is employed in this solver as this method can 
 be easyly vectrized and be parallelized.
 ------------------
 Finite-difference method, curvilinear coodinate system
 Vectorizable and parallelizable on each grid point
 No. of grid points : imax x jmax x kmax including boundaries
 ------------------
 A,B,C:coefficient matrix, wrk1: source term of Poisson equation
 wrk2 : working area, OMEGA : relaxation parameter
 BND:control variable for boundaries and objects ( = 0 or 1)
 P: pressure
********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>
#include "param.h"
#include "../../include/dataCompression.h"

float jacobi(int);
int initmax(int,int,int);
void initmt(int,int);
void initcomm(int,int,int);
void sendp(int,int,int);
void sendp1();
void sendp2();
void sendp3();

double fflop(int,int,int);
double mflops(int,double,double);

static float  p[MIMAX][MJMAX][MKMAX];
static float  a[4][MIMAX][MJMAX][MKMAX],
              b[3][MIMAX][MJMAX][MKMAX],
              c[3][MIMAX][MJMAX][MKMAX];
static float  bnd[MIMAX][MJMAX][MKMAX];
static float  wrk1[MIMAX][MJMAX][MKMAX],
              wrk2[MIMAX][MJMAX][MKMAX];
static float omega;
static int npe,id;

static int ndx,ndy,ndz;
static int imax,jmax,kmax;
static int ndims=3,iop[3];
static int npx[2],npy[2],npz[2];
MPI_Comm     mpi_comm_cart;
MPI_Datatype ijvec,ikvec,jkvec;

//todo
static float cr = 0; //compression rate
static int cr_num = 0;
static double duration = 0;
//modify CT
int CT = 0;

int
main(int argc,char *argv[])
{
	if(argc > 1) CT = atoi(argv[1]);

  int    i,j,k,nn;
  int    mx,my,mz,it;
  float  gosa;
  double cpu,cpu0,cpu1,flop,target;

  target= 60.0;
  omega= 0.8;
  mx= MX0-1;
  my= MY0-1;
  mz= MZ0-1;
  ndx= NDX0;
  ndy= NDY0;
  ndz= NDZ0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npe);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  initcomm(ndx,ndy,ndz);
  it= initmax(mx,my,mz);

  /*
   *    Initializing matrixes
   */
  initmt(mx,it);

  if(id==0){
    printf("Sequential version array size\n");
    printf(" mimax = %d mjmax = %d mkmax = %d\n",MX0,MY0,MZ0);
    printf("Parallel version array size\n");
    printf(" mimax = %d mjmax = %d mkmax = %d\n",MIMAX,MJMAX,MKMAX);
    printf("imax = %d jmax = %d kmax =%d\n",imax,jmax,kmax);
    printf("I-decomp = %d J-decomp = %d K-decomp =%d\n",ndx,ndy,ndz);
  }

  nn= 3;
  if(id==0){
    printf(" Start rehearsal measurement process.\n");
    printf(" Measure the performance in %d times.\n\n",nn);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  cpu0= MPI_Wtime();
  gosa= jacobi(nn);
  cpu1= MPI_Wtime() - cpu0;

  MPI_Allreduce(&cpu1,
                &cpu,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD);

  flop= fflop(mz,my,mx);
  if(id == 0){
    printf(" MFLOPS: %f time(s): %f %e\n\n",
           mflops(nn,cpu,flop),cpu,gosa);
  }
  nn= (int)(target/(cpu/3.0));
  nn= 12; // for simgrid

  if(id == 0){
    printf(" Now, start the actual measurement process.\n");
    printf(" The loop will be excuted in %d times\n",nn);
    printf(" This will take about one minute.\n");
    printf(" Wait for a while\n\n");
  }

  /*
   *    Start measuring
   */
  MPI_Barrier(MPI_COMM_WORLD);
  cpu0 = MPI_Wtime();
  gosa = jacobi(nn);
  cpu1 = MPI_Wtime() - cpu0;

  MPI_Allreduce(&cpu1,
                &cpu,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD);

  if(id == 0){
    printf("cpu : %f sec.\n", cpu);
    printf("Loop executed for %d times\n",nn);
    printf("Gosa : %e \n",gosa);
    printf("MFLOPS measured : %f\n",mflops(nn,cpu,flop));
    printf("Score based on Pentium III 600MHz : %f\n",
           mflops(nn,cpu,flop)/82.84);
    //todo
    printf("Compression rate: %f \n", 1/(cr/cr_num));
    // printf("Execution time: %f \n", duration/nn);
  }

  MPI_Finalize();
  
  return (0);
}

double
fflop(int mx,int my, int mz)
{
  return((double)(mz-2)*(double)(my-2)*(double)(mx-2)*34.0);
}

double
mflops(int nn,double cpu,double flop)
{
  return(flop/cpu*1.e-6*(double)nn);
}

void
initmt(int mx,int it)
{
  int i,j,k;

  for(i=0 ; i<MIMAX ; ++i)
    for(j=0 ; j<MJMAX ; ++j)
      for(k=0 ; k<MKMAX ; ++k){
        a[0][i][j][k]=0.0;
        a[1][i][j][k]=0.0;
        a[2][i][j][k]=0.0;
        a[3][i][j][k]=0.0;
        b[0][i][j][k]=0.0;
        b[1][i][j][k]=0.0;
        b[2][i][j][k]=0.0;
        c[0][i][j][k]=0.0;
        c[1][i][j][k]=0.0;
        c[2][i][j][k]=0.0;
        p[i][j][k]=0.0;
        wrk1[i][j][k]=0.0;
        wrk2[i][j][k]=0.0;
        bnd[i][j][k]=0.0;
      }

  for(i=0 ; i<imax ; ++i)
    for(j=0 ; j<jmax ; ++j)
      for(k=0 ; k<kmax ; ++k){
        a[0][i][j][k]=1.0;
        a[1][i][j][k]=1.0;
        a[2][i][j][k]=1.0;
        a[3][i][j][k]=1.0/6.0;
        b[0][i][j][k]=0.0;
        b[1][i][j][k]=0.0;
        b[2][i][j][k]=0.0;
        c[0][i][j][k]=1.0;
        c[1][i][j][k]=1.0;
        c[2][i][j][k]=1.0;
        p[i][j][k]=(float)((i+it)*(i+it))/(float)((mx-1)*(mx-1));
        wrk1[i][j][k]=0.0;
        wrk2[i][j][k]=0.0;
        bnd[i][j][k]=1.0;
      }
}

float
jacobi(int nn)
{
  int i,j,k,n;
  float gosa,wgosa,s0,ss;

  for(n=0 ; n<nn ; ++n){
    gosa = 0.0;
    wgosa= 0.0;

    for(i=1 ; i<imax-1 ; ++i)
      for(j=1 ; j<jmax-1 ; ++j)
        for(k=1 ; k<kmax-1 ; ++k){
          s0 = a[0][i][j][k] * p[i+1][j  ][k  ]
             + a[1][i][j][k] * p[i  ][j+1][k  ]
             + a[2][i][j][k] * p[i  ][j  ][k+1]
             + b[0][i][j][k] * ( p[i+1][j+1][k  ] - p[i+1][j-1][k  ]
                               - p[i-1][j+1][k  ] + p[i-1][j-1][k  ] )
             + b[1][i][j][k] * ( p[i  ][j+1][k+1] - p[i  ][j-1][k+1]
                               - p[i  ][j+1][k-1] + p[i  ][j-1][k-1] )
             + b[2][i][j][k] * ( p[i+1][j  ][k+1] - p[i-1][j  ][k+1]
                               - p[i+1][j  ][k-1] + p[i-1][j  ][k-1] )
             + c[0][i][j][k] * p[i-1][j  ][k  ]
             + c[1][i][j][k] * p[i  ][j-1][k  ]
             + c[2][i][j][k] * p[i  ][j  ][k-1]
             + wrk1[i][j][k];

          ss = ( s0 * a[3][i][j][k] - p[i][j][k] ) * bnd[i][j][k];
          wgosa += ss*ss;

          wrk2[i][j][k] = p[i][j][k] + omega * ss;
        }

    for(i=1 ; i<imax-1 ; ++i)
      for(j=1 ; j<jmax-1 ; ++j)
        for(k=1 ; k<kmax-1 ; ++k)
          p[i][j][k] = wrk2[i][j][k];

    double start_time = MPI_Wtime();
    sendp(ndx,ndy,ndz);
    double end_time = MPI_Wtime();
    //printf("execution time: %f\n", end_time-start_time);
    duration += end_time-start_time;

    MPI_Allreduce(&wgosa,
                  &gosa,
                  1,
                  MPI_FLOAT,
                  MPI_SUM,
                  MPI_COMM_WORLD);
  } /* end n loop */

  return(gosa);
}


void
initcomm(int ndx,int ndy,int ndz)
{
  int  i,j,k,tmp;
  int  ipd[3],idm[3],ir;
  MPI_Comm  icomm;

  if(ndx*ndy*ndz != npe){
    if(id==0){
      printf("Invalid number of PE\n");
      printf("Please check partitioning pattern or number of PE\n");
    }
    MPI_Finalize();
    exit(0);
  }

  icomm= MPI_COMM_WORLD;

  idm[0]= ndx;
  idm[1]= ndy;
  idm[2]= ndz;

  ipd[0]= 0;
  ipd[1]= 0;
  ipd[2]= 0;
  ir= 0;


  MPI_Cart_create(icomm,
                  ndims,
                  idm,
                  ipd,
                  ir,
                  &mpi_comm_cart);
  MPI_Cart_get(mpi_comm_cart,
               ndims,
               idm,
               ipd,
               iop);

  if(ndz > 1){
    MPI_Cart_shift(mpi_comm_cart,
                   2,
                   1,
                   &npz[0],
                   &npz[1]);
  }                     
  if(ndy > 1){
    MPI_Cart_shift(mpi_comm_cart,
                   1,
                   1,
                   &npy[0],
                   &npy[1]);
  }                     
  if(ndx > 1){
    MPI_Cart_shift(mpi_comm_cart,
                   0,
                   1,
                   &npx[0],
                   &npx[1]);
  }                     

}

int
initmax(int mx,int my,int mz)
{
  int  i,tmp,it;
  int  mx1[NDX0+1],my1[NDY0+1],mz1[NDZ0+1];
  int  mx2[NDX0+1],my2[NDY0+1],mz2[NDZ0+1];

  tmp= mx/ndx;
  mx1[0]= 0;
  for(i=1;i<=ndx;i++){
    if(i <= mx%ndx)
      mx1[i]= mx1[i-1] + tmp + 1;
    else
      mx1[i]= mx1[i-1] + tmp;
  }
  tmp= my/ndy;
  my1[0]= 0;
  for(i=1;i<=ndy;i++){
    if(i <= my%ndy)
      my1[i]= my1[i-1] + tmp + 1;
    else
      my1[i]= my1[i-1] + tmp;
  }
  tmp= mz/ndz;
  mz1[0]= 0;
  for(i=1;i<=ndz;i++){
    if(i <= mz%ndz)
      mz1[i]= mz1[i-1] + tmp + 1;
    else
      mz1[i]= mz1[i-1] + tmp;
  }

  for(i=0 ; i<ndx ; i++){
    mx2[i] = mx1[i+1] - mx1[i];
    if(i != 0)     mx2[i] = mx2[i] + 1;
    if(i != ndx-1) mx2[i] = mx2[i] + 1;
  }
  for(i=0 ; i<ndy ; i++){
    my2[i] = my1[i+1] - my1[i];
    if(i != 0)     my2[i] = my2[i] + 1;
    if(i != ndy-1) my2[i] = my2[i] + 1;
  }
  for(i=0 ; i<ndz ; i++){
    mz2[i] = mz1[i+1] - mz1[i];
    if(i != 0)     mz2[i] = mz2[i] + 1;
    if(i != ndz-1) mz2[i] = mz2[i] + 1;
  }

  imax = mx2[iop[0]];
  jmax = my2[iop[1]];
  kmax = mz2[iop[2]];

  if(iop[0] == 0)
    it= mx1[iop[0]];
  else
    it= mx1[iop[0]] - 1;

  if(ndx > 1){
    MPI_Type_vector(jmax,
                    kmax,
                    MKMAX,
                    MPI_FLOAT,
                    &jkvec);
    MPI_Type_commit(&jkvec);
  }                    
  if(ndy > 1){
    MPI_Type_vector(imax,
                    kmax,
                    MJMAX*MKMAX,
                    MPI_FLOAT,
                    &ikvec);
    MPI_Type_commit(&ikvec);
  }                    
  if(ndz > 1){
    MPI_Type_vector(imax*jmax,
                    1,
                    MKMAX,
                    MPI_FLOAT,
                    &ijvec);
    MPI_Type_commit(&ijvec);
  }                    

  return(it);
}

void
sendp(int ndx,int ndy,int ndz)
{
  if(ndz > 1)
    sendp3();

  if(ndy > 1)
    sendp2();

  if(ndx > 1)
    sendp1();
}

void
sendp3()
{
  MPI_Status   st[4];
  MPI_Request  req[4];

  //todo
  MPI_Status   st_bitmask[16];
  MPI_Request  req_bitmask[16];

  MPI_Status   st_plus[12];
  MPI_Request  req_plus[12];

  MPI_Status   st_bitwise[8];
  MPI_Request  req_bitwise[8];  

  MPI_Status   st_sz[4];
  MPI_Request  req_sz[4];  

  if(CT == 7)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npz[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npz[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 3, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 3, kmax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], imax*jmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], imax*jmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    int type[2] = {0, 0};
    float medium[2];
    medium[0] = med_dataset_float(data_small[0], imax*jmax, &type[0]);
    medium[1] = med_dataset_float(data_small[1], imax*jmax, &type[1]);
    char float_arr0[32+1];
    char float_arr1[32+1];
    floattostr(&medium[0], float_arr0);
    floattostr(&medium[1], float_arr1);
    char mask0[1+8+8];
    char mask1[1+8+8];
    strncpy(mask0, float_arr0, 1+8+8);
    strncpy(mask1, float_arr1, 1+8+8);

    myCompress_bitwise_mask(data_small[0], imax*jmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0], type[0], mask0);
    myCompress_bitwise_mask(data_small[1], imax*jmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1], type[1], mask1);    

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npz[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npz[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    cr += data_bytes_send[0]*8.0/(imax*jmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(imax*jmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];
    float medium_recv[2];
    int type_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npz[1], 2, mpi_comm_cart, req_bitmask);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npz[1], 3, mpi_comm_cart, req_bitmask+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npz[0], 4, mpi_comm_cart, req_bitmask+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npz[0], 5, mpi_comm_cart, req_bitmask+3);  
    MPI_Irecv(&medium_recv[0], 1, MPI_FLOAT, npz[1], 6, mpi_comm_cart, req_bitmask+4);
    MPI_Irecv(&medium_recv[1], 1, MPI_FLOAT, npz[0], 7, mpi_comm_cart, req_bitmask+5);
    MPI_Irecv(&type_recv[0], 1, MPI_INT, npz[1], 8, mpi_comm_cart, req_bitmask+6);
    MPI_Irecv(&type_recv[1], 1, MPI_INT, npz[0], 9, mpi_comm_cart, req_bitmask+7);    
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npz[0], 2, mpi_comm_cart, req_bitmask+8); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npz[0], 3, mpi_comm_cart, req_bitmask+9); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npz[1], 4, mpi_comm_cart, req_bitmask+10); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npz[1], 5, mpi_comm_cart, req_bitmask+11); 
    MPI_Isend(&medium[0], 1, MPI_FLOAT, npz[0], 6, mpi_comm_cart, req_bitmask+12); 
    MPI_Isend(&medium[1], 1, MPI_FLOAT, npz[1], 7, mpi_comm_cart, req_bitmask+13); 
    MPI_Isend(&type[0], 1, MPI_INT, npz[0], 8, mpi_comm_cart, req_bitmask+14); 
    MPI_Isend(&type[1], 1, MPI_INT, npz[1], 9, mpi_comm_cart, req_bitmask+15);     
    MPI_Waitall(16, req_bitmask, st_bitmask);    

    char float_arr0_recv[32+1];
    char float_arr1_recv[32+1];
    floattostr(&medium_recv[0], float_arr0_recv);
    floattostr(&medium_recv[1], float_arr1_recv);
    char mask0_recv[1+8+8];
    char mask1_recv[1+8+8];
    strncpy(mask0_recv, float_arr0_recv, 1+8+8);
    strncpy(mask1_recv, float_arr1_recv, 1+8+8);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise_mask(data_bits_recv[0], data_bytes_recv[0], imax*jmax, type_recv[0], mask0_recv);
    decompressed_data[1] = myDecompress_bitwise_mask(data_bits_recv[1], data_bytes_recv[1], imax*jmax, type_recv[1], mask1_recv);

    int pointer = 0;
    for(int a=0; a<imax; a++)
    {
      for(int b=0; b<jmax; b++)
      {
        p[a][b][kmax-1] = decompressed_data[0][pointer] + data_min_recv[0];
        p[a][b][0] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  }

  if(CT == 6)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npz[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npz[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 3, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 3, kmax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], imax*jmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], imax*jmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    myCompress_bitwise_np(data_small[0], imax*jmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0]);
    myCompress_bitwise_np(data_small[1], imax*jmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1]);

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npz[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npz[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    // printf("npz %d %d \n", npz[0], npz[1]);  
    // printf("send %d %d \n", data_bytes_send[0], data_bytes_send[1]);
    // printf("recv %d %d \n", data_bytes_recv[0], data_bytes_recv[1]);

    cr += data_bytes_send[0]*8.0/(imax*jmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(imax*jmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npz[1], 2, mpi_comm_cart, req_bitwise);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npz[1], 3, mpi_comm_cart, req_bitwise+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npz[0], 4, mpi_comm_cart, req_bitwise+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npz[0], 5, mpi_comm_cart, req_bitwise+3);  
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npz[0], 2, mpi_comm_cart, req_bitwise+4); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npz[0], 3, mpi_comm_cart, req_bitwise+5); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npz[1], 4, mpi_comm_cart, req_bitwise+6); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npz[1], 5, mpi_comm_cart, req_bitwise+7); 
    MPI_Waitall(8, req_bitwise, st_bitwise);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise_np(data_bits_recv[0], data_bytes_recv[0], imax*jmax);
    decompressed_data[1] = myDecompress_bitwise_np(data_bits_recv[1], data_bytes_recv[1], imax*jmax);
    int pointer = 0;
    for(int a=0; a<imax; a++)
    {
      for(int b=0; b<jmax; b++)
      {
        p[a][b][kmax-1] = decompressed_data[0][pointer] + data_min_recv[0];
        p[a][b][0] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  }

  if(CT == 5)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npz[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npz[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 3, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 3, kmax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], imax*jmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], imax*jmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    myCompress_bitwise(data_small[0], imax*jmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0]);
    myCompress_bitwise(data_small[1], imax*jmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1]);

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npz[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npz[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    // printf("npz %d %d \n", npz[0], npz[1]);  
    // printf("send %d %d \n", data_bytes_send[0], data_bytes_send[1]);
    // printf("recv %d %d \n", data_bytes_recv[0], data_bytes_recv[1]);

    cr += data_bytes_send[0]*8.0/(imax*jmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(imax*jmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npz[1], 2, mpi_comm_cart, req_bitwise);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npz[1], 3, mpi_comm_cart, req_bitwise+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npz[0], 4, mpi_comm_cart, req_bitwise+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npz[0], 5, mpi_comm_cart, req_bitwise+3);  
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npz[0], 2, mpi_comm_cart, req_bitwise+4); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npz[0], 3, mpi_comm_cart, req_bitwise+5); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npz[1], 4, mpi_comm_cart, req_bitwise+6); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npz[1], 5, mpi_comm_cart, req_bitwise+7); 
    MPI_Waitall(8, req_bitwise, st_bitwise);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise(data_bits_recv[0], data_bytes_recv[0], imax*jmax);
    decompressed_data[1] = myDecompress_bitwise(data_bits_recv[1], data_bytes_recv[1], imax*jmax);
    int pointer = 0;
    for(int a=0; a<imax; a++)
    {
      for(int b=0; b<jmax; b++)
      {
        p[a][b][kmax-1] = decompressed_data[0][pointer] + data_min_recv[0];
        p[a][b][0] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  }

  if(CT == 4)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npz[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npz[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 3, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 3, kmax-2, imax, jmax, kmax);   

    unsigned char* data_bits_send[2] = {NULL, NULL};

    char binfile0[64];
    sprintf(binfile0, "dataset/r%dp3data0.dat", id); 
    char binfile1[64];
    sprintf(binfile1, "dataset/r%dp3data1.dat", id); 
    writetobinary_float(binfile0, data[0], imax*jmax); //.txt --> .dat
    writetobinary_float(binfile1, data[1], imax*jmax); //.txt --> .dat    
    char sz_comp_cmd0[64];
    char sz_comp_cmd1[64];
    sprintf(sz_comp_cmd0, "%s%g%sdataset/r%dp3data0%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, id, sz_comp_cmd_suffix2, imax*jmax);
    sprintf(sz_comp_cmd1, "%s%g%sdataset/r%dp3data1%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, id, sz_comp_cmd_suffix2, imax*jmax);
    //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
    int iret_comp0 = system(sz_comp_cmd0); //.dat --> .dat.sz
    int iret_comp1 = system(sz_comp_cmd1); //.dat --> .dat.sz
    char binfile_sz0[64];
    sprintf(binfile_sz0, "dataset/r%dp3data0.dat.sz", id);
    char binfile_sz1[64];
    sprintf(binfile_sz1, "dataset/r%dp3data1.dat.sz", id);
    data_bits_send[0] = readfrombinary_char(binfile_sz0, &data_bytes_send[0]);  
    data_bits_send[1] = readfrombinary_char(binfile_sz1, &data_bytes_send[1]);   

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npz[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npz[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    cr += data_bytes_send[0]*8.0/(imax*jmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(imax*jmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npz[1], 2, mpi_comm_cart, req_sz);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npz[0], 3, mpi_comm_cart, req_sz+1);  
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npz[0], 2, mpi_comm_cart, req_sz+2); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npz[1], 3, mpi_comm_cart, req_sz+3); 
    MPI_Waitall(4, req_sz, st_sz);    

    char binfile_zs0[64];
    sprintf(binfile_zs0, "dataset/r%dp3data0.dat.zs", id);
    char binfile_zs1[64];
    sprintf(binfile_zs1, "dataset/r%dp3data1.dat.zs", id);
    writetobinary_char(binfile_zs0, data_bits_recv[0], data_bytes_recv[0]); //.dat.zs
    writetobinary_char(binfile_zs1, data_bits_recv[1], data_bytes_recv[1]); //.dat.zs
    char sz_decomp_cmd0[64];
    char sz_decomp_cmd1[64];
    sprintf(sz_decomp_cmd0, "%sdataset/r%dp3data0%s%d", sz_decomp_cmd_prefix, id, sz_decomp_cmd_suffix, imax*jmax);
    sprintf(sz_decomp_cmd1, "%sdataset/r%dp3data1%s%d", sz_decomp_cmd_prefix, id, sz_decomp_cmd_suffix, imax*jmax);
    //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
    int iret_decomp0 = system(sz_decomp_cmd0); //.dat.zs --> .dat.zs.out
    int iret_decomp1 = system(sz_decomp_cmd1); //.dat.zs --> .dat.zs.out
    char binfile_out0[64];
    sprintf(binfile_out0, "dataset/r%dp3data0.dat.zs.out", id);
    char binfile_out1[64];
    sprintf(binfile_out1, "dataset/r%dp3data1.dat.zs.out", id);
    char txtfile0[64];
    sprintf(txtfile0, "dataset/r%dp3data0.dat.zs.out.txt", id);  
    char txtfile1[64];
    sprintf(txtfile1, "dataset/r%dp3data1.dat.zs.out.txt", id);  

    float* decompressed_data[2];
    decompressed_data[0] = readfrombinary_writetotxt_float(binfile_out0, txtfile0, imax*jmax);
    decompressed_data[1] = readfrombinary_writetotxt_float(binfile_out1, txtfile1, imax*jmax);

    int pointer = 0;
    for(int a=0; a<imax; a++)
    {
      for(int b=0; b<jmax; b++)
      {
        p[a][b][kmax-1] = decompressed_data[0][pointer];
        p[a][b][0] = decompressed_data[1][pointer];
        pointer++;
      }
    }
  }  

  //revised
  if(CT == 1)
  {
    int array_float_len_send[2] = {0, 0};
    int array_float_len_recv[2] = {0, 0};

    float* array_float_send[2] = {NULL, NULL};
    char* array_char_send[2] = {NULL, NULL};
    int* array_char_displacement_send[2] = {NULL, NULL};

    MPI_Irecv(&array_float_len_recv[0], 1, MPI_INT, npz[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&array_float_len_recv[1], 1, MPI_INT, npz[0], 1, mpi_comm_cart, req+1);
    float* data_0 = transform_3d_array_to_1d_array(p, 3, 1, imax, jmax, kmax);
    float* data_1 = transform_3d_array_to_1d_array(p, 3, kmax-2, imax, jmax, kmax);

    array_float_len_send[0] = myCompress(data_0, &array_float_send[0], &array_char_send[0], &array_char_displacement_send[0], imax*jmax);
    array_float_len_send[1] = myCompress(data_1, &array_float_send[1], &array_char_send[1], &array_char_displacement_send[1], imax*jmax);
    MPI_Isend(&array_float_len_send[0], 1, MPI_INT, npz[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&array_float_len_send[1], 1, MPI_INT, npz[1], 1, mpi_comm_cart, req+3); 
    //printf("send1 %d %d \n", array_float_len_send[0], array_float_len_send[1]);
    MPI_Waitall(4, req, st);

    //printf("recv1 %d %d \n", array_float_len_recv[0], array_float_len_recv[1]);

    float* array_float_recv[2]; 
    char* array_char_recv[2]; 
    int* array_char_displacement_recv[2]; 

    int num_p_recv_0 = array_float_len_recv[0], num_c_recv_0 = imax*jmax - array_float_len_recv[0];
    int num_p_recv_1 = array_float_len_recv[1], num_c_recv_1 = imax*jmax - array_float_len_recv[1];
    array_float_recv[0] = (float*) malloc(sizeof(float)*num_p_recv_0);
    array_char_recv[0] = (char*) malloc(sizeof(char)*num_c_recv_0);
    array_char_displacement_recv[0] = (int*) malloc(sizeof(int)*num_c_recv_0);
    array_float_recv[1] = (float*) malloc(sizeof(float)*num_p_recv_1);
    array_char_recv[1] = (char*) malloc(sizeof(char)*num_c_recv_1);
    array_char_displacement_recv[1] = (int*) malloc(sizeof(int)*num_c_recv_1);
    MPI_Irecv(array_float_recv[0], num_p_recv_0, MPI_FLOAT, npz[1], 2, mpi_comm_cart, req_plus);
    MPI_Irecv(array_char_recv[0], num_c_recv_0, MPI_CHAR, npz[1], 3, mpi_comm_cart, req_plus+1);
    MPI_Irecv(array_char_displacement_recv[0], num_c_recv_0, MPI_INT, npz[1], 4, mpi_comm_cart, req_plus+2);    
    MPI_Irecv(array_float_recv[1], num_p_recv_1, MPI_FLOAT, npz[0], 5, mpi_comm_cart, req_plus+3);
    MPI_Irecv(array_char_recv[1], num_c_recv_1, MPI_CHAR, npz[0], 6, mpi_comm_cart, req_plus+4);
    MPI_Irecv(array_char_displacement_recv[1], num_c_recv_1, MPI_INT, npz[0], 7, mpi_comm_cart, req_plus+5);  
    //mycompress
    int num_p_send_0 = array_float_len_send[0], num_c_send_0 = imax*jmax - array_float_len_send[0];
    int num_p_send_1 = array_float_len_send[1], num_c_send_1 = imax*jmax - array_float_len_send[1];
    MPI_Isend(array_float_send[0], num_p_send_0, MPI_FLOAT, npz[0], 2, mpi_comm_cart, req_plus+6); 
    MPI_Isend(array_char_send[0], num_c_send_0, MPI_CHAR, npz[0], 3, mpi_comm_cart, req_plus+7); 
    MPI_Isend(array_char_displacement_send[0], num_c_send_0, MPI_INT, npz[0], 4, mpi_comm_cart, req_plus+8); 
    MPI_Isend(array_float_send[1], num_p_send_1, MPI_FLOAT, npz[1], 5, mpi_comm_cart, req_plus+9); 
    MPI_Isend(array_char_send[1], num_c_send_1, MPI_CHAR, npz[1], 6, mpi_comm_cart, req_plus+10); 
    MPI_Isend(array_char_displacement_send[1], num_c_send_1, MPI_INT, npz[1], 7, mpi_comm_cart, req_plus+11); 
    MPI_Waitall(12, req_plus, st_plus);

    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + 1.0*((float)num_p_send_0/(imax*jmax));
    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_1/(imax*jmax)) + 1.0*((float)num_p_send_1/(imax*jmax));
    cr_num += 2;

    //calculate bitwise compress ratio
    // printf("%d %d\n", num_p_send_0, num_c_send_0);
    // printf("%d %d\n", num_p_send_1, num_c_send_1);
    // printf("%f \n", calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0));
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax));
    // printf("%f \n", (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax)));
    // printf("%f \n", cr);
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_1/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[1], num_p_send_1)*((float)num_p_send_1/(imax*jmax));
    // printf("%f \n", cr);
    // cr_num += 2;
    // printf("%d \n", cr_num);

    float* decompressed_data_0 = myDecompress(array_float_recv[0], array_char_recv[0], array_char_displacement_recv[0], imax*jmax);
    int pointer_0 = 0;
    // for(int a=0; a<imax; a++)
    // {
    //   for(int b=0; b<jmax; b++)
    //   {
    //     p[a][b][kmax-1] = decompressed_data_0[pointer_0++];
    //   }
    // }
    float* decompressed_data_1 = myDecompress(array_float_recv[1], array_char_recv[1], array_char_displacement_recv[1], imax*jmax);
    int pointer_1 = 0;
    for(int a=0; a<imax; a++)
    {
      for(int b=0; b<jmax; b++)
      {
        p[a][b][kmax-1] = decompressed_data_0[pointer_0++];
        p[a][b][0] = decompressed_data_1[pointer_1++];
      }
    }
  }

  if (CT == 0)
  {
    MPI_Irecv(&p[0][0][kmax-1],
              1,
              ijvec,
              npz[1],
              1,
              mpi_comm_cart,
              req);
    MPI_Irecv(&p[0][0][0],
              1,
              ijvec,
              npz[0],
              2,
              mpi_comm_cart,
              req+1);
  }

  //todo
  // if(CT == 1)
  // {
  //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 3, 1, imax, jmax, kmax);
  //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 3, kmax-2, imax, jmax, kmax);
  // }
  // else if(CT == 2)
  // {
  //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 3, 1, imax, jmax, kmax);
  //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 3, kmax-2, imax, jmax, kmax); 
  // } 
  // else if(CT == 3)    
  // {
  //   cr += calcCompressionRatio_himeno_nolossy_area(p, 3, 1, imax, jmax, kmax);
  //   cr += calcCompressionRatio_himeno_nolossy_area(p, 3, kmax-2, imax, jmax, kmax);     
  // }
  // else if(CT == 4)    
  // {
  //   cr += calcCompressionRatio_himeno_sz(p, 3, 1, imax, jmax, kmax);
  //   cr += calcCompressionRatio_himeno_sz(p, 3, kmax-2, imax, jmax, kmax);     
  // }
  // cr_num += 2;

  if(CT == 0)
  {
    MPI_Isend(&p[0][0][1],
              1,
              ijvec,
              npz[0],
              1,
              mpi_comm_cart,
              req+2);
    MPI_Isend(&p[0][0][kmax-2],
              1,
              ijvec,
              npz[1],
              2,
              mpi_comm_cart,
              req+3);
    MPI_Waitall(4,
                req,
                st);

    //printf("npz %d %d \n", npz[0], npz[1]);  
  }
}

void
sendp2()
{
  MPI_Status  st[4];
  MPI_Request req[4];

  //todo
  MPI_Status   st_bitmask[16];
  MPI_Request  req_bitmask[16];

  MPI_Status   st_plus[12];
  MPI_Request  req_plus[12];

  MPI_Status   st_bitwise[8];
  MPI_Request  req_bitwise[8];  

  MPI_Status   st_sz[4];
  MPI_Request  req_sz[4];    

  if(CT == 7)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npy[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npy[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 2, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 2, jmax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], imax*kmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], imax*kmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    int type[2] = {0, 0};
    float medium[2];
    medium[0] = med_dataset_float(data_small[0], imax*kmax, &type[0]);
    medium[1] = med_dataset_float(data_small[1], imax*kmax, &type[1]);
    char float_arr0[32+1];
    char float_arr1[32+1];
    floattostr(&medium[0], float_arr0);
    floattostr(&medium[1], float_arr1);
    char mask0[1+8+8];
    char mask1[1+8+8];
    strncpy(mask0, float_arr0, 1+8+8);
    strncpy(mask1, float_arr1, 1+8+8);

    myCompress_bitwise_mask(data_small[0], imax*kmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0], type[0], mask0);
    myCompress_bitwise_mask(data_small[1], imax*kmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1], type[1], mask1);    

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npy[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npy[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    cr += data_bytes_send[0]*8.0/(imax*kmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(imax*kmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];
    float medium_recv[2];
    int type_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npy[1], 2, mpi_comm_cart, req_bitmask);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npy[1], 3, mpi_comm_cart, req_bitmask+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npy[0], 4, mpi_comm_cart, req_bitmask+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npy[0], 5, mpi_comm_cart, req_bitmask+3);  
    MPI_Irecv(&medium_recv[0], 1, MPI_FLOAT, npy[1], 6, mpi_comm_cart, req_bitmask+4);
    MPI_Irecv(&medium_recv[1], 1, MPI_FLOAT, npy[0], 7, mpi_comm_cart, req_bitmask+5);
    MPI_Irecv(&type_recv[0], 1, MPI_INT, npy[1], 8, mpi_comm_cart, req_bitmask+6);
    MPI_Irecv(&type_recv[1], 1, MPI_INT, npy[0], 9, mpi_comm_cart, req_bitmask+7);    
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npy[0], 2, mpi_comm_cart, req_bitmask+8); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npy[0], 3, mpi_comm_cart, req_bitmask+9); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npy[1], 4, mpi_comm_cart, req_bitmask+10); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npy[1], 5, mpi_comm_cart, req_bitmask+11); 
    MPI_Isend(&medium[0], 1, MPI_FLOAT, npy[0], 6, mpi_comm_cart, req_bitmask+12); 
    MPI_Isend(&medium[1], 1, MPI_FLOAT, npy[1], 7, mpi_comm_cart, req_bitmask+13); 
    MPI_Isend(&type[0], 1, MPI_INT, npy[0], 8, mpi_comm_cart, req_bitmask+14); 
    MPI_Isend(&type[1], 1, MPI_INT, npy[1], 9, mpi_comm_cart, req_bitmask+15);     
    MPI_Waitall(16, req_bitmask, st_bitmask);    

    char float_arr0_recv[32+1];
    char float_arr1_recv[32+1];
    floattostr(&medium_recv[0], float_arr0_recv);
    floattostr(&medium_recv[1], float_arr1_recv);
    char mask0_recv[1+8+8];
    char mask1_recv[1+8+8];
    strncpy(mask0_recv, float_arr0_recv, 1+8+8);
    strncpy(mask1_recv, float_arr1_recv, 1+8+8);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise_mask(data_bits_recv[0], data_bytes_recv[0], imax*kmax, type_recv[0], mask0_recv);
    decompressed_data[1] = myDecompress_bitwise_mask(data_bits_recv[1], data_bytes_recv[1], imax*kmax, type_recv[1], mask1_recv);

    int pointer = 0;
    for(int a=0; a<imax; a++)
    {
      for(int b=0; b<kmax; b++)
      {
        p[a][jmax-1][b] = decompressed_data[0][pointer] + data_min_recv[0];
        p[a][0][b] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  }

  if(CT == 6)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npy[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npy[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 2, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 2, jmax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], imax*kmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], imax*kmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    myCompress_bitwise_np(data_small[0], imax*kmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0]);
    myCompress_bitwise_np(data_small[1], imax*kmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1]);

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npy[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npy[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    // printf("npy %d %d \n", npy[0], npy[1]);  
    // printf("send %d %d \n", data_bytes_send[0], data_bytes_send[1]);
    // printf("recv %d %d \n", data_bytes_recv[0], data_bytes_recv[1]);

    cr += data_bytes_send[0]*8.0/(imax*kmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(imax*kmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npy[1], 2, mpi_comm_cart, req_bitwise);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npy[1], 3, mpi_comm_cart, req_bitwise+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npy[0], 4, mpi_comm_cart, req_bitwise+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npy[0], 5, mpi_comm_cart, req_bitwise+3);  
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npy[0], 2, mpi_comm_cart, req_bitwise+4); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npy[0], 3, mpi_comm_cart, req_bitwise+5); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npy[1], 4, mpi_comm_cart, req_bitwise+6); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npy[1], 5, mpi_comm_cart, req_bitwise+7); 
    MPI_Waitall(8, req_bitwise, st_bitwise);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise_np(data_bits_recv[0], data_bytes_recv[0], imax*kmax);
    decompressed_data[1] = myDecompress_bitwise_np(data_bits_recv[1], data_bytes_recv[1], imax*kmax);
    int pointer = 0;
    for(int a=0; a<imax; a++)
    {
      for(int b=0; b<kmax; b++)
      {
        p[a][jmax-1][b] = decompressed_data[0][pointer] + data_min_recv[0];
        p[a][0][b] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  } 

  if(CT == 5)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npy[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npy[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 2, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 2, jmax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], imax*kmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], imax*kmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    myCompress_bitwise(data_small[0], imax*kmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0]);
    myCompress_bitwise(data_small[1], imax*kmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1]);

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npy[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npy[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    // printf("npy %d %d \n", npy[0], npy[1]);  
    // printf("send %d %d \n", data_bytes_send[0], data_bytes_send[1]);
    // printf("recv %d %d \n", data_bytes_recv[0], data_bytes_recv[1]);

    cr += data_bytes_send[0]*8.0/(imax*kmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(imax*kmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npy[1], 2, mpi_comm_cart, req_bitwise);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npy[1], 3, mpi_comm_cart, req_bitwise+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npy[0], 4, mpi_comm_cart, req_bitwise+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npy[0], 5, mpi_comm_cart, req_bitwise+3);  
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npy[0], 2, mpi_comm_cart, req_bitwise+4); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npy[0], 3, mpi_comm_cart, req_bitwise+5); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npy[1], 4, mpi_comm_cart, req_bitwise+6); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npy[1], 5, mpi_comm_cart, req_bitwise+7); 
    MPI_Waitall(8, req_bitwise, st_bitwise);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise(data_bits_recv[0], data_bytes_recv[0], imax*kmax);
    decompressed_data[1] = myDecompress_bitwise(data_bits_recv[1], data_bytes_recv[1], imax*kmax);
    int pointer = 0;
    for(int a=0; a<imax; a++)
    {
      for(int b=0; b<kmax; b++)
      {
        p[a][jmax-1][b] = decompressed_data[0][pointer] + data_min_recv[0];
        p[a][0][b] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  } 

  if(CT == 4)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npy[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npy[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 2, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 2, jmax-2, imax, jmax, kmax);   

    unsigned char* data_bits_send[2] = {NULL, NULL};

    char binfile0[64];
    sprintf(binfile0, "dataset/r%dp2data0.dat", id); 
    char binfile1[64];
    sprintf(binfile1, "dataset/r%dp2data1.dat", id);     
    writetobinary_float(binfile0, data[0], imax*kmax); //.txt --> .dat
    writetobinary_float(binfile1, data[1], imax*kmax); //.txt --> .dat    
    char sz_comp_cmd0[64];
    char sz_comp_cmd1[64];
    sprintf(sz_comp_cmd0, "%s%g%sdataset/r%dp2data0%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, id, sz_comp_cmd_suffix2, imax*kmax);
    sprintf(sz_comp_cmd1, "%s%g%sdataset/r%dp2data1%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, id, sz_comp_cmd_suffix2, imax*kmax);
    //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
    int iret_comp0 = system(sz_comp_cmd0); //.dat --> .dat.sz
    int iret_comp1 = system(sz_comp_cmd1); //.dat --> .dat.sz
    char binfile_sz0[64];
    sprintf(binfile_sz0, "dataset/r%dp2data0.dat.sz", id);
    char binfile_sz1[64];
    sprintf(binfile_sz1, "dataset/r%dp2data1.dat.sz", id);    
    data_bits_send[0] = readfrombinary_char(binfile_sz0, &data_bytes_send[0]);  
    data_bits_send[1] = readfrombinary_char(binfile_sz1, &data_bytes_send[1]);   

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npy[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npy[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    cr += data_bytes_send[0]*8.0/(imax*kmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(imax*kmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npy[1], 2, mpi_comm_cart, req_sz);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npy[0], 3, mpi_comm_cart, req_sz+1);  
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npy[0], 2, mpi_comm_cart, req_sz+2); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npy[1], 3, mpi_comm_cart, req_sz+3); 
    MPI_Waitall(4, req_sz, st_sz);    

    char binfile_zs0[64];
    sprintf(binfile_zs0, "dataset/r%dp2data0.dat.zs", id);
    char binfile_zs1[64];
    sprintf(binfile_zs1, "dataset/r%dp2data1.dat.zs", id);    
    writetobinary_char(binfile_zs0, data_bits_recv[0], data_bytes_recv[0]); //.dat.zs
    writetobinary_char(binfile_zs1, data_bits_recv[1], data_bytes_recv[1]); //.dat.zs
    char sz_decomp_cmd0[64];
    char sz_decomp_cmd1[64];
    sprintf(sz_decomp_cmd0, "%sdataset/r%dp2data0%s%d", sz_decomp_cmd_prefix, id, sz_decomp_cmd_suffix, imax*kmax);
    sprintf(sz_decomp_cmd1, "%sdataset/r%dp2data1%s%d", sz_decomp_cmd_prefix, id, sz_decomp_cmd_suffix, imax*kmax);
    //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
    int iret_decomp0 = system(sz_decomp_cmd0); //.dat.zs --> .dat.zs.out
    int iret_decomp1 = system(sz_decomp_cmd1); //.dat.zs --> .dat.zs.out
    char binfile_out0[64];
    sprintf(binfile_out0, "dataset/r%dp2data0.dat.zs.out", id);
    char binfile_out1[64];
    sprintf(binfile_out1, "dataset/r%dp2data1.dat.zs.out", id);
    char txtfile0[64];
    sprintf(txtfile0, "dataset/r%dp2data0.dat.zs.out.txt", id);  
    char txtfile1[64];
    sprintf(txtfile1, "dataset/r%dp2data1.dat.zs.out.txt", id);      

    float* decompressed_data[2];
    decompressed_data[0] = readfrombinary_writetotxt_float(binfile_out0, txtfile0, imax*kmax);
    decompressed_data[1] = readfrombinary_writetotxt_float(binfile_out1, txtfile1, imax*kmax);

    int pointer = 0;
    for(int a=0; a<imax; a++)
    {
      for(int b=0; b<kmax; b++)
      {
        p[a][jmax-1][b] = decompressed_data[0][pointer];
        p[a][0][b] = decompressed_data[1][pointer];
        pointer++;
      }
    }
  }    

  if(CT == 1)
  {
    int array_float_len_send[2] = {0, 0};
    int array_float_len_recv[2] = {0, 0};

    float* array_float_send[2] = {NULL, NULL};
    char* array_char_send[2] = {NULL, NULL};
    int* array_char_displacement_send[2] = {NULL, NULL};

    MPI_Irecv(&array_float_len_recv[0], 1, MPI_INT, npy[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&array_float_len_recv[1], 1, MPI_INT, npy[0], 1, mpi_comm_cart, req+1);
    float* data_0 = transform_3d_array_to_1d_array(p, 2, 1, imax, jmax, kmax);
    float* data_1 = transform_3d_array_to_1d_array(p, 2, jmax-2, imax, jmax, kmax);
    array_float_len_send[0] = myCompress(data_0, &array_float_send[0], &array_char_send[0], &array_char_displacement_send[0], imax*kmax);
    array_float_len_send[1] = myCompress(data_1, &array_float_send[1], &array_char_send[1], &array_char_displacement_send[1], imax*kmax);
    MPI_Isend(&array_float_len_send[0], 1, MPI_INT, npy[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&array_float_len_send[1], 1, MPI_INT, npy[1], 1, mpi_comm_cart, req+3); 
    //printf("send1 %d %d \n", array_float_len_send[0], array_float_len_send[1]);
    MPI_Waitall(4, req, st);

    //printf("recv1 %d %d \n", array_float_len_recv[0], array_float_len_recv[1]);

    float* array_float_recv[2]; 
    char* array_char_recv[2]; 
    int* array_char_displacement_recv[2]; 

    int num_p_recv_0 = array_float_len_recv[0], num_c_recv_0 = imax*kmax - array_float_len_recv[0];
    int num_p_recv_1 = array_float_len_recv[1], num_c_recv_1 = imax*kmax - array_float_len_recv[1];
    array_float_recv[0] = (float*) malloc(sizeof(float)*num_p_recv_0);
    array_char_recv[0] = (char*) malloc(sizeof(char)*num_c_recv_0);
    array_char_displacement_recv[0] = (int*) malloc(sizeof(int)*num_c_recv_0);
    array_float_recv[1] = (float*) malloc(sizeof(float)*num_p_recv_1);
    array_char_recv[1] = (char*) malloc(sizeof(char)*num_c_recv_1);
    array_char_displacement_recv[1] = (int*) malloc(sizeof(int)*num_c_recv_1);
    MPI_Irecv(array_float_recv[0], num_p_recv_0, MPI_FLOAT, npy[1], 2, mpi_comm_cart, req_plus);
    MPI_Irecv(array_char_recv[0], num_c_recv_0, MPI_CHAR, npy[1], 3, mpi_comm_cart, req_plus+1);
    MPI_Irecv(array_char_displacement_recv[0], num_c_recv_0, MPI_INT, npy[1], 4, mpi_comm_cart, req_plus+2);    
    MPI_Irecv(array_float_recv[1], num_p_recv_1, MPI_FLOAT, npy[0], 5, mpi_comm_cart, req_plus+3);
    MPI_Irecv(array_char_recv[1], num_c_recv_1, MPI_CHAR, npy[0], 6, mpi_comm_cart, req_plus+4);
    MPI_Irecv(array_char_displacement_recv[1], num_c_recv_1, MPI_INT, npy[0], 7, mpi_comm_cart, req_plus+5);  
    //mycompress
    int num_p_send_0 = array_float_len_send[0], num_c_send_0 = imax*kmax - array_float_len_send[0];
    int num_p_send_1 = array_float_len_send[1], num_c_send_1 = imax*kmax - array_float_len_send[1];
    MPI_Isend(array_float_send[0], num_p_send_0, MPI_FLOAT, npy[0], 2, mpi_comm_cart, req_plus+6); 
    MPI_Isend(array_char_send[0], num_c_send_0, MPI_CHAR, npy[0], 3, mpi_comm_cart, req_plus+7); 
    MPI_Isend(array_char_displacement_send[0], num_c_send_0, MPI_INT, npy[0], 4, mpi_comm_cart, req_plus+8); 
    MPI_Isend(array_float_send[1], num_p_send_1, MPI_FLOAT, npy[1], 5, mpi_comm_cart, req_plus+9); 
    MPI_Isend(array_char_send[1], num_c_send_1, MPI_CHAR, npy[1], 6, mpi_comm_cart, req_plus+10); 
    MPI_Isend(array_char_displacement_send[1], num_c_send_1, MPI_INT, npy[1], 7, mpi_comm_cart, req_plus+11); 
    MPI_Waitall(12, req_plus, st_plus);

    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*kmax)) + 1.0*((float)num_p_send_0/(imax*kmax));
    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_1/(imax*kmax)) + 1.0*((float)num_p_send_1/(imax*kmax));
    cr_num += 2;

    //calculate bitwise compress ratio
    // printf("%d %d\n", num_p_send_0, num_c_send_0);
    // printf("%d %d\n", num_p_send_1, num_c_send_1);
    // printf("%f \n", calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0));
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax));
    // printf("%f \n", (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax)));
    // printf("%f \n", cr);
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_1/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[1], num_p_send_1)*((float)num_p_send_1/(imax*jmax));
    // printf("%f \n", cr);
    // cr_num += 2;
    // printf("%d \n", cr_num);

    float* decompressed_data_0 = myDecompress(array_float_recv[0], array_char_recv[0], array_char_displacement_recv[0], imax*kmax);
    int pointer_0 = 0;
    // for(int a=0; a<imax; a++)
    // {
    //   for(int b=0; b<jmax; b++)
    //   {
    //     p[a][b][kmax-1] = decompressed_data_0[pointer_0++];
    //   }
    // }
    float* decompressed_data_1 = myDecompress(array_float_recv[1], array_char_recv[1], array_char_displacement_recv[1], imax*kmax);
    int pointer_1 = 0;
    for(int a=0; a<imax; a++)
    {
      for(int b=0; b<kmax; b++)
      {
        p[a][jmax-1][b] = decompressed_data_0[pointer_0++];
        p[a][0][b] = decompressed_data_1[pointer_1++];
      }
    }
  }  

  if(CT == 0)
  {
    MPI_Irecv(&p[0][jmax-1][0],
              1,
              ikvec,
              npy[1],
              1,
              mpi_comm_cart,
              req);
    MPI_Irecv(&p[0][0][0],
              1,
              ikvec,
              npy[0],
              2,
              mpi_comm_cart,
              req+1);
    //todo
    // if (CT == 1)
    // {
    //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 2, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 2, jmax-2, imax, jmax, kmax); 
    // }
    // else if(CT == 2)
    // {
    //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 2, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 2, jmax-2, imax, jmax, kmax); 
    // } 
    // else if(CT == 3)    
    // {
    //   cr += calcCompressionRatio_himeno_nolossy_area(p, 2, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_nolossy_area(p, 2, jmax-2, imax, jmax, kmax);     
    // }   
    // else if(CT == 4)    
    // {
    //   cr += calcCompressionRatio_himeno_sz(p, 2, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_sz(p, 2, jmax-2, imax, jmax, kmax);     
    // }      
    // cr_num += 2;         
    MPI_Isend(&p[0][1][0],
              1,
              ikvec,
              npy[0],
              1,
              mpi_comm_cart,
              req+2);
    MPI_Isend(&p[0][jmax-2][0],
              1,
              ikvec,
              npy[1],
              2,
              mpi_comm_cart,
              req+3);

    MPI_Waitall(4,
                req,
                st);
  }


}


void
sendp1()
{
  MPI_Status  st[4];
  MPI_Request req[4];

  //todo
  MPI_Status   st_bitmask[16];
  MPI_Request  req_bitmask[16];

  MPI_Status   st_plus[12];
  MPI_Request  req_plus[12];

  MPI_Status   st_bitwise[8];
  MPI_Request  req_bitwise[8];  

  MPI_Status   st_sz[4];
  MPI_Request  req_sz[4];   

  if(CT == 7)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npx[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npx[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 1, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 1, imax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], jmax*kmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], jmax*kmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    int type[2] = {0, 0};
    float medium[2];
    medium[0] = med_dataset_float(data_small[0], jmax*kmax, &type[0]);
    medium[1] = med_dataset_float(data_small[1], jmax*kmax, &type[1]);
    char float_arr0[32+1];
    char float_arr1[32+1];
    floattostr(&medium[0], float_arr0);
    floattostr(&medium[1], float_arr1);
    char mask0[1+8+8];
    char mask1[1+8+8];
    strncpy(mask0, float_arr0, 1+8+8);
    strncpy(mask1, float_arr1, 1+8+8);

    myCompress_bitwise_mask(data_small[0], jmax*kmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0], type[0], mask0);
    myCompress_bitwise_mask(data_small[1], jmax*kmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1], type[1], mask1);    

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npx[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npx[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    cr += data_bytes_send[0]*8.0/(jmax*kmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(jmax*kmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];
    float medium_recv[2];
    int type_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npx[1], 2, mpi_comm_cart, req_bitmask);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npx[1], 3, mpi_comm_cart, req_bitmask+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npx[0], 4, mpi_comm_cart, req_bitmask+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npx[0], 5, mpi_comm_cart, req_bitmask+3);  
    MPI_Irecv(&medium_recv[0], 1, MPI_FLOAT, npx[1], 6, mpi_comm_cart, req_bitmask+4);
    MPI_Irecv(&medium_recv[1], 1, MPI_FLOAT, npx[0], 7, mpi_comm_cart, req_bitmask+5);
    MPI_Irecv(&type_recv[0], 1, MPI_INT, npx[1], 8, mpi_comm_cart, req_bitmask+6);
    MPI_Irecv(&type_recv[1], 1, MPI_INT, npx[0], 9, mpi_comm_cart, req_bitmask+7);    
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npx[0], 2, mpi_comm_cart, req_bitmask+8); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npx[0], 3, mpi_comm_cart, req_bitmask+9); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npx[1], 4, mpi_comm_cart, req_bitmask+10); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npx[1], 5, mpi_comm_cart, req_bitmask+11); 
    MPI_Isend(&medium[0], 1, MPI_FLOAT, npx[0], 6, mpi_comm_cart, req_bitmask+12); 
    MPI_Isend(&medium[1], 1, MPI_FLOAT, npx[1], 7, mpi_comm_cart, req_bitmask+13); 
    MPI_Isend(&type[0], 1, MPI_INT, npx[0], 8, mpi_comm_cart, req_bitmask+14); 
    MPI_Isend(&type[1], 1, MPI_INT, npx[1], 9, mpi_comm_cart, req_bitmask+15);     
    MPI_Waitall(16, req_bitmask, st_bitmask);    

    char float_arr0_recv[32+1];
    char float_arr1_recv[32+1];
    floattostr(&medium_recv[0], float_arr0_recv);
    floattostr(&medium_recv[1], float_arr1_recv);
    char mask0_recv[1+8+8];
    char mask1_recv[1+8+8];
    strncpy(mask0_recv, float_arr0_recv, 1+8+8);
    strncpy(mask1_recv, float_arr1_recv, 1+8+8);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise_mask(data_bits_recv[0], data_bytes_recv[0], jmax*kmax, type_recv[0], mask0_recv);
    decompressed_data[1] = myDecompress_bitwise_mask(data_bits_recv[1], data_bytes_recv[1], jmax*kmax, type_recv[1], mask1_recv);

    int pointer = 0;
    for(int a=0; a<jmax; a++)
    {
      for(int b=0; b<kmax; b++)
      {
        p[imax-1][a][b] = decompressed_data[0][pointer] + data_min_recv[0];
        p[0][a][b] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  }

  if(CT == 6)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npx[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npx[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 1, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 1, imax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], jmax*kmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], jmax*kmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    myCompress_bitwise_np(data_small[0], jmax*kmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0]);
    myCompress_bitwise_np(data_small[1], jmax*kmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1]);

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npx[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npx[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    // printf("npx %d %d \n", npx[0], npx[1]);  
    // printf("send %d %d \n", data_bytes_send[0], data_bytes_send[1]);
    // printf("recv %d %d \n", data_bytes_recv[0], data_bytes_recv[1]);

    cr += data_bytes_send[0]*8.0/(jmax*kmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(jmax*kmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npx[1], 2, mpi_comm_cart, req_bitwise);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npx[1], 3, mpi_comm_cart, req_bitwise+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npx[0], 4, mpi_comm_cart, req_bitwise+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npx[0], 5, mpi_comm_cart, req_bitwise+3);  
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npx[0], 2, mpi_comm_cart, req_bitwise+4); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npx[0], 3, mpi_comm_cart, req_bitwise+5); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npx[1], 4, mpi_comm_cart, req_bitwise+6); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npx[1], 5, mpi_comm_cart, req_bitwise+7); 
    MPI_Waitall(8, req_bitwise, st_bitwise);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise_np(data_bits_recv[0], data_bytes_recv[0], jmax*kmax);
    decompressed_data[1] = myDecompress_bitwise_np(data_bits_recv[1], data_bytes_recv[1], jmax*kmax);
    int pointer = 0;
    for(int a=0; a<jmax; a++)
    {
      for(int b=0; b<kmax; b++)
      {
        p[imax-1][a][b] = decompressed_data[0][pointer] + data_min_recv[0];
        p[0][a][b] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  }  

  if(CT == 5)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npx[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npx[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 1, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 1, imax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], jmax*kmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], jmax*kmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    myCompress_bitwise(data_small[0], jmax*kmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0]);
    myCompress_bitwise(data_small[1], jmax*kmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1]);

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npx[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npx[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    // printf("npx %d %d \n", npx[0], npx[1]);  
    // printf("send %d %d \n", data_bytes_send[0], data_bytes_send[1]);
    // printf("recv %d %d \n", data_bytes_recv[0], data_bytes_recv[1]);

    cr += data_bytes_send[0]*8.0/(jmax*kmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(jmax*kmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npx[1], 2, mpi_comm_cart, req_bitwise);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npx[1], 3, mpi_comm_cart, req_bitwise+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npx[0], 4, mpi_comm_cart, req_bitwise+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npx[0], 5, mpi_comm_cart, req_bitwise+3);  
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npx[0], 2, mpi_comm_cart, req_bitwise+4); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npx[0], 3, mpi_comm_cart, req_bitwise+5); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npx[1], 4, mpi_comm_cart, req_bitwise+6); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npx[1], 5, mpi_comm_cart, req_bitwise+7); 
    MPI_Waitall(8, req_bitwise, st_bitwise);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise(data_bits_recv[0], data_bytes_recv[0], jmax*kmax);
    decompressed_data[1] = myDecompress_bitwise(data_bits_recv[1], data_bytes_recv[1], jmax*kmax);
    int pointer = 0;
    for(int a=0; a<jmax; a++)
    {
      for(int b=0; b<kmax; b++)
      {
        p[imax-1][a][b] = decompressed_data[0][pointer] + data_min_recv[0];
        p[0][a][b] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  }  

  if(CT == 4)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npx[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npx[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 1, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 1, imax-2, imax, jmax, kmax);   

    unsigned char* data_bits_send[2] = {NULL, NULL};

    char binfile0[64];
    sprintf(binfile0, "dataset/r%dp1data0.dat", id); 
    char binfile1[64];
    sprintf(binfile1, "dataset/r%dp1data1.dat", id);     
    writetobinary_float(binfile0, data[0], jmax*kmax); //.txt --> .dat
    writetobinary_float(binfile1, data[1], jmax*kmax); //.txt --> .dat    
    char sz_comp_cmd0[64];
    char sz_comp_cmd1[64];
    sprintf(sz_comp_cmd0, "%s%g%sdataset/r%dp1data0%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, id, sz_comp_cmd_suffix2, jmax*kmax);
    sprintf(sz_comp_cmd1, "%s%g%sdataset/r%dp1data1%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, id, sz_comp_cmd_suffix2, jmax*kmax);
    //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
    int iret_comp0 = system(sz_comp_cmd0); //.dat --> .dat.sz
    int iret_comp1 = system(sz_comp_cmd1); //.dat --> .dat.sz
    char binfile_sz0[64];
    sprintf(binfile_sz0, "dataset/r%dp1data0.dat.sz", id);
    char binfile_sz1[64];
    sprintf(binfile_sz1, "dataset/r%dp1data1.dat.sz", id);    
    data_bits_send[0] = readfrombinary_char(binfile_sz0, &data_bytes_send[0]);  
    data_bits_send[1] = readfrombinary_char(binfile_sz1, &data_bytes_send[1]);   

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npx[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npx[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    cr += data_bytes_send[0]*8.0/(jmax*kmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(jmax*kmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npx[1], 2, mpi_comm_cart, req_sz);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npx[0], 3, mpi_comm_cart, req_sz+1);  
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npx[0], 2, mpi_comm_cart, req_sz+2); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npx[1], 3, mpi_comm_cart, req_sz+3); 
    MPI_Waitall(4, req_sz, st_sz);    

    char binfile_zs0[64];
    sprintf(binfile_zs0, "dataset/r%dp1data0.dat.zs", id);
    char binfile_zs1[64];
    sprintf(binfile_zs1, "dataset/r%dp1data1.dat.zs", id);    
    writetobinary_char(binfile_zs0, data_bits_recv[0], data_bytes_recv[0]); //.dat.zs
    writetobinary_char(binfile_zs1, data_bits_recv[1], data_bytes_recv[1]); //.dat.zs
    char sz_decomp_cmd0[64];
    char sz_decomp_cmd1[64];
    sprintf(sz_decomp_cmd0, "%sdataset/r%dp1data0%s%d", sz_decomp_cmd_prefix, id, sz_decomp_cmd_suffix, jmax*kmax);
    sprintf(sz_decomp_cmd1, "%sdataset/r%dp1data1%s%d", sz_decomp_cmd_prefix, id, sz_decomp_cmd_suffix, jmax*kmax);
    //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
    int iret_decomp0 = system(sz_decomp_cmd0); //.dat.zs --> .dat.zs.out
    int iret_decomp1 = system(sz_decomp_cmd1); //.dat.zs --> .dat.zs.out 
    char binfile_out0[64];
    sprintf(binfile_out0, "dataset/r%dp1data0.dat.zs.out", id);
    char binfile_out1[64];
    sprintf(binfile_out1, "dataset/r%dp1data1.dat.zs.out", id);
    char txtfile0[64];
    sprintf(txtfile0, "dataset/r%dp1data0.dat.zs.out.txt", id);  
    char txtfile1[64];
    sprintf(txtfile1, "dataset/r%dp1data1.dat.zs.out.txt", id);       

    float* decompressed_data[2];
    decompressed_data[0] = readfrombinary_writetotxt_float(binfile_out0, txtfile0, jmax*kmax);
    decompressed_data[1] = readfrombinary_writetotxt_float(binfile_out1, txtfile1, jmax*kmax);

    int pointer = 0;
    for(int a=0; a<jmax; a++)
    {
      for(int b=0; b<kmax; b++)
      {
        p[imax-1][a][b] = decompressed_data[0][pointer];
        p[0][a][b] = decompressed_data[1][pointer];
        pointer++;
      }
    }
  }   

  if(CT == 1)
  {
    int array_float_len_send[2] = {0, 0};
    int array_float_len_recv[2] = {0, 0};

    float* array_float_send[2] = {NULL, NULL};
    char* array_char_send[2] = {NULL, NULL};
    int* array_char_displacement_send[2] = {NULL, NULL};

    MPI_Irecv(&array_float_len_recv[0], 1, MPI_INT, npx[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&array_float_len_recv[1], 1, MPI_INT, npx[0], 1, mpi_comm_cart, req+1);
    float* data_0 = transform_3d_array_to_1d_array(p, 1, 1, imax, jmax, kmax);
    float* data_1 = transform_3d_array_to_1d_array(p, 1, imax-2, imax, jmax, kmax);
    array_float_len_send[0] = myCompress(data_0, &array_float_send[0], &array_char_send[0], &array_char_displacement_send[0], jmax*kmax);
    array_float_len_send[1] = myCompress(data_1, &array_float_send[1], &array_char_send[1], &array_char_displacement_send[1], jmax*kmax);
    MPI_Isend(&array_float_len_send[0], 1, MPI_INT, npx[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&array_float_len_send[1], 1, MPI_INT, npx[1], 1, mpi_comm_cart, req+3); 
    //printf("send1 %d %d \n", array_float_len_send[0], array_float_len_send[1]);
    MPI_Waitall(4, req, st);

    //printf("recv1 %d %d \n", array_float_len_recv[0], array_float_len_recv[1]);

    float* array_float_recv[2]; 
    char* array_char_recv[2]; 
    int* array_char_displacement_recv[2]; 

    int num_p_recv_0 = array_float_len_recv[0], num_c_recv_0 = jmax*kmax - array_float_len_recv[0];
    int num_p_recv_1 = array_float_len_recv[1], num_c_recv_1 = jmax*kmax - array_float_len_recv[1];
    array_float_recv[0] = (float*) malloc(sizeof(float)*num_p_recv_0);
    array_char_recv[0] = (char*) malloc(sizeof(char)*num_c_recv_0);
    array_char_displacement_recv[0] = (int*) malloc(sizeof(int)*num_c_recv_0);
    array_float_recv[1] = (float*) malloc(sizeof(float)*num_p_recv_1);
    array_char_recv[1] = (char*) malloc(sizeof(char)*num_c_recv_1);
    array_char_displacement_recv[1] = (int*) malloc(sizeof(int)*num_c_recv_1);
    MPI_Irecv(array_float_recv[0], num_p_recv_0, MPI_FLOAT, npx[1], 2, mpi_comm_cart, req_plus);
    MPI_Irecv(array_char_recv[0], num_c_recv_0, MPI_CHAR, npx[1], 3, mpi_comm_cart, req_plus+1);
    MPI_Irecv(array_char_displacement_recv[0], num_c_recv_0, MPI_INT, npx[1], 4, mpi_comm_cart, req_plus+2);    
    MPI_Irecv(array_float_recv[1], num_p_recv_1, MPI_FLOAT, npx[0], 5, mpi_comm_cart, req_plus+3);
    MPI_Irecv(array_char_recv[1], num_c_recv_1, MPI_CHAR, npx[0], 6, mpi_comm_cart, req_plus+4);
    MPI_Irecv(array_char_displacement_recv[1], num_c_recv_1, MPI_INT, npx[0], 7, mpi_comm_cart, req_plus+5);  
    //mycompress
    int num_p_send_0 = array_float_len_send[0], num_c_send_0 = jmax*kmax - array_float_len_send[0];
    int num_p_send_1 = array_float_len_send[1], num_c_send_1 = jmax*kmax - array_float_len_send[1];
    MPI_Isend(array_float_send[0], num_p_send_0, MPI_FLOAT, npx[0], 2, mpi_comm_cart, req_plus+6); 
    MPI_Isend(array_char_send[0], num_c_send_0, MPI_CHAR, npx[0], 3, mpi_comm_cart, req_plus+7); 
    MPI_Isend(array_char_displacement_send[0], num_c_send_0, MPI_INT, npx[0], 4, mpi_comm_cart, req_plus+8); 
    MPI_Isend(array_float_send[1], num_p_send_1, MPI_FLOAT, npx[1], 5, mpi_comm_cart, req_plus+9); 
    MPI_Isend(array_char_send[1], num_c_send_1, MPI_CHAR, npx[1], 6, mpi_comm_cart, req_plus+10); 
    MPI_Isend(array_char_displacement_send[1], num_c_send_1, MPI_INT, npx[1], 7, mpi_comm_cart, req_plus+11); 
    MPI_Waitall(12, req_plus, st_plus);

    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_0/(jmax*kmax)) + 1.0*((float)num_p_send_0/(jmax*kmax));
    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_1/(jmax*kmax)) + 1.0*((float)num_p_send_1/(jmax*kmax));
    cr_num += 2;

    //calculate bitwise compress ratio
    // printf("%d %d\n", num_p_send_0, num_c_send_0);
    // printf("%d %d\n", num_p_send_1, num_c_send_1);
    // printf("%f \n", calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0));
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax));
    // printf("%f \n", (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax)));
    // printf("%f \n", cr);
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_1/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[1], num_p_send_1)*((float)num_p_send_1/(imax*jmax));
    // printf("%f \n", cr);
    // cr_num += 2;
    // printf("%d \n", cr_num);

    float* decompressed_data_0 = myDecompress(array_float_recv[0], array_char_recv[0], array_char_displacement_recv[0], jmax*kmax);
    int pointer_0 = 0;
    // for(int a=0; a<imax; a++)
    // {
    //   for(int b=0; b<jmax; b++)
    //   {
    //     p[a][b][kmax-1] = decompressed_data_0[pointer_0++];
    //   }
    // }
    float* decompressed_data_1 = myDecompress(array_float_recv[1], array_char_recv[1], array_char_displacement_recv[1], jmax*kmax);
    int pointer_1 = 0;
    for(int a=0; a<jmax; a++)
    {
      for(int b=0; b<kmax; b++)
      {
        p[imax-1][a][b] = decompressed_data_0[pointer_0++];
        p[0][a][b] = decompressed_data_1[pointer_1++];
      }
    }
  } 

  if(CT == 0)
  {
    MPI_Irecv(&p[imax-1][0][0],
              1,
              jkvec,
              npx[1],
              1,
              mpi_comm_cart,
              req);
    MPI_Irecv(&p[0][0][0],
              1,
              jkvec,
              npx[0],
              2,
              mpi_comm_cart,
              req+1);
    //todo
    // if (CT == 1)
    // {
    //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 1, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 1, imax-2, imax, jmax, kmax); 
    // }
    // else if(CT == 2)
    // {
    //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 1, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 1, imax-2, imax, jmax, kmax); 
    // } 
    // else if(CT == 3)    
    // {
    //   cr += calcCompressionRatio_himeno_nolossy_area(p, 1, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_nolossy_area(p, 1, imax-2, imax, jmax, kmax);     
    // }    
    // else if(CT == 4)    
    // {
    //   cr += calcCompressionRatio_himeno_sz(p, 1, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_sz(p, 1, imax-2, imax, jmax, kmax);     
    // }   
    // cr_num += 2;               
    MPI_Isend(&p[1][0][0],
              1,
              jkvec,
              npx[0],
              1,
              mpi_comm_cart,
              req+2);
    MPI_Isend(&p[imax-2][0][0],
              1,
              jkvec,
              npx[1],
              2,
              mpi_comm_cart,
              req+3);

    MPI_Waitall(4,
                req,
                st);
    }


}
            
