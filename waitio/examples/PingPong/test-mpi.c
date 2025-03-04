/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil -*- */
#include <mpi.h>
#include "waitio.h"
#include "waitio_mpi.h"

int main (int argc, char *argv[]) {
  int  data[2], *buf = data;
  waitio_group_t grp1;
  int ret;
  int wrank;
  
  if((ret = MPI_Init( &argc, &argv)) != MPI_SUCCESS) {
      fprintf(stderr, "MPI_Init failed code %d\n", ret);
      exit(ret);
  }

  waitio_create_universe (&grp1);
  if(grp1 == NULL) {
      fprintf(stderr, "waitio_create_universe failed code %d\n", ret);
      MPI_Finalize();
      exit(ret);
  }

  waitio_group_rank(grp1, &wrank);
  if((wrank%2) == 1) {
      if((ret = waitio_mpi_recv(buf, 4, WAITIO_MPI_CHAR, wrank-1, 0, grp1)) != 0) {
          fprintf(stderr, "%d waitio_mpi_recv error %d\n", wrank, ret);
      }
      else {
          fprintf(stderr, "rank%d: waitio_mpi_recv from rank%d\n", wrank,  wrank-1);
      }
      if((ret = waitio_mpi_send(buf, 4, WAITIO_MPI_CHAR, wrank-1, 0, grp1)) != 0) {
          fprintf(stderr, "%d waitio_send error %d\n", wrank, ret);
      }
      else {
          fprintf(stderr, "rank%d: waitio_mpi_send to rank%d\n", wrank,  wrank-1);
      }
  }
  else {
      if((ret = waitio_mpi_send(buf, 4, WAITIO_MPI_CHAR, wrank+1, 0, grp1)) != 0) {
          fprintf(stderr, "%d waitio_send error %d\n", wrank, ret);
      }
      else {
          fprintf(stderr, "rank%d: waitio_mpi_send to rank%d\n", wrank,  wrank+1);
      }
      if((ret = waitio_mpi_recv(buf, 4, WAITIO_MPI_CHAR, wrank+1, 0, grp1)) != 0) {
          fprintf(stderr, "%d waitio_recv error %d\n", wrank, ret);
      }
      else {
          fprintf(stderr, "rank%d: waitio_mpi_recv from rank%d\n", wrank,  wrank+1);
      }
  }
  waitio_finalize();
  MPI_Finalize();
}
