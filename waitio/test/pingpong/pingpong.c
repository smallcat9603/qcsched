/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil -*- */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <mpi.h>

int main (int argc, char *argv[]) {
  int  data[2], *buf = data;
  MPI_Status status;
  int ret;
  int rank;
  
  if((ret = MPI_Init( &argc, &argv)) != MPI_SUCCESS) {
      fprintf(stderr, "MPI_Init failed code %d\n", ret);
      exit(ret);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if((rank%2) == 1) {
      if((ret = MPI_Recv(buf, 4,MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, &status)) != 0) {
          fprintf(stderr, "%d MPI_Recv error %d\n", rank, ret);
      }
      else {
          fprintf(stderr, "rank%d: MPI_Recv from rank%d\n", rank,  rank-1);
      }
      if((ret = MPI_Send(buf, 4, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD)) != 0) {
          fprintf(stderr, "%d MPI_Send error %d\n", rank, ret);
      }
      else {
          fprintf(stderr, "rank%d: MPI_Send to rank%d\n", rank,  rank-1);
      }
  }
  else {
      if((ret = MPI_Send(buf, 4, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD)) != 0) {
          fprintf(stderr, "%d MPI_Send error %d\n", rank, ret);
      }
      else {
          fprintf(stderr, "rank%d: MPI_Send to rank%d\n", rank,  rank+1);
      }
      if((ret = MPI_Recv(buf, 4, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD, &status)) != 0) {
          fprintf(stderr, "%d MPI_Recv error %d\n", rank, ret);
      }
      else {
          fprintf(stderr, "rank%d: MPI_Recv from rank%d\n", rank,  rank+1);
      }
  }
  MPI_Finalize();
}
