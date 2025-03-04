/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil -*- */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main (int argc, char *argv[]) {
    int ret, rank;
    int  data[2], *buf = data;
    MPI_Request req;
    MPI_Status sta;
  
    if((ret = MPI_Init( &argc, &argv)) != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Init failed code %d\n", ret);
        exit(ret);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if((rank%2) == 1) {
        MPI_Irecv((char *)buf, 4, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, &req);
        if((ret = MPI_Wait(&req, &sta)) != 0) {
            fprintf(stderr, "%d MPI_Irecv error %d\n", rank, ret);
        }
        else {
            fprintf(stderr, "rank%d: MPI_Irecv from rank%d\n", rank,  rank-1);
        }
        MPI_Isend((char *)buf, 4, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, &req);
        if((ret = MPI_Wait(&req, &sta)) != 0) {
            fprintf(stderr, "%d MPI_Isend error %d\n", rank, ret);
        }
        else {
            fprintf(stderr, "rank%d: MPI_Isend to rank%d\n", rank,  rank-1);
        }
    }
    else {
        MPI_Isend((char *)buf, 4, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD, &req);
        if((ret = MPI_Wait(&req, &sta)) != 0) {
            fprintf(stderr, "%d MPI_Isend error %d\n", rank, ret);
        }
        else {
            fprintf(stderr, "rank%d: MPI_Isend to rank%d\n", rank,  rank+1);
        }
        MPI_Irecv((char *)buf, 4, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD, &req);
        if((ret = MPI_Wait(&req, &sta)) != 0) {
            fprintf(stderr, "%d MPI_Irecv error %d\n", rank, ret);
        }
        else {
            fprintf(stderr, "rank%d: MPI_Irecv from rank%d\n", rank,  rank+1);
        }
    }
    MPI_Finalize();
}
