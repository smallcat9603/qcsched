/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil -*- */
#include <mpi.h>
#include "waitio.h"

int truef(int pbid, int n) { return 1; }

int main (int argc, char *argv[]) {
    int ret, wrank;
    waitio_filter_func_t func[4]= {truef, truef, NULL, NULL};
    int  array[4] = {1, 2, 0, 0};
    int  data[2], *buf = data;
    waitio_req_t req;
    waitio_group_t grp1;
  
    if((ret = MPI_Init( &argc, &argv)) != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Init failed code %d\n", ret);
        exit(ret);
    }
    if((ret = waitio_init(10)) != 0) {  
        fprintf(stderr, "waitio_init failed code %d\n", ret);
        MPI_Abort(MPI_COMM_WORLD, ret);
        exit(ret);
    }

    grp1 = waitio_create_group(0, func, array);
    if(grp1 == NULL) {
        fprintf(stderr, "waitio_create_group failed code %d\n", ret);
        MPI_Finalize();
        exit(ret);
    }

    waitio_group_rank(grp1, &wrank);
    if((wrank%2) == 1) {
        waitio_irecv(grp1, wrank-1, (char *)buf, 4, 0, &req);
        if((ret = waitio_wait(&req)) != 0) {
            fprintf(stderr, "%d waitio_irecv error %d\n", wrank, ret);
        }
        else {
            fprintf(stderr, "rank%d: waitio_irecv from rank%d\n", wrank,  wrank-1);
        }
        waitio_isend(grp1, wrank-1, (char *)buf, 4, 0, &req);
        if((ret = waitio_wait(&req)) != 0) {
            fprintf(stderr, "%d waitio_isend error %d\n", wrank, ret);
        }
        else {
            fprintf(stderr, "rank%d: waitio_isend to rank%d\n", wrank,  wrank-1);
        }
    }
    else {
        waitio_isend(grp1, wrank+1, (char *)buf, 4, 0, &req);
        if((ret = waitio_wait(&req)) != 0) {
            fprintf(stderr, "%d waitio_isend error %d\n", wrank, ret);
        }
        else {
            fprintf(stderr, "rank%d: waitio_isend to rank%d\n", wrank,  wrank+1);
        }
        waitio_irecv(grp1, wrank+1, (char *)&buf, 4, 0, &req);
        if((ret = waitio_wait(&req)) != 0) {
            fprintf(stderr, "%d waitio_irecv error %d\n", wrank, ret);
        }
        else {
            fprintf(stderr, "rank%d: waitio_irecv from rank%d\n", wrank,  wrank+1);
        }
    }
    waitio_finalize();
    MPI_Finalize();
}
