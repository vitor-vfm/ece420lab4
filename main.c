#define LAB4_EXTEND

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "Lab4_IO.h"
#include "timer.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

#define THRESHOLD 0.0001

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    double start_time, end_time;

    int rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    struct node *nodehead;
    int nodecount;
    int *num_in_links, *num_out_links;
    double *r, *r_pre;
    int i, j;
    double damp_const;
    int iterationcount = 0;

    if (get_node_stat(&nodecount, &num_in_links, &num_out_links)) return 254;

    // Calculate the result
    if (node_init(&nodehead, num_in_links, num_out_links, 0, nodecount)) return 254;

    GET_TIME(start_time);

    const unsigned chunk_size = nodecount / comm_size;
    const unsigned start = rank * chunk_size;
    const unsigned end = start + chunk_size <= nodecount ? start + chunk_size : nodecount;

    r = malloc(nodecount * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));
    double *local_buf = calloc(sizeof(double), chunk_size);
    for ( i = 0; i < nodecount; ++i)
        r[i] = 1.0 / nodecount;
    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;

    // CORE CALCULATION
    do{
        ++iterationcount;
        vec_cp(r, r_pre, nodecount);

        // Calculate for chunk
        for (i = start; i < end; ++i){
            local_buf[i-start] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j)
                local_buf[i-start] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            local_buf[i-start] *= DAMPING_FACTOR;
            local_buf[i-start] += damp_const;
        }

        // send the chunk to other processes
        MPI_Allgather(local_buf, end-start, MPI_DOUBLE, r, end-start, MPI_DOUBLE, MPI_COMM_WORLD);

    } while(rel_error(r, r_pre, nodecount) >= EPSILON);

    GET_TIME(end_time);

    free(local_buf);

    if (rank == 0) {
        Lab4_saveoutput(r, nodecount, end_time - start_time);
    }

    // post processing
    node_destroy(nodehead, nodecount);
    free(num_in_links);
    free(num_out_links);

    MPI_Finalize();
}
