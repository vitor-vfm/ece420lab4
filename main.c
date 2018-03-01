#define LAB4_EXTEND

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Lab4_IO.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

#define THRESHOLD 0.0001

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

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

    const size_t chunk_size = nodecount / comm_size;
    const size_t start = rank * chunk_size;
    const size_t end = start + chunk_size <= nodecount ? start + chunk_size : nodecount;

    r = malloc(nodecount * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));
    for ( i = 0; i < nodecount; ++i)
        r[i] = 1.0 / nodecount;
    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;

    // CORE CALCULATION
    do{
        ++iterationcount;
        vec_cp(r, r_pre, nodecount);

        // Calculate for chunk
        for (i = start; i < end; ++i){
            r[i] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j)
                r[i] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            r[i] *= DAMPING_FACTOR;
            r[i] += damp_const;
        }

        // send the chunk to other processes

        // read the chunks of other processes

    } while(rel_error(r, r_pre, nodecount) >= EPSILON);

    // post processing
    node_destroy(nodehead, nodecount);
    free(num_in_links);
    free(num_out_links);

    MPI_Finalize();
}
