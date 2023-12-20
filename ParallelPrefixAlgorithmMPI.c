#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv) {
    int num_procs = 15;
    int rank;
    int proc_n = num_procs + 1;
    int n = 100; 
    int k = log2(proc_n);
    int num_leaves = proc_n / 2;
    int m = num_leaves;
    int block = (n % num_leaves ? n / num_leaves + 1 : n / num_leaves);
    int block_dim = num_procs * block;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int L[m], R[m], S[m], T[m];
    
    int *local_x = (int *)malloc(block * sizeof(int));
    int *local_y = (int *)malloc(block * sizeof(int));
    int *local_y2 = (int *)malloc(block * sizeof(int));
    int *x = (int *)malloc(block_dim * sizeof(int)); 
    int *y = (int *)malloc(block_dim * sizeof(int));
    int *proc_ranks;
    MPI_Group world_group, leaves_group;
    MPI_Comm leaves_comm;
    proc_ranks = (int *)malloc(num_leaves * sizeof(int));
    for (int i = 0; i < num_leaves; i++) {
        proc_ranks[i] = num_leaves - 1 + i;
    }
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_incl(world_group, num_leaves, proc_ranks, &leaves_group);
    MPI_Comm_create(MPI_COMM_WORLD, leaves_group, &leaves_comm);

    if(rank == num_leaves - 1) {
        for (int i = 0; i < block_dim; ++i) {
            if (i >= n) {
                x[i] = 0;
            }
            else {
                x[i] = 1;
            }
            y[i] = 0;
        }
    }

    for (int i = 0; i < block; ++i) {
        local_x[i] = 0;
    }

    if (MPI_COMM_NULL != leaves_comm) {
        MPI_Scatter(x, block, MPI_INT, local_x, block, MPI_INT, 0, leaves_comm);
    }

    // Penjanje - 1. korak
    int sum = 0;
    for (int i = 0; i < block; ++i) {
        sum += local_x[i];
    }
    if (rank >= m - 1) {
        MPI_Send(&sum, 1, MPI_INT, (rank - 1) / 2, 0, MPI_COMM_WORLD);
    }

    // Penjanje - ostali koraci
    for (int j = 2; j <= k - 1; ++j) {
        m = m / 2;
        for (int i = 1; i <= m; ++i) {
            if (rank == i + m - 2) {
                MPI_Recv(&L[rank], 1, MPI_INT, 2 * rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&R[rank], 1, MPI_INT, 2 * rank + 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                S[rank] = L[rank] + R[rank];
                MPI_Send(&S[rank], 1, MPI_INT, (rank - 1) / 2, 0, MPI_COMM_WORLD);
            }
        }
    }

    // Spuštanje - k-ti korak, početak spuštanja
    m = m / 2;
    if (rank == 0) {
        MPI_Recv(&L[rank], 1, MPI_INT, 2 * rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&R[rank], 1, MPI_INT, 2 * rank + 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        S[rank] = L[rank] + R[rank];
        MPI_Send(&L[rank], 1, MPI_INT, 2 * rank + 2, 0, MPI_COMM_WORLD);
    }

    // Spuštanje - ostali koraci
    for (int j = 1; j <= k - 2; ++j) {
        m = 2 * m;
        if (rank == m - 1) {
            MPI_Send(&L[rank], 1, MPI_INT, 2 * rank + 2, 0, MPI_COMM_WORLD);
        }
        else {
            for (int i = 2; i <= m; ++i) {
                if (rank == i + m - 2) {
                    MPI_Recv(&T[rank], 1, MPI_INT, (rank - 1) / 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(&T[rank], 1, MPI_INT, 2 * rank + 1, 0, MPI_COMM_WORLD);
                    S[rank] = T[rank] + L[rank];
                    MPI_Send(&S[rank], 1, MPI_INT, 2 * rank + 2, 0, MPI_COMM_WORLD);
                }
            }
        }
    }

    m = 2 * m;
    for (int i = 2; i <= m; ++i) {
        if (rank == i + m - 2) {
            MPI_Recv(&T[rank], 1, MPI_INT, (rank - 1) / 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    int prefix_sum = 0;
    for (int i = 0; i < block; ++i) {
        prefix_sum += local_x[i];
        local_y2[i] = prefix_sum;
    }

    if (rank != num_leaves - 1) {
        for (int i = 0; i < block; ++i) {
            local_y[i] = local_y2[i] + T[rank];
        }
    }
    else {
        for (int i = 0; i < block; ++i) {
            local_y[i] = local_y2[i];
        }
    }

    if (MPI_COMM_NULL != leaves_comm) {
        MPI_Gather(local_y, block, MPI_INT, y, block, MPI_INT, 0, leaves_comm);

        if (rank == num_leaves - 1) {
            printf("y = ");
            for (int i = 0; i < n; i++) {
                printf("%d ", y[i]);
            }
            printf("\n");
        }
    }

    MPI_Finalize();

    if (rank == num_leaves - 1) {
        free(x);
        free(y);
    }
    
    free(proc_ranks);
    free(local_x);
    free(local_y);
    free(local_y2);

    return 0;
}