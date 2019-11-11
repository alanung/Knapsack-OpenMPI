/* Knapsack calculation based on that of */
/* https://www.tutorialspoint.com/cplusplus-program-to-solve-knapsack-problem-using-dynamic-programming */

#include <stdio.h>
#include <sys/time.h>
#include <stdint.h>
#include <mpi.h>

long int knapSack(long int C, long int w[], long int v[], int n);

uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t) 1000000 + tv.tv_usec;
}

int main(int argc, char *argv[]) {
    long int C;    /* capacity of backpack */
    int n;    /* number of items */
    int i;    /* loop counter */

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        scanf("%ld", &C);
        scanf("%d", &n);
    }

    // broadcast problem parameters to all ranks
    MPI_Bcast(&C, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    long int v[n], w[n];        /* value, weight */

    if (rank == 0) {
        for (i = 0; i < n; i++) {
            scanf("%ld %ld", &v[i], &w[i]);
        }
    }

    MPI_Bcast(w, n, MPI_LONG, 0, MPI_COMM_WORLD);
//    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(v, n, MPI_LONG, 0, MPI_COMM_WORLD);

    uint64_t start = GetTimeStamp();
    long int ks = knapSack(C, w, v, n);

    if (rank == 0) {
        printf("knapsack occupancy %ld\n", ks);
        printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
    }

    MPI_Finalize();

    return 0;
}

/* PLACE YOUR CHANGES BELOW HERE */
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

long int max(long int x, long int y) {
    return (x > y) ? x : y;
}

long int knapSack(long int C, long int w[], long int v[], int n) {
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // block width (# columns allocated to each rank)
    long int width = (long int) ceil((double) (C + 1) / nprocs);
//    if (rank == 0)
//        printf("width: %ld\n", width);
    // column offset for this rank
    long int offset = rank * width;
//    printf("rank: %d, offset: %ld\n", rank, offset);

    // allocate memory to store row of DP memorisation table
    long int *K = calloc(width * nprocs, sizeof *K);
    assert(K);
    // allocate memory as buffer for new block results
    long int *temp = malloc(width * sizeof *temp);
    assert(temp);

//    if (rank == 0){
//        printf("                          ");
//        for (long int i = 0; i < width * nprocs; i++) {
//            if (i == C)
//                printf("[");
//            printf("%3ld ", i);
//            if (i == C)
//                printf("]");
//        }
//        printf("\n");
//    }

//    omp_set_num_threads(2);
    for (long int item = 0; item < n; item++) {
        // determine new block values from previously broadcast results

        for (long int col = offset; col < offset + width; col++) {
            // if knapsack capacity is 0, no value possible
            if (col == 0)
                temp[col - offset] = 0;
                // else if the item can fit, compare including and excluding the item
            else if (w[item] <= col)
                temp[col - offset] = max(v[item] + K[col - w[item]], K[col]);
                // else (if the item cannot fit), the max value is unchanged
            else
                temp[col - offset] = K[col];
        }
        MPI_Allgather(temp, (int) width, MPI_LONG, K, (int) width, MPI_LONG, MPI_COMM_WORLD);
//        if (rank == 0) {
//            printf("[weight: %3ld, value: %3ld] ", w[item], v[item]);
//            for (long int i = 0; i < width * nprocs; i++) {
//                if (i == C)
//                    printf("[");
//                printf("%3ld ", K[i]);
//                if (i == C)
//                    printf("]");
//            }
//            printf("\n");
//            fflush(stdout);
//        }
    }

    long int result = K[C];
    free(K);
    return result;
}