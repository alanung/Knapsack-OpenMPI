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

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    long int v[n], w[n];        /* value, weight */

    if (rank == 0) {
        for (i = 0; i < n; i++) {
            scanf("%ld %ld", &v[i], &w[i]);
        }

    }

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
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>

long int max(long int x, long int y) {
    return (x > y) ? x : y;
}


long int knapSack(long int C, long int w[], long int v[], int n) {
    // broadcast parameters from rank 0 to all ranks
    MPI_Bcast(&C, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(w, n, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(v, n, MPI_LONG, 0, MPI_COMM_WORLD);


    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // allocate memory for DP memorisation table
    long int **K = malloc((2) * sizeof *K);
    assert(K);
    // fix C
    long int C_fix = ((C + 1) / nprocs + 1) * nprocs;
//    printf("C_fix: %ld\n", C_fix);

    for (long int i = 0; i < 2; i++) {
        K[i] = malloc((C_fix) * sizeof *K[i]);
        assert(K[i]);
    }

    // initialise all cells to default value -1
    for (long int row = 0; row <2; row++)
        for (long int col = 0; col < C_fix; col++)
            K[row][col] = -1;

    long int chunk_size = C_fix / nprocs;
    long int start_col = rank * chunk_size;

//    printf("chunk size: %ld\n", chunk_size);


    while (true) {
        for (int row = 0; row < n + 1; row++) {
            long int *chunk_result = malloc((chunk_size) * sizeof(*chunk_result));

            for (long int i = 0; i < chunk_size; i++) {
                chunk_result[i] = -1;
            }

            for (long int i = 0; i < chunk_size; i++) {
                // initialise all values in the first row and col to 0
                long int col = i + start_col;
//                printf("rank: %d, col %ld\n", rank, col);
                if (row == 0 || col == 0) {
                    K[row][col] = 0;
                } else if (w[row - 1] <= col) {
                    K[row][col] = max(v[row - 1] + K[row - 1][col - w[row - 1]],  // (include the item)
                                      K[row - 1][col]);  // (exclude the item)
                } else {
                    K[row][col] = K[row - 1][col];
                }
                chunk_result[i] = K[row][col];
            }
            if (rank == 3) {
                printf("chunk_result: ");
                for (int j = 0; j < chunk_size; ++j) {
                    printf("%ld ", chunk_result[j]);
                }
//                for (long int i = 0; i < C_fix; i++) {
//                    printf("%ld ", K[row][i]);
//                }
                printf("\n");
            }

            long int *row_result = malloc((C_fix) * sizeof(*row_result));
            for (long int i = 0; i < C_fix; i++) {
                row_result[i] = -1;
            }
            if (rank == 0) {
                printf("\nrow result before: ");
                for (long int i = 0; i < C_fix; i++) {
                    printf("%ld ", row_result[i]);
                }
            }

            MPI_Allgather(chunk_result, 2, MPI_LONG, row_result, 8, MPI_LONG, MPI_COMM_WORLD);
            if (rank == 0) {
                printf("\nrow result after: ");
                for (long int i = 0; i < C_fix; i++) {
                    printf("%ld ", row_result[i]);
                }
            }
//            if (rank == 0) {
//                printf("row result: ");
//            }
//            for (long int i = 0; i < C_fix; i++) {
//                // if not from it self
//                if (i / nprocs != rank) {
//                    K[row][i] = row_result[i];
//                }
//                if (rank == 0) {
//                    printf("%ld ", row_result[i]);
//                }
//            }
//            if (rank == 0) {
//                printf("\n");
//            }
//            if (rank == 0) {
//                for (long int i = 0; i < C_fix; i++) {
//                    // if not from it self
//                    printf("%ld ", K[row][i]);
//                }
//                printf("\n");
//            }
        }
        return K[n][C];
    }
}





