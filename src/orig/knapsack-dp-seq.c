/* Knapsack calculation based on that of */
/* https://www.tutorialspoint.com/cplusplus-program-to-solve-knapsack-problem-using-dynamic-programming */

#include <stdio.h>
#include <sys/time.h>
#include <stdint.h>
#include <mpi.h>

long int knapSack(long int C, long int w[], long int v[], int n);

uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(int argc, char *argv[]) {
    long int C;    /* capacity of backpack */
    int n;    /* number of items */
    int i;    /* loop counter */

    MPI_Init (&argc, &argv);
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        scanf ("%ld", &C);
        scanf ("%d", &n);
    }

    long int v[n], w[n];        /* value, weight */

    if (rank == 0) {
        for (i = 0; i < n; i++) {
            scanf ("%ld %ld", &v[i], &w[i]);
        }

    }

    uint64_t start = GetTimeStamp ();
    long int ks = knapSack(C, w, v, n);

    if (rank == 0) {
        printf ("knapsack occupancy %ld\n", ks);
        printf ("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
    }

    MPI_Finalize ();

    return 0;
}

/* PLACE YOUR CHANGES BELOW HERE */
#include <strings.h>
#include <assert.h>

long int max(long int x, long int y) {
   return (x > y) ? x : y;
}

long int knapSack(long int C, long int w[], long int v[], int n) {
    int i;
    long int wt;
    long int K[n+1][C+1];
    for (i = 0; i <= n; i++) {
        for (wt = 0; wt <= C; wt++) {
            if (i == 0 || wt == 0)
                K[i][wt] = 0;
            else if (w[i-1] <= wt)
                K[i][wt] = max(v[i-1] + K[i-1][wt - w[i-1]], K[i-1][wt]);
            else
                K[i][wt] = K[i-1][wt];
        }
    }
    return K[n][C];
}
