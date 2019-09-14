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
#include <assert.h>
#include <stdbool.h>
#include <strings.h>

#define TAG_TASK_REQUEST 0
#define TAG_TASK_REFUSAL 1

typedef struct {
    int requester;
    int task;
} TaskRequest;

long int max(long int x, long int y) {
   return (x > y) ? x : y;
}

long int knapSack(long int C, long int w[], long int v[], int n) {
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    printf("rank: %d\n", rank);

    MPI_Receive(...);


    // int i;
    // long int wt;
    // long int K[n+1][C+1];
    // for (i = 0; i <= n; i++) {
    //     for (wt = 0; wt <= C; wt++) {
    //         if (i == 0 || wt == 0)
    //             K[i][wt] = 0;
    //         else if (w[i-1] <= wt)
    //             K[i][wt] = max(v[i-1] + K[i-1][wt - w[i-1]], K[i-1][wt]);
    //         else
    //             K[i][wt] = K[i-1][wt];
    //     }
    // }
    // return K[n][C];
    return 0;
}

/**
 * Request to be allocated a task (column of the DP memoisation table).
 */
void task_request(int requester, int task) {
    TaskRequest request = {requester, task};

    // create MPI struct type corresponding to TaskRequest
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Datatype MPI_TaskRequest;
    MPI_Aint offsets[2];
    MPI_Type_create_struct(2, blocklengths, offsets, types, &MPI_TaskRequest);
    MPI_Type_commit(&MPI_TaskRequest);

    // send the message to the next rank
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Bsend(&request, 1, MPI_TaskRequest, (requester + 1) % nprocs,
              TAG_TASK_REQUEST, MPI_COMM_WORLD);
}

/**
 * Notify a rank that the column it has requested has already been requested.
 */
void task_refuse(TaskRequest request) {
    MPI_Bsend(&request.task, 1, MPI_INT, request.requester, TAG_TASK_REFUSAL,
              MPI_COMM_WORLD);
}

void process_messages() {
    while (true) {

    }
}