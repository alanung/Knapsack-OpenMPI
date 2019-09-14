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
#include <assert.h>
#include <stdbool.h>
#include <strings.h>
#include <unistd.h>

#define TAG_TASK_REQUEST 0
#define TAG_TASK_REFUSAL 1

typedef struct {
    int requester;
    int task;
} TaskRequest;

typedef struct {
    int task;
} TaskRefusal;

typedef struct {
    int row;
    int col;
    int val;
} CellResult;

typedef struct {
    int row;
    int col;
} TaskStatus;

#define MPI_TASK_STATUS 0
#define MPI_CELL_RESULT 1
#define MPI_TASK_REFUSAL 2
#define MPI_TASK_REQUEST 3

MPI_Datatype structs_register(int type) {
    if (type == MPI_TASK_STATUS) {
        TaskStatus task_status;
        MPI_Datatype mpi_task_status;
        MPI_Datatype old_types[2] = {MPI_INT, MPI_INT};
        MPI_Aint indices[2];
        int block_length[2] = {1, 1};
        MPI_Address(&task_status, &indices[0]);
        MPI_Address(&task_status.col, &indices[1]);
        indices[1] -= indices[0];
        indices[0] = 0;
        MPI_Type_struct(2, block_length, indices, old_types, &mpi_task_status);
        MPI_Type_commit(&mpi_task_status);
        return mpi_task_status;
    } else if (type == MPI_CELL_RESULT) {
        CellResult cell_result;
        MPI_Datatype mpi_cell_result;
        MPI_Datatype old_types[3] = {MPI_INT, MPI_INT, MPI_INT};
        MPI_Aint indices[3];
        int block_length[3] = {1, 1, 1};
        MPI_Address(&cell_result, &indices[0]);
        MPI_Address(&cell_result.col, &indices[1]);
        MPI_Address(&cell_result.val, &indices[2]);
        indices[1] -= indices[0];
        indices[2] -= indices[0];
        indices[0] = 0;
        MPI_Type_struct(3, block_length, indices, old_types, &mpi_cell_result);
        MPI_Type_commit(&mpi_cell_result);
        return mpi_cell_result;
    } else if (type == MPI_TASK_REQUEST) {
        TaskRequest task_request;
        MPI_Datatype mpi_task_request;
        MPI_Datatype old_types[2] = {MPI_INT, MPI_INT};
        MPI_Aint indices[2];
        int block_length[2] = {1, 1};
        MPI_Address(&task_request, &indices[0]);
        MPI_Address(&task_request.task, &indices[1]);
        indices[1] -= indices[0];
        indices[0] = 0;
        MPI_Type_struct(2, block_length, indices, old_types, &mpi_task_request);
        MPI_Type_commit(&mpi_task_request);
        return mpi_task_request;
    } else if (type == MPI_TASK_REFUSAL) {
        TaskRefusal task_refusal;
        MPI_Datatype mpi_task_refusal;
        MPI_Datatype old_types[1] = {MPI_INT};
        MPI_Aint indices[1];
        int block_length[1] = {1};
        MPI_Address(&task_refusal, &indices[0]);
        indices[0] = 0;
        MPI_Type_struct(1, block_length, indices, old_types, &mpi_task_refusal);
        MPI_Type_commit(&mpi_task_refusal);
        return mpi_task_refusal;
    }
}

void task_request(int requester, int task, MPI_Datatype mpi_task_status);

void task_refuse(TaskRequest request);

void receiveAllMsg(int rank, int nprocs, MPI_Datatype mpi_task_status);

long int max(long int x, long int y) {
    return (x > y) ? x : y;
}

long int knapSack(long int C, long int w[], long int v[], int n) {
    int rank, nprocs;
    //register mpi struct types
    MPI_Datatype mpi_task_status = structs_register(MPI_TASK_STATUS);
    MPI_Datatype mpi_cell_result = structs_register(MPI_CELL_RESULT);
    MPI_Datatype mpi_task_request = structs_register(MPI_TASK_REQUEST);
    MPI_Datatype mpi_task_refusal = structs_register(MPI_TASK_REFUSAL);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    task_request(rank, 888, mpi_task_status);

    receiveAllMsg(rank, nprocs, mpi_task_status);


    printf("[*] P-%d started.\n", rank);
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
void task_request(int requester, int task, MPI_Datatype mpi_task_status) {
    TaskStatus task_status;
    for (int i = 0; i < 5; i++) {
        task_status.row = requester * 10 + i;
        task_status.col = requester * 100 + i;
        // send the message to the next rank
        int nprocs;
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Bsend(&task_status, 1, mpi_task_status, (requester + 1) % nprocs,
                  TAG_TASK_REQUEST, MPI_COMM_WORLD);
    }

}

/**
 * Notify a rank that the column it has requested has already been requested.
 */
void task_refuse(TaskRequest request) {
    TaskRefusal refusal;
    MPI_Bsend(&request.task, 1, MPI_INT, request.requester, TAG_TASK_REFUSAL,
              MPI_COMM_WORLD);
}

/**
 * Listen for a request message.
 */
//void listen_for_requests(int rank, int nprocs, MPI_Datatype mpi_task_status) {
void receiveAllMsg(int rank, int nprocs, MPI_Datatype mpi_task_status) {
    TaskStatus task_status;
//    sleep(3);
    MPI_Request request;
    MPI_Status status;

    int flag;
    do {
        MPI_Irecv(&task_status, 1, mpi_task_status, MPI_ANY_SOURCE, TAG_TASK_REQUEST,
                  MPI_COMM_WORLD, &request);
        MPI_Test(&request, &flag, &status);
        printf("[%d received] task_status.row: %d,task_status.col: %d\n", rank, task_status.row, task_status.col);
        if (!flag) {
            break;
        }
    } while (flag);


}


/**
 * Listen for a refusal message.
 */
void listen_for_refusals() {
    TaskRefusal refusal;
    while (true) {
        MPI_Recv(&refusal, 1, MPI_INT, MPI_ANY_SOURCE, TAG_TASK_REQUEST,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("[refusal received] task: %d\n", refusal.task);
    }
}
