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

typedef struct {
    int requester;
    long int task;
} TaskRequest;

typedef struct {
    int task;
} TaskRefusal;

typedef struct {
    int row;
    long int col;
    long int val;
    int handling_by;
} CellResult;

typedef struct {
    int row;
    long int col;
    int handling_by;
} TaskStatus;

MPI_Datatype structs_register(int type);

void print_info(int rank, long int C, long int *K, long int *v, long int *w, int n);
void knapsack_master(long int C, int n, long int **K);
void knapsack_worker(long int C, long int *w, long int *v, int n, long int **K);

#define TAG_TASK_REQUEST 0
#define TAG_TASK_REFUSAL 1
#define TAG_RESULT 2
#define TAG_TASK_STATUS 2

#define MPI_TASK_STATUS 0
#define MPI_CELL_RESULT 1
#define MPI_TASK_REFUSAL 2
#define MPI_TASK_REQUEST 3

long int max(long int x, long int y) {
    return (x > y) ? x : y;
}

long int knapSack(long int C, long int w[], long int v[], int n) {
    // broadcast parameters from rank 0 to all ranks
    MPI_Bcast(&C, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(w, n, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(v, n, MPI_LONG, 0, MPI_COMM_WORLD);

    // register MPI struct types
    MPI_Datatype mpi_task_status = structs_register(MPI_TASK_STATUS);
    MPI_Datatype mpi_cell_result = structs_register(MPI_CELL_RESULT);
    MPI_Datatype mpi_task_request = structs_register(MPI_TASK_REQUEST);
    MPI_Datatype mpi_task_refusal = structs_register(MPI_TASK_REFUSAL);

    // declare structs to hold different types of messages
    TaskStatus task_status;
    CellResult cell_result;
    TaskRefusal task_refusal;
    TaskRequest task_request;

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // allocate memory for DP memoisation table
    long int **K = malloc((n + 1) * sizeof *K);
    assert(K);
    for (long int i = 0; i < n + 1; i++) {
        K[i] = malloc((C + 1) * sizeof(*K[i]));
        assert(K[i]);

    }

    if (rank == 0)
        knapsack_master(C, n, K);
    else
        knapsack_worker(C, w, v, n, K);

    // BELOW NOT NEEDED

    // initialise all cells to default value INT_MIN
    for (long int row = 0; row < n + 1; row++)
        for (long int col = 0; col < C + 1; col++)
            K[row][col] = -1;


    long int col = rank;

    if (rank == 0) {

    }

    while (true) {



//        if (rank == 1) {
//            printf("[*] %d is handling %ld, start from %d\n", rank, col, row);
//            print_info(rank, C, K, v, w, n);
//
//        }
//        print_info(rank, C, task_status_arr, K, v, w, n);
        for (int row = 0; row < n + 1; row++) {

            // skip this cell if it has already been calculated
            if (K[row * (C + 1) + col] != -1) continue;

            // initialise all values in the first row and col to 0
            if (row == 0 || col == 0) {
                K[row * (C + 1) + col] = 0;

                // if the considered item (row - 1) can fit in the current capacity (col)...
            } else if (w[row - 1] <= col) {
                // and the cell dependencies are available, calculate the value of this cell
                if (K[(row - 1) * (C + 1) + col - w[row - 1]] != -1 && K[(row - 1) * (C + 1) + col] != -1)
                    K[row * (C + 1) + col] = max(
                            v[row - 1] + K[(row - 1) * (C + 1) + col - w[row - 1]],  // (include the item)
                            K[(row - 1) * (C + 1) + col]);  // (exclude the item)
                    // otherwise, if the cell dependencies are not available, stop calculating this column
                else
                    break;

                // otherwise, if it cannot fit...
            } else {
                // take the value from the previous row if available
                if (K[(row - 1) * (C + 1) + col] != -1)
                    K[row * (C + 1) + col] = K[(row - 1) * (C + 1) + col];
                    // otherwise, stop calculating this column
                else
                    break;
            }

            // broadcast cell result
//            broadcast_cell_result(rank, nprocs, col, row, K, mpi_cell_result);
            CellResult tmp_cell_result[nprocs];
            cell_result.val = K[row * (C + 1) + col];
            cell_result.col = col;
            cell_result.row = row;
            cell_result.handling_by = rank;


            MPI_Allgather(&cell_result, 1, mpi_cell_result, tmp_cell_result,
                          1, mpi_cell_result, MPI_COMM_WORLD);

            for (int i = 0; i < nprocs; ++i) {
                K[tmp_cell_result[i].row * (C + 1) + tmp_cell_result[i].col] = tmp_cell_result[i].val;
            }
        }

        if (K[n * (C + 1) + C] != -1) {
            return K[n * (C + 1) + C];
        }

        col += nprocs;
        col %= C + 1;
    }
}

void knapsack_master(long int C, int n, long int **K) {
    MPI_Datatype mpi_cell_result = structs_register(MPI_CELL_RESULT);
    CellResult result;
    for (int i = 0; i < (C + 1) * (n + 1); i++) {
        MPI_Recv(&result, 1, mpi_cell_result, MPI_ANY_SOURCE, TAG_RESULT,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        K[result.row][result.col] = result.val;
    }
}

void knapsack_worker(long int C, long int *w, long int *v, int n, long int **K) {
    MPI_Datatype mpi_cell_result = structs_register(MPI_CELL_RESULT);
    CellResult result;
    for () {
        // all cells in row 0 and col 0 have value 0
        if (row == 0 || col == 0) {
            K[row * (C + 1) + col] = 0;

        // if the considered item (row - 1) can fit in the current capacity (col), compare including and excluding
        } else if (w[row - 1] <= col) {
            // TODO: GET VALUES FROM MASTER
            K[row * (C + 1) + col] = max(K[(row - 1) * (C + 1) + col],  // (exclude the item)
                                         v[row - 1] + K[(row - 1) * (C + 1) + col - w[row - 1]]);  // (include the item)

        // otherwise, if it cannot fit, take the value from the previous row
        } else {
            // TODO: GET VALUE FROM MASTER
            K[row * (C + 1) + col] = K[(row - 1) * (C + 1) + col];
        }
        cell_result.val = K[row][col];
        cell_result.col = col;
        cell_result.row = row;
        cell_result.handling_by = rank;
        MPI_Send(&result, 1, mpi_cell_result, 1, TAG_RESULT, MPI_COMM_WORLD);
    }
}

MPI_Datatype structs_register(int type) {
    if (type == MPI_TASK_STATUS) {
//        typedef struct {
//            int row;
//            int handling_by;
//            int current_pos;
//        } TaskStatus;
        TaskStatus task_status;
        MPI_Datatype mpi_task_status;
        MPI_Datatype old_types[3] = {MPI_INT, MPI_LONG, MPI_INT};
        MPI_Aint indices[3];
        int block_length[3] = {1, 1, 1};
        MPI_Get_address(&task_status, &indices[0]);
        MPI_Get_address(&task_status.col, &indices[1]);
        MPI_Get_address(&task_status.handling_by, &indices[2]);
        indices[2] -= indices[0];
        indices[1] -= indices[0];
        indices[0] = 0;
        MPI_Type_create_struct(3, block_length, indices, old_types, &mpi_task_status);
        MPI_Type_commit(&mpi_task_status);
        return mpi_task_status;
    } else if (type == MPI_CELL_RESULT) {
//        typedef struct {
//            int row;
//            long int col;
//            long int val;
//            int handling_by;
//        } CellResult;
        CellResult cell_result;
        MPI_Datatype mpi_cell_result;
        MPI_Datatype old_types[4] = {MPI_INT, MPI_LONG, MPI_LONG, MPI_INT};
        MPI_Aint indices[4];
        int block_length[4] = {1, 1, 1, 1};
        MPI_Get_address(&cell_result, &indices[0]);
        MPI_Get_address(&cell_result.col, &indices[1]);
        MPI_Get_address(&cell_result.val, &indices[2]);
        MPI_Get_address(&cell_result.handling_by, &indices[3]);
        indices[1] -= indices[0];
        indices[2] -= indices[0];
        indices[3] -= indices[0];
        indices[0] = 0;
        MPI_Type_create_struct(4, block_length, indices, old_types, &mpi_cell_result);
        MPI_Type_commit(&mpi_cell_result);
        return mpi_cell_result;
    } else if (type == MPI_TASK_REQUEST) {
        TaskRequest task_request;
        MPI_Datatype mpi_task_request;
        MPI_Datatype old_types[2] = {MPI_INT, MPI_LONG};
        MPI_Aint indices[2];
        int block_length[2] = {1, 1};
        MPI_Get_address(&task_request, &indices[0]);
        MPI_Get_address(&task_request.task, &indices[1]);
        indices[1] -= indices[0];
        indices[0] = 0;
        MPI_Type_create_struct(2, block_length, indices, old_types, &mpi_task_request);
        MPI_Type_commit(&mpi_task_request);
        return mpi_task_request;
    } else if (type == MPI_TASK_REFUSAL) {
        TaskRefusal task_refusal;
        MPI_Datatype mpi_task_refusal;
        MPI_Datatype old_types[1] = {MPI_INT};
        MPI_Aint indices[1];
        int block_length[1] = {1};
        MPI_Get_address(&task_refusal, &indices[0]);
        indices[0] = 0;
        MPI_Type_create_struct(1, block_length, indices, old_types, &mpi_task_refusal);
        MPI_Type_commit(&mpi_task_refusal);
        return mpi_task_refusal;
    }
}

void print_info(int rank, long int C, long int *K, long int *v, long int *w, int n) {
    printf("==========\n");


    printf("\n----------\n");

    printf("\n");

    for (int i = 0; i < n + 1; ++i) {
        if (i != 0) {
            printf("%3ld %3ld: ", v[i - 1], w[i - 1]);

        } else {
            printf("%3d %3d: ", 0, 0);

        }
        for (long int j = 0; j < C + 1; ++j) {
            printf("%4ld ", K[i * (C + 1) + j]);
        }
        printf("\n");
    }

}





