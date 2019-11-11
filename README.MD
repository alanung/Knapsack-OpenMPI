# Knapsack-OpenMPI
## Requirements
+ Open MPI
+ GCC
+ MPICC
## Compile
```shell script
cd $PROJECT_PATH
mkdir bin -p
mpicc src/orig/knapsack-dp-seq.c -o bin/knapsack-dp-seq
mpicc src/knapsack-dp-mpi.c -o bin/knapsack-dp-mpi -std=c99 -lm
gcc util/genknap.c -lm -o bin/genknap
```
## Usage
Generate graph,
```shell script
./bin/genknap 1000 1000 100 20 17 > simple.in
```
Use sequential algorithm,
```shell script
user@machine $ ./bin/knapsack-dp-seq < simple.in            
knapsack occupancy 16614
Time: 8005 us

```
Use parallel algorithm,
```shell script
user@machine $ ./bin/knapsack-dp-mpi < simple.in            
knapsack occupancy 16614
Time: 5632 us

```

Generate large graph,
```shell script
./bin/genknap 1000 20000 100 20 17 > simple.in
```
Use sequential algorithm,
```shell script
user@machine $ ./bin/knapsack-dp-seq < simple.in            
[1]    4264 segmentation fault (core dumped)  ./bin/knapsack-dp-seq < simple.in
```
Use parallel algorithm,
```shell script
user@machine $ ./bin/knapsack-dp-mpi < simple.in             
knapsack occupancy 60634
Time: 116925 us
```