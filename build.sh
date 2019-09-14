echo -e "[*] Compiling knapsack-dp-seq.c to bin/knapsack ... \c"
mpicc src/orig/knapsack-dp-seq.c -o bin/knapsack-dp-seq
echo -e "done."

echo -e "[*] Compiling knapsack-dp-mpi.c to bin/knapsack ... \c"
mpicc src/knapsack-dp-mpi.c -o bin/knapsack-dp-mpi
echo -e "done."

echo -e "[*] Compiling generator.c to bin/generator ... \c"
gcc util/generator.c -lm -o bin/generator
echo -e "done."