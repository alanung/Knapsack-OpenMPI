echo -e "[*] Compiling knapsack-dp-seq.c to bin/knapsack-dp-seq ... \c"
mpicc src/orig/knapsack-dp-seq.c -o bin/knapsack-dp-seq
echo -e "done."

echo -e "[*] Compiling knapsack-dp-mpi.c to bin/knapsack-dp-mpi ... \c"
mpicc src/knapsack-dp-mpi.c -o bin/knapsack-dp-mpi -std=c99 -lm
echo -e "done."

echo -e "[*] Compiling generator.c to bin/generator ... \c"
gcc util/genknap.c -lm -o bin/genknap
echo -e "done."
