echo -e "[*] Compiling knapsack.c to bin/knapsack ... \c"
mpicc src/orig/knapsack.c -o bin/knapsack
echo -e "done."

echo -e "[*] Compiling generator.c to bin/generator ... \c"
gcc util/generator.c -lm -o bin/generator
echo -e "done."