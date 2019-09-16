test() {
  ./build.sh >/dev/null
  echo "The knapsack-dp-mpi:"
  ./bin/genknap $1 $2 $3 $4 13 >tmp.in
  mpirun ./bin/knapsack-dp-mpi <tmp.bin
  rm tmp.in
}

for ((i = $1; i < $2 + 1; i += $3)); do
  echo "    Testing $i ... "
  test $i $4 $5 $6
  echo -e "    done.\r\n\r\n"
done
