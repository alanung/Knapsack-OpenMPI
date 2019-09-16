test() {
  ./build.sh > /dev/null
  ./bin/genknap $1 $2 $3 $4 101010 >in_sadfasdf
  #  cat in_sadfasdf
  echo "The knapsack-dp-seq:"
  cat in_sadfasdf | mpirun ./bin/knapsack-dp-seq

  echo "The knapsack-dp-mpi:"
  cat in_sadfasdf | mpirun ./bin/knapsack-dp-mpi -np 4
  rm in_sadfasdf
}

capacity=0

for ((i = $1; i < $2 + 1; i += $3)); do
  echo "    Testing $i ... "
  test $i $4 $5 $6
  echo -e "    done.\r\n\r\n"
done
