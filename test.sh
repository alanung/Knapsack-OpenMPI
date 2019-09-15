test() {
  ./build.sh
  ./bin/genknap $1 $2 $3 $4 101010 > in_sadfasdf
  cat in_sadfasdf
  echo "The knapsack-dp-seq:"
  cat in_sadfasdf | mpirun ./bin/knapsack-dp-seq

  echo "The knapsack-dp-mpi:"
  cat in_sadfasdf | mpirun ./bin/knapsack-dp-mpi -np 4
  rm in_sadfasdf
}

if [ ! $1 ] || [ ! $2 ] || [ ! $3 ]; then
  echo -e "Usage:\r\n  ./tesh.sh \$C \$n \$mean\r\n"
  echo -e "Example:\r\n  tesh.sh 10 7 8\r\n"

else
  test $1 $2 $3 $4
fi
