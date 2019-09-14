test() {

  ./build.sh
  Problem=$(./bin/generator $1 $2 $3)
  echo -e $Problem
  echo $Problem | ./bin/knapsack
}

if [ ! $1 ] || [ ! $2 ] || [ ! $3 ] ; then
  echo -e "Usage:\r\n  ./tesh.sh \$C \$n \$mean\r\n"
  echo -e "Example:\r\n  tesh.sh 10 7 8\r\n"

else
  test $1 $2 $3
fi
