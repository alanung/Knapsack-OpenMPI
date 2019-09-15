//
// genknap.c
//
// Generate a random (0-1) knapsack problem on stdout.
//
// Inspired by Lachlan Andrew's earlier program, but more
// controllable.  Poorer quality pseudo-random numbers,
// but good enough for our purposes.
//
// Originally written by Les Kitchen <comp90025@po.ljk.id.au> for
// the subject COMP90025 Parallel and Multicore Computing
// Semester 2, 2019.
// Copyright rests with the University of Melbourne, as work-for-hire.
// Students may use this code for their study in COMP90025.


// Aside from using // comments, this is pretty straight C.
// You might need -std=c99 on older compilers.

#include <stdio.h>
#include <stdlib.h>

void
usage( char *prog, size_t status )
{
  fprintf(stderr, "Usage: %s <capacity:long int> <numitems:int> <maxitmval:long int> <maxitmwgt:long int> <seed:unsigned int>\n", prog);
  exit(status);
}

int debug = 0;

int main( int argc, char **argv )
{
  // Overall capacity of the knapsack.
  long int capacity;
  // Number of items in pool
  // (I'd make it size_t, but here int for consistency with previous code.)
  int numitems;
  // Start seed for srand()
  unsigned int seed;
  // Item weight range and value range.
  // For this version min is fixed, max is argument, range is computed.
  // For this version lower is fixed and upper is argument, range is computed.
  long int minitmval = 1, maxitmval, itmvalrange;
  long int minitmwgt = 1, maxitmwgt, itmwgtrange;

  // Five commandline args.  (See usage().)

  // Overall argv check.
  if ( argc != 6 ) {
    usage( argv[0], 10);
  }

  // Argument: knapsack capacity 
  {
    int s = sscanf( argv[1], "%ld", &capacity );
    if ( s != 1 || capacity < 0) {
      fprintf(stderr, "%s: Bad sscanf %d on capacity %s giving %ld\n", argv[0], s, argv[1], capacity);
      usage( argv[0], 20);
    }
  }

  // Argument: number of items in pool
  {
    int s = sscanf( argv[2], "%d", &numitems );
    if ( s != 1 || numitems < 0) {
      fprintf(stderr, "%s: Bad sscanf %d on number of items %s giving %d\n", argv[0], s, argv[2], numitems);
      usage( argv[0], 30);
    }
  }

  // Argument: maximum item value, and immediately derived range
  {
    int s = sscanf( argv[3], "%ld", &maxitmval );
    if ( s != 1 || maxitmval < minitmval) {
      fprintf(stderr, "%s: Bad sscanf %d on maxitmval %s giving %ld\n", argv[0], s, argv[3], maxitmval);
      usage( argv[0], 40);
    }
    itmvalrange = maxitmval - minitmval + 1;
  }

  // Argument: maximum item weight, and immediately derived range
  {
    int s = sscanf( argv[4], "%ld", &maxitmwgt );
    if ( s != 1 || maxitmwgt < minitmwgt) {
      fprintf(stderr, "%s: Bad sscanf %d on maxitmwgt %s giving %ld\n", argv[0], s, argv[4], maxitmwgt);
      usage( argv[0], 50);
    }
    itmwgtrange = maxitmwgt - minitmwgt + 1;
  }


  // Argument: RNG seed
  {
    int s = sscanf( argv[5], "%u", &seed );
    if ( s != 1 ) {
      fprintf(stderr, "%s: Bad sscanf %d on seed %s giving %u\n", argv[0], s, argv[5], seed);
      usage( argv[0], 60);
    }
  }

  if ( debug ) {
    fprintf(stderr, "%s arg summary: %ld %d %ld %ld %u\n",
            argv[0],
            capacity,
            numitems,
            maxitmval,
            maxitmwgt,
            seed);
  }

  // Set up RNG.  Using using older version for better compatibility.
  srand( seed );

  // Output generated knapsack problem in expected format.
  // Note, currently no checking for output errors, like
  // no-space-on-device.  We have blind faith.
  fprintf(stdout, "%ld\n%d\n", capacity, numitems);
  
  // No strong reason for running the loop backwards, except that
  // test against zero would be infinitesimally quicker.
  for ( int k = numitems;  k > 0;  --k ) {
    // Note, division by floated 1+RAND_MAX to get range [0.0,1.0)
    fprintf(stdout,
            "%ld %ld\n",
            (long int) (minitmval + itmvalrange * ((double) rand() / (1.0+RAND_MAX))),
            (long int) (minitmwgt + itmwgtrange * ((double) rand() / (1.0+RAND_MAX)))
            );
  }

  return 0;

}
