#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int
main (int argc, char *argv[])
{
    if (argc < 4) {
        fprintf (stderr, "usage: %s C n over-subscription\n", argv[0]);
        exit (1);
    }
    
    long int C = atol (argv[1]);        /* capacity */
    int n = atoi (argv[2]);             /* number of items */
    float over_sub = atof (argv[3]);    /* (expected sum of weights)/capacity */
    
    long int mean = (int) (over_sub * C / n);

    printf ("%ld\n%d\n", C, n);
    int i;
    for (i = 0; i < n; i++) {
        printf ("%ld %ld\n",
                1 + (long int)(floor(-log (rand() / RAND_MAX) * mean)),
                1 + rand() % mean);
    }
}
