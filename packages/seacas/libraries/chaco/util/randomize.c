/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include "chaco_random.h"

/* Randomly permute elements of an array. */
void      randomize(int *array, int n)
                		/* array of integer values */
            			/* number of values */
{
    double    value;		/* random value */
    int       index;		/* array index to swap with */
    int       temp;		/* holds value being swapped */
    int       i;		/* loop counter */
    double    drandom(void);

    for (i = 1; i <= n; i++) {
	value = drandom();
	index = n * value + 1;
	temp = array[i];
	array[i] = array[index];
	array[index] = temp;
    }
}

double    drandom(void)
{
  return rand_rect_port();
}

void      setrandom(long int seed)
{
    init_rand_port(seed);
}

