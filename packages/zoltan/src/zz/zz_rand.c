// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_rand.h"
#include "zz_const.h"


/****************************************************************************/
 /* Linear congruential number generator, with right shift to account for
 * the lack of randomness (even/odd/even/odd pattern) in low order bit.
 *
 * Needed because different random number implementations on different 
 * machines produced different answers!  This generator provides a portable, 
 * fast, algorithm with adequate random number generation. 
 *
 * Number generated should be the same on all platforms that have
 * 32 bit integers.
 */

static unsigned int zidum = ZOLTAN_RAND_INIT;

unsigned int Zoltan_Seed()
{
/* Function that returns the current value of the Zoltan seed. */
  return zidum;
}


unsigned Zoltan_Rand(unsigned int *myidum) {
/* 
 * If myidum is non-NULL, use *myidum as the generator value.  This feature
 * allows synchronization of the RNG across processors.
 * If myidum is NULL, use zidum.
 */
unsigned int *idum;

  if (myidum) 
    idum = myidum;
  else
    idum = &zidum;
  *idum = ((1664525U * *idum) + 1013904223U) >> 1;
  return (*idum);
}



void Zoltan_Srand (unsigned int seed, unsigned int *myidum) {
/* 
 * If myidum is non-NULL, set *myidum to the seed.  
 * If myidum is NULL, set zidum.
 */
unsigned int *idum;

  if (myidum) 
    idum = myidum;
  else
    idum = &zidum;
  *idum = seed;
}



void Zoltan_Srand_Sync(
  unsigned int seed, 
  unsigned int *myidum,
  MPI_Comm comm
)
{
/* Synchronize the random number seed across processor within a communicator.
 * Proc 0's seed is the broadcast seed used by all procs. 
 * If myidum is non-NULL, set *myidum to the seed.  
 * If myidum is NULL, set zidum.
 */

unsigned int *idum;

  if (myidum) 
    idum = myidum;
  else
    idum = &zidum;
  *idum = seed;
  MPI_Bcast(idum, 1, MPI_UNSIGNED, 0, comm);
}


unsigned int Zoltan_Rand_InRange (unsigned int *myidum, unsigned int n)
{
  double denom = ZOLTAN_RAND_MAX + 1.0;

  return (int) ((double) n * (double) Zoltan_Rand(myidum) / denom);
}
    
/* Randomly permute an array of ints. */
void Zoltan_Rand_Perm_Int (int *data, int n, unsigned int *myidum)
{
int i, number, temp;
double denom = ZOLTAN_RAND_MAX + 1.0;
/* Scaling of random number to appropriate range is done as recommended
 * in Numerical Recipes in C.
 */

  for (i = n; i > 0; i--) {
    number       = (int) ((double) i * (double) Zoltan_Rand(myidum) / denom);
    temp         = data[number];
    data[number] = data[i-1];
    data[i-1]    = temp;
  }
}

/* Randomly permute an array of ZOLTAN_GNO_TYPEs. TODO64 - still a good permutation for 8 byte ints? */

void Zoltan_Rand_Perm_Gno (ZOLTAN_GNO_TYPE *data, ZOLTAN_GNO_TYPE n, unsigned int *myidum)
{
ZOLTAN_GNO_TYPE i, number, temp;
double denom = ZOLTAN_RAND_MAX + 1.0;
/* Scaling of random number to appropriate range is done as recommended
 * in Numerical Recipes in C.
 */

  for (i = n; i > 0; i--) {
    number       = (ZOLTAN_GNO_TYPE) ((double) i * (double) Zoltan_Rand(myidum) / denom);
    temp         = data[number];
    data[number] = data[i-1];
    data[i-1]    = temp;
  }
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
