/* Listing 2
 *    rand_por[t].c
 *  see
 *    L'Ecuyer - Comm. of the ACM, Oct. 1990, vol. 33.
 *    Numerical Recipes in C, 2nd edition, pp. 278-86 (NO CODE IS FROM HERE)
 *    L'Ecuyer and Cote, ACM Transactions on Mathematical
 *       Software, March 1991
 *    Russian peasant algorithm -- Knuth, vol. II, pp. 442-43
 *  Copyright (c) 2005-2017 National Technology & Engineering Solutions
 *  of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 *  NTESS, the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *      * Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *
 *      * Redistributions in binary form must reproduce the above
 *        copyright notice, this list of conditions and the following
 *        disclaimer in the documentation and/or other materials provided
 *        with the distribution.
 *
 *      * Neither the name of NTESS nor the names of its
 *        contributors may be used to endorse or promote products derived
 *        from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>

#include "chaco_random.h"

#define MOD 2147483647L /* modulus for generator  */
#define MULT 41358L     /* multiplier             */
                        /* modulus = mult*quotient  + remainder  */
#define Q 51924L        /* int(modulus / multiplier) */
#define R 10855L        /* remainder                 */
#define MAX_VALUE (MOD - 1)

#define EXP_VAL 1285562981L /* value for 10,000th draw */

#define IMPOSSIBLE_RAND (-1)
#define STARTUP_RANDS                                                                              \
  16 /* throw away this number of                                                                  \
   initial random numbers */

static long rand_num = IMPOSSIBLE_RAND;

/* initialize random number generator with seed */
long init_rand_port(long seed)
{
  extern long rand_num;
  int         i;

  if (seed < 1 || seed > MAX_VALUE) { /* if seed out of range */
    seed = get_init_rand_port();      /* get seed */
  }

  rand_num = seed;
  for (i = 0; i < STARTUP_RANDS; i++) { /* and throw away */
    rand_num = genr_rand_port(rand_num);
  } /* some initial
                                        ones */

  return seed;
}

/* get a long initial seed for gererator
  assumes that rand returns a short integer */
long get_init_rand_port(void)
{
  long seed;

  srand((unsigned int)time(NULL)); /* initialize system generator */
  do {
    seed = ((long)rand()) * rand();
    seed += ((long)rand()) * rand();
  } while (seed > MAX_VALUE);

  assert(seed > 0);

  return seed;
}

/* generate the next value in sequence from generator
    uses approximate factoring
    residue = (a * x) mod modulus
           = a*x - [(a*x)/modulus]*modulus
    where
        [(a*x)/modulus] = integer less than or equal to (a*x)/modulus
    approximate factoring avoids overflow associated with a*x and
        uses equivalence of above with
    residue = a * (x - q * k) - r* k + (k-k1) * modulus
    where
        modulus = a * q + r
        q = [modulus/a]
        k = [x/q]  (= [ax/aq])
        k1 = [a*x/modulus]
    assumes
        a, m > 0
        0 < init_rand < modulus
        a * a <= modulus
        [a*x/a*q]-[a*x/modulus] <= 1
            (for only one addition of modulus below) */
long genr_rand_port(long init_rand)
{
  long k, residue;

  k       = init_rand / Q;
  residue = MULT * (init_rand - Q * k) - R * k;
  if (residue < 0) {
    residue += MOD;
  }

  assert(residue >= 1 && residue <= MAX_VALUE);
  return residue;
}

/* get a random number */
long rand_port(void)
{
  extern long rand_num;
  if (rand_num == IMPOSSIBLE_RAND) {
    /* if not initialized, do it now */
    rand_num = 1;
    init_rand_port(rand_num);
  }

  rand_num = genr_rand_port(rand_num);

  return rand_num;
}

/* generates a value on (0,1) with mean of .5
  range of values is [1/(MAX_VALUE+1), MAX_VALUE/(MAX_VALUE+1)]
  to get [0,1], use (double)(rand_port()-1)/(double)(MAX_VALUE-1) */
double rand_rect_port(void) { return (double)rand_port() / (double)(MAX_VALUE + 1); }

  /* skip ahead in recursion
    residue = (a^skip * init) mod modulus
    Use Russian peasant algorithm  */

  /* calculate residue = (a * x) mod modulus for arbitrary a and x
    without overflow assume 0 < a < modulus and 0 < x < modulus
    use Russian peasant algorithm followed by approximate factoring */

#if defined(TESTING)
/* Test the generator */
#include <stdio.h>
int main(void)
{
  long seed;
  int  i;
  seed = init_rand_port(1);
  printf("Seed for random number generator is %ld\n", seed);
  i = STARTUP_RANDS; /* threw away STARTUP_RANDS */
  do {
    rand_port();
    i++;
  } while (i < 9999);

  printf("On draw 10000, random number should be %ld\n", EXP_VAL);
  printf("On draw %d, random number is %ld\n", i + 1, rand_port());
}
#endif /* TESTING */
/* End of File */
