/******************************************************************
 *                                                                *
 * File          : llnlmath.c                                     *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL    *
 * Version of    : 1 September 1994                               *
 *----------------------------------------------------------------*
 * This is the implementation file for a C math library.          *
 *                                                                *
 ******************************************************************/


#include <stdio.h>
#include <math.h>
#include "llnlmath.h"
#include "llnltyps.h"

namespace CVODE {

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)


float UnitRoundoff(void)
{
  float u;
  volatile float one_plus_u;

  u = ONE;
  one_plus_u = ONE + u;
  while (one_plus_u != ONE) {
    u /=  TWO;
    one_plus_u = ONE + u;
  }
  u *=  TWO;

  return(u);
}


float RPowerI(float base, int exponent)
{
  int i, expt;
  float prod;

  prod = ONE;
  expt = ABS(exponent);
  for(i=1; i <= expt; ++i ) prod *= base;
  if (exponent < 0) prod = ONE/prod;
  return(prod);
}


float RPowerR(float base, float exponent)
{

  if (base <= ZERO) return(ZERO);

  return((float)pow((double)base,(double)exponent));
}


float RSqrt(float x)
{
  if (x <= ZERO) return(ZERO);

  return((float) sqrt((double) x));
}

} // namespace CVODE

