/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif



/****************************************************************************/
/*  Random Number generator due to Knuth found in Numerical Recipes in C
   (2nd edition) by Press, Vetterling, Teukolsky, Flannery (Page 284.)
   Currently, for use in the hypergraph library only. Needed because different
   random number implimentations on different machines produced different
   answers! Hence this is a single portable, fast, algorithm with adequate
   random number generation. NOTE: this assumes 32 bit ints, previously
   these variables were unsigned long (as was the return value) which
   gave problems on stratus (which assumed 64 bit longs.) */

static unsigned int idum = 123456789U;

unsigned Zoltan_HG_Rand (void) {
   return idum = (1664525U * idum) + 1013904223U;
   }



void Zoltan_HG_Srand (unsigned int seed) {
   idum = seed;
   }


/* Randomly permute an array of ints. */
void Zoltan_HG_Rand_Perm_Int (int *data, int n)
{
int i, number, temp;

  for (i = n; i > 0; i--) {
     number       = Zoltan_HG_Rand() % i;
     temp         = data[number];
     data[number] = data[i-1];
     data[i-1]    = temp;
     }
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
