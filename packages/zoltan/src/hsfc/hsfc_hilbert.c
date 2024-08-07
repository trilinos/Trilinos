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

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include "hsfc_hilbert_const.h"  /* contains state tables and documentation */
#include "hsfc.h"
#include "zz_const.h"

/* see maintenance note at the end of this file for information about extending
the precision of these routines. */

/* Given a 1-d coordinate in [0,1], returns it as the Hilbert key) */
double Zoltan_HSFC_InvHilbert1d (ZZ *zz, double *coord)
   {
   char *yo = "Zoltan_HSFC_InvHilbert1d";

   /* sanity check for input arguments */
   if (coord[0] < 0.0)
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Spatial Coordinates out of range.");

   return *coord;
   }



/* Given x,y coordinates in [0,1]x[0,1], returns the Hilbert key [0,1] */
double Zoltan_HSFC_InvHilbert2d (ZZ *zz, double *coord)
   {
   static const unsigned *d[]={idata2d,  idata2d  +4, idata2d  +8, idata2d  +12};
   static const unsigned *s[]={istate2d, istate2d +4, istate2d +8, istate2d +12};

   int level;
   unsigned int key[2], c[2], temp, state;
   const int MAXLEVEL = 28; /* 56 bits of significance, 28 per dimension */
   char *yo = "Zoltan_HSFC_InvHilbert2d";

   /* sanity check for input arguments */
   if ((coord[0] < 0.0) || (coord[0] > 1.0) || (coord[1] < 0.0) 
    || (coord[1] > 1.0))
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Spatial Coordinates out of range.");

   /* convert x,y coordinates to integers in range [0, IMAX] */
   c[0] = (unsigned int) (coord[0] * (double) IMAX);               /* x */
   c[1] = (unsigned int) (coord[1] * (double) IMAX);               /* y */

   /* use state tables to convert nested quadrant's coordinates level by level */
   key[0] = key[1] = 0;
   state = 0;
   for (level = 0; level < MAXLEVEL; level++) {
      temp = ((c[0] >> (30-level)) & 2)    /* extract 2 bits at current level */
           | ((c[1] >> (31-level)) & 1);

      /* treat key[] as long shift register, shift in converted coordinate */
      key[0] = (key[0] << 2) | (key[1] >> 30);
      key[1] = (key[1] << 2) | *(d[state] + temp);

      state = *(s[state] + temp);
      }

   /* convert 2 part Hilbert key to double and return */
   return ldexp ((double) key[0], -24)  +  ldexp ((double) key[1], -56);
   }



/* Given x,y,z coordinates in [0,1]x[0,1]x[0,1], returns Hilbert key in [0,1] */
double Zoltan_HSFC_InvHilbert3d (ZZ *zz, double *coord)
   {
   static const unsigned int *d[] =
     {idata3d,      idata3d +8,   idata3d +16,  idata3d +24,
      idata3d +32,  idata3d +40,  idata3d +48,  idata3d +56,
      idata3d +64,  idata3d +72,  idata3d +80,  idata3d +88,
      idata3d +96,  idata3d +104, idata3d +112, idata3d +120,
      idata3d +128, idata3d +136, idata3d +144, idata3d +152,
      idata3d +160, idata3d +168, idata3d +176, idata3d +184};

   static const unsigned int *s[] =
     {istate3d,      istate3d +8,   istate3d +16,  istate3d +24,
      istate3d +32,  istate3d +40,  istate3d +48,  istate3d +56,
      istate3d +64,  istate3d +72,  istate3d +80,  istate3d +88,
      istate3d +96,  istate3d +104, istate3d +112, istate3d +120,
      istate3d +128, istate3d +136, istate3d +144, istate3d +152,
      istate3d +160, istate3d +168, istate3d +176, istate3d +184};

   int level;
   unsigned int key[2], c[3], temp, state;
   const int MAXLEVEL = 19; /* 56 bits of significance, 18+ per dimension */
   char *yo = "Zoltan_HSFC_InvHilbert3d";

   /* sanity check for input arguments */
   if ((coord[0] < 0.0)  || (coord[0] > 1.0) || (coord[1] < 0.0)
     || (coord[1] > 1.0) || (coord[2] < 0.0) || (coord[2] > 1.0))
        ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Spatial Coordinates out of range.");

   /* convert x,y,z coordinates to integers in range [0,IMAX] */
   c[0] = (unsigned int) (coord[0] * (double) IMAX);     /* x */
   c[1] = (unsigned int) (coord[1] * (double) IMAX);     /* y */
   c[2] = (unsigned int) (coord[2] * (double) IMAX);     /* z */

   /* use state tables to convert nested quadrant's coordinates level by level */
   key[0] = key[1] = 0;
   state = 0;
   for (level = 0; level < MAXLEVEL; level++) {
      temp = ((c[0] >> (29-level)) & 4)  /* extract 3 bits at current level */
           | ((c[1] >> (30-level)) & 2)
           | ((c[2] >> (31-level)) & 1);

      /* treat key[] as long shift register, shift in converted coordinate */
      key[0] = (key[0] << 3) |  (key[1] >> 29);
      key[1] = (key[1] << 3) | *(d[state] + temp);

      state = *(s[state] + temp);
      }

   /* convert 2 part Hilbert key to double and return */
   return ldexp ((double) key[0], -25)  +  ldexp ((double) key[1], -57);
   }


   
/* Note: the following code has been tested and is fine.  It was necessary
during the testing for the new box assign algorithm.  Since it is potentially
useful, I am leaving this code intact, but ifdef'ing it out because of the
SQA coverage requirement. */
   
#ifdef RTHRTH
/* Given the Hilbert key, returns it as the coordinate in [0,1] */
void Zoltan_HSFC_Hilbert1d (ZZ *zz, double *coord, double key)
   {
   char *yo = "Zoltan_HSFC_Hilbert1d";

   /* sanity check for input argument */
   if ((key < 0.0) || (key > 1.0))
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Hilbert coordinate out of range.");

   *coord = key;
   }



/* Given the Hilbert key, returns the 2-d coordinates in [0,1]x[0,1] */
void Zoltan_HSFC_Hilbert2d (ZZ *zz, double *coord, double key)
   {
   static const unsigned *d[] = {data2d,  data2d  +4, data2d  +8, data2d  +12};
   static const unsigned *s[] = {state2d, state2d +4, state2d +8, state2d +12};
   int level, state;
   unsigned int c[2], ikey[2], temp;
   double t;
   static const int MAXLEVEL = 28;  /* only 56 significant bits, 28 per dimension */
   char *yo = "Zoltan_HSFC_Hilbert2d";

   /* sanity check for input argument */
   if ((key < 0.0) || (key > 1.0))
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Hilbert coordinate out of range.");

   ikey[1] = (unsigned int) (modf (key * (double) IMAX, &t) * (double) IMAX);
   ikey[0] = (unsigned int) t;

   /* use state tables to convert nested quadrant's coordinates level by level */
   c[0] = c[1] = 0;
   state = 0;
   for (level = 0; level < MAXLEVEL; level++) {
      temp  = *(d[state] + (ikey[0] >> 30));      /* 2 bits: xy */
      state = *(s[state] + (ikey[0] >> 30));

      c[0] |= ((temp & 2) << (30-level));         /* insert x bit */
      c[1] |= ((temp & 1) << (31-level));         /* insert y bit */

      /* treating ikey[] as shift register, shift next 2 bits in place */
      ikey[0] = (ikey[0] << 2) | (ikey[1] >> 30);
      ikey[1] =  ikey[1] << 2;
      }

   /* convert integer coordinates to doubles */
   coord[0] = (double) c[0] / (double) IMAX;      /* x in [0,1] */
   coord[1] = (double) c[1] / (double) IMAX;      /* y in [0,1] */
   }



/* Given the Hilbert key, returns the 3-d coordinates in [0,1]x[0,1]x[0,1] */
void Zoltan_HSFC_Hilbert3d (ZZ *zz, double *coord, double key)
   {
   static const unsigned int *d[] =
     {data3d,      data3d +8,   data3d +16,  data3d +24,
      data3d +32,  data3d +40,  data3d +48,  data3d +56,
      data3d +64,  data3d +72,  data3d +80,  data3d +88,
      data3d +96,  data3d +104, data3d +112, data3d +120,
      data3d +128, data3d +136, data3d +144, data3d +152,
      data3d +160, data3d +168, data3d +176, data3d +184};

   static const unsigned int *s[] =
     {state3d,      state3d +8,   state3d +16,  state3d +24,
      state3d +32,  state3d +40,  state3d +48,  state3d +56,
      state3d +64,  state3d +72,  state3d +80,  state3d +88,
      state3d +96,  state3d +104, state3d +112, state3d +120,
      state3d +128, state3d +136, state3d +144, state3d +152,
      state3d +160, state3d +168, state3d +176, state3d +184};

   int level, state;
   unsigned int c[3], ikey[2], temp;
   double t;
   static const int MAXLEVEL = 19;   /* 56 significant bits, 18+ per dimension */
   char *yo = "Zoltan_HSFC_Hilbert3d";

   /* sanity check for input argument */
   if ((key < 0.0) || (key > 1.0))
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Hilbert coordinate out of range.");

   ikey[1] = (unsigned int) (modf (key * (double) IMAX, &t) * (double) IMAX);
   ikey[0] = (unsigned int) t;

   /* use state tables to convert nested quadrant's coordinates level by level */
   c[0] = c[1] = c[2] = 0;
   state = 0;
   for (level = 0; level < MAXLEVEL; level++) {
      temp  = *(d[state] + (ikey[0] >> 29));     /* get next 3 bits: xyz */
      state = *(s[state] + (ikey[0] >> 29));

      c[0] |= ((temp & 4) << (29-level));        /* insert x bit */
      c[1] |= ((temp & 2) << (30-level));        /* insert y bit */
      c[2] |= ((temp & 1) << (31-level));        /* insert z bit */

      /* treating ikey[] as shift register, shift next 3 bits in place */
      ikey[0] = (ikey[0] << 3) | (ikey[1] >> 29);
      ikey[1] =  ikey[1] << 3;
      }

   /* convert coordinates to doubles */
   coord[0] = (double) c[0] / (double) IMAX;     /* x in [0,1] */
   coord[1] = (double) c[1] / (double) IMAX;     /* y in [0,1] */
   coord[2] = (double) c[2] / (double) IMAX;     /* z in [0,1] */
   }
#endif

/* MAINTENANCE NOTE:  Per the design review 04/15/03, this section discusses
how to increase the precision of the routines in this file. The table driven
logic can be extended to arbitrary precision.  Currently, these routines only
take or return the precision of a double.

First, consider how to extend the Zoltan_HSFC_InvHilbertxx() routines to return
the full precision for the 2 or 3 dimensional input:
The statement
      c[0] = (unsigned int) (coord[0] * (double) IMAX); (etc. for each dimension)
returns about 32 bits of usable information per dimension.
Change the declaration:
   const int MAXLEVEL = 28;
to
   const int MAXLEVEL = 32;
Then the key array contains 64 bits in 2d or 96 bits in 3d. Change its
declaration to either 2 or 3 ints as appropriate.
The easiest thing is to modify the routine to return the key array itself to
preserve the maximum precision.  This requires changing the function definition:
   double Zoltan_HSFC_InvHilbert2d (ZZ *zz, double *coord)
to
   int *Zoltan_HSFC_InvHilbert2d (ZZ *zz, double *coord).
and changing the return statement:
   return ldexp ((double) key[0], -24)  +  ldexp ((double) key[1], -56);
to
   return key;


Second, consider how to extend the Zoltan_HSFC_Hilbertxx() routines. If the key
array is used as the return argument above, then these Hilbertxx routines should
be modified to accept the key array as the Hilbert coordinate rather than the
current double.  This requires changing the function definition:
   void Zoltan_HSFC_Hilbert2d (ZZ *zz, double *coord, double key)
to
   void Zoltan_HSFC_Hilbert2d (ZZ *zz, double *coord, int* ikey)
Change the declaration
   static const int MAXLEVEL = 28; (2d or 19 in 3d)
to
   static const int MAXLEVEL = 32;
Eliminate the lines that convert the orginal double arguement of key to ikey.
*/   


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

