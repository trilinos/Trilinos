/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * This software is distributed under the GNU Lesser General Public License. *
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

#include <stdio.h>
#include "comm.h"

/* Given a list of processors & message data, sort them by proc ID */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int       Zoltan_Comm_Sort_Ints(
int      *vals_sort,		/* values to be sorted */
int      *vals_other,		/* other array to be reordered w/ sort */
int       nvals)		/* length of these two arrays */
{
/* This routine will ensure that the ordering */
/* produced by the invert_map routines is deterministic.  This should */
/* make bugs more reproducible.  This is accomplished by sorting */
/* the message lists by processor ID. */

    int       temp;		/* swapping value */
    int       lo, hi;		/* counters from bottom and top of array */
    int       pivot;		/* value to partition with */

    if (nvals <= 1) return(ZOLTAN_OK);

    /* Partition */
    lo = nvals/2;
    if (lo == nvals - 1) --lo;
    pivot = vals_sort[lo];

    lo = -1;
    hi = nvals;

    while (lo < hi) {
	do {
	    hi--;
        } while (vals_sort[hi] > pivot);

	do {
	    lo++;
	} while (vals_sort[lo] < pivot);

        if (lo < hi) {	/* Swap low and high items */
	    temp = vals_sort[lo];
	    vals_sort[lo] = vals_sort[hi];
	    vals_sort[hi] = temp;

	    temp = vals_other[lo];
	    vals_other[lo] = vals_other[hi];
	    vals_other[hi] = temp;
	}
    } 

    /* Recurse */
    if (hi + 1 > 1) Zoltan_Comm_Sort_Ints(vals_sort, vals_other, hi + 1);
    if (nvals - hi - 1 > 1)
	Zoltan_Comm_Sort_Ints(&vals_sort[hi + 1], &vals_other[hi + 1], nvals - hi - 1);


    return(ZOLTAN_OK);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
