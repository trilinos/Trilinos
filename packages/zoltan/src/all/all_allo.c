/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/
#ifndef lint
static char *cvs_all_allo_c =
  "$Id$";
#endif

#include <stdio.h>
#include <stdlib.h>
#include "all_allo_const.h"

#ifdef __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif


/******************************************************************************
 *
 *                    Dynamic Allocation of Multidimensional Arrays
 *-----------------------------------------------------------------------------
 *
 * Example Usage:
 *
 *     typedef  struct
 *       {      int     bus1;
 *              int     bus2;
 *              int     dest;
 *      }       POINT;
 *
 *      POINT    **points, corner;
 *
 *      points = (POINT **) LB_array_alloc (2, x, y, sizeof(POINT));
 *                               ^ ^ ^
 *                               | | |
 *         number of dimensions--+ | |
 *                                 | |
 *          first dimension max----+ |
 *                                   |
 *         second dimension max------+
 *
 *         (points may be now be used as if it were declared
 *          POINT points[x][y])
 *
 *          This particular version is limited to dimensions of 3 or less.
 *
 *      corner = points[2][3]; (refer to the structure as you would any array)
 *
 *      free (points); (frees the entire structure in one fell swoop)
 *
 *****************************************************************************/
/******************************************************************************
*       The following section is a commented section containing
*       an example main code:
*******************************************************************************
*double *array_alloc();
*main()
*{
*  int ***temp;
*   int *temp2;
*   int i, j, k;
*   int il, jl, kl;
*
*   malloc_debug(2);
*   il = 2;
*   jl = 3;
*   kl = 3;
*   temp = (int ***) array_alloc(3,il,jl,kl,sizeof(int));
*   for (i=0; i<il; i++) {
*      for (j=0; j<jl; j++) {
*         for (k=0; k<kl; k++) temp[i][j][k] = 1;
*      }
*   }
*
*   temp2 = (int *) malloc(10*sizeof(int));
*   for (i=0; i<10; i++) temp2[i] = 0;
*
*   for (i=0; i<il; i++) {
*      for (j=0; j<jl; j++) {
*         for (k=0; k<kl; k++) (void) printf(" %d\n", temp[i][j][k]);
*      }
*   }
*   malloc_verify();
*}
******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef __STDC__

double *LB_array_alloc(int numdim, ...)

#else

double *LB_array_alloc(va_alist)
va_dcl

#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

{
  int i, j;
  struct dim {
    long index;  /* Number of elements in the dimension  */
    long total;  /* Total number of elements             */
    long size;   /* Size of a single element in bytes    */
    long off;    /* offset from beginning of array       */
  } dim[4];      /* Info about each dimension            */

#ifndef __STDC__
  int numdim;           /* Number of dimensions                 */
#endif

  long    total;        /* Total size of the array              */
  double *dfield;       /* ptr to avoid lint complaints         */
  char   *field;        /* The multi-dimensional array          */
  char  **ptr;          /* Pointer offset                       */
  char   *data;         /* Data offset                          */
  va_list va;           /* Current pointer in the argument list */

#ifdef __STDC__
  va_start(va, numdim);
#else
  va_start(va);
  numdim = va_arg(va, int);
#endif

  if (numdim <= 0) {
    fprintf(stderr, "LB_array_alloc ERROR: number of dimensions, %d, is <=0\n",
            numdim);
    exit(-1);
  }
  else if (numdim > 4) {
    fprintf(stderr, "LB_array_alloc ERROR: number of dimensions, %d, is > 4\n",
            numdim);
    exit(-1);
  }

  dim[0].index = va_arg(va, int);

  if (dim[0].index <= 0) {
#ifdef DEBUG
    fprintf(stderr, "WARNING, LB_array_alloc called with first "
            "dimension <= 0, %d\n\twill return the nil pointer\n", 
            dim[0].index);
#endif
    return((double *) NULL);
  }

  dim[0].total = dim[0].index;
  dim[0].size  = sizeof(void *);
  dim[0].off   = 0;
  for (i = 1; i < numdim; i++) {
    dim[i].index = va_arg(va, int);
    if (dim[i].index <= 0) {
      fprintf(stderr, "WARNING: LB_array_alloc called with dimension %d <= 0, "
              "%d\n", i+1, dim[i].index);
      fprintf(stderr, "\twill return the nil pointer\n");
      return((double *) NULL);
    }
    dim[i].total = dim[i-1].total * dim[i].index;
    dim[i].size  = sizeof(void *);
    dim[i].off   = dim[i-1].off + dim[i-1].total * dim[i-1].size;
  }

  dim[numdim-1].size = va_arg(va, int);
  va_end(va);

  /* Round up the last offset value so data is properly aligned. */

  dim[numdim-1].off = dim[numdim-1].size *
    ((dim[numdim-1].off+dim[numdim-1].size-1)/dim[numdim-1].size);

  total = dim[numdim-1].off + dim[numdim-1].total * dim[numdim-1].size;

  dfield = (double *) LB_smalloc((int) total);
  field  = (char *) dfield;

  for (i = 0; i < numdim - 1; i++) {
    ptr  = (char **) (field + dim[i].off);
    data = (char *) (field + dim[i+1].off);
    for (j = 0; j < dim[i].total; j++) {
      ptr[j] = data + j * dim[i+1].size * dim[i+1].index;
    }
  }

  return dfield;

} /* LB_array_alloc */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Safe version of malloc.  Does not initialize memory .*/

/* Modified by Scott Hutchinson (1421) 20 January 1993 */

double *LB_smalloc(int n)

{
  double *pntr;           /* return value */

  if (n < 0) {
    fprintf(stderr, "LB_smalloc ERROR: Non-positive argument. (%d)\n", n);
    exit(-1);
  }
  else if (n == 0)
    pntr = NULL;
  else
    pntr = (double *) malloc((unsigned) n);

  if (pntr == NULL && n != 0) {
    fprintf(stderr, "LB_smalloc : Out of space - number of bytes "
            "requested = %d\n", n);
    exit(0);
  }

  return pntr;

} /* LB_smalloc */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_safe_free (void **ptr)
{
/*
 *  This version of free calls the system's free function
 *  with maximum error checking. It also doesn't call free if ptr is
 *  the NULL pointer.
 */
 
  if (*ptr != NULL) {
 
    free(*ptr);
 
    /*
     *  Set the value of ptr to NULL, so that further references
     *  to it will be flagged.
     */
 
    *ptr = NULL;
  }
}  /* LB_safe_free */

/*****************************************************************************/
/*                      END of all_allo.c                                     */
/*****************************************************************************/
