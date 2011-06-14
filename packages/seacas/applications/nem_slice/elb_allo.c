/*
 * Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include "elb_allo_const.h"

static void *smalloc (size_t n);

/******************************************************************************
 *
 *                    Dynamic Allocation of Multidimensional Arrays
 *-----------------------------------------------------------------------------
 *
 * Example Usage:
 *
 *     typedef	struct
 *       {	int	bus1;
 *              int	bus2;
 *              int	dest;
 *      }       POINT;
 *
 *      POINT    **points, corner;
 *
 *      points = (POINT **) array_alloc (2, x, y, sizeof(POINT));
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
/*****************************************************************************/
/*****************************************************************************/
void *array_alloc (int numdim , ...)
{
   int	i, j;
   struct dim {
      size_t	index;	/* Number of elements in the dimension	*/
      size_t	total;	/* Total number of elements 		*/
      size_t	size;	/* Size of a single element in bytes	*/
      size_t	off;	/* offset from beginning of array	*/
   }	dim[3];		/* Info about each dimension 		*/
   size_t total;	/* Total size of the array		*/
   void *dfield;	/* ptr to avoid lint complaints		*/
   char *field;		/* The multi-dimensional array		*/
   char **ptr;		/* Pointer offset			*/
   char *data;		/* Data offset				*/
   va_list va;		/* Current pointer in the argument list	*/

   va_start(va, numdim);

   if (numdim <= 0) {
     fprintf(stderr,
             "array_alloc ERROR: number of dimensions, %d, is <=0\n", numdim);
        return NULL;
   } else if (numdim > 3) {
     fprintf(stderr,
             "array_alloc ERROR: number of dimensions, %d, is > 3\n", numdim);
     return NULL;
   }

   dim[0].index = va_arg(va, int);

   if (dim[0].index <= 0)
     return NULL;

   dim[0].total = dim[0].index;
   dim[0].size = sizeof(void *);
   dim[0].off = 0;
   for (i=1; i<numdim; i++) {
      dim[i].index = va_arg(va, int);
      if (dim[i].index <= 0)
	return NULL;

      dim[i].total = dim[i-1].total * dim[i].index;
      dim[i].size = sizeof(void *);
      dim[i].off = dim[i-1].off + dim[i-1].total * dim[i-1].size;
   }
   dim[numdim-1].size = va_arg(va, int);
   va_end(va);

   /* Round up the last offset value so data is properly aligned. */
   dim[numdim-1].off = dim[numdim-1].size *
       ((dim[numdim-1].off+dim[numdim-1].size-1)/dim[numdim-1].size);

   total = dim[numdim-1].off + dim[numdim-1].total * dim[numdim-1].size;

   dfield = (void *) smalloc(total);
   field = (char *) dfield;

   for (i=0; i<numdim - 1; i++) {
      ptr = (char **) (field + dim[i].off);
      data = (char *) (field + dim[i+1].off);
      for (j=0; j<dim[i].total; j++) {
	 ptr[j] = data + j * dim[i+1].size * dim[i+1].index;
      }
   }
   return(dfield);
}

/* Safe version of malloc.  Does not initialize memory .*/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void *smalloc (size_t n)
{
	void *pntr;		/* return value */

	if (n == 0)
	  pntr = NULL;
	else
	  pntr = malloc(n);

	if(pntr == NULL && n != 0) {
	  fprintf(stderr, "smalloc: Out of space - number of bytes "
		  "requested = %ld\n", n);
	  exit(0);
	}
   	return (pntr);
}

/*****************************************************************************/
/*                       END of elb_allo.c                                   */
/*****************************************************************************/
