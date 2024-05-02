/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for smalloc_ret
#include <stddef.h>  // for size_t
#include <stdio.h>   // for NULL

/* Dynamically allocate a 2 dimensional array. */
/* Return instead of dying if out of space. */

void *array_alloc_2D_ret(size_t dim1, size_t dim2, size_t size)
/* size of first dimension */
/* size of second dimension */
/* size of array elements */
{
  size_t total;       /* Total size of the array */
  size_t aligned_dim; /* dim1 or dim1+1 to ensure data alignment */
  size_t offset;      /* offset of array elements */
  char  *field;       /* The multi-dimensional array */
  char **ptr;         /* Pointer offset */
  char  *data;        /* Data offset */
  size_t j;           /* loop counter */

  aligned_dim = (dim1 % 2) ? dim1 + 1 : dim1;
  offset      = aligned_dim * sizeof(void *);
  total       = offset + dim1 * dim2 * size;
  field       = smalloc_ret(total);

  if (field != NULL) {
    ptr  = (char **)field;
    data = field;
    data += offset;
    for (j = 0; j < dim1; j++) {
      ptr[j] = data + j * size * dim2;
    }
  }

  return (field);
}
