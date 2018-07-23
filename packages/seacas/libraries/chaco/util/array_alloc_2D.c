/*
 * Copyright (c) 2005-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
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
 *     * Neither the name of NTESS nor the names of its
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
  char * field;       /* The multi-dimensional array */
  char **ptr;         /* Pointer offset */
  char * data;        /* Data offset */
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
