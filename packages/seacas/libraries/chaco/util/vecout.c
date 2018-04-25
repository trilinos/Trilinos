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

#include "defs.h"
#include <math.h>
#include <stdio.h>

/* Print vertically range of double vector. */
void vecout(double *vec, int beg, int end, char *tag, char *file_name)
{
  FILE *file;
  int   i;

  if (file_name != NULL) {
    file = fopen(file_name, "w");
  }
  else {
    file = stdout;
  }

  fprintf(file, "%s:\n", tag);
  for (i = beg; i <= end; i++) {
    if (fabs(vec[i]) >= 1.0e-16) {
      fprintf(file, "%2d.   %24.16f\n", i, vec[i]);
    }
    else {
      fprintf(file, "%2d.         %g \n", i, vec[i]);
    }
  }
  if (file_name != NULL) {
    fclose(file);
  }
}

/* Scale the eigenvector such that the first element is non-negative.
 * This is an attempt to reduce machine-to-machine variance which can
 * result from calculated eigenvectors being mirror-images of each
 * other due to small roundoff..
 */
void vecnorm(double *vec, int beg, int end)
{
  int i;
  if (vec[beg] < 0.0) {
    for (i = beg; i <= end; i++) {
      vec[i] *= -1.0;
    }
  }
}
