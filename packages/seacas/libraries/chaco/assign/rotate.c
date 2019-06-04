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

#include <math.h> // for cos, sin

void rotate2d(double **yvecs,   /* ptr to list of y-vectors (rotated) */
              int      nmyvtxs, /* length of yvecs */
              double   theta    /* angle to rotate by */
)
{
  double temp1; /* hold values for a while */
  double c, s;  /* cosine and sine of theta */
  int    i;     /* loop counter */

  s = sin(theta);
  c = cos(theta);

  for (i = 1; i <= nmyvtxs; i++) {
    temp1       = yvecs[1][i];
    yvecs[1][i] = c * temp1 + s * yvecs[2][i];
    yvecs[2][i] = -s * temp1 + c * yvecs[2][i];
  }
}

void rotate3d(double **yvecs,                         /* ptr to list of y-vectors (to be rotated) */
              int      nmyvtxs,                       /* length of yvecs */
              double theta, double phi, double gamma2 /* rotational parameters */
)
{
  double temp1, temp2;   /* hold values for a while */
  double ctheta, stheta; /* cosine and sine of theta */
  double cphi, sphi;     /* cosine and sine of phi */
  double cgamma, sgamma; /* cosine and sine of gamma */
  double onemcg;         /* 1.0 - cosine(gamma) */
  double a1, a2, a3;     /* rotation matrix entries */
  double b1, b2, b3;     /* rotation matrix entries */
  double c1, c2, c3;     /* rotation matrix entries */
  int    i;              /* loop counter */

  stheta = sin(theta);
  ctheta = cos(theta);
  sphi   = sin(phi);
  cphi   = cos(phi);
  sgamma = sin(gamma2);
  cgamma = cos(gamma2);

  onemcg = 1.0 - cgamma;

  a1 = cgamma + cphi * ctheta * onemcg * cphi * ctheta;
  a2 = sgamma * sphi + cphi * stheta * onemcg * cphi * ctheta;
  a3 = -sgamma * cphi * stheta + sphi * onemcg * cphi * ctheta;

  b1 = -sgamma * sphi + cphi * ctheta * onemcg * cphi * stheta;
  b2 = cgamma + cphi * stheta * onemcg * cphi * stheta;
  b3 = sgamma * cphi * ctheta + sphi * onemcg * cphi * stheta;

  c1 = sgamma * cphi * stheta + cphi * ctheta * onemcg * sphi;
  c2 = -sgamma * cphi * ctheta + cphi * stheta * onemcg * sphi;
  c3 = cgamma + sphi * onemcg * sphi;

  for (i = 1; i <= nmyvtxs; i++) {
    temp1 = yvecs[1][i];
    temp2 = yvecs[2][i];

    yvecs[1][i] = a1 * temp1 + b1 * temp2 + c1 * yvecs[3][i];
    yvecs[2][i] = a2 * temp1 + b2 * temp2 + c2 * yvecs[3][i];
    yvecs[3][i] = a3 * temp1 + b3 * temp2 + c3 * yvecs[3][i];
  }
}
