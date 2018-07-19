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

#include <math.h>

static double s, c;     /* sign and cosine of angle */
static double s2, cos2; /* squares of sign and cosine of angle */

double func2d(double coeffs[5], /* five different 4-way products */
              double theta      /* angular parameter */
              )

/* Returns value of penalty function at given angle. */
{
  double val; /* functional value */

  s    = sin(theta);
  c    = cos(theta);
  cos2 = c * c;
  s2   = s * s;

  val = (cos2 * cos2 + s2 * s2) * (coeffs[0] + coeffs[4]);
  val += 12 * cos2 * s2 * coeffs[2];
  val += 4 * (s2 * s * c - cos2 * c * s) * (coeffs[3] - coeffs[1]);

  return (val);
}

double grad2d(double coeffs[5], /* five different 4-way products */
              double theta      /* angular parameter */
              )

/* Returns 1st derivative of penalty function at given angle. */
{
  double val; /* functional value */

  s    = sin(theta);
  c    = cos(theta);
  cos2 = c * c;
  s2   = s * s;

  val = 4 * (cos2 * cos2 + s2 * s2) * (coeffs[1] - coeffs[3]);
  val += 24 * cos2 * s2 * (coeffs[3] - coeffs[1]);
  val += 4 * (s2 * s * c - cos2 * c * s) * (coeffs[0] + coeffs[4] - 6 * coeffs[2]);

  return (val);
}

double hess2d(double coeffs[5] /* five different 4-way products */
              )

/* Returns 2nd derivative of penalty function at given angle. */
{
  double val; /* functional value */

  val = -4 * (cos2 * cos2 + s2 * s2) * (coeffs[0] + coeffs[4] - 6 * coeffs[2]);
  val += 24 * s2 * cos2 * (coeffs[0] + coeffs[4] - 6 * coeffs[2]);
  val += 64 * (s2 * s * c - cos2 * c * s) * (coeffs[1] - coeffs[3]);

  return (val);
}
