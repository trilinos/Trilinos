/*
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
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

/*******************************************************************
 *	STAT.C: Statistics routines:
 *
 *	newsample(n):	Add a new sample to the mean/average totals.
 *	mean():	        Returns the mean of the samples.
 *	deviation():	Returns the standard deviation of the sample.
 *	reset_mean():	Resets everything to 0.
 *
 ********************************************************************/

#include <math.h>

void newsample(int);
double mean(void);
double deviation(void);
double variance(void);
void reset_mean(void);

static long Numnums = 0;
static double Mean = 0.0;
static double StdDev = 0.0;

void newsample (int n)
{
  double TMean;

  /* See Knuth, TAOCP vol 2, 3rd edition, page 232 */
  TMean = Mean;
  Numnums++;
  Mean = TMean + (n - TMean) / Numnums;

  if (Numnums > 1)
    StdDev += (n - TMean) * (n - Mean);
}

double mean (void)
{
  return Mean;
}

double variance (void)
{
  return (Numnums > 1) ? StdDev/(Numnums-1) : 0.0;
}

double deviation (void)
{
  return sqrt(variance());
}

void reset_mean (void)
{
  Numnums = 0;
  Mean = 0;
  StdDev = 0;
}
