/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#include "blas_dh.h"

#undef __FUNC__
#define __FUNC__ "matvec_euclid_seq"
void
matvec_euclid_seq (int n, int *rp, int *cval, double *aval, double *x,
		   double *y)
{
  START_FUNC_DH int i, j;
  int from, to, col;
  double sum;

  if (np_dh > 1)
    SET_V_ERROR ("only for sequential case!\n");

  {
    for (i = 0; i < n; ++i)
      {
	sum = 0.0;
	from = rp[i];
	to = rp[i + 1];
	for (j = from; j < to; ++j)
	  {
	    col = cval[j];
	    sum += (aval[j] * x[col]);
	  }
	y[i] = sum;
      }
  }
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "Axpy"
void
Axpy (int n, double alpha, double *x, double *y)
{
  START_FUNC_DH int i;

  for (i = 0; i < n; ++i)
    {
      y[i] = alpha * x[i] + y[i];
    }
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "CopyVec"
void
CopyVec (int n, double *xIN, double *yOUT)
{
  START_FUNC_DH int i;

  for (i = 0; i < n; ++i)
    {
      yOUT[i] = xIN[i];
    }
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "ScaleVec"
void
ScaleVec (int n, double alpha, double *x)
{
  START_FUNC_DH int i;

  for (i = 0; i < n; ++i)
    {
      x[i] *= alpha;
    }
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "InnerProd"
double
InnerProd (int n, double *x, double *y)
{
  START_FUNC_DH double result, local_result = 0.0;

  int i;

  for (i = 0; i < n; ++i)
    {
      local_result += x[i] * y[i];
    }

  if (np_dh > 1)
    {
      MPI_Allreduce (&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, comm_dh);
    }
  else
    {
      result = local_result;
    }

END_FUNC_VAL (result)}

#undef __FUNC__
#define __FUNC__ "Norm2"
double
Norm2 (int n, double *x)
{
  START_FUNC_DH double result, local_result = 0.0;
  int i;

  for (i = 0; i < n; ++i)
    {
      local_result += (x[i] * x[i]);
    }

  if (np_dh > 1)
    {
      MPI_Allreduce (&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, comm_dh);
    }
  else
    {
      result = local_result;
    }
  result = sqrt (result);
END_FUNC_VAL (result)}
