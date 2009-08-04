/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
