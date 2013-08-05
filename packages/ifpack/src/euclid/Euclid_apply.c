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

#include "Euclid_dh.h"
#include "Mat_dh.h"
#include "Factor_dh.h"
#include "Parser_dh.h"
#include "TimeLog_dh.h"
#include "SubdomainGraph_dh.h"


static void scale_rhs_private (Euclid_dh ctx, double *rhs);
static void permute_vec_n2o_private (Euclid_dh ctx, double *xIN,
				     double *xOUT);
static void permute_vec_o2n_private (Euclid_dh ctx, double *xIN,
				     double *xOUT);

#undef __FUNC__
#define __FUNC__ "Euclid_dhApply"
void
Euclid_dhApply (Euclid_dh ctx, double *rhs, double *lhs)
{
  START_FUNC_DH double *rhs_, *lhs_;
  double t1, t2;

  t1 = MPI_Wtime ();

  /* default settings; for everything except PILU */
  ctx->from = 0;
  ctx->to = ctx->m;

  /* case 1: no preconditioning */
  if (!strcmp (ctx->algo_ilu, "none") || !strcmp (ctx->algo_par, "none"))
    {
      int i, m = ctx->m;
      for (i = 0; i < m; ++i)
	lhs[i] = rhs[i];
      goto END_OF_FUNCTION;
    }

  /*----------------------------------------------------------------
   * permute and scale rhs vector
   *----------------------------------------------------------------*/
  /* permute rhs vector */
  if (ctx->sg != NULL)
    {

      permute_vec_n2o_private (ctx, rhs, lhs);
      CHECK_V_ERROR;
      rhs_ = lhs;
      lhs_ = ctx->work2;
    }
  else
    {
      rhs_ = rhs;
      lhs_ = lhs;
    }

  /* scale rhs vector */
  if (ctx->isScaled)
    {

      scale_rhs_private (ctx, rhs_);
      CHECK_V_ERROR;
    }

  /* note: rhs_ is permuted, scaled; the input, "rhs" vector has
     not been disturbed.
   */

  /*----------------------------------------------------------------
   * big switch to choose the appropriate triangular solve
   *----------------------------------------------------------------*/

  /* sequential and mpi block jacobi cases */
  if (np_dh == 1 || !strcmp (ctx->algo_par, "bj"))
    {
      Factor_dhSolveSeq (rhs_, lhs_, ctx);
      CHECK_V_ERROR;
    }


  /* pilu case */
  else
    {
      Factor_dhSolve (rhs_, lhs_, ctx);
      CHECK_V_ERROR;
    }

  /*----------------------------------------------------------------
   * unpermute lhs vector
   * (note: don't need to unscale, because we were clever)
   *----------------------------------------------------------------*/
  if (ctx->sg != NULL)
    {
      permute_vec_o2n_private (ctx, lhs_, lhs);
      CHECK_V_ERROR;
    }

END_OF_FUNCTION:;

  t2 = MPI_Wtime ();
  /* collective timing for triangular solves */
  ctx->timing[TRI_SOLVE_T] += (t2 - t1);

  /* collective timing for setup+krylov+triSolves
     (intent is to time linear solve, but this is
     at best probelematical!)
   */
  ctx->timing[TOTAL_SOLVE_TEMP_T] = t2 - ctx->timing[SOLVE_START_T];

  /* total triangular solve count */
  ctx->its += 1;
  ctx->itsTotal += 1;

END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "scale_rhs_private"
void
scale_rhs_private (Euclid_dh ctx, double *rhs)
{
  START_FUNC_DH int i, m = ctx->m;
  REAL_DH *scale = ctx->scale;

  /* if matrix was scaled, must scale the rhs */
  if (scale != NULL)
    {
      for (i = 0; i < m; ++i)
	{
	  rhs[i] *= scale[i];
	}
    }
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "permute_vec_o2n_private"
void
permute_vec_o2n_private (Euclid_dh ctx, double *xIN, double *xOUT)
{
  START_FUNC_DH int i, m = ctx->m;
  int *o2n = ctx->sg->o2n_col;
  for (i = 0; i < m; ++i)
    xOUT[i] = xIN[o2n[i]];
END_FUNC_DH}


#undef __FUNC__
#define __FUNC__ "permute_vec_n2o_private"
void
permute_vec_n2o_private (Euclid_dh ctx, double *xIN, double *xOUT)
{
  START_FUNC_DH int i, m = ctx->m;
  int *n2o = ctx->sg->n2o_row;
  for (i = 0; i < m; ++i)
    xOUT[i] = xIN[n2o[i]];
END_FUNC_DH}
