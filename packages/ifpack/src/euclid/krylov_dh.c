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
#include "krylov_dh.h"
#include "Mem_dh.h"
#include "Parser_dh.h"
#include "Mat_dh.h"


#undef __FUNC__
#define __FUNC__ "bicgstab_euclid"
void
bicgstab_euclid (Mat_dh A, Euclid_dh ctx, double *x, double *b, int *itsOUT)
{
  START_FUNC_DH int its, m = ctx->m;
  bool monitor;
  int maxIts = ctx->maxIts;
  double atol = ctx->atol, rtol = ctx->rtol;

  /* scalars */
  double alpha, alpha_1,
    beta_1,
    widget, widget_1, rho_1, rho_2, s_norm, eps, exit_a, b_iprod, r_iprod;

  /* vectors */
  double *t, *s, *s_hat, *v, *p, *p_hat, *r, *r_hat;

  monitor = Parser_dhHasSwitch (parser_dh, "-monitor");

  /* allocate working space */
  t = (double *) MALLOC_DH (m * sizeof (double));
  s = (double *) MALLOC_DH (m * sizeof (double));
  s_hat = (double *) MALLOC_DH (m * sizeof (double));
  v = (double *) MALLOC_DH (m * sizeof (double));
  p = (double *) MALLOC_DH (m * sizeof (double));
  p_hat = (double *) MALLOC_DH (m * sizeof (double));
  r = (double *) MALLOC_DH (m * sizeof (double));
  r_hat = (double *) MALLOC_DH (m * sizeof (double));

  /* r = b - Ax */
  Mat_dhMatVec (A, x, s);	/* s = Ax */
  CopyVec (m, b, r);		/* r = b */
  Axpy (m, -1.0, s, r);		/* r = b-Ax */
  CopyVec (m, r, r_hat);	/* r_hat = r */

  /* compute stopping criteria */
  b_iprod = InnerProd (m, b, b);
  CHECK_V_ERROR;
  exit_a = atol * atol * b_iprod;
  CHECK_V_ERROR;		/* absolute stopping criteria */
  eps = rtol * rtol * b_iprod;	/* relative stoping criteria (residual reduction) */

  its = 0;
  while (1)
    {
      ++its;
      rho_1 = InnerProd (m, r_hat, r);
      if (rho_1 == 0)
	{
	  SET_V_ERROR ("(r_hat . r) = 0; method fails");
	}

      if (its == 1)
	{
	  CopyVec (m, r, p);	/* p = r_0 */
	  CHECK_V_ERROR;
	}
      else
	{
	  beta_1 = (rho_1 / rho_2) * (alpha_1 / widget_1);

	  /* p_i = r_(i-1) + beta_(i-1)*( p_(i-1) - w_(i-1)*v_(i-1) ) */
	  Axpy (m, -widget_1, v, p);
	  CHECK_V_ERROR;
	  ScaleVec (m, beta_1, p);
	  CHECK_V_ERROR;
	  Axpy (m, 1.0, r, p);
	  CHECK_V_ERROR;
	}

      /* solve M*p_hat = p_i */
      Euclid_dhApply (ctx, p, p_hat);
      CHECK_V_ERROR;

      /* v_i = A*p_hat */
      Mat_dhMatVec (A, p_hat, v);
      CHECK_V_ERROR;

      /* alpha_i = rho_(i-1) / (r_hat^T . v_i ) */
      {
	double tmp = InnerProd (m, r_hat, v);
	CHECK_V_ERROR;
	alpha = rho_1 / tmp;
      }

      /* s = r_(i-1) - alpha_i*v_i */
      CopyVec (m, r, s);
      CHECK_V_ERROR;
      Axpy (m, -alpha, v, s);
      CHECK_V_ERROR;

      /* check norm of s; if small enough:
       * set x_i = x_(i-1) + alpha_i*p_i and stop.
       * (Actually, we use the square of the norm)
       */
      s_norm = InnerProd (m, s, s);
      if (s_norm < exit_a)
	{
	  SET_INFO ("reached absolute stopping criteria");
	  break;
	}

      /* solve M*s_hat = s */
      Euclid_dhApply (ctx, s, s_hat);
      CHECK_V_ERROR;

      /* t = A*s_hat */
      Mat_dhMatVec (A, s_hat, t);
      CHECK_V_ERROR;

      /* w_i = (t . s)/(t . t) */
      {
	double tmp1, tmp2;
	tmp1 = InnerProd (m, t, s);
	CHECK_V_ERROR;
	tmp2 = InnerProd (m, t, t);
	CHECK_V_ERROR;
	widget = tmp1 / tmp2;
      }

      /* x_i = x_(i-1) + alpha_i*p_hat + w_i*s_hat */
      Axpy (m, alpha, p_hat, x);
      CHECK_V_ERROR;
      Axpy (m, widget, s_hat, x);
      CHECK_V_ERROR;

      /* r_i = s - w_i*t */
      CopyVec (m, s, r);
      CHECK_V_ERROR;
      Axpy (m, -widget, t, r);
      CHECK_V_ERROR;

      /* check convergence; continue if necessary;
       * for continuation it is necessary thea w != 0.
       */
      r_iprod = InnerProd (m, r, r);
      CHECK_V_ERROR;
      if (r_iprod < eps)
	{
	  SET_INFO ("stipulated residual reduction achieved");
	  break;
	}

      /* monitor convergence */
      if (monitor && myid_dh == 0)
	{
	  fprintf (stderr, "[it = %i] %e\n", its, sqrt (r_iprod / b_iprod));
	}

      /* prepare for next iteration */
      rho_2 = rho_1;
      widget_1 = widget;
      alpha_1 = alpha;

      if (its >= maxIts)
	{
	  its = -its;
	  break;
	}
    }

  *itsOUT = its;

  FREE_DH (t);
  FREE_DH (s);
  FREE_DH (s_hat);
  FREE_DH (v);
  FREE_DH (p);
  FREE_DH (p_hat);
  FREE_DH (r);
  FREE_DH (r_hat);
END_FUNC_DH}

#undef __FUNC__
#define __FUNC__ "cg_euclid"
void
cg_euclid (Mat_dh A, Euclid_dh ctx, double *x, double *b, int *itsOUT)
{
  START_FUNC_DH int its, m = A->m;
  double *p, *r, *s;
  double alpha, beta, gamma, gamma_old, eps, bi_prod, i_prod;
  bool monitor;
  int maxIts = ctx->maxIts;
  /* double atol = ctx->atol */
  double rtol = ctx->rtol;

  monitor = Parser_dhHasSwitch (parser_dh, "-monitor");

  /* compute square of absolute stopping threshold  */
  /* bi_prod = <b,b> */
  bi_prod = InnerProd (m, b, b);
  CHECK_V_ERROR;
  eps = (rtol * rtol) * bi_prod;

  p = (double *) MALLOC_DH (m * sizeof (double));
  s = (double *) MALLOC_DH (m * sizeof (double));
  r = (double *) MALLOC_DH (m * sizeof (double));

  /* r = b - Ax */
  Mat_dhMatVec (A, x, r);	/* r = Ax */
  CHECK_V_ERROR;
  ScaleVec (m, -1.0, r);	/* r = b */
  CHECK_V_ERROR;
  Axpy (m, 1.0, b, r);		/* r = r + b */
  CHECK_V_ERROR;

  /* solve Mp = r */
  Euclid_dhApply (ctx, r, p);
  CHECK_V_ERROR;

  /* gamma = <r,p> */
  gamma = InnerProd (m, r, p);
  CHECK_V_ERROR;

  its = 0;
  while (1)
    {
      ++its;

      /* s = A*p */
      Mat_dhMatVec (A, p, s);
      CHECK_V_ERROR;

      /* alpha = gamma / <s,p> */
      {
	double tmp = InnerProd (m, s, p);
	CHECK_V_ERROR;
	alpha = gamma / tmp;
	gamma_old = gamma;
      }

      /* x = x + alpha*p */
      Axpy (m, alpha, p, x);
      CHECK_V_ERROR;

      /* r = r - alpha*s */
      Axpy (m, -alpha, s, r);
      CHECK_V_ERROR;

      /* solve Ms = r */
      Euclid_dhApply (ctx, r, s);
      CHECK_V_ERROR;

      /* gamma = <r,s> */
      gamma = InnerProd (m, r, s);
      CHECK_V_ERROR;

      /* set i_prod for convergence test */
      i_prod = InnerProd (m, r, r);
      CHECK_V_ERROR;

      if (monitor && myid_dh == 0)
	{
	  fprintf (stderr, "iter = %i  rel. resid. norm: %e\n", its,
		   sqrt (i_prod / bi_prod));
	}

      /* check for convergence */
      if (i_prod < eps)
	break;

      /* beta = gamma / gamma_old */
      beta = gamma / gamma_old;

      /* p = s + beta p */
      ScaleVec (m, beta, p);
      CHECK_V_ERROR;
      Axpy (m, 1.0, s, p);
      CHECK_V_ERROR;

      if (its >= maxIts)
	{
	  its = -its;
	  break;
	}
    }

  *itsOUT = its;

  FREE_DH (p);
  FREE_DH (s);
  FREE_DH (r);
END_FUNC_DH}
