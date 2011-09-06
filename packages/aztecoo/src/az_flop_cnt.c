/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/

#include <stdio.h>
#include <math.h>

#include "az_aztec.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_calc_solve_flops(int options[], int total_its, double total_time,
                           int gn, double gnnz, int data_org[],
                           int proc_config[])

/*******************************************************************************

  Function which determines the Mflop/s rate of the iterative solver.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     double, MFlop/s rate.
  ============

  Parameter list:
  ===============

  options:         Determines specific solution method and other parameters.

  total_its:       Number of iterations performed by the iterative solver.

  total_time:      Total time required to perform total_its (seconds).

  gn:              Number of unknowns in the system (length of vectors).

  gnnz:            Number of nonzeros in the coefficient matrix.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  double inner_flops, daxpy_flops, matvec_flops, scale_flops;
  double Meg = 1.0e+06, iter_flops, precond_flops, total_flops, lu_flops;
  double blk_size, proc_blk_flops;
  int    krysub, solver_flag, conv_flag, N, N_Blk;

  /**************************** execution begins ******************************/

  solver_flag = options[AZ_solver];
  conv_flag   = options[AZ_conv];
  krysub      = options[AZ_kspace];
  N           = data_org[AZ_N_internal] + data_org[AZ_N_border];
  N_Blk       = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

  if (total_its < 0) total_its = -total_its;

  /* calculate the flops required for the basic linear algebra operations */

  daxpy_flops  = inner_flops = 2.0*((double) gn);
  matvec_flops = 2.0 * gnnz - ((double) gn);

  if (options[AZ_scaling] == AZ_BJacobi) {
    if (N_Blk == 0) blk_size = 1.0;
    else            blk_size = (double) N / (double) N_Blk;

    proc_blk_flops = ((double) N_Blk) * (2.0 * blk_size * blk_size * blk_size -
                                         1.5 * blk_size * blk_size +
                                         5.0 / 2.0 * blk_size);

    lu_flops         = (double) proc_config[AZ_N_procs] * proc_blk_flops;
    scale_flops   = 2.0 * gnnz + 2.0 * gnnz * blk_size + lu_flops;
  }
  else if (options[AZ_scaling] == AZ_row_sum) {
    scale_flops = 2.0 * gnnz + (double) gn;
  }
  else if (options[AZ_scaling] == AZ_none) {
    scale_flops = 0.0;
  }
  else if (options[AZ_scaling] == AZ_Jacobi) {
    scale_flops =  gnnz;
  }
  else {
    (void) AZ_printf_out("\t\tFlops not available for options[AZ_scaling] = %d\n",
                  options[AZ_scaling]);
    return -1.0;
  }

  /*
   * Calculate the flops required for the iterative solver minus the
   * preconditioning.
   */

  iter_flops = AZ_calc_iter_flops(solver_flag, inner_flops, daxpy_flops,
                                  matvec_flops, total_its, gnnz,
                                  (double) krysub);
  if (iter_flops < 0.0) return -1.0;

  if (conv_flag == 3)
    iter_flops += (double) total_its * (double) gn;

  /* calculate the flops required by the preconditioner */

  precond_flops = AZ_calc_precond_flops(solver_flag, options, 
                                        daxpy_flops, matvec_flops, total_its,
                                        gn, gnnz, data_org, proc_config);
  if (precond_flops < 0.0) return -1.0;

  /* I'm commenting this out because I can't figure out why it's here and it
     appears to be an error. SAH, 2/15/96

     if (solver_flag == AZ_gmres) precond_flops /= (double) krysub;
     */

  total_flops = scale_flops + iter_flops + precond_flops;

  if (proc_config[AZ_node] == 0) {
    (void) AZ_printf_out("\t\tscale_flops: %e\titer_flops: %e\n", scale_flops,
                  iter_flops);
    (void) AZ_printf_out("\t\tprecond_flops: %e\ttotal_flops: %e\n\n", precond_flops,
                  total_flops);
  }

  if (iter_flops < 0.0 || precond_flops < 0.0)
    return 0.0;
  else {
    if (Meg * total_time == 0.0)
      return 0.0;
    else
      return (total_flops / (Meg*total_time));
  }

} /* AZ_calc_solve_flops */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_calc_iter_flops(int solver_flag, double inner_flops,
                          double daxpy_flops, double matvec_flops,
                          int total_its, double gnnz, double K)

/*******************************************************************************

  Function which determines the number of flops for all iterations of the
  given solver.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     double, number of flops.
  ============

  Parameter list:
  ===============

  solver_flag:     Determines the specific iterative method used.

  inner_flops:     Number of flops in an inner product.

  daxpy_flops:     Number of flops in a daxpy operation.

  matvec_flops:    Number of flops in a matrix-vector multiply.

  total_its:       Number of iterations performed by the iterative solver.

  gnnz:            Number of nonzeros in the coefficient matrix.

  K:               Krylov subspace size (GMRES only).

*******************************************************************************/

{

  /**************************** execution begins ******************************/

  switch (solver_flag) {
  case AZ_cg:
    return (gnnz + 2.0*inner_flops + matvec_flops +
      total_its*(3.0*(inner_flops + daxpy_flops) + matvec_flops + 8.0));

  case AZ_analyze:
  case AZ_fixed_pt:
    return ( total_its*(inner_flops + 2.*daxpy_flops + matvec_flops ));

  case AZ_GMRESR:  /* this is not really right */
    return (gnnz + inner_flops + 2.0*matvec_flops +
      total_its*(matvec_flops + 1.5*inner_flops + 14.0) +
      0.5*total_its*K*(2.0*daxpy_flops + inner_flops + 7.0) +
      0.5*K*K*total_its + total_its/K*(2.0*matvec_flops + inner_flops + 3.0));

  case AZ_gmres:
    return (gnnz + inner_flops + 2.0*matvec_flops +
      total_its*(matvec_flops + 1.5*inner_flops + 14.0) +
      0.5*total_its*K*(2.0*daxpy_flops + inner_flops + 7.0) +
      0.5*K*K*total_its + total_its/K*(2.0*matvec_flops + inner_flops + 3.0));

  case AZ_cgs:
    return (gnnz + 1.5*inner_flops + 2.0*matvec_flops + daxpy_flops +
      total_its*(5.5*inner_flops + 2*matvec_flops + 3.0*daxpy_flops + 3.0));

  case AZ_tfqmr:
    return (gnnz + 2.5*inner_flops + 3.0*matvec_flops + daxpy_flops +
      total_its*(8.5*inner_flops + 2.0*matvec_flops + 4.0*daxpy_flops + 18.0));

  case AZ_bicgstab:
    return (gnnz + 1.5*inner_flops + 2.0*matvec_flops + daxpy_flops +
      total_its*(7.0*inner_flops + 2.0*matvec_flops + 2.0*daxpy_flops + 6.0));

  case AZ_lu:
    (void) AZ_printf_err(
                   "\t\tWARNING: Flop count not implemented for lu solver\n");
  return -1.0;

  default:
    (void) AZ_printf_out( "\t\tFlops not available for options[AZ_solver] = %d\n",
                   solver_flag);
  return -1.0;
  }

} /* AZ_calc_iter_flops */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_calc_precond_flops(int solver_flag, int options[], 
                             double daxpy_flops, double matvec_flops,
                             int total_its, int gn, double gnnz, int data_org[],
                             int proc_config[])

/*******************************************************************************

  Function which determines the number of flops for preconditioning done over
  all iterations.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     double, number of flops in preconditioning.
  ============

  Parameter list:
  ===============

  solver_flag:     Determines the specific iterative method used.

  options:         Determines specific solution method and other parameters.

  daxpy_flops:     Number of flops in a daxpy operation.

  matvec_flops:    Number of flops in a matrix-vector multiply.

  total_its:       Number of iterations performed by the iterative solver.

  gn:              Number of unknowns in the system (length of vectors).

  gnnz:            Number of nonzeros in the coefficient matrix.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see file Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int    power;
  double flops, fpower;
  double blk_size, proc_blk_flops, lu_flops;
  int    N_Blk, N, Num_Proc;

  /**************************** execution begins ******************************/

  N_Blk    = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
  N        = data_org[AZ_N_internal] + data_org[AZ_N_border];
  Num_Proc = proc_config[AZ_N_procs];
  power    = options[AZ_poly_ord];
  fpower   = fabs((double) power);

  switch(options[AZ_precond]) {
  case AZ_none:                 /* no preconditioning */
    return 0.0;

  case AZ_Jacobi:                       /* Jacobi preconditioning */
    if (N_Blk == 0) blk_size = 1.0;
    else            blk_size = (double) N / (double) N_Blk;

    proc_blk_flops = ((double) N_Blk) * (2.0 * blk_size * blk_size * blk_size -
                                         1.5 * blk_size * blk_size +
                                         5.0 / 2.0 * blk_size);
    lu_flops = (double) Num_Proc * proc_blk_flops;

    flops = ((double) total_its) * 4.0 * blk_size * (double) gn + lu_flops;
    break;

  case AZ_Neumann:              /* Neumann series polynomial precond. */
    if (!power) return 0.0;     /* diagonal scaling */

    flops  = (double) gn;
    flops += fpower * matvec_flops;
    flops += 3.0 * fpower * ((double) gn);

    if (solver_flag != AZ_cg && solver_flag != AZ_gmres)
      flops *= 2.0 * (double) total_its;
    else
      flops *= (double) total_its;

    flops += gnnz;
    break;

  case AZ_ls:                   /* Least squares polynomial precond. */
    if (!power) return 0.0;     /* diagonal scaling */

    flops  = (double) gn;
    flops += fpower * matvec_flops;
    flops += fpower * daxpy_flops;

    if (solver_flag != AZ_cg && solver_flag != AZ_gmres)
      flops *= 2.0 * (double) total_its;
    else
      flops *= (double) total_its;

    flops += gnnz;
    break;

  case AZ_sym_GS:               /* Symmetric GS. */
    if (!power) return 0.0;     /* diagonal scaling */

    flops  = (double) gn;
    flops += 2.0 * fpower * matvec_flops;
    flops += 2.0 * fpower * daxpy_flops;

    if (solver_flag != AZ_cg && solver_flag != AZ_gmres)
      flops *= 2.0 * (double) total_its;
    else
      flops *= (double) total_its;

    flops += gnnz;
    break;

  default:
    (void) AZ_printf_out( "\t\tFlops not available for options[AZ_precond] = %d\n",
                   options[AZ_precond]);
  return -1.0;
  }

  return flops;

} /* AZ_calc_precond_flops */
