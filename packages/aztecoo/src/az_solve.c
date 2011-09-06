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

/*
 * All of these routines form the interface between the code drivers and the
 * user's choice of algorithm and PDE problem.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"

#define AZ_oldprecond         0
#define AZ_oldsubdomain_solve 1
#define AZ_oldreorder         2
#define AZ_oldpre_calc        3
#define AZ_SAVE_SIZE          5

int az_iterate_id = 0;
int az_iterate_id_increment = 10000;
int az_iterate_recursion_level = 0;

void AZ_iterate(double x[], double b[], int options[], double params[],
                double status[], int proc_config[], AZ_MATRIX *Amat,
                AZ_PRECOND *precond, struct AZ_SCALING *scaling)

/*******************************************************************************
This is the new Aztec interface. This routine calls AZ_oldsolve() passing
in for example Amat->indx for the indx[] parameter in AZ_oldsolve().

NOTE: Users can still invoke AZ_solve() in the old Aztec way. AZ_solve
      also calls AZ_oldsolve(). However, matrix-free and coarse grid
       capabilities are not available via AZ_solve().


  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============


  x:               Contains the result vector x.

  b:               Contains the vector b.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  Amat:            Structure used to represent the matrix (see az_aztec.h
                   and Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  precond:         Structure used to represent the preconditionner
                   (see az_aztec.h ad Aztec User's Guide).

  scaling:         Structure used to represent the  scaling
                   (see az_aztec.h ad Aztec User's Guide).

******************************************************************************/

{
  double  total_time, start_t;
  int     prec_allocated = 0;
  struct AZ_SCALING *scale2;

  /* az_iterate_id becomes nonzero only if AZ_iterate is being called
     recursively -- i.e., if az_iterate_recursion_level has already
     been incremented. */
  az_iterate_id -= az_iterate_id_increment*az_iterate_recursion_level;

  /* now increment az_iterate_recursion_level */
  ++az_iterate_recursion_level;

  if (scaling == NULL)
    scale2 = AZ_scaling_create();
  else scale2 = scaling;

  AZ__MPI_comm_space_ok();
  if (Amat->mat_create_called != 1) {
    if (proc_config[AZ_node] == 0) {
      AZ_printf_out("AZ_iterate: AZ_matrix_create(int) should be called to\n");
      AZ_printf_out("            create matrix object (Amat) to be solved!\n");
    }
    exit(1);
  }
  if (precond == NULL) {
    if (options[AZ_precond] == AZ_user_precond) {
      if (proc_config[AZ_node] == 0) {
        AZ_printf_out("AZ_iterate: Can not use NULL for precond argument when\n");
        AZ_printf_out("            options[AZ_precond] == AZ_user_precond.\n");
      }
      exit(1);
    }
    precond = AZ_precond_create(Amat, AZ_precondition, NULL);
    prec_allocated = 1;
  }
  if (precond->prec_create_called != 1) {
    if (proc_config[AZ_node] == 0) {
      AZ_printf_out("AZ_iterate: AZ_precond_create should be called to\n   ");
      AZ_printf_out("       create preconditioning object!\n");
    }
    exit(1);
  }
  if (precond->Pmat->mat_create_called != 1) {
    if (proc_config[AZ_node] == 0) {
      AZ_printf_out("AZ_iterate: AZ_matrix_create(int) should be called to\n   ");
      AZ_printf_out("       create preconditioning matrix object (precond->Pmat)!\n");
    }
    exit(1);
  }
  if (Amat->matvec == NULL) {
    if (proc_config[AZ_node] == 0) {
      AZ_printf_out("AZ_iterate: Matrix vector product needs to be set via ");
      AZ_printf_out("AZ_set_MSR(...),\n             AZ_set_VBR(...), or ");
      AZ_printf_out("AZ_set_MATFREE(...).\n");
    }
    exit(1);
  }

  Amat->data_org[AZ_name] += az_iterate_id;
  precond->Pmat->data_org[AZ_name] += az_iterate_id;

  scale2->mat_name = Amat->data_org[AZ_name];

  AZ_iterate_setup(options, params, proc_config, Amat, precond);
  AZ_sync(proc_config);
  start_t = AZ_second();
  AZ_oldsolve(x, b, options, params, status, proc_config, Amat, precond,scale2);
  total_time = AZ_gmax_double(AZ_second() - start_t, proc_config);

  status[AZ_solve_time] = total_time;
  if ((options[AZ_output] != AZ_none) && (options[AZ_output] != AZ_warnings)) {
    if (proc_config[AZ_node] == 0) {
      (void) AZ_printf_out("\n\n\t\tSolution time: %f (sec.)\n", total_time);
      (void) AZ_printf_out("\t\ttotal iterations: %d\n", (int) status[AZ_its]);
    }
  }

  if (options[AZ_diagnostics]==AZ_all) 
    AZ_flop_rates(Amat->data_org,Amat->indx,Amat->bpntr, Amat->bindx,
		  options, status, total_time, proc_config);

  AZ_iterate_finish(options, Amat, precond);

  precond->Pmat->data_org[AZ_name] -= az_iterate_id;
  Amat->data_org[AZ_name] -= az_iterate_id;

  /* decrement az_iterate_recursion_level in preparation for leaving
     this AZ_iterate scope. */
  --az_iterate_recursion_level;

  /* if az_iterate_recursion_level is nonzero at this point, AZ_iterate is
     being called recursively and so we'll alter az_iterate_id by the
     opposite amount that it was altered above when we entered this function.
  */
  az_iterate_id += az_iterate_id_increment*az_iterate_recursion_level;

  if (prec_allocated)  AZ_precond_destroy(&precond);
  if (scaling == NULL) AZ_scaling_destroy(&scale2);
}

/* AZ_iterate*/


void AZ_solve(double x[], double b[], int options[], double params[],
              int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
              double val[], int data_org[], double status[], int proc_config[])

     /*******************************************************************************
  In order to assure backward compatibility  Aztec's previous AZ_solve()
  has been  renamed to AZ_oldsolve(), and this routine defines the 3 new
  parameters used by  AZ_oldsolve and called AZ_oldsolve.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               On input, contains the initial guess. On output contains the                   
                   solution to the linear system.

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  val:             Array containing the nonzero entries of the matrix (see Aztec User's Guide).

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/


{
  AZ_MATRIX *Amat;
  AZ_PRECOND *precond;
  double  total_time, start_t;
  struct AZ_SCALING *scale2;

  scale2 = AZ_scaling_create();
  AZ__MPI_comm_space_ok();
  Amat    = AZ_matrix_create(data_org[AZ_N_internal]+data_org[AZ_N_border]);
  precond = AZ_precond_create(Amat, AZ_precondition, NULL);

  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX)
    AZ_set_MSR(Amat, bindx, val, data_org, 0, NULL, AZ_LOCAL);
  else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
    AZ_set_VBR(Amat, rpntr, cpntr, bpntr, indx, bindx, val,
               data_org, 0, NULL, AZ_LOCAL);
  else {
    AZ_printf_err("Unknown matrix type (%d)\n",data_org[AZ_matrix_type]);
    AZ_printf_err("Matrix-free is now available via AZ_iterate()\n");
    exit(1);
  }

  if (options[AZ_precond] == AZ_user_precond) {
    AZ_printf_err("Unknown preconditioning options[AZ_precond] =  (%d)\n",options[AZ_precond]);
    AZ_printf_err("User preconditioning is now available via AZ_iterate()\n");
    exit(1);
  }

  options[AZ_recursion_level] = 0;

  if (options[AZ_pre_calc] != AZ_reuse) {
    (void) AZ_manage_memory(0,AZ_EVERYBODY_BUT_CLEAR,(Amat->data_org)[AZ_name],"kvecs",(int *) 0);
  }
  (void) AZ_manage_memory(0, AZ_CLEAR, AZ_SYS+az_iterate_id, (char *) 0, (int *) 0);

  /* output solver, scaling, and preconditioning options */

  AZ_print_call_iter_solve(options, params, proc_config[AZ_node], 0, Amat, precond);


  AZ_sync(proc_config);
  start_t = AZ_second();

  AZ_oldsolve(x, b, options, params, status, proc_config, Amat, precond,scale2);

  total_time = AZ_gmax_double(AZ_second() - start_t, proc_config);

  status[AZ_solve_time] = total_time;
  if ((options[AZ_output] != AZ_none) && (options[AZ_output] != AZ_warnings)) {
    if (proc_config[AZ_node] == 0) {
      (void) AZ_printf_out("\n\n\t\tSolution time: %f (sec.)\n", total_time);
      (void) AZ_printf_out("\t\ttotal iterations: %d\n", (int) status[AZ_its]);
    }
  }

  if (options[AZ_diagnostics]==AZ_all) 
    AZ_flop_rates(data_org,indx,bpntr, bindx, options, status, total_time,
		  proc_config);

  if (options[AZ_keep_info] == 0)
    (void) AZ_manage_memory(0,AZ_CLEAR,(Amat->data_org)[AZ_name],(char *) 0,
                            (int *) 0);

  (void) AZ_manage_memory(0, AZ_CLEAR, AZ_SYS+az_iterate_id, (char *) 0, (int *) 0);

  AZ_precond_destroy(&precond);
  AZ_matrix_destroy(&Amat);
  AZ_scaling_destroy(&scale2);

}
/* AZ_solve*/

#include <string.h>
#ifdef USING_ML
#include "ml_include.h"
#endif

void AZ_oldsolve(double x[], double b[], int options[], double params[],
                 double status[], int proc_config[], AZ_MATRIX *Amat,
                 AZ_PRECOND *precond, struct AZ_SCALING *scaling)


     /*******************************************************************************
  Aztec's previous AZ_solve() is renamed to AZ_oldsolve() with 3 new
  parameters appended to it: Amat, precond, scaling.
  This routine is never called directly by an application. It is only
  used internally by Aztec.

  Solve the system of equations given in the VBR format using an iterative
  method specified by 'options[AZ_solver]'. Store the result in 'x'.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  val:             Array containing the nonzero entries of the matrix (see
                   Aztec User's Guide).

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  Amat:            Structure used to represent the matrix (see az_aztec.h
                   and Aztec User's Guide).

  precond:         Structure used to represent the preconditionner
                   (see az_aztec.h ad Aztec User's Guide).

  scaling:         Structure used to represent the  scaling
                   (see az_aztec.h ad Aztec User's Guide).


*******************************************************************************/

{

  /* local variables */

  int     i, j;
  int t1, NNN;
  double *newparams, largest, dtemp;
  int     *data_org;
  char tag[80];
  int   save_old_values[7], changed = 0;
  char  tstr[15];
  struct AZ_CONVERGE_STRUCT *conv_info;
  int size1;
#ifdef HAVE_AZLU
  int     itemp2, itemp3, nonzeros;
  int     N, N_Blk;
  struct  context context;
  struct aztec_choices aztec_choices;
  int    nz_used;
  int bandwidth, extra_factor_nonzeros;
  AZ_MATRIX Amat2;
  int *bindx2, N_nz_factors;
  double *val2;
#endif
#ifdef ML
  int largest_index;
  double global_largest;
  int size2, *ibuf, allocated, row_length, *already_printed, *ibuf2, row2_length;
  double *tempv, *tempy, *tempz, *dbuf;
  char boundary;
  ML *ml;
#endif


  /**************************** execution begins ******************************/
  newparams = params;
  /* If not using AztecOO convergence test, we must create one */
  if (options[AZ_conv]!=AZTECOO_conv_test)
    conv_info = AZ_converge_create();
  else
    conv_info = Amat->conv_info;
  data_org  = Amat->data_org;
  conv_info->scaling = scaling;
  AZ__MPI_comm_space_ok();
  status[AZ_Aztec_version] = -1.;
  save_old_values[6] = options[AZ_ignore_scaling];
  if ( options[AZ_conv] == AZ_expected_values) {
    options[AZ_ignore_scaling] = AZ_TRUE;
    NNN = data_org[AZ_N_internal] + data_org[AZ_N_border];
    sprintf(tag, "some weights %d %d %d", data_org[AZ_name],
            options[AZ_recursion_level],NNN);
    newparams = (double *) AZ_manage_memory((AZ_PARAMS_SIZE+NNN)*
                                            sizeof(double), AZ_ALLOC,
                                            data_org[AZ_name],tag,&i);
    if ((options[AZ_pre_calc] == AZ_reuse)& (options[AZ_scaling] != AZ_none))
      AZ_scale_f(AZ_SCALE_SOL, Amat, options, b, x, proc_config, scaling);

    AZ_abs_matvec_mult(x,&(newparams[AZ_weights]), Amat, proc_config);
    if ((options[AZ_pre_calc] == AZ_reuse)& (options[AZ_scaling] != AZ_none)){
      AZ_scale_f(AZ_INVSCALE_SOL, Amat, options, b, x, proc_config, scaling);
      AZ_scale_f(AZ_INVSCALE_RHS, Amat, options, &(newparams[AZ_weights]),
                 x, proc_config, scaling);
    }
    largest = 0.;
    for (i = 0; i < NNN; i++)
      if (newparams[i+AZ_weights] > largest)
        largest = newparams[AZ_weights+i];
    largest *= 100.;
    for (i = 0; i < NNN; i++)
      if (newparams[i+AZ_weights] == 0.) newparams[AZ_weights+i]=largest;
    for (i = 0; i < AZ_weights; i++) newparams[i] = params[i];
  }


  /* grab the version string and rip out the ending dots and put it into status */
  AZ_version(tstr);
  t1 = 0; j =0;
  for (i = 0; i < (int) strlen(tstr); i++)  {
    if (tstr[i] == '.') { t1++; if (t1 == 1) tstr[j++] = tstr[i];}
    else tstr[j++] = tstr[i];
  }
  tstr[j] = '\0';
  sscanf(&(tstr[6]),"%lf", &(status[AZ_Aztec_version]));
  if (!AZ_oldsolve_setup(x, b, options, newparams, status, proc_config,
                         Amat, precond, save_old_values,scaling)) return;

  /* solve the system */

  AZ_flush_out();
  switch (options[AZ_solver]) {

  case AZ_cg:

    /* conjugate gradient */

    AZ_pcg_f(b,x,&(newparams[AZ_weights]),options,params,proc_config,status,Amat,precond,conv_info);

    break;

  case AZ_cg_condnum:

    /* conjugate gradient with condition number estimate */

    AZ_pcg_f_condnum(b,x,&(newparams[AZ_weights]),options,params,proc_config,status,Amat,precond,conv_info);

    break;

  case AZ_gmres:

    /* GMRES */

    AZ_pgmres(b, x, &(newparams[AZ_weights]), options, params, proc_config,
              status, Amat, precond, conv_info);

    break;

  case AZ_gmres_condnum:

    /* GMRES with condition number estimate */

    AZ_pgmres_condnum(b, x, &(newparams[AZ_weights]), options, params, proc_config,
		      status, Amat, precond, conv_info);

    break;

  case AZ_fixed_pt:
    AZ_fix_pt(b, x, &(newparams[AZ_weights]), options, params, proc_config,
              status, Amat, precond, conv_info);

    break;

  case AZ_analyze:


    for (i=0;i < Amat->data_org[AZ_N_internal]+Amat->data_org[AZ_N_border]; i++)
      b[i] = 0.0;
    AZ_random_vector(x, data_org, proc_config);

    dtemp=AZ_gdot(Amat->data_org[AZ_N_internal]+Amat->data_org[AZ_N_border],
                  x, x, proc_config);
    dtemp = sqrt(dtemp);
    for (j=0;j<Amat->data_org[AZ_N_internal]+Amat->data_org[AZ_N_border];j++)
      x[j] = x[j]/dtemp;

    params[AZ_temp] = 1.;
    params[AZ_tol] = 1.e-40;
    i = options[AZ_max_iter];
    while (i > 0) {
      if (options[AZ_max_iter] > 10) options[AZ_max_iter] = 10;
      AZ_fix_pt(b, x, &(newparams[AZ_weights]), options, params, proc_config,
                status, Amat, precond, conv_info);

      dtemp=AZ_gdot(Amat->data_org[AZ_N_internal]+Amat->data_org[AZ_N_border],
                    x, x, proc_config);
      dtemp = sqrt(dtemp);
      for (j=0;j<Amat->data_org[AZ_N_internal]+Amat->data_org[AZ_N_border];j++)
        x[j] = x[j]/dtemp;
      if (options[AZ_extreme] == AZ_high) {
        if (dtemp < 2.0) params[AZ_temp] *= 100.;
      }
      else {
        if (dtemp > 1.0) { params[AZ_temp] /= (1.1*pow(dtemp,.1) ); changed++; }
        else if ( changed == 0 ) params[AZ_temp] *= 2.;
        else if (changed < 2) { params[AZ_temp] *= .7; changed += 3; }
      }
      i -= options[AZ_max_iter];
      options[AZ_max_iter] = i;
    }

    /*  find the largest point */
    size1 = Amat->data_org[AZ_N_internal]+Amat->data_org[AZ_N_border];
    largest = -1.;
#ifdef ML
    largest_index = -1;
#endif
    for (i=0;i < size1; i++) {
      if ( fabs(x[i]) > largest ) {
        largest = fabs(x[i]);
#ifdef ML
        largest_index = i;
#endif
      }
    }


#ifdef ML
    global_largest = AZ_gmax_double(largest, proc_config);

    size2 = size1 + Amat->data_org[AZ_N_external];
    tempv = (double *) AZ_allocate(size2*sizeof(double));
    tempy = (double *) AZ_allocate(size2*sizeof(double));
    ibuf =  (int    *) AZ_allocate(size2*sizeof(int   ));
    dbuf =  (double *) AZ_allocate(size2*sizeof(double));
    AZ_exchange_bdry(x, Amat->data_org, proc_config);

    /* Print out the interpolation stencil at that point */

    ml = (ML *) AZ_get_precond_data(precond);
    row_length = 0;
    if (global_largest == largest) {
      allocated = size2/2;
      AZ_printf_out("processor %d (with %d rows) has largest value (%e,%d) in the lambda-vec\n",
             proc_config[AZ_node], size1, largest,largest_index);
      ML_get_matrix_row(&(ml->Pmat[0]), 1, &largest_index, &allocated, &ibuf, &dbuf,
                        &row_length, 0);
      AZ_printf_out("interpolation operator at that point\n");
      for (i = 0; i < row_length; i++) AZ_printf_out("\t(%d, %e)\n",ibuf[i],dbuf[i]);
      AZ_printf_out("Interpolation columns corresponding to above c-points\n");
    }
    AZ_flush_out();

    /* Print out interpolation column for each coarse grid    */
    /* point that appears in the above interpolation stencil. */

    row_length = AZ_gsum_int(row_length, proc_config);
    for (i = 0; i < row_length; i++) {
      for (j = 0; j < ml->Pmat[0].invec_leng; j++) tempv[j] = 0.0;
      for (j = 0; j < size1; j++) tempy[j] = 0.0;
      if (global_largest == largest) tempv[ibuf[i]] = 1.;
      else ibuf[i] = 0;
      ibuf[i] = AZ_gsum_int(ibuf[i], proc_config);

      ML_Operator_Apply(&(ml->Pmat[0]),ml->Pmat[0].invec_leng,tempv,size1,tempy);
      for (j = 0; j < size1; j++)
        if (tempy[j] != 0.0) AZ_printf_out("%d: \t (%d,%d %e)\n",proc_config[AZ_node],ibuf[i],j,tempy[j]);
    }
    AZ_flush_out();

    if (global_largest == largest) {

      /* Print out the matrix stencil at the largest point */

      already_printed = (int *) tempy;
      for (i = 0; i < size2; i++) already_printed[i] = 0;

      AZ_printf_out("neighbors (A, lvec)\n");
      ML_get_matrix_row(&(ml->Amat[1]), 1, &largest_index, &allocated, &ibuf, &dbuf,
                        &row_length, 0);
      for (i = 0; i < row_length; i++) {
        boundary = ' ';
        if ((ibuf[i] < size1) && (Amat->val[ibuf[i]] == 1.0)) boundary = 'b';
        AZ_printf_out("%d: \t%d (%e , %e) %c\n",proc_config[AZ_node],ibuf[i], dbuf[i],
               x[ibuf[i]], boundary);
        already_printed[ibuf[i]] = 1;
      }

      /* Print out the distance two values */

      AZ_printf_out("\ndistance 2 neighbors\n");
      ibuf2 = (int *) tempv;
      for (i = 0; i < row_length; i++) {
        ML_get_matrix_row(&(ml->Amat[1]), 1, &(ibuf[i]), &allocated, &ibuf2, &dbuf,
                          &row2_length, 0);
        for (j = 0; j < row2_length; j++) {
          boundary = ' ';
          if ((ibuf2[i] < size1) && (Amat->val[ibuf2[i]] == 1.0)) boundary = 'b';
          if ( already_printed[ ibuf2[j] ] == 0) {
            AZ_printf_out("%d: \t(%d, %e) %c via %d\n",proc_config[AZ_node],ibuf2[j],
                   x[ibuf2[j]],boundary,i);
            already_printed[ibuf2[i]] = 1;
          }
        }
      }
    }
    AZ_flush_out();
#endif



    break;

  case AZ_GMRESR:

    /* GMRESR */

    AZ_pgmresr(b, x, &(newparams[AZ_weights]), options, params, proc_config,
               status, Amat, precond, conv_info);

    break;

  case AZ_cgs:

    /* conjugate gradient squared */

    AZ_pcgs(b, x, &(newparams[AZ_weights]), options, params, proc_config, status,
            Amat, precond, conv_info);

    break;

  case AZ_tfqmr:

    /* transpose-free quasi minimum residual */

    AZ_pqmrs(b, x, &(newparams[AZ_weights]), options, params, proc_config,status,
             Amat, precond, conv_info);
    break;

  case AZ_bicgstab:

    /* stabilized bi-conjugate gradient */
    AZ_pbicgstab(b, x,&(newparams[AZ_weights]), options, params,proc_config,status,
                 Amat, precond, conv_info);
    break;

  case AZ_symmlq:

    /* conjugate gradient squared */

#ifdef eigen
    AZ_psymmlq(b, x, &(newparams[AZ_weights]), options, params, proc_config, status,
               Amat, precond, conv_info);
#else
    AZ_printf_out("symmlq not implemented in this version\n");
#endif
    break;

  case AZ_lu:
#ifdef HAVE_AZLU

    data_org = Amat->data_org;
    N        = data_org[AZ_N_internal] + data_org[AZ_N_border];
    N_Blk    = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

    /* direct sparse LU */

    for (i = 0; i < N ; i++) x[i] = b[i];

    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) nonzeros = (Amat->bindx)[N];
    else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
      nonzeros = (Amat->indx)[(Amat->bpntr)[N_Blk]];
    else {
      AZ_matfree_Nnzs(Amat);
      nonzeros = Amat->N_nz;
    }


    /* save certain options so that we can restore the user's data */

    itemp2 = options[AZ_precond];
    itemp3 = options[AZ_overlap];
    t1     = options[AZ_subdomain_solve];
    options[AZ_overlap] = AZ_none;
    options[AZ_precond] = AZ_lu;
    options[AZ_subdomain_solve] = AZ_lu;

    /* Set up the subdomain data structure */

    context.A_overlapped           = &Amat2;
    context.aztec_choices          = &aztec_choices;
    context.aztec_choices->options = options;
    context.aztec_choices->params  = newparams;
    context.tag                    = tag;
    sprintf(tag,"lu solve %d",data_org[AZ_name]);

    /* Compute the matrix bandwidth and the number of additional */
    /* nonzeros that might be required by the factorization      */

    i = AZ_compute_max_nz_per_row(Amat, N, N_Blk, &bandwidth);

    AZ_space_for_factors(0.0, nonzeros, N, &extra_factor_nonzeros,options,
                         bandwidth,0);

    /* Adjust the space for the factors based */
    /* on how much memory we can allocate.    */

    N_nz_factors = nonzeros + extra_factor_nonzeros;
    AZ_hold_space(&context,N);
    N_nz_factors = AZ_adjust_N_nz_to_fit_memory(N_nz_factors,2,1);
    AZ_free_space_holder(&context);

    /* Allocate and fill the factor arrays with MSR matrix info */

    bindx2 = (int *) AZ_allocate(N_nz_factors*sizeof(int));
    val2   = (double *) AZ_allocate(N_nz_factors*sizeof(double));
    if (val2 == NULL) AZ_perror("Not enough space for factorization\n");
    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
      for (i = 0 ; i < nonzeros ; i++ ) {
        bindx2[i] = (Amat->bindx)[i];
        val2[i] = (Amat->val)[i];
      }
    }
    else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
      AZ_vb2msr(N_Blk,Amat->val,Amat->indx,Amat->bindx,Amat->rpntr,Amat->cpntr,
                Amat->bpntr, val2, bindx2);
    }
    else AZ_matfree_2_msr(Amat,val2,bindx2,nonzeros);

    Amat2.bindx = bindx2;
    Amat2.val   = val2;
    Amat2.data_org = data_org;

    AZ_factor_subdomain(&context, N, N_nz_factors, &nz_used);

    if ((options[AZ_output] != AZ_none ) && (options[AZ_output] != AZ_warnings)){
      AZ_printf_out("\n********************************************************************\n");
      AZ_printf_out("*****  Condition number estimate for preconditioner = %.4e\n",
             AZ_condest(N, &context));
      AZ_printf_out("********************************************************************\n");
    }

    AZ_solve_subdomain(x, N, &context);
    AZ_free(val2);
    AZ_free(bindx2);
    options[AZ_precond] = itemp2;
    options[AZ_overlap] = itemp3;
    options[AZ_subdomain_solve] = t1;

    status[AZ_its]      = (double ) 1.0;
    status[AZ_why]      = (double ) AZ_normal;
    status[AZ_r]        = (double ) 0.0;
    status[AZ_rec_r]    = (double ) 0.0;
    status[AZ_scaled_r] = (double ) 0.0;
#else
    AZ_printf_err("AZ_lu unavailable: configure with --enable-aztecoo-azlu to make available\n");
    exit(1);
#endif
    break;

  default:
    (void) AZ_printf_err("ERROR: options[AZ_solver] has improper value(%d)\n",
                   options[AZ_solver]);
    exit(-1);
  }
  AZ_flush_out();

  /* Must delete conv_info if we created it */
  if (options[AZ_conv]!=AZTECOO_conv_test)
    AZ_converge_destroy(&conv_info);

  AZ_oldsolve_finish(x, b, options, proc_config, Amat, save_old_values,scaling);
  options[AZ_ignore_scaling] = save_old_values[6];
} /* AZ_oldsolve */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_print_call_iter_solve(int options[], double params[], int az_proc,
                              int recur, AZ_MATRIX * Amat, AZ_PRECOND *precond)

     /*******************************************************************************

  Print out the type of solver called, scaling and preconditioning information.

  Author:          SNL
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  az_proc:         Current processor.

  recur:           Indicates the level of recursion that this call corresponds to.

*******************************************************************************/

{
  char prefix[40];
  int i, super_special = 0, str_leng = 0;

  /**************************** execution begins ******************************/

  if ( (options[AZ_output] == AZ_last) || (options[AZ_output] == AZ_none) ||
       (options[AZ_output] == AZ_warnings) ||  (az_proc != 0) ) return;

  prefix[str_leng++] = '\t'; prefix[str_leng++] = '\t'; prefix[str_leng++] = '*';
  prefix[str_leng++] = '*'; prefix[str_leng++] = '*'; prefix[str_leng++] = '*';
  prefix[str_leng++] = '*'; prefix[str_leng++] = ' ';
  for (i = 0 ; i < recur ; i++ ) {
    prefix[str_leng++] = ' ';
    prefix[str_leng++] = ' ';
  }
  prefix[str_leng++] = '\0';

  if ( recur == 0 )
    (void) AZ_printf_out("\n\t\t****************************************"
                  "***************\n");

  /* First print problem description (if available) */


  if (Amat != NULL) {
    if (Amat->print_string != NULL) {
      (void) AZ_printf_out(prefix);
      (void) AZ_printf_out("Problem: ");
      (void) AZ_printf_out("%s\n",Amat->print_string);
    }
  }

  (void) AZ_printf_out(prefix);

  /* next, print out chosen solver */

  switch (options[AZ_solver]) {

  case AZ_cg:
    (void) AZ_printf_out("Preconditioned CG"); break;
  case AZ_cg_condnum:
    (void) AZ_printf_out("Preconditioned CG (with condnum)"); break;
  case AZ_gmres:
    (void) AZ_printf_out("Preconditioned GMRES"); break;
  case AZ_gmres_condnum:
    (void) AZ_printf_out("Preconditioned GMRES (with condnum)"); break;
  case AZ_analyze:
    (void) AZ_printf_out("Preconditioned analysis"); break;
  case AZ_GMRESR:
    (void) AZ_printf_out("Preconditioned GMRESR"); break;
  case AZ_fixed_pt:
    (void) AZ_printf_out("Preconditioned fixed-point iter."); break;
  case AZ_cgs:
    (void) AZ_printf_out("Preconditioned CGS"); break;
  case AZ_tfqmr:
    (void) AZ_printf_out("Preconditioned TFQMR"); break;
  case AZ_bicgstab:
    (void) AZ_printf_out("Preconditioned BICGSTAB"); break;
  case AZ_symmlq:
    (void) AZ_printf_out("Preconditioned SYMMLQ-like"); break;
  case AZ_lu:
    (void) AZ_printf_out("LU");
  }

  (void) AZ_printf_out(" solution\n");

  /* next output preconditioning options */

  (void) AZ_printf_out(prefix);

  if ((precond != NULL) && (precond->prec_function != AZ_precondition)) {
    if (precond->print_string == NULL) (void) AZ_printf_out("user ");
    else AZ_printf_out("%s",precond->print_string);
  }
  else {
    switch (options[AZ_precond]) {
    case AZ_none:
      (void) AZ_printf_out("No preconditioning"); break;
    case AZ_Neumann:
      (void) AZ_printf_out("Order %d Neumann series polynomial", options[AZ_poly_ord]);
      break;
    case AZ_ls:
      (void) AZ_printf_out("Order %d least-squares polynomial", options[AZ_poly_ord]);
      break;
    case AZ_Jacobi:
      (void) AZ_printf_out("%d step block Jacobi", options[AZ_poly_ord]);
      break;
    case AZ_smoother:
      (void) AZ_printf_out("%d step loc avg smoother", options[AZ_poly_ord]);
      break;
    case AZ_sym_GS:
      (void) AZ_printf_out("%d step symmetric Gauss-Seidel", options[AZ_poly_ord]);
      break;
    case AZ_dom_decomp:
      if (options[AZ_subdomain_solve] == AZ_bilu)
        AZ_printf_out("BILU(%d) domain decomp. with", options[AZ_graph_fill]);
      /* Begin Aztec 2.1 mheroux mod */
      else if (options[AZ_subdomain_solve] == AZ_bilu_ifp) {
        AZ_printf_out("IFPACK BILU(%d) ( ATHRESH = %.3e, RTHRESH = %.3e)\n ",
               options[AZ_graph_fill],params[AZ_athresh], params[AZ_rthresh]);
        AZ_printf_out(prefix); AZ_printf_out("with");
      }
      /* End Aztec 2.1 mheroux mod */
      else if (options[AZ_subdomain_solve] == AZ_ilut) {
        AZ_printf_out("ILUT( fill-in = %.3e, drop = %.3e)\n ",
               params[AZ_ilut_fill], params[AZ_drop]);
        AZ_printf_out(prefix); AZ_printf_out("with");
      }
      else if (options[AZ_subdomain_solve] == AZ_ilu)
        AZ_printf_out("ILU(%d) domain decomp. with", options[AZ_graph_fill]);
      else if (options[AZ_subdomain_solve] == AZ_rilu)
        AZ_printf_out("RILU(%d,%.2f) domain decomp. with",options[AZ_graph_fill],
               params[AZ_omega]);
      else if (options[AZ_subdomain_solve] == AZ_lu)
        AZ_printf_out("LU domain decomp. with");
      else if (options[AZ_subdomain_solve] == AZ_icc)
        AZ_printf_out("icc(%d) domain decomp. with",options[AZ_graph_fill]);
      else if (options[AZ_subdomain_solve] < AZ_SOLVER_PARAMS)
        (void) AZ_printf_out("iterative subdomain solve with");
      else {
        (void) AZ_printf_out("Unknown subdomain solver (%d)\n",
                      options[AZ_subdomain_solve]);
        exit(1);
      }
      if      (options[AZ_overlap] == AZ_none) AZ_printf_out("out overlap");
      else if (options[AZ_overlap] == AZ_diag) AZ_printf_out(" diagonal overlap");
      else if (options[AZ_type_overlap] == AZ_symmetric) AZ_printf_out(" symmetric");
      if (options[AZ_overlap] >= 1) AZ_printf_out(" overlap = %d ",options[AZ_overlap]);
      super_special = 1;
      break;
    case AZ_icc:
      (void) AZ_printf_out("incomplete Choleski decomposition");
      super_special = 1;
      break;
    case AZ_user_precond:
      (void) AZ_printf_out("user ");
    default:
      if (options[AZ_precond] < AZ_SOLVER_PARAMS)
        (void) AZ_printf_out("iterative preconditioner");
    }
  }

  (void) AZ_printf_out("\n");

  (void) AZ_printf_out(prefix);

  /* lastly, print out the scaling information */

  switch (options[AZ_scaling]) {

  case AZ_none:
    (void) AZ_printf_out("No"); break;
  case AZ_Jacobi:
    (void) AZ_printf_out("block Jacobi"); break;
  case AZ_BJacobi:
    (void) AZ_printf_out("block Jacobi"); break;
  case AZ_row_sum:
    (void) AZ_printf_out("left row-sum"); break;
  case AZ_sym_diag:
    (void) AZ_printf_out("symmetric diagonal"); break;
  case AZ_sym_row_sum:
    (void) AZ_printf_out("symmetric row sum"); break;
  case AZ_equil:
    (void) AZ_printf_out("equilibrated");
  }

  (void) AZ_printf_out(" scaling\n");

  if (super_special==1) {
    (void) AZ_printf_out("%sNOTE: convergence VARIES when the total number "
                  "of\n",prefix);
    (void) AZ_printf_out("%s      processors is changed.\n",prefix);
  }

  if ( recur == 0 )
    (void) AZ_printf_out("\t\t****************************************"
                  "***************\n");

} /* AZ_print_call_iter_solver */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_output_matrix(double val[], int indx[], int bindx[], int rpntr[],
                      int cpntr[], int bpntr[], int proc_config[],
                      int data_org[])

     /*******************************************************************************

  Routine to perform full matrix dump on each processor.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  a:               Array containing the nonzero entries of the matrix.

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  int  iblk_row, i, j, ib1, ib2, n1, iblk, jblk, m1, ipoint, jpoint;
  int  ival = 0;
  int  k,num_nonzeros;
  int  num_total_nodes, N_external_nodes;
  int  Proc;
  char str[5];
  char nstr[40];

  /********** execution begins **********/

  Proc               = proc_config[AZ_node];
  N_external_nodes = data_org[AZ_N_external];

  if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
    num_total_nodes    = data_org[AZ_N_int_blk]+data_org[AZ_N_bord_blk];
    N_external_nodes = data_org[AZ_N_ext_blk];

    /* Print out the VBR indexing information for the matrix */

    AZ_print_sync_start(Proc, AZ_TRUE, proc_config);

    (void) AZ_printf_out("\n----- Proc: %d indx -----\n\n", Proc);
    for (iblk = 0; iblk < num_total_nodes; iblk++) {
      for (i = *(bpntr+iblk); i < *(bpntr+iblk+1); i++)
        (void) AZ_printf_out("%d ", *(indx+i));

      if (iblk == num_total_nodes - 1)
        (void) AZ_printf_out("%d\n",*(indx+i));
      else
        (void) AZ_printf_out("\n");
    }

    (void) AZ_printf_out("\n----- Proc: %d bindx -----\n\n", Proc);
    for (iblk = 0; iblk < num_total_nodes; iblk++) {
      for (i = *(bpntr+iblk); i < *(bpntr+iblk+1); i++)
        (void) AZ_printf_out("%d ", *(bindx+i));
      (void) AZ_printf_out("\n");
    }

    (void) AZ_printf_out("\n----- Proc: %d rpntr -----\n\n", Proc);
    for (i = 0; i < num_total_nodes + 1; i++)
      (void) AZ_printf_out("%d ", *(rpntr+i));
    (void) AZ_printf_out("\n");

    (void) AZ_printf_out("\n----- Proc: %d cpntr -----\n\n", Proc);
    for (i = 0; i < num_total_nodes + N_external_nodes + 1; i++)
      (void) AZ_printf_out("%d ", *(cpntr+i));
    (void) AZ_printf_out("\n");

    (void) AZ_printf_out("\n----- Proc: %d bpntr -----\n\n", Proc);
    for (i = 0; i < num_total_nodes + 1; i++)
      (void) AZ_printf_out("%d ", *(bpntr+i));
    (void) AZ_printf_out("\n");

    AZ_print_sync_end(proc_config, AZ_TRUE);

    AZ_print_sync_start(Proc, AZ_TRUE, proc_config);
    (void) AZ_printf_out("AZ_solve debug output - full matrix dump: Processor %d\n",
                  Proc);

    /* loop over block rows */

    for (iblk_row = 0; iblk_row < num_total_nodes; iblk_row++) {

      /* number of rows in the current row block */

      m1 = rpntr[iblk_row+1] - rpntr[iblk_row];

      /* starting index of current row block */

      ival = indx[bpntr[iblk_row]];

      /* loop over all the blocks in the current block-row */

      for (j = bpntr[iblk_row]; j < bpntr[iblk_row+1]; j++) {
        jblk = bindx[j];

        /* the starting point column index of the current block */

        ib1 = cpntr[jblk];

        /* ending point column index of the current block */

        ib2 = cpntr[jblk+1];

        /* number of columns in the current block */

        n1 = ib2 - ib1;

        (void) AZ_printf_out("\nProc: %d Block Row: %d Block Column: %d "
                      "Row Pointer: %d Column Pointer: %d\n", Proc, iblk_row,
                      jblk, rpntr[iblk_row], rpntr[jblk]);
        (void) AZ_printf_out("----------------------------------------"
                      "----------------------------------------\n");

        for (ipoint = 0; ipoint < m1; ipoint++) {
          for (jpoint = 0; jpoint < n1; jpoint++)
            (void) AZ_printf_out("a[%d]: %e ", ival+jpoint*m1+ipoint,
                          val[ival+jpoint*m1+ipoint]);
          (void) AZ_printf_out("\n");
        }

        ival += m1*n1;
      }
    }

    AZ_print_sync_end(proc_config, AZ_TRUE);
  }

  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
    num_total_nodes = data_org[AZ_N_internal]+data_org[AZ_N_border];

    N_external_nodes = data_org[AZ_N_external];
    i                  = num_total_nodes + N_external_nodes;

    if      (i < 10    ) sprintf(str, "%%1d");
    else if (i < 100   ) sprintf(str, "%%2d");
    else if (i < 1000  ) sprintf(str, "%%3d");
    else if (i < 10000 ) sprintf(str, "%%4d");
    else if (i < 100000) sprintf(str," %%5d");
    else sprintf(str, "%%d");
    sprintf(nstr, "a(,%s)=%%8.1e ", str);

    AZ_print_sync_start(Proc, AZ_TRUE, proc_config);

    (void) AZ_printf_out("\n----- Proc: %d -----\n\n", Proc);

    num_nonzeros = bindx[num_total_nodes];
    (void) AZ_printf_out("val:  ");
    for (i = 0; i < num_nonzeros; i++) {
      (void) AZ_printf_out("%9.1e", val[i]);
      if ((i%8) == 7) (void) AZ_printf_out("\n    ");
    }

    (void) AZ_printf_out("\nbindx:");
    for (i = 0; i < num_nonzeros; i++) {
      (void) AZ_printf_out("%9d", bindx[i]);
      if ((i%8) == 7) (void) AZ_printf_out("\n    ");
    }
    (void) AZ_printf_out("\n");

    for (i = 0; i < num_total_nodes; i++ ) {
      (void) AZ_printf_out("\nrow");
      (void) AZ_printf_out(str, i);
      (void) AZ_printf_out(":");
      (void) AZ_printf_out(nstr, i, val[i]);
      k = 0;
      for (j = bindx[i]; j < bindx[i+1]; j++ ) {
        (void) AZ_printf_out(nstr, bindx[j], val[j]);
        k++;
        if (((k%4) == 3) && (j != bindx[i+1]-1))
          (void) AZ_printf_out("\n      ");
      }
    }

    (void) AZ_printf_out("\n");
    AZ_print_sync_end( proc_config, AZ_TRUE);
  }

} /* AZ_output_matrix */

void AZ_flop_rates(int data_org[],int indx[],int bpntr[], int bindx[],
                   int options[], double status[], double total_time,
                   int proc_config[])
{
  int N, N_Blk, gn;
  double MFlops, gnnz;
  /*
   * calculate the solver MFlop/s
   */

#ifdef AZ_FLOP_CNTS
#endif

  N_Blk        = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
  N            = data_org[AZ_N_internal] + data_org[AZ_N_border];
  if ( (options[AZ_output] != AZ_none) && (options[AZ_output] != AZ_warnings) &&
       (data_org[AZ_matrix_type] != AZ_USER_MATRIX)) {
    gn = AZ_gsum_int(N, proc_config);
    if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
      gnnz = AZ_gsum_double((double) (indx[bpntr[N_Blk]]-indx[0]),
                            proc_config);
    else if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX)
      gnnz = AZ_gsum_double((double) (bindx[N]), proc_config);
    else return;

    if (proc_config[AZ_node] == 0) {
      total_time += 1.0e-10;
      MFlops = AZ_calc_solve_flops(options, (int) status[AZ_its],total_time,gn,
                                   gnnz, data_org, proc_config);
      if (MFlops > 0.0) {
        (void) AZ_printf_out("\t\tSolver MFlop rate: %f MFlops/sec.\n\t\t", MFlops);
        (void) AZ_printf_out("Solver processor MFlop rate: %f MFlops/sec.\n\n",
                      MFlops/proc_config[AZ_N_procs]);
      }
    }
  }

#ifdef TIME_VB

  /* calculate the individual kernel times and Mflops */

  AZ_time_kernals(gn, gnnz, val, indx, bindx, rpntr,
                  rpntr, bpntr, x, b, (int) status[AZ_NUMBER_OF_ITS],
                  options, data_org, proc_config);

#endif
}
/**************************************************************/
/**************************************************************/
/**************************************************************/
void AZ_mk_identifier(double *params, int *options,
                      int *data_org, char *tag) {
  /*
 * Make a unique identifier that is linked to preconditioner/
 * solver/scaling/problem size options that the user has set.
 */

  double dtemp;
  int    itemp;
  char   code1,code2,code3;


  dtemp = (params[AZ_ilut_fill] + 3.1415)*(params[AZ_drop] + 2.712)*
    (1.1 + ((double) options[AZ_graph_fill]));
  if (dtemp > 0) dtemp = pow(dtemp,.01);


  itemp = 4*(options[AZ_overlap] + 1) + 2*options[AZ_reorder] +
    options[AZ_type_overlap];

  if (itemp < 94) code1 = '!' + itemp;
  else code1 = itemp%255;

  if ( options[AZ_precond] == AZ_dom_decomp )
    itemp = options[AZ_subdomain_solve];
  else itemp = options[AZ_precond];

  if (itemp < 94) code2 = '!' + itemp;
  else code2 = itemp%255;


  itemp = options[AZ_scaling];
  if (itemp < 94) code3 = '!' + itemp;
  else code3 = itemp%255;

  sprintf(tag,"P%d %c%c%c %.4f", data_org[AZ_N_internal]+data_org[AZ_N_border],
          code1, code2,code3, dtemp);
}

/************************************************************************/
/************************************************************************/
/************************************************************************/

void AZ_mk_context(int options[], double params[], int data_org[],
                   AZ_PRECOND *precond, int proc_config[])
{
  /* Make a context to be associated with the preconditioning data structure
   * and stick it in the field: precond->context. This context should be
   * uniquely based on the name of the matrix which the preconditioner operates
   * on as well as the type of preconditioning options that have been chosen.
   *
   * Note: if this context already exists, this routine will simply set
   * precond->context to point to it.
   *
   *
   */
  char tag[80];
  int  istatus;

  AZ_mk_identifier(params,options,data_org, tag);

  precond->context = (struct context *) AZ_manage_memory(sizeof(struct context),
                                                         AZ_ALLOC,
                                                         data_org[AZ_name],
                                                         tag,&istatus);

  if (istatus == AZ_NEW_ADDRESS) {
    AZ_zero_out_context(precond->context);
    if ((options[AZ_pre_calc] == AZ_reuse) && (proc_config[AZ_node] == 0)){
      AZ_printf_err("Error:\tDid not find previous factorization (");
      AZ_printf_err( "requested \n\tby setting options[AZ_pre_calc] to ");
      AZ_printf_err( "AZ_reuse).\n\tTo find this factorization, the ");
      AZ_printf_err( "following\n\tparameters must match the previous");
      AZ_printf_err(" factorization:");
      AZ_printf_err( "\n\t\t 1) Total number of unknowns.");
      AZ_printf_err( "\n\t\t 2) options[AZ_overlap]");
      AZ_printf_err( "\n\t\t 3) options[AZ_scaling]");
      AZ_printf_err( "\n\t\t 4) options[AZ_precond]");
      AZ_printf_err( "\n\t\t 5) options[AZ_reorder]");
      AZ_printf_err( "\n\t\t 6) options[AZ_type_overlap]");
      AZ_printf_err( "\n\t\t 7) options[AZ_subdomain_solve]");
      AZ_printf_err( "\n\t\t 8) options[AZ_graph_fill]");
      AZ_printf_err( "\n\t\t 9) params[AZ_ilut_fill]");
      AZ_printf_err( "\n\t\t10) params[AZ_drop]");
      AZ_printf_err( "\n\t\t11) data_org[AZ_name]\n");
      AZ_printf_out("XXX%sXXX %d %d\n",tag,data_org[AZ_name],(int) sizeof(struct context));
      (void) AZ_manage_memory(0, -43, AZ_SYS+az_iterate_id, (char *) 0, (int *) 0);
    }
    if (options[AZ_pre_calc] == AZ_reuse) exit(1);

    tag[0] = 'T';
    precond->context->tag = (char *) AZ_manage_memory(sizeof(char)*80,
                                                      AZ_ALLOC,
                                                      data_org[AZ_name],
                                                      tag,&istatus);
    tag[0] = 'P';
    sprintf(precond->context->tag,"%s",tag);
  }
}
/************************************************************************/
/************************************************************************/
/************************************************************************/

void AZ_rm_context(int options[], double params[], int data_org[])
{
  /* Remove context associated with the preconditioning data structure.
   * See 'AZ_mk_context' for more details.
   *
   *
   */
  char tag[80];
  int  istatus;


  AZ_mk_identifier(params,options,data_org, tag);

  tag[0] = 'T';
  AZ_manage_memory(sizeof(char)*80, AZ_SELECTIVE_CLEAR,
                   data_org[AZ_name], tag,&istatus);
  tag[0] = 'P';
  AZ_manage_memory(sizeof(struct context), AZ_SELECTIVE_CLEAR,
                   data_org[AZ_name],tag,&istatus);
}

int AZ_oldsolve_setup(double x[], double b[], int options[], double params[],
                      double status[], int proc_config[], AZ_MATRIX *Amat,
                      AZ_PRECOND *precond, int save_old_values[],
                      struct AZ_SCALING *scaling)


     /*******************************************************************************
  Aztec's previous AZ_solve() is renamed to AZ_oldsolve() with 3 new
  parameters appended to it: Amat, precond, scaling.
  This routine is never called directly by an application. It is only
  used internally by Aztec.

  Solve the system of equations given in the VBR format using an iterative
  method specified by 'options[AZ_solver]'. Store the result in 'x'.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  b:               Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  val:             Array containing the nonzero entries of the matrix (see
                   Aztec User's Guide).

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  Amat:            Structure used to represent the matrix (see az_aztec.h
                   and Aztec User's Guide).

  precond:         Structure used to represent the preconditionner
                   (see az_aztec.h ad Aztec User's Guide).

  scaling:         Structure used to represent the  scaling
                   (see az_aztec.h ad Aztec User's Guide).


*******************************************************************************/

{

  /* local variables */

  int     i;
  int     Proc, Num_Proc;
  int     *data_org;
  int     prec_context, *toptions;
  int *tdata_org;
  AZ_PRECOND *tprecond;
  double *tparams;
  char tag[80];


  /**************************** execution begins ******************************/

  if (options[AZ_recursion_level] == 0) precond->timing[0] = 0.0;
#ifdef AZTEC_MPI
  if ( proc_config[AZ_Comm_Set] != AZ_Done_by_User) {
    AZ_printf_out("Error: Communicator not set. Use AZ_set_comm()\n");
    AZ_printf_out("       (e.g. AZ_set_comm(proc_config,MPI_COMM_WORLD)).\n");
    exit(1);
  }
#endif

  /* save a few variables that Aztec will mess with so that */
  /* we can restore them before returning to the user.      */


  save_old_values[AZ_oldprecond] = options[AZ_precond];
  save_old_values[AZ_oldsubdomain_solve] = options[AZ_subdomain_solve];
  save_old_values[AZ_oldreorder] = options[AZ_reorder];
  save_old_values[AZ_oldpre_calc] = options[AZ_pre_calc];
  i  = 1;
  if (options[AZ_precond] == AZ_lu) sprintf(tag,"AZ_lu");
  else if (options[AZ_precond] == AZ_ilu ) sprintf(tag,"AZ_ilu");
  else if (options[AZ_precond] == AZ_ilut) sprintf(tag,"AZ_ilut");
  else if (options[AZ_precond] == AZ_bilu) sprintf(tag,"AZ_bilu");
  /* Begin Aztec 2.1 mheroux mod */
  else if (options[AZ_precond] == AZ_bilu_ifp) sprintf(tag,"AZ_bilu_ifp");
  /* End Aztec 2.1 mheroux mod */
  else if (options[AZ_precond] == AZ_rilu) sprintf(tag,"AZ_rilu");
  else i = 0;
  if (i == 1) {
    AZ_printf_out("To use this preconditioner, you must now set\n");
    AZ_printf_out("      options[AZ_precond] = AZ_dom_decomp;\n");
    AZ_printf_out("      options[AZ_subdomain_solve] = %s;\n",tag);
    exit(1);
  }


  data_org = Amat->data_org;

  Proc         = proc_config[AZ_node];
  Num_Proc     = proc_config[AZ_N_procs];

  data_org[AZ_internal_use] = 1;

  /* check for inconsistent user input options */

  if (!AZ_check_options(options, Proc, data_org, Num_Proc, params,
                        Amat, precond)) {

    status[AZ_its]      = (double )  0.0;
    status[AZ_why]      = (double )  AZ_param;
    status[AZ_r]        = (double ) -1.0;
    status[AZ_rec_r]    = (double ) -1.0;
    status[AZ_scaled_r] = (double ) -1.0;
    return(0);
  }

#ifdef AZ_ENABLE_CAPTURE_MATRIX
  /* Test to see if we should capture matrix, rhs and partitioning info
     in an ASCII data file.  If the file "AZ_write_matrix_now" exists in
     the current working directory, then the files
     - AZ_capture_matrix.dat
     - AZ_capture_rhs.dat
     - AZ_capture_partition.dat (VBR only)

     will be appended with the current matrix in (i,j,val) format, the
     current RHS and the current partition information.  The existence
     of "AZ_write_matrix_now" is check each time.  Thus, capturing can
     be turned on and off at will during the run of a simulation.
  */

  if (options[AZ_output] != AZ_none) {
    AZ_capture_matrix( Amat, proc_config, data_org, b);
  }
#endif


  /* If desired, print out the matrix and indexing arrays */

  if ((data_org[AZ_matrix_type] == AZ_MSR_MATRIX) ||
      (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)){
    if (options[AZ_output] == AZ_all)
      AZ_print_out((int *) NULL, (int *) NULL, (int *) NULL, (int *) NULL,
                   Amat->val, Amat->indx, Amat->bindx, Amat->rpntr,
                   Amat->cpntr, Amat->bpntr, proc_config, AZ_input_form,
                   data_org[AZ_matrix_type], data_org[AZ_N_int_blk] +
                   data_org[AZ_N_bord_blk],data_org[AZ_N_ext_blk],0);
  }

  toptions  = options;
  tparams   = params;
  tprecond  = precond;
  tdata_org = precond->Pmat->data_org;

  prec_context = 1;
  if (options[AZ_precond] == AZ_multilevel)
    prec_context = AZ_multilevel;

  while ( prec_context ) {

    /* adjust print frequency according to user input */

    i = toptions[AZ_output];

    if  (i == AZ_warnings) toptions[AZ_print_freq] =toptions[AZ_max_iter] + 10;
    else if (i == AZ_none) toptions[AZ_print_freq] =toptions[AZ_max_iter] + 10;
    else if (i == AZ_all ) toptions[AZ_print_freq] = 1;
    else if (i == AZ_last || i == AZ_summary) toptions[AZ_print_freq] =toptions[AZ_max_iter] + 1;
    else                   toptions[AZ_print_freq] =toptions[AZ_output];
    if ((i != AZ_none) && (i != AZ_warnings) &&
        (proc_config[AZ_node] == 0)) AZ_printf_out("\n");

    /*  Setup data structure to record domain decomposition factorization */
    /*  information that can be used in future solves.                    */
    /*  Note: If this information already exists, we recover the previous */
    /*        information.                                               */

    AZ_mk_context(toptions, tparams, tdata_org, tprecond, proc_config);

    /* check if multiple preconditioners are being passed in */
    /* if this is the case we must afix a subdomain context  */
    /* to each one.                                          */

    if ( prec_context == AZ_multilevel) {
      if (tprecond->next_prec == NULL) prec_context = 0;
      else {
        tprecond = tprecond->next_prec;
        toptions = tprecond->options;
        tparams  = tprecond->params;
        tdata_org= tprecond->Pmat->data_org;
      }
    }
    else prec_context = 0;
  }

  /* scale matrix, rhs and initial guess if required */

  AZ_scale_f(AZ_SCALE_MAT_RHS_SOL, Amat, options, b, x, proc_config, scaling);
  return(1);
}



void AZ_oldsolve_finish(double x[], double b[], int options[],
                        int proc_config[], AZ_MATRIX *Amat, int save_old_values[],
                        struct AZ_SCALING *scaling)
{


  AZ_scale_f(AZ_INVSCALE_SOL, Amat, options, b, x, proc_config, scaling);
  AZ_scale_f(AZ_INVSCALE_RHS, Amat, options, b, x, proc_config, scaling);

  options[AZ_precond]         = save_old_values[AZ_oldprecond];
  options[AZ_subdomain_solve] = save_old_values[AZ_oldsubdomain_solve];
  options[AZ_reorder]         = save_old_values[AZ_oldreorder];
  options[AZ_pre_calc]        = save_old_values[AZ_oldpre_calc];

  AZ_flush_out();

} /* AZ_oldsolve_finish  */

void AZ_iterate_setup(int options[], double params[], int proc_config[],
                      AZ_MATRIX *Amat, AZ_PRECOND *precond)

     /*******************************************************************************
This is the new Aztec interface. This routine calls AZ_oldsolve() passing
in for example Amat->indx for the indx[] parameter in AZ_oldsolve().

NOTE: User's can still invoke AZ_solve() in the old Aztec way. AZ_solve
      also calls AZ_oldsolve(). However, matrix-free and coarse grid
       capabilities are not available via AZ_solve().


  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============


  x:               Contains the result vector x.

  b:               Contains the vector b.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  status:          On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  Amat:            Structure used to represent the matrix (see az_aztec.h
                   and Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  precond:         Structure used to represent the preconditionner
                   (see az_aztec.h ad Aztec User's Guide).

  scaling:         Structure used to represent the  scaling
                   (see az_aztec.h ad Aztec User's Guide).

******************************************************************************/

{
  if (Amat->matrix_type == AZ_MSR_MATRIX) {
    Amat->matvec = AZ_MSR_matvec_mult;
    Amat->data_org[AZ_N_int_blk] = Amat->data_org[AZ_N_internal];
    Amat->data_org[AZ_N_bord_blk] = Amat->data_org[AZ_N_border];
  }

  else if (Amat->matrix_type == AZ_VBR_MATRIX)
    Amat->matvec = AZ_VBR_matvec_mult;

  else if (Amat->matrix_type == AZ_USER_MATRIX){
    Amat->data_org[AZ_N_int_blk] = Amat->data_org[AZ_N_internal];
    Amat->data_org[AZ_N_bord_blk] = Amat->data_org[AZ_N_border];
  }


  Amat->data_org[AZ_matrix_type] = Amat->matrix_type;

  /*
    if (options[AZ_precond] != AZ_user_precond)
    precond->prec_function = AZ_precondition;
  */

  options[AZ_recursion_level] = 0;

  if (options[AZ_pre_calc] != AZ_reuse) {
    (void) AZ_manage_memory(0,AZ_EVERYBODY_BUT_CLEAR,
                            (Amat->data_org)[AZ_name],
                            "kvecs",
                            (int *) 0);
    (void) AZ_manage_memory(0,AZ_EVERYBODY_BUT_CLEAR,
                           (precond->Pmat->data_org)[AZ_name],
                            "kvecs", (int *) 0);
  }
  (void) AZ_manage_memory(0, AZ_CLEAR, AZ_SYS+az_iterate_id,
                          (char *) 0, (int *) 0);

  /* output solver, scaling, and preconditioning options */

  AZ_print_call_iter_solve(options, params, proc_config[AZ_node], 0, Amat, precond);
}

void AZ_iterate_finish(int options[], AZ_MATRIX *Amat, AZ_PRECOND *precond)
{
  if (options[AZ_keep_info] == 0) {

#ifdef IFPACK
    /* Delete IFPACK object */
    if (options[AZ_subdomain_solve] == AZ_bilu_ifp)
      ifp_freebiluk(precond->context->precon);

#endif
    (void) AZ_manage_memory(0,AZ_CLEAR,
                            (Amat->data_org)[AZ_name],
                            (char *) 0,(int *) 0);
    (void) AZ_manage_memory(0,AZ_CLEAR,
                           (precond->Pmat->data_org)[AZ_name],
                            (char *) 0,(int *) 0);
  }

  (void) AZ_manage_memory(0, AZ_CLEAR, AZ_SYS+az_iterate_id, (char *) 0, (int *) 0);
}

int AZ_initialize(double x[], double b[], int options[],
                  double params[], double status[], int proc_config[], AZ_MATRIX *Amat,
                  AZ_PRECOND *precond, int save_old_values[], struct AZ_SCALING *scaling)
{
  AZ_iterate_setup(options, params, proc_config, Amat, precond);
  return(AZ_oldsolve_setup(x, b, options, params, status, proc_config, Amat,
                           precond, save_old_values, scaling));

}

void AZ_finalize(double x[], double b[], int options[], int
                 proc_config[], AZ_MATRIX *Amat, AZ_PRECOND *precond, int save_old_values[],
                 struct AZ_SCALING *scaling)
{
  AZ_oldsolve_finish(x, b, options, proc_config, Amat, save_old_values,scaling);
  AZ_iterate_finish(options, Amat, precond);
}
