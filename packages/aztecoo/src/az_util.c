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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "az_aztec.h"
#include "az_blas_wrappers.h"

/*
 * File containing utility functions for solvers.  Note: Some of the fem high
 * level solvers such as the nonlinear solver call these routines as well.
 */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_compute_residual(double b[], double x[], double r[], 
			 int proc_config[], AZ_MATRIX *Amat)

/*******************************************************************************

  Compute the residual r = b - Ax.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val:             Array containing the nonzero entries of the matrix (see
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  b:               Right hand side of linear system.

  x:               On input, contains the initial guess. On output contains the
                   solution to the linear system.

  r:               On output, residual vector.

  options:         Determines specific solution method and other parameters.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  params:          Drop tolerance and convergence tolerance info.

  Amat:            Structure used to represent the matrix (see az_aztec.h
                   and Aztec User's Guide).



*******************************************************************************/

{

  /* local variables */

  register int i;
  int N;

  /**************************** execution begins ******************************/

  N = Amat->data_org[AZ_N_internal] + Amat->data_org[AZ_N_border];

  Amat->matvec(x, r, Amat, proc_config);


  for(i = 0; i < N; i++) r[i] = b[i] - r[i];

} /* AZ_compute_residual */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gmax_vec(int N, double vec[], int proc_config[])

/*******************************************************************************

  Routine to return the maximum element of a distributed vector "vec".

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     double, maximum value in vector 'vec'
  ============

  Parameter list:
  ===============

  N:               Length of vector 'vec'.

  vec:             Vector of length 'N'.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  register int i;
  double       rmax = 0.0;

  /**************************** execution begins ******************************/

  for (i = 0; i < N; i++) rmax = AZ_MAX(rmax, vec[i]);
  rmax = AZ_gmax_double(rmax, proc_config);

  return rmax;

} /* AZ_gmax_vec */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gdot(int N, double r[], double z[], int proc_config[])

/*******************************************************************************

  Routine to perform dot product of r and z with unit stride. This routine call
  the BLAS routine ddot to do the local vector dot product and then uses the
  global summation routine AZ_gsum_double to obtain the reguired global result.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     double, dot product of vectors 'r' and 'z'
  ============

  Parameter list:
  ===============

  N:               Length of vector 'vec'.

  r, z:            Vectors of length 'N'.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_Procs] is the number of processors.

*******************************************************************************/

{

  static int one = 1;
  int        add_N;

  add_N = N;

  return AZ_gsum_double(DDOT_F77(&add_N, r, &one, z, &one), proc_config);

} /* dot */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

#ifdef TIME_VB

void AZ_time_kernals(int gN, double gnzeros, double val[], int indx[],
                     int bindx[], int rpntr[], int cpntr[], int bpntr[],
                     double x[], double y[], int ntimes, int options[],
                     int data_org[], int proc_config[], double params[],
                     AZ_MATRIX *Amat)

/*******************************************************************************

  Solve the system of equations given in the VBR format using an iterative
  method specified by 'options[AZ_solver]'. Store the result in 'x'.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  gN:              Global order of the linear system of equations.

  gnzeros:         Global number of nonzeros in val.

  val:             Array containing the nonzero entries of the matrix (see
                   Aztec User's Guide).

  indx,
  bindx,
  rpntr,
  cpntr,
  bpntr:           Arrays used for DMSR and DVBR sparse matrix storage (see
                   file Aztec User's Guide).

  x, y:            Vectors of length gN.

  ntimes:          Number of times to perform each operation.

  options:         Determines specific solution method and other parameters.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  register int i, j;
  double       start_t, total_t;
  double       sprx_comm_t, sprx_comp_t, dot_comm_t, dot_comp_t;
  double       sprx_t, sprx_overlap_cop_t;
  double       daxpy_comp_t, dzpax_comp_t;
  double       Djunk;
  double      *ptr_vec1, *ptr_vec2, *ptr_vec3, *z1, *z2, *z3, alpha = 1.0;
  int          one = 1, bpntr_index;
  double       Mflop, Mflops_comp, Mflops_node_comp;
  double       Mflops_comm, Mflops_node_comm;
  double       sparax_overlap_border_comp_t, sparax_overlap_internal_comp_t;
  double       read_nonoverlap_t, read_nonoverlap_t_max, read_t_max;
  double       time_overlap_max;
  double       gather_t, write_t, read_t, time;
  double       max, min, avg;
  double       overall_t, start_overall_t;

  int          Num_Proc, Proc;
  int          Num_Internal_Blks;
  int          iout;

  char        *message_recv_add[AZ_MAX_NEIGHBORS];
  char        *message_send_add[AZ_MAX_NEIGHBORS];
  int          message_recv_length[AZ_MAX_NEIGHBORS];
  int          message_send_length[AZ_MAX_NEIGHBORS];
  AZ_MATRIX   Amat2;

  /*
    message_send_add:
    message_send_add[i] points to the beginning of the list of
    values to be sent to the ith neighbor (i.e. data_org[AZ_neighbors+i] or
    sometimes locally defined as proc_num_neighbor[i]). That is,
    *(message_send_add[i] + j) is the jth value to be sent to the
    ith neighbor.

    message_send_length:
    message_send_length[i] is the number of bytes to be sent to
    the ith neighbor (i.e. data_org[AZ_neighbors+i] or sometimes
    locally defined as proc_num_neighbor[i]).

    message_recv_add:
    message_recv_add[i] points to the beginning of the list of
    locations which are to receive values sent by the ith neighbor (i.e.
    data_org[AZ_neighbors+i] or sometimes locally defined as
    proc_num_neighbor[i]). That is, *(message_recv_add[i] + j) is the
    location where the jth value sent from the ith neighbor will be stored.
    message_recv_length:
    message_recv_length[i] is the number of bytes to be sent to
    the ith neighbor (i.e. data_org[AZ_neighbors+i] or sometimes
    locally defined as proc_num_neighbor[i]).
    */

  int          temp1, temp2;

  /**************************** execution begins ******************************/

  Proc              = proc_config[AZ_node];
  Num_Proc          = proc_config[AZ_N_procs];
  Num_Internal_Blks = data_org[AZ_N_int_blk];
  iout              = options[AZ_output];

  if (data_org[AZ_matrix_type] == AZ_USER_MATRIX){
      if (az_proc == 0) {
        (void) AZ_printf_err( "ERROR: matrix-vector timing not available for"
                       "matrix-free.\n"
                       "      (data_org[AZ_matrix_type] = %d)\n\n",
                       data_org[AZ_matrix_type]);
      }
      exit (-1);
  }

  if (data_org[AZ_matrix_type] != AZ_VBR_MATRIX) {
    (void) AZ_printf_err( "I'm not sure if the timing stuff works for \n"
                   "nonVBR matrices. For example I don't think that\n"
                   "we have overlapped communication/computations.\n"
                   "However, probably a lot of the stuff works?\n");
    exit(-1);
  }

  if (ntimes < 0)   ntimes = -ntimes;
  if (ntimes < 200) ntimes = 200;

  if (Proc == 0) {
    (void) AZ_printf_out("++++++++++ timing run ++++++++++\n");
    (void) AZ_printf_out("\nnprocs: %d ntimes: %d\n", Num_Proc, ntimes);
    (void) AZ_printf_out("N: %d\tNzeros: %e\t\n\n", gN, gnzeros);
  }

  /* sparax */

  if (iout > 0 && Proc == 0) {
    (void) AZ_printf_out("timing nonoverlapped matrix vector multiply\n");
  }

  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++)
    Amat->matvec(x, r, Amat, proc_config);


  sprx_t = AZ_second() - start_t;

  /* exchange boundary */

  if (iout > 0 && Proc ==0){
    (void) AZ_printf_out("timing unstructured communication\n");
  }

  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++) {
    AZ_exchange_bdry(x, data_org, proc_config);
    AZ_exchange_bdry(y, data_org, proc_config);
  }

  sprx_comm_t = 0.5*(AZ_second() - start_t);
  sprx_comp_t = sprx_t - sprx_comm_t;

  /* overlapped communication */

  if (iout > 0 && Proc == 0) {
    (void) AZ_printf_out("timing overlapped sparse matrix vector multiply\n");
  }
  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++)
    Amat->matvec(x, r, Amat, proc_config);

  time = AZ_second() - start_t;

  gather_t                       = (double) 0.0;
  write_t                        = (double) 0.0;
  read_t                         = (double) 0.0;
  read_nonoverlap_t              = (double) 0.0;
  sparax_overlap_internal_comp_t = (double) 0.0;
  sparax_overlap_border_comp_t   = (double) 0.0;

  if (iout > 0 && Proc == 0) {
    (void) AZ_printf_out("time the individual routines for sparax_overlap\n");
  }

  start_overall_t = AZ_second();

  for (i = 1; i <= ntimes; i++) {

    /* time the individual routines for sparax_overlap */

    start_t = AZ_second();
    AZ_gather_mesg_info(x, data_org, message_recv_add, message_send_add,
                        message_recv_length, message_send_length);
    gather_t += AZ_second() - start_t;

    start_t = AZ_second();
    AZ_write_local_info(data_org, message_recv_add, message_send_add,
                        message_recv_length, message_send_length);
    write_t += AZ_second() - start_t;

    start_t = AZ_second();

    temp1                      = data_org[AZ_N_bord_blks];
    temp2                      = data_org[AZ_N_border];
    data_org[AZ_N_bord_blks] = 0;
    data_org[AZ_N_border] = 0;

    Amat->matvec(x, r, Amat, proc_config);


    data_org[AZ_N_bord_blks] = temp1;
    data_org[AZ_N_border] = temp2;

    sparax_overlap_internal_comp_t += AZ_second() - start_t;

    start_t = AZ_second();
    AZ_read_local_info(data_org, message_recv_add, message_recv_length);
    read_t += AZ_second() - start_t;

    start_t = AZ_second();

    /* compute boundary portion of the sparse matrix - vector product */

    bpntr_index = bpntr[Num_Internal_Blks];

    temp1                      = data_org[AZ_N_int_blks];
    temp2                      = data_org[AZ_N_internal];
    data_org[AZ_N_int_blks]    = data_org[AZ_N_bord_blks];
    data_org[AZ_N_internal]    = data_org[AZ_N_border];
    data_org[AZ_N_bord_blks]   = 0;
    data_org[AZ_N_border]      = 0;

/*
    Amat2.rpntr = &(Amat->rpntr[Num_Internal_Blks]);
    Amat2.cpntr =   Amat->cpntr;
    Amat2.bpntr = &(Amat->bpntr[Num_Internal_Blks]);
    Amat2.bindx = &(Amat->bindx[bpntr_index]);
    Amat2.matvec      = Amat->matvec;
    Amat2.matrix_type = Amat->matrix_type;
    Amat2.data_org = &(Amat->data_org);
    Amat2.val   = &(Amat->val[bpntr_index]);
    Amat2.indx  = &(Amat->indx[bpntr_index]);
*/

    Amat->matvec(x, r, Amat, proc_config);

/*
                   &bindx[bpntr_index], &rpntr[Num_Internal_Blks], cpntr,
                   &bpntr[Num_Internal_Blks], x, &y[rpntr[Num_Internal_Blks]],
                   0, data_org);
*/

    data_org[AZ_N_bord_blks] = data_org[AZ_N_int_blks];
    data_org[AZ_N_border]    = data_org[AZ_N_internal];
    data_org[AZ_N_int_blks]  = temp1;
    data_org[AZ_N_internal]  = temp2;

    sparax_overlap_border_comp_t += AZ_second() - start_t;
  }

  overall_t = AZ_second() - start_overall_t;

  for (i = 1; i <= ntimes; i++) {

    /* time the reads in the nonoverlapped case */

    AZ_gather_mesg_info(x, data_org, message_recv_add, message_send_add,
                        message_recv_length, message_send_length);

    AZ_write_local_info(data_org, message_recv_add, message_send_add,
                        message_recv_length, message_send_length);

    start_t = AZ_second();
    AZ_read_local_info(data_org, message_recv_add, message_recv_length);
    read_nonoverlap_t += AZ_second() - start_t;
  }

  read_t_max            = AZ_gmax_double(read_t, proc_config);
  read_nonoverlap_t_max = AZ_gmax_double(read_nonoverlap_t, proc_config);
  time_overlap_max      = AZ_gmax_double(time, proc_config);

  /* dot */

  if (iout > 0 && Proc == 0) {
    (void) AZ_printf_out("time the individual routines for ddot\n");
  }

  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++) Djunk = AZ_gdot(N, x, y, proc_config);
  total_t = AZ_second() - start_t;

  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++) AZ_gsum_double(Djunk, proc_config);

  dot_comm_t = AZ_gmax_double(AZ_second() - start_t, proc_config);
  dot_comp_t = AZ_gmax_double(total_t - dot_comm_t, proc_config);

  /* daxpy */

  if (iout > 0 && Proc == 0) {
    (void) AZ_printf_out("time the individual routines for daxpy\n");
  }
  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++) {
    Djunk = (double) i;
    DAXPY_F77(&N, &Djunk, x, &one, y, &one);
  }
  daxpy_comp_t = AZ_gmax_double(AZ_second() - start_t, proc_config);

  if (Proc == 0) {              /* calculate and print results */
    (void) AZ_printf_out("\n********** sparax ***********\n\n");
    (void) AZ_printf_out("nonoverlapped\n");
    (void) AZ_printf_out("\t\tmax\t\tavg\t\tmin\n");
  }

  /* dzpax */

  if (iout > 0 && Proc == 0) {
    (void) AZ_printf_out("time the individual routines for dzpax\n");
  }

  z1 = (double *) AZ_allocate(N * sizeof(double));
  z2 = (double *) AZ_allocate(N * sizeof(double));
  z3 = (double *) AZ_allocate(N * sizeof(double));
  for (i = 0; i < N; i++) {
    z2[i] = (double) i;
    z3[i] = (double) 2*i;
  }

  start_t = AZ_second();
  for (i = 1; i <= ntimes; i++) {
    Djunk = 3.14 * (double) i;
    for (j = 0; j < N; j++)
      z1[j] = z2[j] + Djunk*z3[j];
    /*    dzpax_(&N, &Djunk, z3, &one, z2, &one, z1, &one);*/
  }
  dzpax_comp_t = AZ_gmax_double(AZ_second() - start_t, proc_config);

  AZ_free((void *) z1);
  AZ_free((void *) z2);
  AZ_free((void *) z3);

  if (Proc == 0) {              /* calculate and print results */
    (void) AZ_printf_out("\n********** sparax ***********\n\n");
    (void) AZ_printf_out("nonoverlapped\n");
    (void) AZ_printf_out("\t\tmax\t\tavg\t\tmin\n");
  }

  max = AZ_gmax_double(sprx_t, proc_config);
  avg = AZ_gavg_double(sprx_t, proc_config);
  min = AZ_gmin_double(sprx_t, proc_config);
  if (Proc == 0) (void) AZ_printf_out("total_time\t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(sprx_comp_t, proc_config);
  avg = AZ_gavg_double(sprx_comp_t, proc_config);
  min = AZ_gmin_double(sprx_comp_t, proc_config);
  if (Proc == 0) (void) AZ_printf_out("comp_time\t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(sprx_comm_t, proc_config);
  avg = AZ_gavg_double(sprx_comm_t, proc_config);
  min = AZ_gmin_double(sprx_comm_t, proc_config);
  if (Proc == 0) (void) AZ_printf_out("comm_time\t%e\t%e\t%e\n", max, avg, min);

  sprx_comp_t = AZ_gmax_double(sprx_comp_t, proc_config);
  sprx_t      = AZ_gmax_double(sprx_t, proc_config);

  if (Proc == 0) {
    Mflop             = (double) ntimes * (gnzeros + gnzeros) * 1.0e-6;
    Mflops_comp       = Mflop/sprx_comp_t;
    Mflops_node_comp  = Mflops_comp/(double) Num_Proc;
    Mflops_comm       = Mflop/sprx_t;
    Mflops_node_comm  = Mflops_comm/(double) Num_Proc;

    (void) AZ_printf_out("computation Mflops: %e \n", Mflops_comp);
    (void) AZ_printf_out("computation Mflops per node: %e \n", Mflops_node_comp);
    (void) AZ_printf_out("comp & comm Mflops: %e \n", Mflops_comm);
    (void) AZ_printf_out("comp & comm Mflops per node: %e \n\n", Mflops_node_comm);
  }

  /* statistics */

  if (Proc == 0) {
    (void) AZ_printf_out("overlapped\n\t\tmax\t\tavg\t\tmin\n");
  }

  max = AZ_gmax_double(time, proc_config);
  min = AZ_gmin_double(time, proc_config);
  avg = AZ_gavg_double(time, proc_config);
  if (Proc == 0) (void) AZ_printf_out("total_time\t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(gather_t, proc_config);
  min = AZ_gmin_double(gather_t, proc_config);
  avg = AZ_gavg_double(gather_t, proc_config);
  if (Proc == 0) (void) AZ_printf_out("gather_t\t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(write_t, proc_config);
  min = AZ_gmin_double(write_t, proc_config);
  avg = AZ_gavg_double(write_t, proc_config);
  if (Proc == 0) (void) AZ_printf_out("write_t \t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(sparax_overlap_internal_comp_t, proc_config);
  min = Az_gmin_double(sparax_overlap_internal_comp_t, proc_config);
  avg = AZ_gavg_double(sparax_overlap_internal_comp_t, proc_config);
  if (Proc == 0) (void) AZ_printf_out("internal_t\t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(read_t, proc_config);
  min = Az_gmin_double(read_t, proc_config);
  avg = AZ_gavg_double(read_t, proc_config);
  if (Proc == 0) (void) AZ_printf_out("read_t   \t%e\t%e\t%e\n", max, avg, min);

  max = AZ_gmax_double(sparax_overlap_border_comp_t, proc_config);
  min = AZ_gmin_double(sparax_overlap_border_comp_t, proc_config);
  avg = AZ_gavg_double(sparax_overlap_border_comp_t, proc_config);
  if (Proc == 0) (void) AZ_printf_out("border_t\t%e\t%e\t%e\n", max, avg, min);

  if (Proc == 0) {
    Mflops_comm      = Mflop/time_overlap_max;
    Mflops_node_comm = Mflops_comm/(double) Num_Proc;
    (void) AZ_printf_out("comp & comm Mflops: %e \n", Mflops_comm);
    (void) AZ_printf_out("comp & comm Mflops per node: %e \n\n", Mflops_node_comm);

    (void) AZ_printf_out("Ratio of overlapped/nonoverlapped sparax times: %e\n\n",
                  time_overlap_max/sprx_t);
    (void) AZ_printf_out("Ratio of overlapped/nonoverlapped read times: %e\n\n",
                  read_t_max/read_nonoverlap_t_max);
  }

  max = AZ_gmax_double(time, proc_config);
  if (time  ==  max) {
    (void) AZ_printf_out("max time proc\n");
    (void) AZ_printf_out("\nProc:%d\tNum_Neighbors: %d\n", Proc, data_org[AZ_N_neigh]);
    (void) AZ_printf_out("total_time: %e\n", time);
    (void) AZ_printf_out("gather_t  : %e\n", gather_t);
    (void) AZ_printf_out("write_t   : %e\n", write_t);
    (void) AZ_printf_out("internal_t: %e\n", sparax_overlap_internal_comp_t);
    (void) AZ_printf_out("read_t    : %e\n", read_t);
    (void) AZ_printf_out("border_t  : %e\n\n", sparax_overlap_border_comp_t);
  }

  max = AZ_gmax_double(overall_t, proc_config);
  if (overall_t  ==  max) {
    (void) AZ_printf_out("overall max time proc\n");
    (void) AZ_printf_out("\nProc:%d\tNum_Neighbors: %d\n", Proc, data_org[AZ_N_neigh]);
    (void) AZ_printf_out("overall total_time: %e\n", overall_t);
    (void) AZ_printf_out("total_time: %e\n", time);
    (void) AZ_printf_out("gather_t  : %e\n", gather_t);
    (void) AZ_printf_out("write_t   : %e\n", write_t);
    (void) AZ_printf_out("internal_t: %e\n", sparax_overlap_internal_comp_t);
    (void) AZ_printf_out("read_t    : %e\n", read_t);
    (void) AZ_printf_out("border_t  : %e\n\n", sparax_overlap_border_comp_t);
  }
  /*
    if (Proc == 0) {
    (void) AZ_printf_out("cop overlapped\n");
    (void) AZ_printf_out("\t\tmax\t\tavg\t\tmin\n");
    }

    max = AZ_gmax_double(sprx_overlap_cop_t, proc_config);
    min = AZ_gmin_double(sprx_overlap_cop_t, proc_config);
    avg = AZ_gavg_double(sprx_overlap_cop_t, proc_config);
    if (Proc == 0) (void) AZ_printf_out("total_time\t%e\t%e\t%e\n", max, avg, min);

    if (Proc == 0) {
    Mflops_comm = Mflop/max;
    Mflops_node_comm  = Mflops_comm/(double)Num_Proc;
    (void) AZ_printf_out("comp & comm Mflops: %e \n", Mflops_comm);
    (void) AZ_printf_out("comp & comm Mflops per node: %e \n\n", Mflops_node_comm);
    }
    */
  if (Proc == 0) {
    (void) AZ_printf_out("\n********* exchange **********\n\n");
    (void) AZ_printf_out("comm_time = %7.4e sec.\n", sprx_comm_t);

    (void) AZ_printf_out("\n************ dot ************\n\n");
    Mflop = (double) ntimes * 2.0 * ((double) gN * 1.0e-06);

    (void) AZ_printf_out("comp_time = %7.4e sec.\n", dot_comp_t);
    (void) AZ_printf_out("comm_time = %7.4e sec.\n", dot_comm_t);

    Mflops_comp      = Mflop/dot_comp_t;
    Mflops_node_comp = Mflops_comp/(double) Num_Proc;
    Mflops_comm      = Mflop/(dot_comm_t + dot_comp_t);
    Mflops_node_comm = Mflops_comm/(double) Num_Proc;

    (void) AZ_printf_out("computation Mflops: %e \n", Mflops_comp);
    (void) AZ_printf_out("computation Mflops per node: %e \n", Mflops_node_comp);
    (void) AZ_printf_out("comp & comm Mflops: %e \n", Mflops_comm);
    (void) AZ_printf_out("comp & comm Mflops per node: %e \n\n", Mflops_node_comm);

    (void) AZ_printf_out("\n*********** daxpy ***********\n\n");
    (void) AZ_printf_out("comp_time = %7.4e sec.\n", daxpy_comp_t);
    Mflops_comp      = Mflop/daxpy_comp_t;
    Mflops_node_comp = Mflops_comp/(double) Num_Proc;

    (void) AZ_printf_out("computation Mflops: %e \n", Mflops_comp);
    (void) AZ_printf_out("computation Mflops per node: %e \n", Mflops_node_comp);

    (void) AZ_printf_out("\n*********** dzpax ***********\n\n");
    (void) AZ_printf_out("comp_time = %7.4e sec.\n", dzpax_comp_t);
    Mflops_comp  = Mflop/dzpax_comp_t;
    Mflops_node_comp  = Mflops_comp/(double) Num_Proc;

    (void) AZ_printf_out("computation Mflops: %e \n", Mflops_comp);
    (void) AZ_printf_out("computation Mflops per node: %e \n", Mflops_node_comp);
  }

} /* AZ_time_kernals */

#endif

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_find_index(int key, int list[], int length)

/*******************************************************************************

  Find 'key' in 'list' and return the index number.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     int, -1 = key not found, i = list[i] = key
  ============

  Parameter list:
  ===============

  key:             Element to be search for in list.

  list:            List to be searched.

  length:          Length of list.

*******************************************************************************/

{

  /* local variables */

  int start, end;
  int mid;

  /**************************** execution begins ******************************/

  if (length == 0) return -1;

  start = 0;
  end   = length - 1;

  while (end - start > 1) {
    mid = (start + end) / 2;
    if (list[mid] < key) start = mid;
    else end = mid;
  }

  if (list[start] == key) return start;
  if (list[end] == key)   return end;
  return -1;

} /* AZ_find_index */

int AZ_exit(int input)
{
#ifdef AZTEC_MPI
   MPI_Finalize();
#endif
   exit(input);
   return(1);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_free_memory(int label)

/*******************************************************************************

  Clear all internal memory that has been allocated by previous calls to
  AZ_manage_memory(*,*,label,*,*). This memory is typically includes
  preconditioner and scaling information.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  label:           An integer associated with memory that was allocated with
                   AZ_manage_memory(). On ouput, all memory allocated with this
                   integer will be freed.

*******************************************************************************/

{
  (void) AZ_manage_memory((int) NULL, AZ_CLEAR, label, (char *) NULL,
                          (int *) NULL);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double *AZ_manage_memory(unsigned int input_size, int action, int type, 
                         char *name, int *status)

/*******************************************************************************

  AZ_manage_memory() either frees memory that was previously allocated or
  returns a pointer to memory for future use (this pointer can be a newly
  allocated block or a block that was previously allocated).


  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     double *,  pointer to allocated memory.
  ============

  Parameter list:
  ===============

  input_size:      integer variable. On input, size indicates the amount of
                   memory that is needed (NOTE: ignored if action == AZ_CLEAR).

  action:          On input, action indicates whether to allocate or free
                   memory.
                   = AZ_ALLOC: look for a chunk of memory (in the memory
                               management list pointed to by 'head') that has
                               already been allocated with a memory-type 'type',
                               a size in bytes of 'size', and labelled with the
                               string 'name'.

                               If this memory is found return a pointer to this
                               memory and set *status to AZ_OLD_ADDRESS.
                               If this memory is not found allocate the memory,
                               put it in the memory management list, return a
                               pointer to this piece of memory, and set *status
                               to AZ_NEW_ADDRESS.

                   = AZ_REALLOC:look for a chunk of memory (in the memory
                               management list pointed to by 'head') that has
                               already been allocated with a memory-type 'type',
                               and labelled with the string 'name'. Reallocate
                               this item with size 'size' and change the size
                               field.

                   = AZ_CLEAR: free all memory associated with any chunk of
                               memory that has the memory-type 'type'.

                   = AZ_CLEAR_ALL: free all memory associated with any chunk of
                               memory regardless of 'type'.

                   = AZ_SELECTIVE_CLEAR: free all memory associated with any
                               chunk ofmemory that has the memory-type 'type'
                               and the memory-name 'name'.

  type:            On input, a number associated with each piece of memory that
                   is allocated.  This number is used when checking to see if
                   a requested piece of memory is already allocated. This number
                   is also used when freeing up memory (see AZ_CLEAR).

  name:            On input, a character string associated with each piece of
                   memory that is allocated. This string is used when checking
                   to see if a requested piece of memory is already allocated.
                   Additionally, this string is printed when AZ_manage_memory()
                   fails to allocate memory.

  status:          On output,
                   *status = AZ_NEW_ADDRESS indicates that a new piece of memory
                                            has been allocated.
                   *status = AZ_OLD_ADDRESS indicates that a pointer to a
                                            previously allocated chunk of memory
                                            is returned.
                   NOTE: status is not used if action == AZ_CLEAR has been
                         allocated.


*******************************************************************************/

{

  /* local variables */

  struct mem_ptr {
    char   *name;
    double *address;
    int     type;
    int     size;
    char    special;
    struct mem_ptr *next;
  };

  long int size;
  static struct mem_ptr *head = NULL;
  struct mem_ptr        *current, *temp,*prev, *thenext;
  int                   found = 0, i,j, n2, nn;
  /* FIX: long int aligned_str_mem, aligned_j, aligned_size; */
int aligned_str_mem,aligned_j,aligned_size;
double *dtmp;

  /**************************** execution begins ******************************/
  size = (long int) input_size;
  current = head;

  if (action == -43) {

    /* print the list */

    while (current != NULL) {
      (void) AZ_printf_out("(%8d,%8d  %ld) ==> %s\n", current->type, current->size,
                      (long int) current->address, current->name);
      prev    = current;
      current = current->next;
    }
    return (double *) 0;
  }

  else if (action == AZ_ALLOC) {
    *status = AZ_OLD_ADDRESS;
    if (size == 0) return (double *) 0;

    /* first look for entry */
    while ( current != NULL) {
      if ( (current->size == size) && (current->type == type) &&
           (strcmp(current->name,name) == 0) ) {
        found = 1;
        break;
      }
      prev    = current;
      current = current->next;
    }

    if (found == 0) {

      /*
       * Put in a new entry if it has the type 0 we will put it at the front of
       * the list else put it at the end of the list This is done for efficiency
       * reasons as type = 0 is used a lot.
       */

      j = strlen(name) + 1;
      aligned_str_mem = sizeof(struct mem_ptr);
      aligned_j       = j;
      aligned_size    = size;

      aligned_str_mem +=  (sizeof(double) - (aligned_str_mem%sizeof(double)));
      aligned_j       +=  (sizeof(double) - (aligned_j%sizeof(double)));
      aligned_size    +=  (sizeof(double) - (aligned_size%sizeof(double)));

      dtmp = (double *) AZ_allocate((unsigned int) (aligned_str_mem+aligned_j+
                                                aligned_size) );

      if (dtmp== NULL) {
        (void) AZ_printf_err( "Error: Not enough space to allocate\n");
        (void) AZ_printf_err( "       '%s' in memory_manage()\n", name);
        (void) AZ_printf_err( "        Asked for %ld bytes. Perhaps\n", size);
        (void) AZ_printf_err( "        a smaller problem should be run\n");
        exit(-1);
      }
      temp = (struct mem_ptr *) &(dtmp[aligned_size/sizeof(double)]);
      temp->name = (char *)     &(dtmp[(aligned_str_mem+aligned_size)/
				          sizeof(double)]);
      temp->address = dtmp;
      for (i = 0 ; i < j; i++ ) (temp->name)[i] = name[i];



      temp->special = 'N';
      temp->type = type;
      temp->size = size;

      /* put at the head of the list */

      temp->next = head;
      head       = temp;

      *status = AZ_NEW_ADDRESS;
      return temp->address;
    }
    else {
      *status = AZ_OLD_ADDRESS;
      return current->address;
    }
  }

  else if (action == AZ_CLEAR) {
    prev = NULL;
    while (current != NULL) {
      if (current->type == type) {
        if (prev == NULL) head       = current->next;
        else              prev->next = current->next;
        temp = current->next;
        AZ_free(current->address);
        current = temp;
      }
      else {
        prev    = current;
        current = current->next;
      }
    }
    return (double *) 0;
  }

  else if (action == AZ_CLEAR_ALL) {
    prev = NULL;
    while (current != NULL) {
      if (prev == NULL) head       = current->next;
      else              prev->next = current->next;
      temp = current->next;
      AZ_free(current->address);
      current = temp;
    }
    return (double *) 0;
  }

  else if (action == AZ_EVERYBODY_BUT_CLEAR) {
    prev = NULL;

    while (current != NULL) {
      if ( (current->type == type) &&
           (strncmp(current->name,name,5) != 0) ) {
        if (prev == NULL) head       = current->next;
        else              prev->next = current->next;

        temp = current->next;

        AZ_free(current->address);
        current = temp;
      }
      else {
        prev    = current;
        current = current->next;
      }
    }
    return (double *) 0;
  }
  else if (action == AZ_SUBSELECTIVE_CLEAR) {
    prev = NULL;

    while (current != NULL) {
      if ( (current->type == type) &&
           (strncmp(current->name,name,5) == 0) ) {
        if (prev == NULL) head       = current->next;
        else              prev->next = current->next;

        temp = current->next;

        AZ_free(current->address);
        current = temp;
      }
      else {
        prev    = current;
        current = current->next;
      }
    }
    return (double *) 0;
  }
  else if (action == AZ_SELECTIVE_CLEAR) {
    prev = NULL;

    while (current != NULL) {
      if ( (current->type == type) &&
         (strcmp(current->name,name) == 0) ) {
        if (prev == NULL) head       = current->next;
        else              prev->next = current->next;

        temp = current->next;

        AZ_free(current->address);
        current = temp;
      }
      else {
        prev    = current;
        current = current->next;
      }
    }
    return (double *) 0;
  }

  else if ( (action == AZ_REALLOC) || (action == AZ_SPEC_REALLOC)) {
     prev = NULL;

    /* first look for entry */

    while ( current != NULL) {
      if ( (current->type == type) &&
           (strcmp(current->name,name) == 0) ) {
        found = 1;
        break;
      }
      prev    = current;
      current = current->next;
    }
    if (current == NULL) {
      (void) AZ_printf_err( "memory_management error: %s with type %d not",
                     name, type);
      (void) AZ_printf_err( " found during reallocation\n");
      exit(-1);
    }
    if ( (action == AZ_REALLOC) && (current->special != 'N')) {
        *status = AZ_SPECIAL;
        return(current->address);
    }
    else if (current->special != 'N') {
       *status = AZ_SPECIAL;
    }
    else *status = AZ_OLD_ADDRESS;

    if (current->size == size) return(current->address);


    j = strlen(name) + 1;
    aligned_str_mem = sizeof(struct mem_ptr);
    aligned_j       = j;
    aligned_size    = size;

    aligned_str_mem +=  (sizeof(double) - (aligned_str_mem%sizeof(double)));
    aligned_j       +=  (sizeof(double) - (aligned_j%sizeof(double)));
    aligned_size    +=  (sizeof(double) - (aligned_size%sizeof(double)));

    thenext = current->next;
    dtmp    = current->address;
    dtmp    = (double *) AZ_realloc((char *) dtmp,(unsigned int) 
                                    aligned_str_mem+aligned_j+aligned_size);
    if (dtmp == NULL) {
      (void) AZ_printf_err( "Error:Not enough memory for '%s'\n", name);
      (void) AZ_printf_err( "      Asked for %ld bytes. Perhaps\n", size);
      (void) AZ_printf_err( "      a smaller problem should be run\n");
      exit(-1);
    }
    temp = (struct mem_ptr *) &(dtmp[aligned_size/sizeof(double)]);
    temp->name = (char *)     &(dtmp[(aligned_str_mem+aligned_size)/
                                          sizeof(double)]);
    temp->address = dtmp;
    for (i = 0 ; i < j; i++ ) (temp->name)[i] = name[i];
    temp->special = 'N';
    if (action == AZ_SPEC_REALLOC) temp->special = 'S';
    temp->type = type;
    temp->size = size;

    if (prev == NULL) head  = temp;
    else  prev->next = temp;
    temp->next = thenext;

    return temp->address;

  }
  else if (action == AZ_RESET_STRING) {
      prev = NULL;
 
     /* first look for entry */
 
     n2 = strlen(name);
     while ( current != NULL) {
       nn = strlen(current->name);
       if ( (current->type == type) &&
            (strncmp(current->name,name,(unsigned) nn) == 0) &&
            (n2 < nn*2)) {
         found = 1;
         break;
       }
       prev    = current;
       current = current->next;
     }
     if (current == NULL) {
       (void) AZ_printf_err( "memory_management error: %s with type %d not",
                      name, type);
       (void) AZ_printf_err( " found while changing name\n");
       exit(-1);
     }
     sprintf(current->name,"%s",&(name[nn]));
     *status = AZ_OLD_ADDRESS;
  }
  else if (action == AZ_EMPTY ) {
     if (current == NULL) *status = -1;
     else *status = 1;
  }
  else if (action == AZ_LOOKFOR_PRINT) {
    /* first look for entry */

    AZ_printf_out("Looking in system for possible mismatches\n");
    while ( current != NULL) {
      if (current->name[0] == 'P') {
         if (current->type != type) {
            AZ_printf_out("AZ_name/type does not match %d %d\n\n", current->type,type);
            return(NULL);
         }
         i = strcmp(current->name,name);
         if (i != 0) {
            AZ_printf_out("option keys do not match (%s) (%s)\n", current->name,name);
         
            nn = 1;
            while ( (current->name[nn] != ' ' )&&(nn < (int) strlen(current->name))) {
               if (current->name[nn] != name[nn]) break;
               nn++;
            }
            if ((current->name[nn] != ' ') || (name[nn] != ' ')) {
               sscanf(&(current->name[1]),"%d", &i);
               sscanf(&(name[1]),"%d", &j);
               AZ_printf_out("found context with different matrix size (%d vs. %d)\n",
                      i,j);
               return(NULL);
            }
            i = (int) ( name[++nn] - '!');
            j = (int) ( current->name[nn] - '!');
            if (j != i) {
               AZ_printf_out("==> check overlapping and reordering choices\n");
               return(NULL);
            }
            i = (int) ( name[++nn] - '!');
            j = (int) ( current->name[nn] - '!');
            if (j != i) {
              AZ_printf_out("==> check preconditioner and subdomain solver choices\n");
              return(NULL);
            }
            i = (int) ( name[++nn] - '!');
            j = (int) ( current->name[nn] - '!');
            if (j != i) {
              AZ_printf_out("==> check scaling choices\n");
              return(NULL);
            }
           AZ_printf_out("==> check AZ_ilut_fill, AZ_drop & AZ_graph_fill choices\n");
           return(NULL);
         }
      }
      prev    = current;
      current = current->next;
    }
  }
  else {
    (void) AZ_printf_err( "Error: Invalid action(%d) in AZ_manage_memory()\n",
                   action);
    exit(-1);
  }
  return((double *) NULL);

} /* AZ_manage_memory */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_mysleep(int i)

/*******************************************************************************

  A dummy routine that just wastes some time.  We essentially try and mimick the
  unix sleep() function which is not available on the nCUBE.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  i:               On input, indicates how much time to spend in this routine.
                   The larger the value of 'i', the more time is spent in this
                   routine. On the Ncube 2, this usually corresponds to about i
                   seconds.

*******************************************************************************/

{

  /* local variables */

  int    j, k, top, inner;
  double a;

  /**************************** execution begins ******************************/

  a = 10.01;
  top = 700 * i;
  inner = 10 * 10 * 10;

  for (k = 0; k < top; k++)
    for (j = 0; j < inner; j++)
      a = 2.0 / (a - 1);

} /* AZ_mysleep */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_p_error(char *str, int proc)

/*******************************************************************************

  Print the string if we are processor 0.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  str:             String to be printed.

  proc:            Current processor number.

*******************************************************************************/

{
  if (proc == 0) (void) AZ_printf_out( "%s", str);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_get_new_eps(double *epsilon, double recursive, double actual,
                   int options[], int proc_config[])

/*******************************************************************************

  Routine which decides what to do when the computed residual has converged but
  not the true residual.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  epsilon:         Residual tolerance.

  recursive:       On input, the norm of the residual produced by an iterative
                   method via recursion (e.g. r = r_0 + delta_1 + delta_2 +
                   delta_3 ... )

  actual:          On input, the norm of the explicitly computed residual (i.e.
                   r = b - A x ). In the absence of rounding error, recursive
                   and actual should be the same value. However, due to rounding
                   errors these two might differ.

  options:         Aztec options array: used to test output state.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  double difference;

  /**************************** execution begins ******************************/

  difference = fabs(actual - recursive);
  if (difference > *epsilon) return AZ_QUIT;
  else {
    *epsilon = *epsilon - 1.5 * difference;
    while (*epsilon < 0.0) *epsilon += .1*difference;
   }

  if (proc_config[AZ_node] == 0 && options[AZ_output]!=AZ_none)
    (void) AZ_printf_out("\n\t\tTrying to reduce actual residual "
                  "further\n\t\t     (recursive = %e, actual = %e)\n\n",
                  recursive, actual);

  return AZ_CONTINUE;

} /* AZ_get_new_eps */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_terminate_status_print(int situation, int iter, double status[],
                               double rec_residual, double params[],
                               double scaled_r_norm, double actual_residual,
                               int options[], int proc_config[])

/*******************************************************************************

  Routine to output conditions under which iterative solver was terminated if
  other than specified convergence.

  Author:          Scott A. Hutchinson, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  situation:       On input, indicates why iterative method terminated.
                   = AZ_normal   : iterative solver terminates normally (i.e.
                                   convergence criteria met).
                   = AZ_breakdown: iterative solver has broken down.
                   = AZ_loss     : iterative solver has terminated due to a lack
                                   of accuracy in the recursive residual
                                   (caused by rounding errors).
                   = AZ_maxits   : iterative solver has reached the maximum
                                   number of requested iterations without
                                   convergence.
                   = AZ_ill_cond : the upper triangular factor of the
                                   Hessenberg within GMRES is ill-conditioned.
                                   There is a possibility that the matrix
                                   is singular and we have the least-squares
                                   solution.

  iter:            Number of iterations completed.

  status:          !!!!!THIS MAY NOT BE ACCURATE FOR THIS ROUTINE!!!!!
                   On output, indicates termination status:
                    0:  terminated normally.
                   -1:  maximum number of iterations taken without achieving
                        convergence.
                   -2:  Breakdown. The algorithm can not proceed due to
                        numerical difficulties (usually a divide by zero).
                   -3:  Internal residual differs from the computed residual due
                        to a significant loss of precision.

  rec_residual:    On input, the norm of the residual produced by an iterative
                   method via recursion (e.g. r = r_0 + delta_1 + delta_2 +
                   delta_3 ... )

  params:          Drop tolerance and convergence tolerance info.

  scaled_r_norm:   Residual expression (requested by the user via
                   OPTIONS[AZ_conv] ... see user's guide and
                   AZ_compute_global_scalars()) which is used in printing and
                   compared with PARAMS[AZ_TOL] when determining convergence.

  actual_residual: On input, the norm of the explicitly computed residual (i.e.
                   r = b - A x ).  In the absence of rounding error, recursive
                   and actual should be the same value. However, due to rounding
                   errors these two might differ.

  options:         Determines specific solution method and other parameters.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  static int iterations = 0;
  char      *solver_name;
  int        solver_flag, conv_flag;
  double     eps;

  /**************************** execution begins ******************************/

  eps = params[AZ_tol];

  /* set status */

  if (scaled_r_norm < eps) situation = AZ_normal;

  status[AZ_its] = (double ) iter;
  status[AZ_why] = (double ) situation;
  status[AZ_r]   = actual_residual;

  if (actual_residual == -1.0) status[AZ_r] = rec_residual;

  status[AZ_rec_r]    = rec_residual;
  status[AZ_scaled_r] = scaled_r_norm;

  if (situation == AZ_normal) return; /* nothing else to do */

  if (options[AZ_output] ==  AZ_none) return;

  /* initialize */

  solver_flag = options[AZ_solver];
  conv_flag   = options[AZ_conv];
  if (!iterations) iterations = iter;

  switch (solver_flag) {
  case AZ_cg:
    solver_name = (char *) AZ_allocate(3*sizeof(char));
    (void) sprintf(solver_name, "cg");
    break;
  case AZ_cg_condnum:
    solver_name = (char *) AZ_allocate(11*sizeof(char));
    (void) sprintf(solver_name, "cg_condnum");
    break;
  case AZ_fixed_pt:
    solver_name = (char *) AZ_allocate(9*sizeof(char));
    (void) sprintf(solver_name, "fixed-pt");
    break;
  case AZ_GMRESR:
    solver_name = (char *) AZ_allocate(7*sizeof(char));
    (void) sprintf(solver_name, "gmresr");
    break;
  case AZ_analyze:
    solver_name = (char *) AZ_allocate(8*sizeof(char));
    (void) sprintf(solver_name, "analyze");
    break;
  case AZ_gmres:
    solver_name = (char *) AZ_allocate(6*sizeof(char));
    (void) sprintf(solver_name, "gmres");
    break;
  case AZ_gmres_condnum:
    solver_name = (char *) AZ_allocate(14*sizeof(char));
    (void) sprintf(solver_name, "gmres_condnum");
    break;
  case AZ_cgs:
    solver_name = (char *) AZ_allocate(4*sizeof(char));
    (void) sprintf(solver_name, "cgs");
    break;
  case AZ_tfqmr:
    solver_name = (char *) AZ_allocate(7*sizeof(char));
    (void) sprintf(solver_name, "tfqmr");
    break;
  case AZ_bicgstab:
    solver_name = (char *) AZ_allocate(9*sizeof(char));
    (void) sprintf(solver_name, "bicgstab");
    break;
  case AZ_symmlq:
    solver_name = (char *) AZ_allocate(7*sizeof(char));
    (void) sprintf(solver_name, "symmlq");
    break;
  case AZ_lu:
    solver_name = (char *) AZ_allocate(4*sizeof(char));
    (void) sprintf(solver_name, "lu");
    break;
  default:
    (void) AZ_printf_err(
                   "Error: invalid solver flag %d passed to terminate_status\n",
                   solver_flag);
  exit(-1);
  }

  if (proc_config[AZ_node] == 0) {
    (void) AZ_printf_err( "\n\n");
    (void) AZ_printf_err("\t*************************************************"
                   "**************\n\n");

    switch (situation) {
    case AZ_ill_cond:
      (void) AZ_printf_err( "\tWarning: the GMRES Hessenberg matrix is "
                     "ill-conditioned.  This may \n\tindicate that the "
                     "application matrix is singular. In this case, GMRES\n"
                     "\tmay have a least-squares solution.\n");
    break;
    case AZ_breakdown:
      if (!solver_flag) {
        (void) AZ_printf_err( "\tWarning: direction vector is no longer "
                       "A-conjugate \n\tto previous vector but solution has "
                       "not converged.\n");
      }
      else {
        (void) AZ_printf_err( "\tWarning: a breakdown in this "
                       "method\n\thas occurred and solution has not "
                       "converged.\n");
      }
      break;

    case AZ_loss:
      (void) AZ_printf_err( "\tWarning: recursive residual indicates "
                     "convergence\n\tthough the true residual is too large.\n");
      (void) AZ_printf_err( "\n\tSometimes this occurs when storage is ");
      (void) AZ_printf_err( "overwritten (e.g. the\n\tsolution vector was not ");
      (void) AZ_printf_err( "dimensioned large enough to hold\n\texternal ");
      (void) AZ_printf_err( "variables). Other times, this is due to roundoff. ");
      (void) AZ_printf_err( "In\n\tthis case, the solution has either ");
      (void) AZ_printf_err( "converged to the accuracy\n\tof the machine or ");
      (void) AZ_printf_err( "intermediate roundoff errors ");
      (void) AZ_printf_err( "occurred\n\tpreventing full convergence. In the ");
      (void) AZ_printf_err( "latter case, try solving\n\tagain using the new ");
      (void) AZ_printf_err( "solution as an initial guess.\n");
      break;

    case AZ_maxits:
      (void) AZ_printf_err( "\tWarning: maximum number of iterations "
                     "exceeded\n\twithout convergence\n");
    break;

    default:
      (void) AZ_printf_err( "\tError: improper code passed from solver %s\n\n",
                     solver_name);
    (void) AZ_printf_err("\t***********************************************%s",
                   "**********\n\n");
    exit(-1);
    }

    (void) AZ_printf_out("\n\tSolver:\t\t\t%s\n", solver_name);
    (void) AZ_printf_out("\tnumber of iterations:\t%d\n\n", iter);

    if (actual_residual != -1.0)
      (void) AZ_printf_out("\tActual residual = %11.4e",actual_residual);

    (void) AZ_printf_out("\tRecursive residual = %11.4e\n\n",rec_residual);
    (void) AZ_printf_out("\tCalculated Norms\t\t\t\tRequested Norm\n");
    (void) AZ_printf_out("\t--------------------------------------------");
    (void) AZ_printf_out("\t--------------\n\n");

    switch (conv_flag) {
    case AZ_noscaled:
      (void) AZ_printf_out( "\t||r||_2 :\t\t%e\t%e\n", scaled_r_norm,
                     eps);
    break;
    case AZ_r0:
      (void) AZ_printf_out( "\t||r||_2 / ||r0||_2:\t\t%e\t%e\n", scaled_r_norm,
                     eps);
    break;
    case AZ_rhs:
      (void) AZ_printf_out( "\t||r||_2 / ||b||_2:\t\t%e\t%e\n", scaled_r_norm,
                     eps);
    break;
    case AZ_Anorm:
      (void) AZ_printf_out( "\t||r||_2 / ||A||_inf:\t\t%e\t%e\n",
                     scaled_r_norm, eps);
    break;
    case AZ_sol:
      (void) AZ_printf_out( "\t\t||r||_inf\n");
    (void) AZ_printf_out( "\t----------------------------- : %e\t%e\n",
                   scaled_r_norm, eps);
    (void) AZ_printf_out( "\t||A||_inf ||x||_1 + ||b||_inf\n");
    break;
    case AZ_expected_values:
    case AZ_weighted:
      (void) AZ_printf_out( "\t||r||_WRMS:\t\t%e\t%e\n", scaled_r_norm, eps);
    break;
    case AZTECOO_conv_test:
      (void) AZ_printf_out( "\tUser-defined AztecOO_StatusTest\n");
    break;
    default:
      (void) AZ_printf_err( "terminate_status: ERROR: convergence test %d "
                     "not implemented\n", conv_flag);
    exit(-1);
    }

    (void) AZ_printf_err("\n\t*********************************************"
                   "******************\n\n");
  }

  /* free memory */

  if (solver_name != NULL)
    AZ_free((void *) solver_name);

} /* AZ_terminate_status_print */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_breakdown_f(int N, double v[], double w[], double inner,
                   int proc_config[])

/*******************************************************************************

  Determine if the 2 vectors v and w are orthogonal. In particular,

       |<v,w>| <  100 ||v||_2 ||w||_2 DBL_EPSILON

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     int
  ============

  Parameter list:
  ===============

  N:               Length of vectors v and w.

  v, w:            Two vectors of length N to be checked for othogonality.

  inner:           <v,w>

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  double v_norm, w_norm;

  v_norm = AZ_gvector_norm(N, 2, v, proc_config);
  w_norm = AZ_gvector_norm(N, 2, w, proc_config);

  return (fabs(inner) <= 100.0 * v_norm * w_norm * DBL_EPSILON);

} /* breakdown */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_random_vector(double u[], int data_org[], int proc_config[])

/*******************************************************************************

  Set the the vector u to a random vector.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  u:               Vector to be initialized.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  static int seed = 0;
  static int start = 1;
  int        i,N;
  int maxint = 2147483647;

  /*********************** BEGIN EXECUTION *********************************/

  /* Distribute the seeds evenly in [-maxint,maxint].  This guarantees nothing
   * about where in a random number stream we are, but avoids overflow situations
   * in parallel when multiplying by a PID. */

  if (start) {
    seed = (int)(maxint * (1.0 - 2.0*(proc_config[AZ_node]+1)/(proc_config[AZ_N_procs]+1.0)));
    start = 0;
  }
  N    = data_org[AZ_N_internal] + data_org[AZ_N_border];

  for (i = 0; i < N; i++) u[i] = AZ_srandom1(&seed);

} /* AZ_random_vector */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_srandom1(int *seed)

/*******************************************************************************

  Random number generator.

  Author:
  =======

  Return code:     double, random number.
  ============

  Parameter list:
  ===============

  seed:            Random number seed.

*******************************************************************************/

{

  /* local variables */

  int    a = 16807;
  int    m = 2147483647;
  int    q = 127773;
  int    r = 2836;

  int    lo, hi, test;
  double rand_num;

  /**************************** execution begins ******************************/

  hi   = *seed / q;
  lo   = *seed % q;
  test = a * lo - r * hi;

  if (test > 0) *seed = test;
  else *seed = test + m;

  rand_num = (double) *seed / (double) m;

  return rand_num;

} /* AZ_srandom1 */



int allo_count = 0, free_count = 0;

#ifdef MEM_CHECK
/* sophisticated wrappers for allocating memory */

struct widget {                  /* widget is used to maintain a linked   */
   int order;                    /* list of all memory that was allocated */
   int size;
   char *address;
   struct widget *next;
};
struct widget *widget_head =  NULL;  /* points to first element of allocated */
                                     /* memory linked list.                  */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

void AZ_print_it() {
/*
 * Print out the allocated memory linked list
 */

   struct widget *current;

   current = widget_head;
   while (current != NULL) {
      AZ_printf_out("(%d,%d,%u)\n",current->order, current->size, current->address);
      current = current->next;
   }
}
extern void AZ_print_it();

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

char *AZ_allocate(unsigned int isize) {

/* 
 *  Allocate memory and record the event by placing an entry in the 
 *  widget_head list. Also recored the size of this entry as well.
 *
 *  Note: we actually allocate more memory than is requested (7 doubles more).
 *  This additional memory is used to record the 'size' and to mark the
 *  memory with a header and trailer which we can later check to see if
 *  they were overwritten.
 *
 */

    char *ptr, *header_start, *header_end;
    struct widget *widget;
    /* FIX: int *size_ptr, i;
       FIX: unsigned int size; */
    int *size_ptr, i, size;
    double *dptr;

    /* FIX: size = isize; */
    size = (int) isize;

    size = size + 7*sizeof(double);
    widget = (struct widget *) malloc(sizeof(struct widget));
    if (widget == NULL) return(NULL);
    ptr = (char *) malloc(size);
    if (ptr == NULL) {
       free(widget);
       return(NULL);
    }
    allo_count++;

    /* put trash in the space to make sure nobody is expecting zeros */
    for (i = 0 ; i < size/sizeof(char) ; i++ ) 
       ptr[i] = 'f';


    /* record the entry */

    widget->order = allo_count;
if (size == 7*sizeof(double) ) {
AZ_printf_out("allocating 0 space %u (%d)\n",ptr,size);
 i = 0;
 size = 1/i;
 widget = NULL;
}
    widget->size  = size - 7*sizeof(double);
    widget->next  = widget_head;
    widget_head   = widget;
    widget->address = ptr;

    size_ptr = (int *) ptr;
    size_ptr[0] = size - 7*sizeof(double);
    dptr     = (double *) ptr;

    /* mark the header */

    header_start = (char *) &(dptr[1]);

    for (i = 0 ; i < 3*sizeof(double)/sizeof(char) ; i++ )
       header_start[i] = 'x';

    /* mark the trailer */

    header_end = &(ptr[ (size/sizeof(char)) - 1]);
    header_start = (char *) &(dptr[4]);
    header_start = & (header_start[(size-7*sizeof(double))/sizeof(char)]);

    while (header_start <= header_end) {
       *header_start = 'x';
       header_start++;
    }

    return( (char *) &(dptr[4]) );
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

void AZ_free(void *vptr) {
/*
 * Free memory and remove the corresponding entry from the widget_head
 * list. Additionally, check that the size stored in the header is correct
 * and that there has not been any memory overwritten in the header or
 * the trailer.
 *
 */

   struct widget *current, *prev;
   double *dptr;
   int *iptr, size, i;
   char *header_start, *header_end, *ptr;

    ptr = (char *) vptr;
    free_count++;
    if (ptr == NULL) {
       AZ_printf_out("Trying to free a NULL ptr\n");
i = 0;
size = 1/i;
widget_head = NULL;
    }
    else {
       current = widget_head;
       prev    = NULL;
       dptr = (double *) ptr;
       --dptr;
       --dptr;
       --dptr;
       --dptr;
       ptr = (char *) dptr;
       while (current != NULL) {
          if (current->address == ptr) break;
          else { prev = current; current = current->next; }
       }
       if (current == NULL) {
          AZ_printf_out("the pointer %u was not found and thus can not be freed.\n",
                  ptr);
          exit(1);
       }
       else {
           /* check to see if the header is corrupted */
           iptr = (int *) ptr;
           header_start = (char *) &(dptr[1]);

           for (i = 0 ; i < 3*sizeof(double)/sizeof(char) ; i++ ) {
              if (header_start[i] != 'x') {
                 AZ_printf_out("header is corrupted for %u (%d,%d)\n",ptr,
                         current->size,current->order);
                 size =  0;
                 size = 1/size;
              }
           }
           size = iptr[0];

           /* check to see if the sizes are different */

           if (current->size != size) {
              AZ_printf_out("Freeing %u whose size has changed (%d,%d)\n",
                     current->address,current->size,size);
              exit(1);
           }

           /* check to see if the trailer is corrupted */

           header_end = &(ptr[ ((size+7*sizeof(double))/sizeof(char)) - 1]);
           header_start = (char *) &(dptr[4]);
           header_start = &(header_start[size/sizeof(char)]);

           while (header_start <= header_end) {
              if ( *header_start != 'x') {
                 AZ_printf_out("trailer is corrupted for %u (%d,%d)\n",
                         ptr, size,
                         current->order);
                 size =  0;
                 size = 1/size;
              }
              header_start++;
           }

           /* free the space and the widget */

           free(ptr);
           if (widget_head == current) widget_head = current->next;
           else prev->next = current->next;
           free(current);

       }
   }

}

char *AZ_realloc(void *vptr, unsigned int new_size) {

   struct widget *current, *prev;
   int i, *iptr, size, *new_size_ptr;
   char *header_start, *header_end, *ptr;
   char *data1, *data2, *new_ptr, *new_header_start, *new_header_end;
   int newmsize, smaller;
   double *dptr, *new_dptr;

    ptr = (char *) vptr;
    if (ptr == NULL) {
       AZ_printf_out("Trying to realloc a NULL ptr\n");
       exit(1);
    }
    else {
       current = widget_head;
       prev    = NULL;
data1 = ptr;
       dptr = (double *) ptr;
       --dptr;
       --dptr;
       --dptr;
       --dptr;
       ptr = (char *) dptr;
       while (current != NULL) {
          if (current->address == ptr) break;
          else { prev = current; current = current->next; }
       }
       if (current == NULL) {
          AZ_printf_out("the pointer %u was not found and thus can not be realloc.\n",
                  ptr);
          exit(1);
       }
       else {
           /* check to see if the header is corrupted */
           iptr = (int *) ptr;
           header_start = (char *) &(dptr[1]);

           for (i = 0 ; i < 3*sizeof(double)/sizeof(char) ; i++ ) {
              if (header_start[i] != 'x') {
                 AZ_printf_out("realloc header is corrupted for %u (%d,%d)\n",ptr,
                         current->size,current->order);
                 size =  0;
                 size = 1/size;
              }
/* DO WE CHECK THE TRAILER ???? */
           }
           size = iptr[0];


    newmsize = new_size + 7*sizeof(double);
    new_ptr = (char *) malloc(newmsize);
    if (new_ptr == NULL) return(NULL);


    new_size_ptr = (int *) new_ptr;
    new_size_ptr[0] = new_size;
    new_dptr     = (double *) new_ptr;
data2 = (char *) &(new_dptr[4]);

    new_header_start = (char *) &(new_dptr[1]);

    for (i = 0 ; i < 3*sizeof(double)/sizeof(char) ; i++ )
       new_header_start[i] = 'x';

    new_header_end = &(new_ptr[ (newmsize/sizeof(char)) - 1]);
    new_header_start = (char *) &(new_dptr[4]);
    new_header_start= &(new_header_start[new_size/sizeof(char)]);

    while (new_header_start <= new_header_end) {
       *new_header_start = 'x';
       new_header_start++;
    }

    smaller = current->size;
    if (smaller > new_size ) smaller = new_size;
    for (i = 0 ; i < smaller; i++) data2[i] = data1[i];

    free(dptr);
    current->size  = new_size;
    current->address = (char *) new_dptr;
    return( (char *) &(new_dptr[4]));

       }
    }
}
   




void spit_it_out()
{
AZ_printf_out("malloc/free %d %d\n",allo_count,free_count);
if (allo_count != free_count) 
AZ_printf_out("WHOA XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
AZ_print_it();
}

#else
/* Simple wrappers for malloc */

char *AZ_allocate(unsigned int size) {
char *temp;
temp = malloc (size);
if (temp != NULL) allo_count++;
return ( temp );
}

void AZ_free(void *ptr) {
free_count++; free(ptr); }

char *AZ_realloc(void *ptr, unsigned int size) {
   return( realloc(ptr, size) );
}
extern void spit_it_out(void);

void spit_it_out()
{
if (allo_count != free_count) 
   AZ_printf_out("malloc/free %d %d XXXXXXXX\n",allo_count,free_count);
else AZ_printf_out("malloc/free %d %d\n",allo_count,free_count);
}

#endif

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
void AZ_perror(char *string)
{
/* Print out 'string' and exit */

   (void) AZ_printf_err("%s",string);
   exit(-1);
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void AZ_lower_tsolve(double *x, int n, double *val, int *bindx, int *iu, 
                    double *y)
{
/*
 *  Do a lower triangular MSR backsolve (I can't remember whether
 *  this is an MSR C or Fortran matrix).
 *
 */
    double dtmp;
    int    i, j, k, i1, i2;

    --y;
    --iu;
    --x;
    --val;
    --bindx;

    i1 = n;
    for (i = 1; i <= i1; ++i) {
        /*  Loop over the columns, up to the diagonal. */
	dtmp = 0.;
	i2 = iu[i] - 1;
	for (k = bindx[i]; k <= i2; ++k) {
	    j = bindx[k];
	    dtmp += val[k] * y[j];
	}
        /*     compute x(i). l(i,i) = 1.0. */
	y[i] = x[i] - dtmp;
    }
} /* AZ_lower_tsolve */

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void AZ_upper_tsolve(double *x, int n, double *val, int *bindx, int *iu)
{
/*
 *  Do an upper triangular MSR backsolve (I can't remember whether
 *  this is an MSR C or Fortran matrix).
 *
 */
    double dtmp;
    int   i, j, k, i1;

    /* Parameter adjustments */
    --iu;
    --x;
    --val;
    --bindx;

    /* Function Body */
    for (i = n; i >= 1; --i) {
	dtmp = 0.;
	i1 = bindx[i + 1] - 1;
	for (k = iu[i]; k <= i1; ++k) {
	    j = bindx[k];
	    dtmp += val[k] * x[j];
	}
	x[i] = (x[i] - dtmp) * val[i];
    }
} /* AZ_upper_tsolve */

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

int AZ_set_solver_parameters(double *params, int *options, AZ_MATRIX *A, 
	AZ_PRECOND *P, struct AZ_SCALING *S)
{
/*
 * Record all the parameters to AZ_oldsolve() in a manage_memory
 * data field so that they can be recovered later.
 *
 */
  static int name = AZ_SOLVER_PARAMS;
  double *iparams, *istatus;
  int    *ioptions, i, i2, i3;
  char   str[80];
  AZ_MATRIX  *Aptr;
  AZ_PRECOND *Pptr;
  struct AZ_SCALING *Sptr;
  unsigned   int size;

  name--;

  if ((options[AZ_conv] == AZ_weighted) || (options[AZ_conv] == AZ_expected_values)) {
     AZ_printf_err("Sorry weighted norm can not be used with this routine.\n");
     exit(1);
  }


   /* allocate each one of these things with manage memory */

  size = AZ_get_sol_param_size();


  sprintf(str,"sol_param %d",name);
  iparams = AZ_manage_memory(size, AZ_ALLOC, -777, str, &i);
  istatus = &(iparams[AZ_PARAMS_SIZE]);
  Aptr    =  (AZ_MATRIX *) &(istatus[AZ_STATUS_SIZE]);
  i       = (int) ( ((double) sizeof(AZ_MATRIX)) / ((double) sizeof(double) ));
  i++;
  Pptr    =  (AZ_PRECOND *) &(istatus[AZ_STATUS_SIZE+i]);
  i2      = (int) ( ((double) sizeof(AZ_PRECOND)) / ((double) sizeof(double) ));
  i2++;
  ioptions= (int *)  &(istatus[AZ_STATUS_SIZE+i + i2]);
  i3      = (int) ( ((double) AZ_OPTIONS_SIZE*sizeof(int) ) / ((double) sizeof(double) ));
  Sptr    =  (struct AZ_SCALING *) &(istatus[AZ_STATUS_SIZE+i + i2 + i3]);

  *Aptr = *A;
  *Pptr = *P;
  *Sptr = *S;

  for (i = 0 ; i < AZ_PARAMS_SIZE ; i++ ) iparams[i]  = params[i];
  for (i = 0 ; i < AZ_OPTIONS_SIZE; i++ ) ioptions[i] = options[i];
  for (i = 0 ; i < AZ_STATUS_SIZE ; i++ ) istatus[i] = 0.0;

  return(name);
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/
extern unsigned int AZ_get_sol_param_size(void); 
unsigned int AZ_get_sol_param_size() 
{
  /* Utility routine which just gives us the size of all Aztec's */
  /* solution parameters.                                        */


  return((AZ_OPTIONS_SIZE)*sizeof(int) +
         (AZ_PARAMS_SIZE + AZ_STATUS_SIZE + 2)*sizeof(double) +
         sizeof(AZ_MATRIX)  + sizeof(AZ_PRECOND) + sizeof(struct AZ_SCALING) );
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void AZ_recover_sol_params(int instance, int **sub_options, double **sub_params,
	double **sub_status, AZ_MATRIX **sub_matrix, AZ_PRECOND **sub_precond,
        struct AZ_SCALING **sub_scaling)
{
/* This routine will require the previously stored parameters associated
 * with 'instance' (see AZ_set_solver_parameters). These parameters correspond
 * to all the variables need to invoke AZ_oldsolve(). This routine is used
 * to recursively call the iteratve solvers or to use the iterative solvers
 * on a subdomain.
 */

  double     *tparams, *tstatus;
  int        *toptions, i, i2, i3;
  char       str[80];
  AZ_MATRIX  *Aptr;
  AZ_PRECOND *Pptr;
  struct AZ_SCALING *scaling;
  
  unsigned   int size;

  size = AZ_get_sol_param_size();

  sprintf(str,"sol_param %d",instance);
  tparams = AZ_manage_memory(size, AZ_ALLOC, -777, str, &i);
  if (i == AZ_NEW_ADDRESS) {
     AZ_printf_err("Error:\tSolver parameters corresponding to ");
     AZ_printf_err("the internal solver = %d\n\twere not found.\n",
             instance);
     exit(1);
  }
  tstatus = &(tparams[AZ_PARAMS_SIZE]);
  Aptr    =  (AZ_MATRIX *) &(tstatus[AZ_STATUS_SIZE]);
  i       = (int) ( ((double) sizeof(AZ_MATRIX))/((double) sizeof(double)));
  i++;
  Pptr    =  (AZ_PRECOND *) &(tstatus[AZ_STATUS_SIZE+i]);
  i2      = (int) ( ((double) sizeof(AZ_PRECOND))/((double) sizeof(double)));
  i2++;
  toptions= (int *)  &(tstatus[AZ_STATUS_SIZE+i + i2]);
  i3      = (int) ( ((double) AZ_OPTIONS_SIZE*sizeof(int) ) / ((double) sizeof(double) ));
  scaling =  (struct AZ_SCALING *) &(tstatus[AZ_STATUS_SIZE+i + i2 + i3]);

  *sub_options = toptions;
  *sub_params  = tparams;
  *sub_status  = tstatus;
  *sub_matrix  = Aptr;
  *sub_precond = Pptr;
  *sub_scaling = scaling;
}

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

void AZ_clear_solver_parameters(int handle)
{
/*
 * Clear the context assoicated with handle.
 *
 */
  int    i;
  char   str[80];
  unsigned int size;


  size = AZ_get_sol_param_size();
  sprintf(str,"sol_param %d",handle);
  AZ_manage_memory(size, AZ_SELECTIVE_CLEAR, -777, str, &i);
}
void AZ_set_precond_print_string(struct AZ_PREC_STRUCT *precond, const char str[])
{
   if ( precond->print_string != NULL) AZ_free( precond->print_string);
   precond->print_string = (char *) AZ_allocate( (strlen(str)+1)*sizeof(char));
   if (precond->print_string == NULL) {
      AZ_printf_out("AZ_set_precond_print_string: Not enough space to allocate string\n");
      exit(1);
   }
   sprintf(precond->print_string,"%s",str);
}

void AZ_set_matrix_print_string(AZ_MATRIX *Amat, const char str[])
{
   if ( Amat->print_string != NULL) AZ_free( Amat->print_string);
   Amat->print_string = (char *) AZ_allocate( (strlen(str)+1)*sizeof(char));
   if (Amat->print_string == NULL) {
      AZ_printf_out("AZ_set_matrix_print_string: Not enough space to allocate string\n");
      exit(1);
   }
   sprintf(Amat->print_string,"%s",str);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
AZ_PRECOND *AZ_precond_create(AZ_MATRIX *Pmat, void (*prec_fun)(
	double *, int *, int *, double *, struct AZ_MATRIX_STRUCT  *,
               struct AZ_PREC_STRUCT *), void *data)
{
   AZ_PRECOND *precond;

   precond = (AZ_PRECOND *) AZ_allocate(sizeof(AZ_PRECOND));
   if (precond == NULL) {
      AZ_printf_out("Error: Not enough space in AZ_precond_create().\n");
      exit(1);
   }
   precond->Pmat          = Pmat;
   precond->prec_function = prec_fun;
   precond->precond_data  = data;
   precond->prec_create_called = 1;
   precond->options = NULL;
   precond->params  = NULL;
   precond->print_string = NULL;
/*
precond->options = (int    *) AZ_allocate(AZ_OPTIONS_SIZE*sizeof(int));
precond->params  = (double *) AZ_allocate(AZ_PARAMS_SIZE*sizeof(double));
AZ_defaults(precond->options,precond->params);
   if (prec_fun != AZ_precondition) precond->options[AZ_precond] = AZ_user_precond;
*/
   return(precond);
}


AZ_MATRIX *AZ_matrix_create(int local)
{
/*
 * Create an Aztec AZ_MATRIX structure and fill in the noncommunication
 * related fields of data_org[]. 
 * Note: This matrix will not work properly with Aztec's AZ_exchange_bdry()
 *       subroutine. Instead, it is intended that this matrix be used for
 *       matrix-free users and matrices which do not require communication.
 *
 * Parameters
 * ========
 *   local              Number of matrix equations residing on this processor.
 *
 *   additional         local+additional is the required size of a vector, x,
 *                      that will be applied to a matrix when performing a
 *                      matrix-vector product. The first 'local' components of
 *                      'x' must contain the appropriate data. The remaining
 *                      'additional' components are used as workspace inside
 *                      the user's matrix vector product.
 *  matrix_type         Either AZ_MSR_MATRIX, AZ_VBR_MATRIX, or AZ_USER_MATRIX.
 *  local_blks          When matrix_type == AZ_VBR_MATRIX, 'local_blks' 
 *                      indicates how many block equations reside on this node.
 */

   AZ_MATRIX  *Amat;

   Amat     = (AZ_MATRIX *) AZ_allocate(sizeof(AZ_MATRIX));
   if (Amat == NULL) {
      AZ_printf_out("Error: Not enough space in AZ_matrix_create().\n");
      exit(1);
   }
   AZ_matrix_init(Amat, local);

   return(Amat);
} 
void AZ_matrix_init(AZ_MATRIX *Amat, int local)
{
/*
 * Create an Aztec AZ_MATRIX structure and fill in the noncommunication
 * related fields of data_org[]. 
 * Note: This matrix will not work properly with Aztec's AZ_exchange_bdry()
 *       subroutine. Instead, it is intended that this matrix be used for
 *       matrix-free users and matrices which do not require communication.
 *
 * Parameters
 * ========
 *   local              Number of matrix equations residing on this processor.
 *
 *   additional         local+additional is the required size of a vector, x,
 *                      that will be applied to a matrix when performing a
 *                      matrix-vector product. The first 'local' components of
 *                      'x' must contain the appropriate data. The remaining
 *                      'additional' components are used as workspace inside
 *                      the user's matrix vector product.
 *  matrix_type         Either AZ_MSR_MATRIX, AZ_VBR_MATRIX, or AZ_USER_MATRIX.
 *  local_blks          When matrix_type == AZ_VBR_MATRIX, 'local_blks' 
 *                      indicates how many block equations reside on this node.
 */

/*
   AZ_MATRIX  *Amat;

   Amat     = (AZ_MATRIX *) AZ_allocate(sizeof(AZ_MATRIX));
   if (Amat == NULL) {
      AZ_printf_out("Error: Not enough space in AZ_matrix_create().\n");
      exit(1);
   }
*/
   Amat->matrix_type = AZ_none;
   Amat->N_local = local;
   Amat->N_ghost = 0;
   Amat->mat_create_called = 1;
   Amat->must_free_data_org = 0;
   Amat->rpntr       = NULL;
   Amat->cpntr       = NULL;
   Amat->bpntr       = NULL;
   Amat->bindx       = NULL;
   Amat->indx        = NULL;
   Amat->val         = NULL;
   Amat->data_org    = NULL;
   Amat->matvec      = NULL;
   Amat->getrow      = NULL;
   Amat->user_comm   = NULL;
   Amat->matrix_norm = -1.0;
   Amat->aux_ival    = NULL;
   Amat->aux_dval    = NULL;
   Amat->aux_ptr     = NULL;
   Amat->matvec_data = NULL;
   Amat->getrow_data = NULL;
   Amat->aux_matrix  = NULL;
   Amat->N_nz        = -1;
   Amat->max_per_row = -1;
   Amat->largest_band= -1;
   Amat->print_string = NULL;
} 
/* Begin Aztec 2.1 mheroux mod */
void AZ_set_MSR(AZ_MATRIX *Amat, int bindx[], double val[], int
                data_org[], int N_update, int update[], int option)
/* End Aztec 2.1 mheroux mod */
{
   Amat->bindx    = bindx;
   Amat->val      = val;
   Amat->data_org = data_org;
   Amat->N_update = N_update;
   Amat->update = update;
   Amat->matrix_type = AZ_MSR_MATRIX;
   Amat->matvec   = AZ_MSR_matvec_mult;
   Amat->N_nz  = -1;
   Amat->max_per_row = -1;
   Amat->largest_band= -1;
   if (option == AZ_LOCAL)
     Amat->has_global_indices = 0;
   else
     Amat->has_global_indices = 1;
}
/* Begin Aztec 2.1 mheroux mod */
void AZ_set_VBR(AZ_MATRIX *Amat, int rpntr[], int cpntr[], int bpntr[], 
		int indx[], int bindx[], double val[], int data_org[], 
		int N_update, int update[], int option)
/* End Aztec 2.1 mheroux mod */
{
   Amat->rpntr = rpntr;
   Amat->cpntr = cpntr;
   Amat->bpntr = bpntr;
   Amat->indx  = indx;
   Amat->bindx = bindx;
   Amat->val   = val;
   Amat->data_org = data_org;
   Amat->N_update = N_update;
   Amat->update = update;
   Amat->matrix_type = AZ_VBR_MATRIX;
   Amat->matvec   = AZ_VBR_matvec_mult;
   Amat->N_nz  = -1;
   Amat->max_per_row = -1;
   Amat->largest_band= -1;
   if (option == AZ_LOCAL)
     Amat->has_global_indices = 0;
   else
     Amat->has_global_indices = 1;
}
void AZ_matrix_destroy(AZ_MATRIX **Amat)
{
  if ( (*Amat) == NULL) return;
   if ((*Amat)->must_free_data_org == 1) {
      AZ_free((*Amat)->data_org);
      (*Amat)->data_org = NULL;
   }
   if ( (*Amat)->print_string != NULL) AZ_free( (*Amat)->print_string);
   AZ_free(*Amat);
   *Amat = NULL;
}
void AZ_precond_destroy(AZ_PRECOND **precond)
{
  if ( (*precond) == NULL) return;
   if ( (*precond)->print_string != NULL) AZ_free( (*precond)->print_string);
   AZ_free(*precond);
   *precond = NULL;
}

void AZ_set_MATNORMINF(AZ_MATRIX *Amat, void* data, double (*matnorminf)(AZ_MATRIX* Amat))
{
  Amat->matnorminf = matnorminf;
}

void AZ_set_MATFREE_matrix_norm(AZ_MATRIX *Amat, double mat_norm)
{
   Amat->matrix_norm = mat_norm;
}
void AZ_set_MATFREE_name(AZ_MATRIX *Amat, int name)
{
   Amat->data_org[AZ_name] = name;
}
void AZ_set_MATFREE(AZ_MATRIX *Amat, void *data, 
    void (*matvec)(double *, double *, struct AZ_MATRIX_STRUCT *, int *))
{
   static int name = 2071;
   int    *data_org;

   if ( Amat->getrow != NULL) {
      AZ_printf_out("Error: AZ_set_MATFREE must be called before AZ_set_MATFREE_getrow\n");
      exit(1);
   }
   Amat->N_nz  = -1;
   Amat->max_per_row = -1;
   Amat->largest_band= -1;
   if (Amat->data_org == NULL) {
      Amat->must_free_data_org = 1;
      data_org = (int  *) AZ_allocate(AZ_COMMLESS_DATA_ORG_SIZE*sizeof(int));
      if (data_org == NULL) {
         AZ_printf_out("Error: Not enough space in AZ_create_matrix().\n");
         exit(1);
      }
      Amat->data_org = data_org;
   }
   else data_org = Amat->data_org;

   name++; 

   data_org[AZ_N_internal]  = 0;    /* Number of rows without data         */
                                    /* dependencies on other procs         */
   data_org[AZ_N_border  ]  = Amat->N_local;
   data_org[AZ_N_external]  = Amat->N_ghost;
   data_org[AZ_name]        = name;
   data_org[AZ_matrix_type] = AZ_USER_MATRIX;
   Amat->matrix_type = AZ_USER_MATRIX;
   data_org[AZ_N_int_blk]   = data_org[AZ_N_internal];
   data_org[AZ_N_bord_blk]  = data_org[AZ_N_border];
   data_org[AZ_N_ext_blk]   = data_org[AZ_N_external];
   data_org[AZ_N_neigh]     = 0;
   data_org[AZ_total_send]  = 0;
   Amat->matvec      = matvec;
   Amat->matvec_data = data;
}
void AZ_set_MATFREE_getrow(AZ_MATRIX *Amat, void *data,
    int  (*getrow)(int *, double *, int *, struct AZ_MATRIX_STRUCT *, 
		   int , int *, int),
    int  (*user_comm)(double *, AZ_MATRIX *), int N_ghost, int proc_config[])
{
   int    *data_org;

   Amat->N_nz  = -1;
   Amat->max_per_row = -1;
   Amat->largest_band= -1;
   Amat->getrow_data = data;
   Amat->N_ghost     = N_ghost;
   AZ_extract_comm_info(&data_org, user_comm, Amat, proc_config,
   			Amat->N_local, Amat->N_ghost);
   Amat->must_free_data_org = 1;
   Amat->user_comm = user_comm;

   data_org[AZ_N_internal]  = 0;    /* Number of rows without data         */
                                    /* dependencies on other procs         */
   data_org[AZ_N_border  ]  = Amat->N_local;
   data_org[AZ_N_external]  = Amat->N_ghost;
   data_org[AZ_matrix_type] = AZ_USER_MATRIX;
   data_org[AZ_N_int_blk] = data_org[AZ_N_internal];
   data_org[AZ_N_bord_blk]= data_org[AZ_N_border];
   data_org[AZ_N_ext_blk] = data_org[AZ_N_external];
   Amat->matrix_type = AZ_USER_MATRIX;
   Amat->getrow = getrow;
   Amat->user_comm = user_comm;
   if (Amat->data_org != NULL) {
      data_org[AZ_name]        = Amat->data_org[AZ_name];
      AZ_free(Amat->data_org);
   }
   Amat->data_org = data_org;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

AZ_MATRIX *AZ_create_matrix(int local, int additional, int matrix_type,
				    int local_blks, int *idata_org)
{
/*
 * Create an Aztec AZ_MATRIX structure and fill in the noncommunication
 * related fields of data_org[]. 
 * Note: This matrix will not work properly with Aztec's AZ_exchange_bdry()
 *       subroutine. Instead, it is intended that this matrix be used for
 *       matrix-free users and matrices which do not require communication.
 *
 * Parameters
 * ========
 *   local              Number of matrix equations residing on this processor.
 *
 *   additional         local+additional is the required size of a vector, x,
 *                      that will be applied to a matrix when performing a
 *                      matrix-vector product. The first 'local' components of
 *                      'x' must contain the appropriate data. The remaining
 *                      'additional' components are used as workspace inside
 *                      the user's matrix vector product.
 *  matrix_type         Either AZ_MSR_MATRIX, AZ_VBR_MATRIX, or AZ_USER_MATRIX.
 *  local_blks          When matrix_type == AZ_VBR_MATRIX, 'local_blks' 
 *                      indicates how many block equations reside on this node.
 */

   static int name = 1071;
   int        *data_org;
   AZ_MATRIX  *Amat;


   name++;
   Amat     = (AZ_MATRIX *) AZ_allocate(sizeof(AZ_MATRIX));
   if (Amat == NULL) {
      AZ_printf_out("Error: Not enough space in AZ_create_matrix().\n");
      exit(1);
   }
   if (idata_org == AZ_NOT_USING_AZTEC_MATVEC) {
      data_org = (int  *) AZ_allocate(AZ_COMMLESS_DATA_ORG_SIZE*sizeof(int));
      if (data_org == NULL) {
         AZ_printf_out("Error: Not enough space in AZ_create_matrix().\n");
         exit(1);
      }
 /* mheroux change */
      Amat->must_free_data_org = 1; /* AZ_matrix_destroy should delete 
                                       data_org */

      data_org[AZ_N_internal] = 0;    /* Number of rows without data         */
                                      /* dependencies on other procs         */
      data_org[AZ_N_border  ] = local;/* Matrix rows with data dependencies  */
                                      /* on other procs                      */
      data_org[AZ_N_external] = 
                          additional;/* When doing y = A x (i.e. matrix-     */
                                     /* vector products) user's can request  */
                                     /* that additional space be passed in   */
                                     /* with the vector 'x'. This additional */
                                     /* space can be used as work space or   */
                                     /* to hold ghost variables when doing   */
                                     /* parallel computations.               */
                                     /* If no additional space is needed,    */
                                     /* this should be set to zero.          */
         
      data_org[AZ_name]       = name;/* A number used to 'name' the matrix  */
      data_org[AZ_matrix_type] = matrix_type;

      if (matrix_type != AZ_VBR_MATRIX) {
         data_org[AZ_N_int_blk] = data_org[AZ_N_internal];
         data_org[AZ_N_bord_blk]= data_org[AZ_N_border];
         data_org[AZ_N_ext_blk] = data_org[AZ_N_external];
      }
      else {
         data_org[AZ_N_int_blk] = data_org[AZ_N_internal];
         data_org[AZ_N_bord_blk]= local_blks;
         data_org[AZ_N_ext_blk] = data_org[AZ_N_external];
      }
      data_org[AZ_N_neigh] = 0;
      data_org[AZ_total_send] = 0;

      Amat->data_org = data_org;
   }
   else Amat->data_org = idata_org;

   Amat->matrix_type = matrix_type;
   Amat->matrix_norm = -1.0;
   Amat->mat_create_called = 1;
   return(Amat);
} 

void AZ_delete_matrix(AZ_MATRIX *ptr)
{
   AZ_free(ptr->data_org);
   AZ_free(ptr);
}


int AZ_compress_msr(int *ibindx[], double *ival[], int allocated, int needed,
		  int name, struct context *context)
{
/*
 * We can't seem to rely on realloc() on the paragon or the tflop
 * machine. So let's try and do the best job using malloc and 
 * free to reduce the amount of space required by val and bindx.
 *
 */

   int bindx_length, val_length, target_bindx, target_val;
   int extra_size, *bind2, *bindx;
   double *ptr1,  *extra, *val, combined_bindx_val = 0;
   int i, modified_need, modified_target, j;
   char label[200];

   if (needed == 0) return(0);
   bindx_length = allocated*sizeof(int);
   val_length   = allocated*sizeof(double);
   sprintf(label,"val %s",context->tag);
   val   = AZ_manage_memory((unsigned) val_length,AZ_REALLOC,name,label, &j);
   if (j == AZ_SPECIAL) return(1);

   bindx = *ibindx;
   val   = *ival;
   target_bindx = needed*sizeof(int);
   target_val   = needed*sizeof(double);
   modified_target = target_bindx;

   i = target_bindx%sizeof(double);
   modified_target += (sizeof(double) - i);    /* add something so things */
   modified_need = modified_target/sizeof(int);/* fall on double bdries   */

   extra_size  = target_val + modified_target - bindx_length + 4*sizeof(double);
                                                        /* 4sizeof(double) */
                                                        /* added for       */
                                                        /* manage_mem().   */

   if ( extra_size <= 0) {   

      /* There is enough room to store both bindx and val */
      /* in the bindx array. So we copy val into bindx,   */
      /* free val, then allocate a new bindx and val      */
      /* (whose length should be less than the old val),  */
      /* store the bindx/val information in these new     */
      /* arrays and then free the old bindx array.        */

      ptr1 = (double *) &(bindx[modified_need]);
      for (i = 0 ; i < needed; i++) ptr1[i] = val[i];

      sprintf(label,"val %s",context->tag);
      AZ_manage_memory((unsigned)val_length,AZ_SELECTIVE_CLEAR,name,label, &j);
      val   = AZ_manage_memory((unsigned)target_val,AZ_ALLOC,name,label, &j);

      sprintf(label,"bind2xx %s",context->tag);
      bind2 = (int *) AZ_manage_memory((unsigned)target_bindx,AZ_ALLOC,name,label, &j);

      for (i = 0 ; i < needed; i++) val[i]   = ptr1[i];
      for (i = 0 ; i < needed; i++) bind2[i] = bindx[i];

      sprintf(label,"bindx %s",context->tag);
      AZ_manage_memory((unsigned)bindx_length,AZ_SELECTIVE_CLEAR,name,label, &j);
      bindx = bind2;
      sprintf(label,"bind2xx %sbindx %s",context->tag,context->tag);
      AZ_manage_memory((unsigned)target_bindx,AZ_RESET_STRING,name,label, &j);
   }
   else {
      i = extra_size%sizeof(double);
      extra_size += (sizeof(double) - i); /* add something so things */
                                          /* fall on double bdries   */
      extra = (double *) AZ_allocate((unsigned) extra_size);
      if (extra != NULL) {

         /* There is enough room to store both val into    */
         /* bindx and extra. So we copy val into bindx and */
         /* extra, free val, then allocate a new val and   */
         /* fill the new val. Once this is done, we free   */
         /* extra (the length of extra and the savings in  */
         /* deleting the old val should be enough to hold  */
         /* bindx), allocate a new bindx, fill the new     */
         /* bindx and free the old bindx.                  */

           extra_size /= sizeof(double);
           if (extra_size > needed) extra_size = needed;
           for (i = 0 ; i < extra_size; i++) extra[i] = val[i];
           ptr1 = (double *) &(bindx[modified_need]);
           for (i = extra_size ; i < needed; i++) 
              ptr1[i-extra_size] = val[i];

           sprintf(label,"val %s",context->tag);
           AZ_manage_memory((unsigned) val_length,AZ_SELECTIVE_CLEAR,name,label, &j);
           val   = AZ_manage_memory((unsigned)target_val,AZ_ALLOC,name,label, &j);

           for (i = 0 ; i < extra_size; i++) val[i]   = extra[i];
           for (i = extra_size ; i < needed; i++) 
              val[i] = ptr1[i-extra_size];
           AZ_free(extra);
           bind2 = (int *) AZ_allocate(modified_target+ 4*sizeof(double));
           if (bind2 != NULL) {
              AZ_free(bind2);
              sprintf(label,"bind2xx %s",context->tag);
              bind2 = (int *) AZ_manage_memory((unsigned) target_bindx,AZ_ALLOC,name,
                                               label, &j);
              for (i = 0 ; i < needed; i++) bind2[i] = bindx[i];

              sprintf(label,"bindx %s",context->tag);
              AZ_manage_memory((unsigned)bindx_length,AZ_SELECTIVE_CLEAR,name,label, &j);
              bindx = bind2;
              sprintf(label,"bind2xx %sbindx %s",context->tag,context->tag);
              AZ_manage_memory((unsigned)target_bindx,AZ_RESET_STRING,name,label, &j);
           }
           else {
              sprintf(label,"bindx %s",context->tag);
              bindx = (int *) AZ_manage_memory((unsigned)target_bindx, AZ_REALLOC,
                                               name,label, &j);
           }
      }
      else if ( target_val + modified_target <= val_length) {

         /* There is enough room to store both val into    */
         /* bindx in val. So we copy bindx to val and free */
         /* the space associated with bindx.               */
         /* Note: we need to mark this space as not being  */
         /* really allocated (i.e. we should not free      */
         /* bindx).                                        */

         sprintf(label,"val %s",context->tag);
         val = AZ_manage_memory((unsigned)(target_val+target_bindx), AZ_SPEC_REALLOC,
                                name,label, &j);
         combined_bindx_val = 1;
         bind2 = (int *) &(val[needed]);
         for (i = 0 ; i < needed ; i++) bind2[i] = bindx[i];

         sprintf(label,"bindx %s",context->tag);
         AZ_manage_memory((unsigned) bindx_length,AZ_SELECTIVE_CLEAR,name,label, &j);
         bindx = bind2;
      }
      else {
         sprintf(label,"val %s",context->tag);
         val = AZ_manage_memory((unsigned)target_val, AZ_REALLOC,name,label, &j);
         sprintf(label,"bindx %s",context->tag);
         bindx = (int *) AZ_manage_memory((unsigned)target_bindx,AZ_REALLOC,
                                          name,label,&j);
      }
   }
   *ibindx = bindx;
   *ival   = val;
   return(combined_bindx_val);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void AZ_loc_avg(AZ_MATRIX *Amat, double r[], double newr[], int Nfixed,
	int fixed_pts[], int proc_config[])

{
/*
 * Do some residual smoothing based on the connectivity given by the
 * VBR matrix. The new value at a point is given by taking the average
 * of the neighbors and averaging this with the previous value at the point.
 *
 * Parameters
 * ==========
 *
 * Nrows                   On input, number of rows in the matrix
 *
 * N_blk_rows              On input, number of block rows in the matrix
 *
 * cpntr, bindx, bpntr     On input, VBR matrix giving connectivity
 *
 * r                       On input, array to be averaged.
 *
 * newr                    On output, result of averaging applied to 'r'.
 *
 */


   int i,j, ii, block_col, cardinality, offset;
   double scale_factor;
   int *data_org, *cpntr, *bpntr, *bindx, Nrows, N_blk_rows;
   int flag;

   data_org = Amat->data_org;
   Nrows    = data_org[AZ_N_internal] + data_org[AZ_N_border];
   bindx    = Amat->bindx;

   AZ_exchange_bdry(r, data_org, proc_config);
   for (ii = 0; ii < Nrows; ii++)  newr[ii] = 0.0;

   if (Amat->matrix_type == AZ_VBR_MATRIX) {
      N_blk_rows = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];
      cpntr      = Amat->cpntr;
      bpntr      = Amat->bpntr;

      for ( i = 0; i < N_blk_rows; i++) {
         cardinality = bpntr[i+1] - bpntr[i] - 1;
         if (cardinality != 0) {
            scale_factor = .5/((double) cardinality);

            for (j = bpntr[i]; j < bpntr[i+1] ; j++) {
               block_col = bindx[j];
               if (block_col != i) {
                  offset = cpntr[i];
                  block_col = cpntr[block_col];
                  for (ii = 0; ii < cpntr[i+1] - cpntr[i]; ii++)
                     newr[offset+ii] += r[block_col+ii];
               }
            }
            for (ii = cpntr[i]; ii < cpntr[i+1]; ii++)  newr[ii] *=scale_factor;
            for (ii = cpntr[i]; ii < cpntr[i+1]; ii++)  newr[ii] += (.5*r[ii]);
         }
      }
   }
   else if (Amat->matrix_type == AZ_MSR_MATRIX) {

      for (i = 0; i < Nrows; i++) {
         cardinality = bindx[i+1] - bindx[i];
         if (cardinality != 0) {
            scale_factor = .5/((double) cardinality);

            for (j = bindx[i]; j < bindx[i+1] ; j++) 
               newr[i] += r[bindx[j]];

            newr[i] *=scale_factor;
            newr[i] += (.5*r[i]);
         }
      }
   }
   else {
      AZ_printf_out("Smoothing can only be done with MSR or VBR matrices\n");
      exit(1);
   }

   /* fix the boundary */

   flag = 0;
   for (i = 0; i < Nfixed; i++ ) {
      if ((flag==0) && (fabs(r[fixed_pts[i]]) > 1.0e-9)) {
         AZ_printf_out("boundary not zero %e\n",r[fixed_pts[i]]);
         flag = 1;
      }
      newr[fixed_pts[i]] = r[fixed_pts[i]];
   }
}
void AZ_matfree_Nnzs(AZ_MATRIX *Amat)
{
   int *data_org, N, space, *cols, i, flag, Nnz, max_length, length;
   int j, smallest, largest, bandwidth = 0;
   double *val;

   data_org = Amat->data_org;

   Nnz = 0;
   max_length = 0;
   N     = data_org[AZ_N_internal] + data_org[AZ_N_border];
   if ( (Amat->getrow == NULL) && (N != 0) ) {
      AZ_printf_out("Error: Only matrices with getrow() defined via ");
      AZ_printf_out("AZ_set_MATFREE_getrow(...)\n       can compute");
      AZ_printf_out(" their nonzeros information.\n");
      exit(1);
   }
   space = 30;
   cols    = (int    *) AZ_allocate(sizeof(int)*space);
   val    = (double *) AZ_allocate(sizeof(double)*space);
   if (val == NULL) {
      AZ_printf_out("AZ_matfree_Nnzs: Out of space. Requested %d.\n",space);
      exit(1);
   }
  
   for (i = 0; i < N; i++) {
      flag  = 0;
      while (flag == 0) {
         flag = Amat->getrow( cols, val, &length, Amat, 1, &i, space);
         if (flag == 0) {
            AZ_free(val);
            AZ_free(cols);
            space  = (int) ( ((double)space)*1.5 + 3);
            cols   = (int    *) AZ_allocate(sizeof(int)*space);
            val    = (double *) AZ_allocate(sizeof(double)*space);
            if (val  == NULL) {
               AZ_printf_out("AZ_matfree_Nnzs: Out of space. Requested %d.\n",space);
               exit(1);
            }
         }
      }
      Nnz += length;
      if (max_length < length) max_length = length;
      if (length != 0) {
         smallest = cols[0];  largest = cols[0]; 
         for (j = 1; j < length; j++) {
            if (cols[j] < smallest) smallest = cols[j];
            if (cols[j] > largest ) largest = cols[j];
         }
         if (bandwidth <= largest - smallest) {
             bandwidth = largest - smallest + 1;
         }
      }
   }
   Amat->N_nz = Nnz;
   Amat->max_per_row = max_length;
   Amat->largest_band = bandwidth;
   AZ_free(val);
   AZ_free(cols);
}


void AZ_matfree_2_msr(AZ_MATRIX *Amat, double *val, int *bindx, int N_nz)
{
/*
 * Take the 'matrix' defined by Amat->getrow() and 
 * make an MSR copy of it to be stored in 'output'.
 *
 */
   int N, *data_org, i, length, flag, nz_ptr, k, j;

   if ( (Amat->N_nz < 0) || (Amat->max_per_row < 0)) 
      AZ_matfree_Nnzs(Amat);

   data_org = Amat->data_org;
   N     = data_org[AZ_N_internal] + data_org[AZ_N_border];
   if ( (Amat->getrow == NULL) && (N != 0) ) {
      AZ_printf_out("Error: Only matrices with getrow() defined via ");
      AZ_printf_out("AZ_set_MATFREE_getrow(...) can be converted to MSR \n");
      exit(1);
   }
   if (N_nz < Amat->N_nz) {
      AZ_printf_out("AZ_matfree_2_msr: Something is wrong. The number of nonzeros\n");
      AZ_printf_out("    allocated for preconditioner is less than the number of\n");
      AZ_printf_out("    nonzeros in matrix (%d vs. %d)!\n",N_nz, Amat->N_nz);
      exit(1);
   }
  
   nz_ptr = N+1;
   bindx[0] = nz_ptr;
   val[N] = 0;
   for (i = 0; i < N; i++) {
      flag = Amat->getrow(&(bindx[nz_ptr]), &(val[nz_ptr]), &length, Amat, 1, 
                   &i, N_nz);
      if (flag == 0) {
         AZ_printf_out("AZ_matfree_2_msr: Something is wrong. The number of nonzeros");
         AZ_printf_out("in matrix is\n   greater than the number of nonzeros");
         AZ_printf_out("recorded in Amat->N_nz (%d)\n",Amat->N_nz);
         exit(1);
      }
      /* move diagonal to proper location */

      for (k =0; k < length; k++) 
         if (bindx[nz_ptr + k] == i) break;
      if (k == length) val[i] = 0.0; /* no diagonal */
      else {
         val[i] = val[nz_ptr+k];
         /* move over the rest */
         for (j = nz_ptr+k+1; j < nz_ptr+length; j++) {
            bindx[j-1] = bindx[j];
            val[j-1]   = val[j];
         }
         length--;
      }
      nz_ptr += length;
      bindx[i+1] = nz_ptr;
   }
}
int  AZ_MSR_getrow(int columns[], double values[], int row_lengths[], 
	struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows, 
	int requested_rows[], int allocated_space)
{
   int    *bindx, i, j, count = 0, row;
   double *val;

   bindx = Amat->bindx;
   val   = Amat->val;

   for (i = 0; i < N_requested_rows; i++) {
      row            = requested_rows[i];
      row_lengths[i] = bindx[row+1] - bindx[row] + 1;
      if (count+row_lengths[i] > allocated_space) return(0);

      /* diagonal */

      columns[count  ] = row;
      values[count++]  = val[row];

      /* off-diagonals */

      for (j = bindx[row] ; j < bindx[row+1] ; j++) {
         columns[count  ]   = bindx[j];
         values[count++] = val[j];
      }
   }
   return(1);
}
int  AZ_VBR_getrow(int columns[], double values[], int row_lengths[], 
	struct AZ_MATRIX_STRUCT *Amat, int N_requested_rows, 
	int requested_rows[], int allocated_space)
{
   int blk_row, count, N_rows, N_cols, row_offset, i, j, ii, start, k;
   int *cpntr, *bindx, *rpntr, *indx, *bpntr, oldcount, row;
   double *val;
   int    Nblks, N;


    bindx = Amat->bindx;
    val   = Amat->val;
    rpntr = Amat->rpntr;
    cpntr = Amat->cpntr;
    bpntr = Amat->bpntr;
    indx  = Amat->indx;
    Nblks = Amat->data_org[AZ_N_int_blk] +
                Amat->data_org[AZ_N_bord_blk];
    N     = Amat->data_org[AZ_N_internal] +
                Amat->data_org[AZ_N_border];

   count = 0;
   for (k = 0; k < N_requested_rows; k++) {
      oldcount = count;
      row = requested_rows[k];
      blk_row = row * Nblks;
      blk_row = blk_row/N;
      while (rpntr[blk_row] < row) blk_row++;
      while (rpntr[blk_row] > row) blk_row--;

      N_rows = rpntr[blk_row+1] - rpntr[blk_row];
      row_offset = row - rpntr[blk_row];

      for (i = bpntr[blk_row]; i < bpntr[blk_row+1]; i++) {
         ii = bindx[i];
         N_cols = cpntr[ii+1] - cpntr[ii];
         if (count+N_cols > allocated_space) return(0);
         start  = indx[i] + row_offset;

         for (j = 0; j < N_cols; j++) {
            columns[count] = j + cpntr[ii];
            values[count++] = val[start];
            start += N_rows;
         }
      }
      row_lengths[k] = count - oldcount;
   }
   return(1);
}
void AZ__MPI_comm_space_ok() {
   if (sizeof(MPI_AZComm)/sizeof(int) > SIZEOF_MPI_AZCOMM) {
      AZ_printf_out("AZ__MPI_comm_space_ok: Recompile all Aztec files with -DSIZEOF_MPI_AZCOMM=%d\n",(int) (sizeof(MPI_AZComm)/sizeof(int) + 1));
      exit(1);
   }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void AZ_version(char str[])
{
/*
 * Print out the version number of Aztec.
 *
 * Parameters
 * ==========
 *
 * str               On input, must have space for at least 15 characters.
 *                   On output, contains the version number for Aztec.
 *
 */
   sprintf(str,"Aztec ?.?.?.?");
#ifdef AZ_ver2_1_0_3
   sprintf(str,"Aztec 2.1.0.3");
#endif
#ifdef AZ_ver2_1_0_4
   sprintf(str,"Aztec 2.1.0.4");
#endif
#ifdef AZ_ver2_1_0_5
   sprintf(str,"Aztec 2.1.0.5");
#endif
#ifdef AZ_ver2_1_0_6
   sprintf(str,"Aztec 2.1.0.6");
#endif
#ifdef AZ_ver2_1_0_7
   sprintf(str,"Aztec 2.1.0.7");
#endif
#ifdef AZ_ver2_1_0_8
   sprintf(str,"Aztec 2.1.0.8");
#endif
#ifdef AZ_ver2_1_0_9
   sprintf(str,"Aztec 2.1.0.9");
#endif
#ifdef AZ_ver2_1_0_10
   sprintf(str,"Aztec 2.1.0.10");
#endif
#ifdef AZ_ver2_1_1_0
   sprintf(str,"Aztec 2.1.1.0");
#endif
#ifdef AZ_ver2_1_1_1
   sprintf(str,"Aztec 2.1.1.1");
#endif
#ifdef AZ_ver2_2_0_0
   sprintf(str,"Aztec 2.2.0.0");
#endif
#ifdef AZ_ver2_2_0_1
   sprintf(str,"Aztec 2.2.0.1");
#endif
#ifdef AZ_ver2_2_0_2
   sprintf(str,"Aztec 2.2.0.2");
#endif
#ifdef AZ_ver2_2_0_3
   sprintf(str,"Aztec 2.2.0.3");
#endif
#ifdef AZ_ver2_2_0_4
   sprintf(str,"Aztec 2.2.0.4");
#endif
#ifdef AZ_ver2_2_1_0
   sprintf(str,"Aztec 2.2.1.0");
#endif
#ifdef AZ_ver2_2_1_1
   sprintf(str,"Aztec 2.2.1.1");
#endif
}
void AZ_space_for_kvecs(int request, int **kvec_sizes, double ***saveme, 
                        double **ptap, int *options, int *data_org, char *suffix,
			int proc, double **block)
{
   char label[64];
   int  i, j, NN;

   NN = data_org[AZ_N_internal] + data_org[AZ_N_border] + data_org[AZ_N_external];
   if (request == AZ_NEW_ADDRESS)
      AZ_manage_memory(0, AZ_SUBSELECTIVE_CLEAR,data_org[AZ_name], "kvecs", (int *) 0);

   sprintf(label,"kvecs1 %s",suffix);
   *kvec_sizes = (int *) AZ_manage_memory(2*sizeof(int),AZ_ALLOC, data_org[AZ_name],
                                          label, &j);
   if (j == AZ_NEW_ADDRESS) {
      (*kvec_sizes)[AZ_Nkept]  = 0;
      (*kvec_sizes)[AZ_Nspace] = 0;
      if (request == AZ_NEW_ADDRESS) (*kvec_sizes)[AZ_Nspace] =options[AZ_keep_kvecs];
   }
   if (request == AZ_OLD_ADDRESS) {
     if ((*kvec_sizes)[AZ_Nkept] > (*kvec_sizes)[AZ_Nspace]) {
         if (proc == 0) 
            AZ_printf_out("Number of krylov vectors exceeds space for krylov vectors?\n");
         exit(1);
     }
     if ((*kvec_sizes)[AZ_Nspace] == 0) {
         if ((proc == 0) && (options[AZ_output] != AZ_none)) {
            AZ_printf_out("AZ_kvec_space:  No previous krylov vectors available ");
            AZ_printf_out("for subspace solution.\n");
            AZ_printf_out("  - Do you want options[AZ_apply_kvecs] set to AZ_APPLY?\n");
            AZ_printf_out("  - Was options[AZ_keep_info] = 1 on previous solve?\n");
            AZ_printf_out("  - Was options[AZ_keep_kvecs] > 0 on previous solve?\n");
         }
     }
     else if ((*kvec_sizes)[AZ_Nkept] == 0) {
         if ((proc == 0) && (options[AZ_output] != AZ_none)) {
            AZ_printf_out("AZ_kvec_space: Space allocated but no previous Krylov "); 
            AZ_printf_out("vectors were kept.\n");
         }
     }
   }   

   sprintf(label,"kvecs2 %s",suffix);
   *block   = (double *) AZ_manage_memory(((*kvec_sizes)[AZ_Nspace]*(NN+1)+1)*
                                         sizeof(double), AZ_ALLOC, data_org[AZ_name], 
                                         label, &j);
   sprintf(label,"kvecs3 %s",suffix);
   *saveme  = (double **) AZ_manage_memory((*kvec_sizes)[AZ_Nspace]*sizeof(double *), 
                                         AZ_ALLOC, data_org[AZ_name], label, &j);

   for (i = 0; i < (*kvec_sizes)[AZ_Nspace]; i++)
     (*saveme)[i] = &((*block)[i*NN]);
   *ptap = &((*block)[  (*kvec_sizes)[AZ_Nspace]*NN   ]);


}
int AZ_compare_update_vs_soln(int N, double update_norm, double alpha, double p[], double x[],
        double update_reduction, int output_flag, int proc_config[], int *first_time)
{
  double t1, t2;
  int    converged = AZ_TRUE;

  if (update_norm >= 0.0) t1 = update_norm;
  else {
     if (alpha < 0) alpha = -alpha;
     t1 = alpha*sqrt(AZ_gdot(N, p, p, proc_config));
  }
  t2 = sqrt(AZ_gdot(N, x, x, proc_config));

  if (t1 > t2*update_reduction) {
     if ((output_flag != AZ_none) && *first_time && (proc_config[AZ_node] == 0) ) {
        AZ_printf_out("Update too large, convergence deferred: ||update|| = %10.3e ||sol|| = %10.3e\n",
                t1, t2);
      }
      converged   = AZ_FALSE;
      *first_time = AZ_FALSE;
  }
  return converged;
}
struct AZ_SCALING *AZ_scale_matrix_only(AZ_MATRIX *Amat, int options[],
			int proc_config[])
{
   struct AZ_SCALING *scaling;
   static int name_count = 13071;
   int    *data_org, old_name, size, i;
   double *temp;

   data_org    = Amat->data_org;
   scaling     = AZ_scaling_create();
   old_name = data_org[AZ_name];
   data_org[AZ_name] = name_count++;
   scaling->mat_name    = data_org[AZ_name];
   scaling->scaling_opt = options[AZ_scaling];
   
   size = data_org[AZ_N_internal]+data_org[AZ_N_border]+data_org[AZ_N_external];
   temp = (double *) malloc(sizeof(double)*size);
   if (temp == NULL) {
      AZ_printf_out("AZ_scale_matrix_only: Not enough space\n"); exit(1);
   }
   for (i = 0; i < size; i++) temp[i] = 0.0;
   AZ_scale_f(AZ_SCALE_MAT_RHS_SOL, Amat, options, temp, temp, proc_config,
              scaling);
   free(temp);
   data_org[AZ_name] = old_name;
   return(scaling);
}
void AZ_scale_rhs_only(double b[], AZ_MATRIX *Amat, int options[], 
		       int proc_config[], struct AZ_SCALING *scaling)
{
   int old_name, old_scaling;

   old_name = Amat->data_org[AZ_name];
   old_scaling = options[AZ_scaling];
   Amat->data_org[AZ_name] = scaling->mat_name;
   options[AZ_scaling] = scaling->scaling_opt;
   AZ_scale_f(AZ_SCALE_RHS, Amat, options, b, b, proc_config, scaling);

   Amat->data_org[AZ_name] = old_name;
   options[AZ_scaling] = old_scaling;
}
void AZ_scale_sol_only(double x[], AZ_MATRIX *Amat, int options[], 
		       int proc_config[], struct AZ_SCALING *scaling)
{
   int old_name, old_scaling;

   old_name = Amat->data_org[AZ_name];
   old_scaling = options[AZ_scaling];
   Amat->data_org[AZ_name] = scaling->mat_name;
   options[AZ_scaling] = scaling->scaling_opt;
   AZ_scale_f(AZ_SCALE_SOL, Amat, options, x, x, proc_config, scaling);

   Amat->data_org[AZ_name] = old_name;
   options[AZ_scaling] = old_scaling;
}

void AZ_scale_rhs_sol_before_iterate(double x[], double b[], AZ_MATRIX *Amat, 
	int options[], int proc_config[], struct AZ_SCALING *scaling)
{
  AZ_scale_rhs_only(b, Amat, options, proc_config, scaling);
  AZ_scale_sol_only(x, Amat, options, proc_config, scaling);
  options[AZ_scaling] = AZ_none;
}

void AZ_unscale_after_iterate(double x[], double b[], AZ_MATRIX *Amat,
			      int options[], int proc_config[],
			      struct AZ_SCALING *scaling)
{
   int old_name;

   old_name = Amat->data_org[AZ_name];

   Amat->data_org[AZ_name] = scaling->mat_name;
   options[AZ_scaling] = scaling->scaling_opt;
   AZ_scale_f(AZ_INVSCALE_SOL, Amat, options, b, x, proc_config, scaling);
   AZ_scale_f(AZ_INVSCALE_RHS, Amat, options, b, x, proc_config, scaling);
   Amat->data_org[AZ_name] = old_name;
}
void AZ_clean_scaling(struct AZ_SCALING **scaling)
{
   AZ_free_memory((*scaling)->mat_name);
   AZ_scaling_destroy(scaling);
}

void AZ_restore_unreordered_bindx(int bindx[], double val[], int update[],
				  int update_index[], int external[],
				  int extern_index[], int data_org[])
{
  /* Restore the bindx & val arrays to how they were before AZ_transform was */
  /* invoked. NOTE: this routine should only be used when reordering is      */
  /* turned off inside AZ_transfrom. That is, AZ_ALL in AZ_transform() should*/
  /* be changed to AZ_EXTERNS.                                               */


  int N, Nghost, i, *rev_extern_ind, global_id;
  
  N = data_org[AZ_N_internal] + data_org[AZ_N_border];
  Nghost = data_org[AZ_N_external];

  /* first check that we have an MSR matrix */

  if ( data_org[AZ_matrix_type] != AZ_MSR_MATRIX) {
    AZ_printf_err("AZ_restore_unreordered_bindx: Error! Only MSR matrices can be restored.\n");
    exit(1);
  }

  /* now check to make sure that the matrix has not been reordered */

  for (i = 0; i < N; i++) {
    if (update_index[i] != i) {
      AZ_printf_err("AZ_restore_unreordered_bindx: Only unreordered matrices can be restored.\n");
      AZ_printf_err("                              Change AZ_ALL in the file 'az_tools.c'\n");
      AZ_printf_err("                              during the AZ_order_ele() invokation within 'AZ_transform()' to AZ_EXTERNS'.\n");
      exit(1);
    }
  }

  /* now build the reverse index for the externals */

  rev_extern_ind = (int *) AZ_allocate(sizeof(int)*(Nghost+1));
  if (rev_extern_ind == NULL) {
    AZ_printf_err("AZ_restore_unreordered_bindx: Not enough space\n");
    exit(1);
  }

  for (i = 0; i < Nghost; i++) {
    rev_extern_ind[extern_index[i] - N] = i;
  }

  /* copy the external part of the matrix */


  for (i = N+1; i < bindx[N]; i++) {
    if (bindx[i] < N) bindx[i] = update[bindx[i]];
    else {
      global_id = external[rev_extern_ind[bindx[i]-N]];
      bindx[i] = global_id;
    }
  }

  AZ_free(rev_extern_ind);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_global2local(int data_org[], int bindx[], int update[], int update_index[],
                     int externs[], int extern_index[])

/*******************************************************************************

  Given the global column indices for the matrix and a list of elements updated
  on this processor, restore the local indices using information computed
  from a previous AZ_transform().


  On input, the column index bindx[k] is converted to j on output where

          update[j] = bindx[k]
   or
          external[i - N_update] = bindx[k] where extern_index[j] = i

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  data_org:        Az computed from previous call to AZ_transform.

  bindx:           MSR. On input, they refer to global
                   column indices. On output, they refer to local column indices
                   as described above. See User's Guide for more information.

  update:          List (global indices) of elements updated on this node.

  external:        List (global indices) of external elements on this node.

  N_external:      Number of external elements on this processor.


*******************************************************************************/

{

  /* local variables */

  int  i, j, k;
  int *bins,shift;
  int  start,end;
  int N_update, N_external;

  /**************************** execution begins ******************************/

  N_update = data_org[AZ_N_internal] + data_org[AZ_N_border];
  N_external = data_org[AZ_N_external];
  if ( data_org[AZ_matrix_type] != AZ_MSR_MATRIX) {
    AZ_printf_err("AZ_restore_unreordered_bindx: Error! Only MSR matrices can be restored.\n");
    exit(1);
  }

  /* now check to make sure that the matrix has not been reordered */

  for (i = 0; i < N_update; i++) {
    if (update_index[i] != i) {
      AZ_printf_err("AZ_restore_unreordered_bindx: Only unreordered matrices can be restored.\n");
      AZ_printf_err("                              Change AZ_ALL in the file 'az_tools.c'\n");
      AZ_printf_err("                              during the AZ_order_ele() invokation within 'AZ_transform()' to AZ_EXTERNS'.\n");
      exit(1);
    }
  }

  /* set up some bins so that we will be able to use AZ_quick_find() */

  bins = (int *) AZ_allocate((N_update / 4 + 10)*sizeof(int));
  if  (bins == NULL) {
    (void) AZ_printf_err( "ERROR: Not enough temp space\n");
    exit(-1);
  }

  AZ_init_quick_find(update, N_update, &shift, bins);

  /*
   * Compute the location of the first and last column index that is stored in
   * the bindx[].
   */

  start = bindx[0]; end = bindx[bindx[0]-1]; 
  
  for (j = start; j < end; j++) {
    k = AZ_quick_find(bindx[j], update, N_update,shift,bins);

    if (k != -1) bindx[j] = k;
    else {
       k = AZ_find_index(bindx[j], externs,N_external);
       if (k != -1) bindx[j] = extern_index[k];
       else {
        (void) AZ_printf_err( "ERROR: column number not found %d\n",
                       bindx[j]);
        exit(-1);
      }
    }
  }

  AZ_free((char *) bins);

} 

/***********************************************************************/
/***********************************************************************/
/***********************************************************************/

double AZ_condest(int n, struct context *context)
{
  int i;
  double * tmp, condest;

  tmp = (double *) AZ_allocate(n*sizeof(double));

  for (i = 0 ; i < n ; i++ ) tmp[i] = 1.0;

  AZ_solve_subdomain(tmp, n, context);

  condest = 0.0;
  for (i=0; i<n; i++)
  if (fabs(tmp[i]) > condest)
  condest = fabs(tmp[i]);

  AZ_free((void *) tmp);

  return(condest);
}
