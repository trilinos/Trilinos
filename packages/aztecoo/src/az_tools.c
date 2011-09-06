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

#ifndef TRILINOS_NO_CONFIG_H
#include "AztecOO_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_MALLOC_H
#ifndef __APPLE__
#include <malloc.h>
#else
#include <sys/malloc.h>
#endif
#endif

#include <float.h>
#include "az_aztec.h"

/* After disabling the wrappers in az_fortran_wrap.c, this is always set to
   false.  I kept the value around, rather than ripping out code ifdef'ed to
   AZ_using_fortran = AZ_TRUE because it might be useful for the new Fortran
   wrappers that are being developed. JW
*/
/*extern int AZ_using_fortran;*/
int AZ_using_fortran = AZ_FALSE;

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_order_ele(int update_index[], int extern_index[], int *internal,
                  int *border, int N_update, int bpntr[], int bindx[],
                  int extern_proc[], int N_external, int option, int m_type)

/*******************************************************************************

  Find an ordering for the elements in 'update' and 'external' (if option ==
  AZ_EXTERNS, we only order 'external'). We first order 'update' by numbering
  the internal elements first followed by the border elements. On output
  update_index[] contains the numbers 0->N_update-1 ordered such that
  update_index[k] < update_index[j] if update[k] is internal and update[j] is
  border. Next, we order the external elements such that elements that are
  updated by the same processor are contiguous.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  update_index:    On output, update_index[i] gives the local numbering of
                   global point 'update[i]'. See User's Guide.

  extern_index:    On output, extern_index[i] gives the local numbering of
                   global point 'external[i]'. See User's Guide.

  internal:        On output, equals the number of internal elements on this
                   processor.

  border:          On output, equals the number of border elements on this
                   processor.

  N_update:        Number of elements updated on this processor.

  bpntr, bindx:    Nonzero matrix information (see User's Guide).
                   Note: for MSR arrays the user invokes this routine by calling
                   AZ_order_ele(...,bindx,bindx, ....). That is both 'bpntr' and
                   'bindx' point to 'bindx'.

  extern_proc:     extern_proc[i] is updating processor of external[i].

  N_external:      Number of external elements on this processor.

  option:          option = AZ_ALL ==> order internal, border and extern ele.
                   option = AZ_EXTERNS ==> order only external elements.
  m_type:          On input, m_type = AZ_MSR_MATRIX
                     or      m_type = AZ_VBR_MATRIX

*******************************************************************************/

{

  /* local variables */

  /**************************** execution begins ******************************/

  int  i, j, lilstatus, count;
  int  *t_ptr;
  char *yo = "AZ_order_ele: ";

  *internal = 0;
  *border   = 0;

  /*
   * Sort through update elements. Classify them as either internal or
   * border. Give an index for the internal elements. Temporarily give a
   * negative index for border elements.
   */

  if (option == AZ_ALL) {
    if (m_type == AZ_MSR_MATRIX)      t_ptr = bindx;
    else if (m_type == AZ_VBR_MATRIX) t_ptr = bpntr;
    else {
      (void) AZ_printf_err( "%sERROR: Unknown matrix type (%d)\n", yo, m_type);
      exit(1);
    }

    for (i = 0; i < N_update; i++) {
      lilstatus = 0;

      for (j = t_ptr[i]; j < t_ptr[i + 1]; j++) {
        if (bindx[j] >= N_update) {
          update_index[i] = -*border - 1;

          /*
           * We encode the border elements as negative numbers so that we can
           * pull them out later.
           */

          (*border)++;
          lilstatus = 1;
          break;
        }
      }

      if (lilstatus == 0) {
        update_index[i] = *internal;
        (*internal)++;
      }
    }

    /* Use negative indices to number the border elements */

    for (i = 0; i < N_update; i++) {
      if (update_index[i] < 0)
        update_index[i] = *internal - update_index[i] - 1;
    }
  }

  else {
    for (i = 0; i < N_update; i++) update_index[i] = i;
    *internal = 0;
    *border = N_update;
  }

  /*
   * Sift through the external elements. For each newly encountered external
   * point assign it the next index in the sequence. Then look for other
   * external elements who are update by the same node and assign them the next
   * set of index numbers in the sequence (ie. elements updated by the same node
   * have consecutive indices).
   */

  count = N_update;
  for (i = 0; i < N_external; i++) extern_index[i] = -1;

  for (i = 0; i < N_external; i++) {
    if (extern_index[i] == -1) {
      extern_index[i] = count++;

      for (j = i + 1; j < N_external; j++) {
        if (extern_proc[j] == extern_proc[i]) extern_index[j] = count++;
      }
    }
  }

} /* AZ_order_ele */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_invorder_vec(double vector[], int data_org[], int update_index[], 
	int rpntr[], double newp[])
{

/*******************************************************************************

  Reorder a vector (could be right hand side or solution vector) which 
  conforms to a matrix reordered via AZ_transform or AZ_reorder_matrix
  so that it conforms to the unAZ_transformed matrix.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  vector:          On input, a vector of length 'length'. On output,
                   'vector' is reordered to be consistant with 
                   update_index[]. 

  data_org:        Array use to specifiy communication information. See User's
                   Guide.

  update_index:    update_index[i] gives the local numbering of global point
                   'update[i]'.

 ******************************************************************************/

    int i,j, ii,length, current;

    length = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];

    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) {
      for (i = 0 ; i < length ; i++ ) newp[i] = vector[ update_index[i] ] ;
    }
    else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
      current = 0;
      for (i = 0 ; i < length ; i++ ) {
        ii = update_index[i];
        for (j = rpntr[ii] ; j < rpntr[ii+1] ; j++ ) { 
          newp[current++] = vector[j] ;
        }
      }
    }
    else {
      (void) AZ_printf_err("Error: Unknown matrix type (%d) in AZ_reorder_vec\n",
		     data_org[AZ_matrix_type]);
      exit(-1);
    }

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_reorder_vec(double vector[], int data_org[], int update_index[], 
	int rpntr[])
{

/*******************************************************************************

  Reorder the vector (could be right hand side or solution vector) so
  that it conforms to a matrix reordered via AZ_transform or AZ_reorder_matrix.

  IMPORTANT: This routine assumes that update_index[] contains two sequencies of
  numbers that are ordered but intertwined. For example,

  update_index:  4 5 0 6 1 2 3 7

  seq 1 =>    0   1 2 3

  seq 2 =>4 5   6 7

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  vector:          On input, a vector of length 'length'. On output,
                   'vector' is reordered to be consistant with 
                   update_index[]. 

  data_org:        Array use to specifiy communication information. See User's
                   Guide.

  update_index:    update_index[i] gives the local numbering of global point
                   'update[i]'.

 ******************************************************************************/

    int i,ii,length, *temp;

    length = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
    if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) 
      AZ_sortqlists((char *) vector, 0, update_index, length, sizeof(double), length);
    else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
       temp = (int *) AZ_allocate( (length+1)*sizeof(int));
       if (temp == NULL) {
          (void) AZ_printf_err( "Out of memory in AZ_reorder_vec()\n");
          exit(-1);
       }
       for (i = 0 ; i < length ; i++ ) {
         ii = update_index[i];
         temp[i] = rpntr[ii+1]-rpntr[ii];
       }
       AZ_sortqlists((char *)vector,temp,update_index,rpntr[length],
                     sizeof(double),length);
       AZ_free(temp);
    }
    else {
      (void) AZ_printf_err("Error: Unknown matrix type (%d) in AZ_reorder_vec\n",
		     data_org[AZ_matrix_type]);
      exit(-1);
    }

}


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_reorder_matrix(int N_update, int bindx[], double val[],
                       int update_index[], int extern_index[], int indx[],
                       int rnptr[], int bnptr[], int N_external,
                       int cnptr[], int option, int mat_type)

/*******************************************************************************

  Reorder the matrix so that it corresponds to the new ordering given by
  'update_index' and 'extern_index'.

  IMPORTANT: This routine assumes that update_index[] contains two sequencies of
  numbers that are ordered but intertwined. For example,

  update_index:  4 5 0 6 1 2 3 7

  seq 1 =>    0   1 2 3

  seq 2 =>4 5   6 7

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        Number of elements updated on this processor.

  val,
  rnptr,
  bindx,
  indx,
  bnptr,
  cnptr:           Sparse matrix arrays. See User's Guide.
                   On input, the matrix corresponds to the initial ordering
                   (e.g. row i corresponds to global row update[i]).
                   On output, the matrix rows and columns are renumbered to
                   correspond to the ordering given by 'update_index' and
                   'extern_index'. (e.g. global row update[i] now appears
                   as row update_index[i] in the matrix).

  update_index:    update_index[i] gives the local numbering of global point
                   'update[i]'.

  extern_index:    extern_index[i] gives the local numbering of global point
                   'external[i]'.

  N_external:      Number of external elements on this processor.

  option:

  mat_type:        Indicates whether this is an MSR (= AZ_MSR_MATRIX) or a
                   VBR (= AZ_VBR_MATRIX).

*******************************************************************************/

{

  /* local variables */

  int   start, end;
  int   val_length, indx_length;
  int  *temp;
  int   i, j;
  char *yo = "AZ_reorder_matrix: ";

  /**************************** execution begins ******************************/

  if (mat_type == AZ_MSR_MATRIX) {
    start = N_update+1;      /* first nonzero offdiag */
    end   = bindx[N_update]; /* last nonzero          */
  }
  else if (mat_type == AZ_VBR_MATRIX) {
    start = 0;               /* first nonzero */
    end   = bnptr[N_update]; /* last nonzero  */

    /* reorder cnptr[] */

    /* 1) convert the cnptr array to give the blk size */

    AZ_convert_ptrs_to_values(cnptr, N_update + N_external);

    /* 2) order the internal blocks.  NOTE: AZ_sortqlists() can only be used for
     * the internal elements as it expects the list to correspond to 2 ordered
     * sequences that are intermixed.
     */

    AZ_sortqlists((char *) cnptr, 0, update_index, N_update, sizeof(int),
                  N_update);

    /* 3) order the external blocks */

    temp = (int *) AZ_allocate((unsigned)(N_external + 1)*sizeof(int));
    if (temp == NULL) {
      (void) AZ_printf_err(
                     "%sERROR: not enough memory to malloc temporary space\n",
                     yo);
      exit(-1);
    }

    for (i = 0; i < N_external; i++)
      temp[extern_index[i] - N_update] = cnptr[i + N_update];

    for (i = 0; i < N_external; i++) cnptr[i + N_update] = temp[i];
    AZ_free((char *) temp);

    /* 4) reconvert cnptr to give pointer information */

    AZ_convert_values_to_ptrs(cnptr, N_update + N_external, 0);
  }
  else {
    (void) AZ_printf_err( "%sERROR: matrix is not MSR or VBR\n", yo);
    exit(-1);
  }

  /*
   * Change column indices (bindx) to reflect new ordering depending depending
   * on whether or not a point is internal or external.
   */

  for (i = start; i < end; i++) {
    if (bindx[i] < N_update) bindx[i] = update_index[bindx[i]];
    else                     bindx[i] = extern_index[bindx[i] - N_update];
  }

  if (option == AZ_EXTERNS) return;

  /* reorder rows */

  if (mat_type == AZ_MSR_MATRIX) {

    /* We move the rows in four steps:
     *  1) sort the first N_update values of 'val'.
     *  2) sort the first N_update values of 'bindx'.
     *     We do this by first converting the ptrs to values
     *     representing the number of nonzero off diagonals.
     *  3) sort the off diagonal column indices.
     *  4) sort the off diagonal matrix nozeros.
     */

    j = bindx[0];
    AZ_convert_ptrs_to_values(bindx, N_update);

    AZ_sortqlists((char *) &(bindx[N_update + 1]), bindx, update_index,
                  end - N_update - 1, sizeof(int), N_update);
    AZ_sortqlists((char *) &(val[N_update + 1]), bindx, update_index,
                  end - N_update - 1, sizeof(double), N_update);
    AZ_sortqlists((char *) val, 0, update_index, N_update, sizeof(double),
                  N_update);
    AZ_sortqlists((char *) bindx, 0, update_index, N_update,
                  sizeof(int), N_update);
    AZ_convert_values_to_ptrs(bindx, N_update, j);
  }
  else {
    val_length  = indx[bnptr[N_update]];
    indx_length = bnptr[N_update];
    AZ_convert_ptrs_to_values(indx, indx_length);

    temp = (int *) AZ_allocate((unsigned) (N_update+1)*sizeof(int));
    if (temp == NULL) {
      (void) AZ_printf_err( "%sERROR: Not enough temp space in reorder.\n",
                     yo);
      exit(-1);
    }

    /* move val */

    for (i = 0; i < N_update; i++) {
      temp[i] = 0;
      for (j = bnptr[i]; j < bnptr[i + 1]; j++) temp[i] += indx[j];
    }

    AZ_sortqlists((char *) val, temp, update_index, val_length, sizeof(double),
                  N_update);
    AZ_free((char *) temp);

    AZ_convert_ptrs_to_values(bnptr, N_update);
    AZ_convert_ptrs_to_values(rnptr, N_update);

    /* move indx */

    AZ_sortqlists((char *) indx, bnptr, update_index, indx_length, sizeof(int),
                  N_update);

    /* move bindx */

    AZ_sortqlists((char *) bindx, bnptr, update_index, indx_length, sizeof(int),
                  N_update);

    /* move bnptr */

    AZ_sortqlists((char *) bnptr, 0, update_index, N_update, sizeof(int),
                  N_update);

    /* move rnptr */

    AZ_sortqlists((char *) rnptr, 0, update_index, N_update, sizeof(int),
                  N_update);

    AZ_convert_values_to_ptrs(rnptr, N_update, 0);
    AZ_convert_values_to_ptrs(bnptr, N_update, 0);
    AZ_convert_values_to_ptrs(indx, indx_length, 0);
  }

} /* AZ_reorder_matrix */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_set_message_info(int N_external, int extern_index[], int N_update,
                         int external[], int extern_proc[], int update[],
                         int update_index[], int proc_config[], int cnptr[],
                         int *data_org[], int mat_type)

/*******************************************************************************

  Perform initializations so that local communication can occur to update the
  external elements. This initialization includes:

   1) determine the number of neighbors to which we send information as well as
      the number of neighbors from which we receive information.  These two
      should be the same. If they are not, we will set up things so that we send
      0 length messages to processors from which we receive (but have no
      information to send them).
   2) determine the total number of unknowns that we must send out.  Using this
      information we allocate space for data_org[].
   3) Initialize data_org[] (see User's Guide) so that it contains the number of
      messages to be sent/received, the node numbers of the processors to which
      we must send and receive, the length of the messages that we must send and
      that we expect to receive from each of the neighbors, and finally, a list
      of the indices of the elements that will be send to other processors (in
      the order that they will be sent).

      NOTE: Implicitly the neighbors are numbered using the ordering of the
      external elements. In particular, the external elements are ordered such
      that those elements updated by the same processor are contiguous. In this
      way, the ordering of the external elements defines an ordering for the
      neighbors.

 The information stored in 'data_org[]' is as follows:

     data_org[AZ_neighbors]         = node number of first neighbor
     data_org[AZ_neighbors+1]       = node number of second neighbor
     data_org[AZ_neighbors+2]       = node number of third  neighbor
         . . .
     data_org[AZ_rec_length]        = # of values to receive from 1st neighbor
     data_org[AZ_rec_length+1]      = # of values to receive from 2nd neighbor
     data_org[AZ_rec_length+2]      = # of values to receive from 3rd neighbor
     . . .
     data_org[AZ_send_length]       = # of values to send to 1st neighbor
     data_org[AZ_send_length+1]     = # of values to send to 2nd neighbor
     data_org[AZ_send_length+2]     = # of values to send to 3rd neighbor
     . . .
     data_org[AZ_send_list]         = elements to be sent to neighbors (starting
                                      with elements to the first neighbor [in
                                      the order that they appear on that
                                      neighbor], followed by elements to the
                                      second neighbor, etc).

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_external:      Number of external elements on this processor.

  extern_index:    On output, extern_index[i] gives the local numbering of
                   global point 'external[i]'. See User's Guide.

  N_update:        Number of elements updated on this processor.

  external:        List (global indices) of external elements on this node.

  extern_proc:     extern_proc[i] is updating processor of external[i].

  update:          List (global indices) of elements updated on this node.

  update_index:    On output, update_index[i] gives the local numbering of
                   global point 'update[i]'. See User's Guide.

  proc_config:     proc_config[AZ_node] is node number.
                   proc_config[AZ_N_procs] is the number of processors.

  cnptr:           VBR column pointer array (see User's Guide).

  data_org:        Array use to specifiy communication information. See User's
                   Guide.
                   On output, this array contains neighbor and message.

  mat_type:        Indicates whether this is an MSR (= AZ_MSR_MATRIX) or a
                   VBR (= AZ_VBR_MATRIX).

*******************************************************************************/

{

  /* local variables */

  int   i, j, ii, type, start, cflag, partner, newlength, new_i;
  int   total_to_be_sent, found;
  int  *new_external, *new_extern_proc;
  int  *neighbors, *tempneigh;
  int  *recv_list, *send_list;
  int   tj, end, tt, oldii;
  int   num_recv_neighbors;
  int   num_send_neighbors;
  int  *bins, shift;
  int  *send_ptr, *lens;
  int   proc,nprocs;
  int   firstone,current;
  MPI_AZRequest request[AZ_MAX_NEIGHBORS];  /* Message handle */
  unsigned int length;

  char *yo = "AZ_set_message_info: ";

  /*---------------------- execution begins -----------------------------*/

  AZ__MPI_comm_space_ok();
  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  neighbors = (int *) AZ_allocate((unsigned) nprocs*sizeof(int));
  tempneigh = (int *) AZ_allocate((unsigned) nprocs*sizeof(int));

  /* Produce a list of the external updating processors corresponding */
  /* to each external point in the order given by 'extern_index[]'    */

  new_extern_proc = (int *) AZ_allocate((unsigned) (N_external+1)*sizeof(int));
  if (new_extern_proc == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }
  for (i = 0 ; i < nprocs ; i++)  neighbors[i] = 0;
  for (i = 0 ; i < N_external+1; i++)  new_extern_proc[i] = 0;

  for (i = 0; i < N_external; i++)
    new_extern_proc[extern_index[i] - N_update] = extern_proc[i];

  /*
   * Count the number of neighbors from which we receive information to update
   * our external elements. Additionally, fill the array neighbors in the
   * following way:
   *      neighbors[i] = 0   ==>  No external elements are updated by
   *                              processor i.
   *      neighbors[i] = x   ==>  (x-1)/nprocs elements are updated from
   *                              processor i.
   */

  num_recv_neighbors = 0;
  length             = 1;

  for (i = 0; i < N_external; i++) {
    if (neighbors[new_extern_proc[i]] == 0) {
      num_recv_neighbors++;
      neighbors[new_extern_proc[i]] = 1;
    }

    if (mat_type != AZ_VBR_MATRIX)
      neighbors[new_extern_proc[i]] += nprocs;
    else {
      neighbors[new_extern_proc[i]] += (nprocs * (cnptr[i + 1 + N_update]
                                                  - cnptr[i + N_update]));

    }
  }

  /*
   * Make a list of the neighbors that will send information to update our
   * external elements (in the order that we will receive this information).
   */

  recv_list = (int *) AZ_allocate((unsigned) AZ_MAX_NEIGHBORS*sizeof(int));
  if (recv_list == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }

  j = 0;
  recv_list[j++] = new_extern_proc[0];
  for (i = 1; i < N_external; i++) {
    if (new_extern_proc[i - 1] != new_extern_proc[i]) {
      recv_list[j++] = new_extern_proc[i];
    }
  }

  /* sum over all processors all the neighbors arrays */

  AZ_gsum_vec_int(neighbors, tempneigh, proc_config[AZ_N_procs], proc_config);

  /* decode the combined 'neighbors' array from all the processors */

  num_send_neighbors = neighbors[proc] % nprocs;

  /* decode 'neighbors[proc] to deduce total number of elements we must send */

  total_to_be_sent = (neighbors[proc] - num_send_neighbors) / nprocs;

  AZ_free((char *) neighbors);
  AZ_free((char *) tempneigh);

  /* send a 0 length message to each of our recv neighbors */

  send_list = (int *)AZ_allocate((unsigned) (num_send_neighbors+1)*sizeof(int));
  if (send_list == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }
  for (i = 0 ; i < num_send_neighbors+1 ; i++ ) send_list[i] = 0;


  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS +AZ_MSG_TYPE;


  /* first post receives */

  length = 0;
  for (i = 0; i < num_send_neighbors; i++) {
    partner = -1;
    mdwrap_iread((void *) external, length, &partner, &type, request+i);
  }

  /* send messages */

  for (i = 0; i < num_recv_neighbors; i++) {
    mdwrap_write((void *) external, 0, recv_list[i], type, &cflag);
  }

  /*
   * Receive message from each send neighbor to construct 'send_list'.
   */

  length = 0;
  for (i = 0; i < num_send_neighbors; i++) {
    partner = -1;
    mdwrap_wait((void *) external, length, &partner, &type, &cflag, request+i);
    if (partner != -1) send_list[i] = partner;
  }



  /*
   * Compare the two lists. In most cases they should be the same.  However, if
   * they are not add new entries to the recv list that are in the send list
   * (but not already in the recv list).
   */

  for (j = 0; j < num_send_neighbors; j++) {
    found = 0;
    for (i = 0; i < num_recv_neighbors; i++) {
      if (recv_list[i] == send_list[j])
        found = 1;
    }

    if (found == 0) {
      recv_list[num_recv_neighbors] = send_list[j];
      (num_recv_neighbors)++;
    }
  }

  AZ_free((char *) send_list);
  num_send_neighbors = num_recv_neighbors;

  /* create data_org array */

  if (!AZ_using_fortran) {
    *data_org = (int *) AZ_allocate(((unsigned) total_to_be_sent + AZ_send_list)
                               *sizeof(int));
  }
  if (*data_org == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }
  if (!AZ_using_fortran) {
     for (i = 0 ; i < total_to_be_sent + AZ_send_list; i++ ) 
        (*data_org)[i] = 0;
  }

  (*data_org)[AZ_total_send] = total_to_be_sent;
  send_ptr = &((*data_org)[AZ_send_list]);

  /*
   * Create 'new_external' which explicitly put the external elements in the
   * order given by 'extern_index'
   */

  new_external = (int *) AZ_allocate((unsigned) (N_external+1)*sizeof(int));
  if (new_external == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }

  for (i = 0; i < N_external; i++) {
    new_external[extern_index[i] - N_update] = external[i];
  }

  /*
   * Send each processor the global index list of the external elements in the
   * order that I will want to receive them when updating my external elements
   */



  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  lens = (int *) AZ_allocate((num_recv_neighbors+1)*sizeof(int));
  if (lens == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }

  for (i = 0; i < num_recv_neighbors; i++) {
    length     = sizeof(int);
    partner    = recv_list[i];
    mdwrap_iread((void *) &(lens[i]), length, &partner, &type, request+i);
  }


  j = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    start  = j;
    newlength = 0;

    /* go through list of external elements until updating processor changes */

    while ((j < N_external) && (new_extern_proc[j] == recv_list[i])) {

      if (mat_type != AZ_VBR_MATRIX) newlength++;
      else newlength += (cnptr[j + 1 + N_update] - cnptr[j + N_update]);

      j++;
      if (j == N_external) break;
    }

    (*data_org)[AZ_rec_length+i] = newlength;
    (*data_org)[AZ_neighbors+i]  = recv_list[i];

    length = j - start;
    mdwrap_write((void *) &length, sizeof(int), recv_list[i], type, &cflag);
  }

  /* receive from each neighbor the global index list of external ele */

  for (i = 0; i < num_recv_neighbors; i++) {
    length     = sizeof(int);
    partner    = (*data_org)[AZ_neighbors+i];
    mdwrap_wait((void *) &(lens[i]), length, &partner, &type, &cflag,request+i);
    (*data_org)[AZ_send_length + i] = lens[i];
  }
  AZ_free(lens);







  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;


  j = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    length     = (  (*data_org)[AZ_send_length + i] ) * sizeof(int);
    partner    = (*data_org)[AZ_neighbors+i];
    mdwrap_iread((void *) &(send_ptr[j]), length, &partner, &type, request+i);
    j += (*data_org)[AZ_send_length + i];
  }

  j = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    start  = j;
    newlength = 0;

    /* go through list of external elements until updating processor changes */

    while ((j < N_external) && (new_extern_proc[j] == recv_list[i])) {

      if (mat_type != AZ_VBR_MATRIX) newlength++;
      else newlength += (cnptr[j + 1 + N_update] - cnptr[j + N_update]);

      j++;
      if (j == N_external) break;
    }
    mdwrap_write((void *) &(new_external[start]), (j-start)* sizeof(int),
             recv_list[i], type, &cflag);
  }

  /* receive from each neighbor the global index list of external ele */

  j = 0;
  for (i = 0; i < num_recv_neighbors; i++) {
    length     = (  (*data_org)[AZ_send_length + i] ) * sizeof(int);
    partner    = (*data_org)[AZ_neighbors+i];
    mdwrap_wait((void *) &(send_ptr[j]), length, &partner, &type, &cflag,
		 request+i);
    j += (*data_org)[AZ_send_length + i];
  }



  tj = 0;
  if (mat_type == AZ_VBR_MATRIX) {
    (*data_org)[AZ_matrix_type] = AZ_VBR_MATRIX;

    /* move the information to the back of the send_ptr */

    for (i = j - 1; i >= 0; i--) send_ptr[total_to_be_sent+i-j] = send_ptr[i];
    tj = total_to_be_sent - j;
  }
  else (*data_org)[AZ_matrix_type] = AZ_MSR_MATRIX;

  total_to_be_sent = j;

  /* replace global indices by local indices */

  oldii     = ii = 0;
  tt        = tj;
  firstone  = AZ_send_length;
  current   = firstone;

  while (((*data_org)[current] == 0) && (current-firstone < num_recv_neighbors))
    current++;

  bins = (int *) AZ_allocate((unsigned) (N_update / 4 + 10)*sizeof(int));
  if (bins == NULL) {
    (void) AZ_printf_err( "%sERROR: Not enough dynamic space.\n", yo);
    exit(-1);
  }
  for (i = 0 ; i < N_update / 4 + 10 ; i++ ) bins[i] = 0;

  AZ_init_quick_find(update, N_update, &shift, bins);
  for (i = tj; i < total_to_be_sent + tj; i++) {

    /*
     * j = AZ_find_index(send_ptr[i], update, N_update);
     */

    j = AZ_quick_find(send_ptr[i], update, N_update, shift, bins);
    if (j == -1) {
      (void) AZ_printf_err( "%sERROR: A point (%d) requested from processor "
                     "%d\nwas not found among the update ele of proc %d\n", yo,
                     send_ptr[i], (*data_org)[AZ_neighbors+current-firstone],
                     proc);
      exit(-1);
    }

    if (mat_type != AZ_VBR_MATRIX) {
      send_ptr[ii] = update_index[j];
      ii +=  1;
    }
    else {
      start = cnptr[update_index[j]];
      end   = cnptr[update_index[j] + 1] - 1;
      for (new_i = start; new_i <= end; new_i++) send_ptr[ii++] = new_i;

      if (i - tt + 1 == (*data_org)[current]) {
        (*data_org)[current] = ii - oldii;
        current++;

        while ( ((*data_org)[current] == 0) &&
                (current-firstone < num_recv_neighbors)) current++;
        oldii = ii;
        tt    = i + 1;
      }
    }
  }

  AZ_free((char *) bins);
  AZ_free((char *) recv_list);
  AZ_free((char *) new_external);
  AZ_free((char *) new_extern_proc);

  (*data_org)[AZ_N_neigh] = num_send_neighbors;

} /* AZ_set_message_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_convert_values_to_ptrs(int array[], int length, int start)

/*******************************************************************************

  Change 'array[]' so that on exit
     1) array[0] = start
     2) array[i+1] - array[i] = value of array[i] on entry to this routine.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  array:           On entry to this routine and array[0] = start.
                   On output, array[i+1] - array[i] = value of array[i].

  length:          Length of array[].

*******************************************************************************/

{

  /* local variables */

  int i;

  /**************************** execution begins ******************************/

  for (i = 1; i < length; i++) array[i] += array[i - 1];
  for (i = length; i > 0; i--)  array[i]  = array[i - 1] + start;

  array[0] = start;

} /* AZ_convert_values_to_ptrs */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_convert_ptrs_to_values(int array[], int length)

/*******************************************************************************

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  array:           On entry to this routine and array[0] = start.
                   On output, array[i+1] - array[i] = value of array[i].

  length:          Length of array[].

*******************************************************************************/

{
  int i;

  for (i = 0; i < length; i++) array[i] = array[i + 1] - array[i];
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_msr2vbr(double val[], int indx[], int rnptr[], int cnptr[], int bnptr[],
                int bindx[], int msr_bindx[], double msr_val[],
                int total_blk_rows, int total_blk_cols, int blk_space,
                int nz_space, int blk_type)

/*******************************************************************************

  Convert the MSR matrix defined in [msr_val,msr_bindx] to a VBR matrix.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val,
  rnptr,
  bindx,
  indx,
  bnptr,
  cnptr:           Sparse matrix arrays. See User's Guide.
                   On input, the matrix corresponds to the initial ordering
                   (e.g. row i corresponds to global row update[i]).
                   On output, the matrix rows and columns are renumbered to
                   correspond to the ordering given by 'update_index' and
                   'extern_index'. (e.g. global row update[i] now appears
                   as row update_index[i] in the matrix).

  msr_val,
  msr_bindx:       On input, MSR matrix to be converted to VBR.
                   See User's Guide.

  total_blk_rows:  Number of block rows in resulting local VBR matrix.

  total_blk_cols:  Number of block columns in resulting local VBR matrix.

  blk_space:       Length of storage allocated for bindx[] and indx[]. An error
                   message will be printed if we try to write past these arrays.

  nz_space:        Length of storage allocated for val[]. An error message will
                   be printed if we try to write past this array.

  blk_type:        If blk_type > 0, blk_type indicates that all block rows (and
                   colunns) have the same size given by 'blk_type'. If
                   blk_type < 0, the block rows have different sizes.

*******************************************************************************/

{

  /* local variables */

  int therow, thecol;
  int i, j;

  /**************************** execution begins ******************************/

  for (i = 0; i < total_blk_rows; i++) rnptr[i] = cnptr[i];

  AZ_convert_values_to_ptrs(rnptr, total_blk_rows, 0);
  AZ_convert_values_to_ptrs(cnptr, total_blk_cols, 0);

  indx[0] = bnptr[0] = 0;

  /* go through each block row */

  for (i = 0; i < total_blk_rows; i++) {
    bnptr[i + 1] = bnptr[i];

    for (therow = rnptr[i]; therow < rnptr[i + 1]; therow++) {

      /* add the diagonal entry */

      thecol = therow;
      AZ_add_new_ele(cnptr, therow, i, bindx, bnptr, indx, val, therow,
                     msr_val[therow], total_blk_cols, blk_space, nz_space,
                     blk_type);

      /* add off diagonal entries */

      for (j = msr_bindx[therow]; j < msr_bindx[therow + 1]; j++) {
        thecol = msr_bindx[j];
        AZ_add_new_ele(cnptr, thecol, i, bindx, bnptr, indx, val, therow,
                       msr_val[j], total_blk_cols, blk_space, nz_space,
                       blk_type);
      }
    }
  }

} /* AZ_msr2vbr */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_find_block_col(int cnptr[], int column, int max_blocks, int blk_size)

/*******************************************************************************

  Return the local index of the block column witch contains the point column
  given by local index 'column'.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     int, local index of the block column.
  ============

  Parameter list:
  ===============

  cnptr:           cnptr[0] = 0 and cnptr[i+1] - cnptr[i] gives the column
                   dimension of global block column
                     'update[i]' if  i <  N_update
                     'external[k]' if  i >= N_update
                   where k = i - N_update.

  column:          Local column index of the column for which we are trying to
                   find the block column containing it.

  blk_size:        blk_size > 0 ==> all the blocks are the same size so we can
                   use a shortcut in computing the block column index.
                   blk_size = 0 ==> short cut not used.

*******************************************************************************/

{
  int blk_col;

  if (blk_size > 0)
    blk_col = column / blk_size;
  else
    blk_col = AZ_find_closest_not_larger(column, cnptr, max_blocks);

  return blk_col;

} /* find_block_col */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_find_block_in_row(int bindx[], int bnptr[], int blk_row, int blk_col,
                         int indx[], int no_elements, double val[],
                         int blk_space, int nz_space)

/*******************************************************************************

  Search the block row 'blk_row' looking for the block column 'blk_col'. If it
  is not found, create it (and initialize it to all zeros). Return the value
  'index' where indx[index] points to the start of the block (blk_row,blk_col).

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     int, index (see explanation above).
  ============

  Parameter list:
  ===============

  val,
  bindx,
  indx,
  bnptr:           Sparse matrix arrays. See User's Guide.

  blk_row,
  blk_col:         Block indices of the block for which we are looking.

  no_elements:     Number of elements in current block.

  blk_space:       Length of storage allocated for bindx[] and indx[]. An error
                   message will be printed if we try to write past these arrays.

  nz_space:        Length of storage allocated for val[]. An error message will
                   be printed if we try to write past this array.

*******************************************************************************/

{

  /* local variables */

  int   ii, k;

  char *yo = "find_block_in_row: ";

  /**************************** execution begins ******************************/

  /* look in row 'blk_row' for 'blk_col' */

  for (k = bnptr[blk_row]; k < bnptr[blk_row + 1]; k++) {
    if (bindx[k] == blk_col) return k;
  }

  /* block was not found so let us create a new block */

  if (bnptr[blk_row + 1] + 2 >= blk_space) {
    (void) AZ_printf_err( "%sERROR: not enough space for block ptrs (indx)\n",
                   yo);
    exit(-1);
  }

  if (indx[bnptr[blk_row + 1]] + no_elements >= nz_space) {
    (void) AZ_printf_err( "%sERROR: not enough space for nonzeros (val)\n",
                   yo);
    exit(-1);
  }

  /* create the block (blk_row, blk_col) */

  bindx[bnptr[blk_row + 1]]    = blk_col;
  indx[bnptr[blk_row + 1] + 1] = no_elements + indx[bnptr[blk_row + 1]];

  for (ii = 0; ii < no_elements; ii++) val[ii+indx[bnptr[blk_row + 1]]] = 0.0;
  bnptr[blk_row + 1]++;
  return (bnptr[blk_row + 1] - 1);

} /* find_block_in_row */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_add_new_ele(int cnptr[], int col, int blk_row, int bindx[], int bnptr[],
                    int indx[], double val[], int row, double new_ele,
                    int maxcols, int blk_space, int nz_space, int blk_type)

/*******************************************************************************

  Given a new element 'new_ele' (whose real row and column indices are given by
  'row' and 'col') store it in the VBR matrix given by cnptr[], bindx[],
  bnptr[],indx[], and val[].

  If the new element is in a block that already exists in the data structure,
  then we just add the new entry in 'val[]'. However, if the new element is in a
  block which does not already exist, we must create the new block and return a
  pointer to it (find_block_in_row(...)) before we can add the new entry to
  'val[]'.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  val,
  bindx,
  indx,
  bnptr:           Sparse matrix arrays. See User's Guide.

  blk_row:         Block indices of the block for which we are looking.

  row, col:        Point row and column of new entry.

  new_ele:         New value to be placed in matrix.

  maxcols:         Total number of block columns in the local submatrix stored
                   on this processor.

  blk_space:       Length of storage allocated for bindx[] and indx[]. An error
                   message will be printed if we try to write past these arrays.

  nz_space:        Length of storage allocated for val[]. An error message will
                   be printed if we try to write past this array.

  blk_type:        blk_type > 0 ==> all the blocks are the same size so we can
                   use a shortcut in computing the block column index.
                   blk_type = 0 ==> short cut not used.

*******************************************************************************/

{

  /* local variables */

  int  blk_col, no_elements, k, start_location, little_col, little_row, offset;

  /*---------------------- execution begins -----------------------------*/

  /* find block column containing 'col' */

  blk_col = AZ_find_block_col(cnptr, col, maxcols, blk_type);

  /* compute number of elements in block containing new point */

  no_elements = (cnptr[blk_col + 1] - cnptr[blk_col]) *
    (cnptr[blk_row + 1] - cnptr[blk_row]);

  /*
   * Search the block row looking for 'blk_col'. If it does not exist, create it
   * (and initialize it to all zeros). Return a ptr (actually an index into
   * indx[]) to the block corresponding to (blk_row,blk_col)
   */

  k = AZ_find_block_in_row(bindx, bnptr, blk_row, blk_col, indx, no_elements,
                           val, blk_space, nz_space);

  /* compute the location of the new element in val[] */

  start_location = indx[k];
  little_col     = col - cnptr[blk_col];
  little_row     = row - cnptr[blk_row];

  offset = little_col * (cnptr[blk_row + 1] - cnptr[blk_row]) + little_row;
  val[start_location + offset] = new_ele;

} /* add_new_ele */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_check_msr(int bindx[], int N_update, int N_external, int option,
                  int proc_config[])

/*******************************************************************************

  Check the msr matrix:
  1) check that the number of nonzero offdiagonals in each row is non-negative
     and is not too large.
  2) check that the column numbers are nonnegative and not too large.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  bindx:           MSR array (see User's Guide).

  N_update:        Number of elements updated on this processor.

  N_external:      Number of external elements on this processor.

  option:          AZ_LOCAL ==> the matrix uses local indices. Thus, the number
                                of nonzeros in a row and the largest column
                                index should not exceed the total number of
                                elements on this processor.
                   AZ_GLOBAL==> the matrix uses global indices. Thus, the number
                                of nonzeros in a row and the largest column
                                index should not exceed the total number of
                                elements in the simulation.

*******************************************************************************/

{

  /* local variables */

  int   i, largest, num, total_ele = 0;

  char *yo = "AZ_check_msr: ";

  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  /* compute the total number of elements in this simulation. */

  if (option == AZ_GLOBAL)
    total_ele = AZ_gsum_int(N_update, proc_config);

  /*
   * First check that the number of offdiagonal nonzeros is positive also
   * compute the largest number of nonzeros in a row.
   */

  largest = -1;
  for (i = 0; i < N_update; i++) {
    num = bindx[i + 1] - bindx[i];

    if (num > largest) largest = num;

    if (num < 0) {
      (void) AZ_printf_err( "%sERROR on proc %d: Number of ", yo,
                     proc_config[AZ_node]);
      (void) AZ_printf_err(
                     "nonzeros offdiagonals in row %d = (%d - %d) = %d\n", i,
                     bindx[i + 1], bindx[i], bindx[i + 1] - bindx[i]);
                     (void) AZ_printf_err( "is negative inside MSR check?\n");
    }
  }

  if (option == AZ_LOCAL) {
    if (largest > N_update + N_external) {
      (void) AZ_printf_err( "%sERROR on proc %d: Number of ", yo,
                     proc_config[AZ_node]);
      (void) AZ_printf_err( "offdiagonals in row %d exceeds the", largest);
      (void) AZ_printf_err( " number of elements on the processor %d\n",
                                    N_update + N_external);
    }
  }
  else {
    if (largest > total_ele) {
      (void) AZ_printf_err( "%sERROR on proc %d: Number of ", yo,
                     proc_config[AZ_node]);
      (void) AZ_printf_err( "offdiagonals in row %d exceeds the", largest);
      (void) AZ_printf_err( " total number of elements in simulation %d\n",
                                    total_ele);
    }
  }

  largest = AZ_gmax_int(largest, proc_config);
  if (proc_config[AZ_node] == 0) {
    (void) AZ_printf_out( "The max number of nonzeros in a row = %d\n",
                   largest);
  }

  /* compute the largest column index */

  largest = -1;
  for (i = bindx[0]; i < bindx[N_update]; i++) {
    if (bindx[i] < 0) {
      (void) AZ_printf_err( "%sWARNING on proc %d: Negative column (%d)= %d\n",
                     yo, proc_config[AZ_node], i, bindx[i]);
    }

    if (bindx[i] > largest) largest = bindx[i];
  }

  if (option == AZ_LOCAL) {
    if (largest > N_update + N_external) {
      (void) AZ_printf_err( "%sWARNING on proc %d: Column ", yo,
                     proc_config[AZ_node]);
      (void) AZ_printf_err( "referenced (%d) that does not exist\n", largest);
      (void) AZ_printf_err( "    # of elements update on proc = %d\n",
                     N_update);
      (void) AZ_printf_err( "    # of external elements       = %d\n",
                     N_external);
    }
  }
  else {

    /* find largest possible index */

    if (largest > total_ele) {
      (void) AZ_printf_err( "%sWARNING on proc %d: Column ", yo,
                     proc_config[AZ_node]);
      (void) AZ_printf_err( "referenced (%d) that is larger than ", largest);
      (void) AZ_printf_err( "the total number of elements in simulation %d.\n",
                     total_ele);
    }
  }

} /* AZ_check_msr */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_check_vbr(int N_update, int N_external, int option,
                  int bindx[], int bnptr[], int cnptr[], int rnptr[],
                  int proc_config[])

/*******************************************************************************

  Check the vbr matrix:
  1) check that the largest column block size (cnptr) is not too large and
     not <= 0.
  2) check that rnptr[i] = cnptr[i]   for i < N_update.
  3) check that the number of nonzero blocks (bnptr) is not too large
     and not <= 0
  4) check that the block column indices (bindx) are not too large
     and not <= 0.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        Number of elements updated on this processor.

  N_external:      Number of external elements on this processor.

  option:          AZ_LOCAL ==> the matrix uses local indices. Thus, the number
                                of nonzeros in a row and the largest column
                                index should not exceed the total number of
                                elements on this processor.
                   AZ_GLOBAL==> the matrix uses global indices. Thus, the number
                                of nonzeros in a row and the largest column
                                index should not exceed the total number of
                                elements in the simulation.

  rnptr,
  bindx,
  indx,
  bnptr,
  cnptr:           Sparse matrix arrays. See User's Guide.
                   On input, the matrix corresponds to the initial ordering
                   (e.g. row i corresponds to global row update[i]).
                   On output, the matrix rows and columns are renumbered to
                   correspond to the ordering given by 'update_index' and
                   'extern_index'. (e.g. global row update[i] now appears
                   as row update_index[i] in the matrix).

  proc_config:     proc_config[AZ_node] is node number.
                   proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int   i;
  int   largest, num;
  int   total_ele = 0;
  int   proc;

  char *yo = "AZ_check_vbr: ";

  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  proc   = proc_config[AZ_node];

  /* compute the total number of blocks in this simulation. */

  if (option == AZ_GLOBAL)
    total_ele = AZ_gsum_int(N_update, proc_config);

  /* first check that the number of offdiagonal nonzeros is positive */
  /* also compute the largest number of nonzeros in a row */

  largest = -1;
  for (i = 0; i < N_update; i++) {
    num = bnptr[i + 1] - bnptr[i];

    if (num > largest) largest = num;

    if (num < 0) {
      (void) AZ_printf_err( "%sERROR on proc %d: Number of nonzeros ", yo,
                     proc);
      (void) AZ_printf_err( "blocks in row %d = (%d - %d) = %d\n", i,
                     bnptr[i + 1], bnptr[i], bnptr[i + 1] - bnptr[i]);
      (void) AZ_printf_err( "are negative inside AZ_vbr_check()?\n");
    }
  }

  if (option == AZ_LOCAL) {
    if (largest > N_update + N_external) {
      (void) AZ_printf_err( "ERROR on proc %d: Number of blocks ", proc);
      (void) AZ_printf_err( "in a row (%d) exceeds the number of ", largest);
      (void) AZ_printf_err( "blocks on the processor %d\n",
                     N_update + N_external);
    }
  }
  else {
    if (largest > total_ele) {
      (void) AZ_printf_err( "ERROR on proc %d: Number of blocks ", proc);
      (void) AZ_printf_err( "in row %d exceeds the total number ", largest);
      (void) AZ_printf_err( "of blocks in simulation %d\n", total_ele);
    }
  }

  largest = AZ_gmax_int(largest, proc_config);
  if (proc == 0) {
    (void) AZ_printf_err( "The max number of nonzero blocks in a row = %d\n",
                   largest);
  }

  /* compute the largest column index */

  largest = -1;
  for (i = 0; i < bnptr[N_update]; i++) {
    if (bindx[i] < 0) {
      (void) AZ_printf_err( "Warning on proc %d: Negative ", proc);
      (void) AZ_printf_err( "column (%d)= %d\n", i, bindx[i]);
    }

    if (bindx[i] > largest) largest = bindx[i];
  }

  if (option == AZ_LOCAL) {
    if (largest > N_update + N_external) {
      (void) AZ_printf_err( "Warning on proc %d: Column referenced ", proc);
      (void) AZ_printf_err( "(%d) that does not exist\n", largest);
      (void) AZ_printf_err( "    # of blocks update on proc = %d\n",
                     N_update);
      (void) AZ_printf_err( "    # of external blocks = %d\n", N_external);
    }
  }
  else {
    if (largest > total_ele) {
      (void) AZ_printf_err( "Warning on proc %d: Column referenced ", proc);
      (void) AZ_printf_err( "(%d) that is larger than the total ", largest);
      (void) AZ_printf_err( "number of blocks in simulation %d\n", total_ele);
    }
  }

  largest = AZ_gmax_int(largest, proc_config);
  if (proc == 0) {
    (void) AZ_printf_err( "The largest block column index is = %d\n", largest);
  }

  /* check rnptr */

  for (i = 0; i <= N_update; i++) {
    if (rnptr[i] != cnptr[i]) {
      (void) AZ_printf_err(
                     "ERROR on proc %d: rnptr(%d) != cnptr(%d) (%d vs %d)\n",
                     proc, i, i, rnptr[i], cnptr[i]);
    }
  }

  /* check cnptr */

  largest = -1;
  for (i = 0; i < N_update + N_external; i++) {
    num = cnptr[i + 1] - cnptr[i];
    if (num > largest) largest = num;

    if (num <= 0) {
      (void) AZ_printf_err( "ERROR on proc %d: Block Size for ", proc);
      (void) AZ_printf_err( "column %d = (%d - %d) = %d\n", i, cnptr[i + 1],
                     cnptr[i], cnptr[i + 1] - cnptr[i]);
    }
  }

  largest = AZ_gmax_int(largest, proc_config);
  if (proc == 0) {
    (void) AZ_printf_err( "The largest column block size is = %d\n", largest);
  }

} /* AZ_check_vbr */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_find_closest_not_larger(int key, int list[], int length)

/*******************************************************************************

  Find the closest number to 'key' in 'list' which is not greater than key and
  return the index number.

  On exit, AZ_find_index() returns: i => list[i] = key or list[i] is closest
  number smaller than key.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     int, see explanation above.
  ============

  Parameter list:
  ===============

  key:             Element to be search for in list.

  list:            List (assumed to be in ascending order) to be searched.

  length:          Length of list.

*******************************************************************************/

{

  /* local variables */

  int mid, start, end;

  /**************************** execution begins ******************************/

  if (length == 0) return -1;

  start = 0;
  end   = length - 1;

  while (end - start > 1) {
    mid = (start + end) / 2;
    if (list[mid] > key) end = mid;
    else start = mid;
  }

  if (list[end] > key) return start;
  else return end;

} /* AZ_find_closest_not_larger */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_find_procs_for_externs(int N_update, int update[],
                               int external[], int N_external,
                               int proc_config[], int **extern_proc)

/*******************************************************************************

  Determine which processors update our external elements. This is done in the
  following way:

  1) the external elements are split into subpieces (to avoid using too much
     storage or overflowing message buffers) on each processor.
  2) the current subpieces for all the processors are concatentated to
     form one large list of external elements with the name of the
     processor who contributed each subpiece.
  3) The subpieces list is searched and a neighbor list is created.This list
     corresponds to the processors who contributed subpieces containing
     an external point that is updated by this processor.
  4) A message is sent out to the neighbors who in turn send all of
     their external elements back to this processor.
  5) We recieve the external elements from the neighbors, determine which
     ones we update, and send them back to the neighbors.
  6) Finally, 'extern_proc' is set for those external elements which were found.
  7) Steps 1-6 are repeated until all the external elements are discovered.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        Number of elements updated on this processor.

  update:          List (global indices) of elements updated on this node.

  external:        List (global indices) of external elements on this node.

  N_external:      Number of external elements on this processor.

  proc_config:     proc_config[AZ_node] is node number.
                   proc_config[AZ_N_procs] is the number of processors.

  extern_proc:     extern_proc[i] is updating processor of external[i].

*******************************************************************************/

{

  /* local variables */

  int  *extern_mess, *neigh_list, *tempneigh;
  int  i, j, k = 0, ii, jj, kk;
  int  num_send_neighbors, num_recv_neighbors;
  int  type, cflag, partner, length;
  int  num_sublist_elements, all_sublists_length;
  int  found = 0, flag = 0;
  int  stride = 0;
  int  first_extern;
  int  buf_size;
  int  nprocs, proc;
  MPI_AZRequest request[AZ_MAX_NEIGHBORS];  /* Message handle */
  MPI_AZRequest send_request[2*AZ_MAX_NEIGHBORS];
  int num_sent, exch_neighbor;
  unsigned int exch_length = 0;
  int to_recv;
int oldk, iii, *ext_copy, **exch_buffers, exch_count, 
    *newbuffer, send_counter;


  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  if (!AZ_using_fortran) {
    *extern_proc = (int *) AZ_allocate((unsigned) (N_external+1)*sizeof(int));
    if ( *extern_proc == NULL) {
       (void) AZ_printf_err( "%d: Not enough memory for extern_proc\n",proc);
       exit(-1);
     }
  }

  for (i = 0; i < N_external; i++) (*extern_proc)[i] = -1;

  neigh_list = (int *) AZ_allocate((unsigned) nprocs*sizeof(int));
  tempneigh  = (int *) AZ_allocate((unsigned) nprocs*sizeof(int));

  if ( (neigh_list == NULL) || (tempneigh == NULL)) {
       (void) AZ_printf_err( " Not enough memory to find externals\n");
       exit(-1);
  }

  /*
   * 'extern_mess' must hold the concatentated sublists as well as all the
   * externals on a processor.
   */

  if (N_external > AZ_MAX_MESSAGE_SIZE) i = AZ_MAX_MESSAGE_SIZE;
  else                                  i = N_external;

  buf_size = AZ_gmax_int(i + 1, proc_config);
  buf_size = AZ_MAX(buf_size, (AZ_TEST_ELE+1) * nprocs);

  extern_mess = (int *) AZ_allocate((unsigned) buf_size*sizeof(int));
  if (extern_mess == NULL) {
    (void) AZ_printf_err( " Not enough memory to find the external procs\n");
    exit(-1);
  }
  for (i = 0 ; i < buf_size ; i++) extern_mess[i] = 0;

  first_extern = 0;

  while (flag == 0) {

    /*
     * Initialize  extern_mess to a sublist of external (of size
     * AZ_TEST_ELE). Additionally, encode the processor number at the head of
     * this list as -(proc+1)
     */

    extern_mess[0]       = -proc - 1;
    num_sublist_elements = AZ_MIN(AZ_TEST_ELE, N_external - first_extern);

    if (num_sublist_elements != 0)
      stride = (N_external - first_extern) / num_sublist_elements;

    for (i = 0; i < num_sublist_elements; i++)
      extern_mess[i + 1] = external[i * stride + first_extern];

    all_sublists_length = num_sublist_elements + 1;

    /*
     * Append together all the sublists found on all processors.  After this
     * operation each processor has a list of all the external elements
     * contained in all the sublists.
     */

    AZ_gappend(extern_mess, &all_sublists_length, buf_size, proc_config);

    /*
     * Figure out which external elements we update by searching all the
     * elements in the concatinated sublists.  Store the processor numbers in
     * the front of 'extern_mess' corresponding to the external elements which
     * we update.
     */

    for (i = 0; i < nprocs; i++) neigh_list[i] = 0;
    num_send_neighbors = 0;

    for (i = 0; i < all_sublists_length; i++) {

      /*
       * Is extern_mess[i] in the current sequence or is it the head of a new
       * subsequence.
       */

      if (extern_mess[i] >= 0) {
        if (found == -1) {      /* no point in current subsequence */

          /* updated by this processor */

          found = AZ_find_index(extern_mess[i], update, N_update);
          if (found != -1) {

            /* this point updated by this processor */

            neigh_list[-k - 1]                = 1;
            extern_mess[num_send_neighbors++] = k;
          }
        }
      }
      else {
        found = -1;
        k     = extern_mess[i]; /* k = processor of current sublist */
      }
    }

    /* find the total number of neighbors from who we receive */

    AZ_gsum_vec_int(neigh_list, tempneigh, proc_config[AZ_N_procs],proc_config);
    num_recv_neighbors = neigh_list[proc];

    if ((num_sublist_elements != 0) && (num_recv_neighbors == 0)) {
      (void) AZ_printf_err( "A column referenced on one processor was not\n");
      (void) AZ_printf_err(
                     "found on any other processors. This means that\n");
      (void) AZ_printf_err( "no processor is assigned to the row equal to \n");
      (void) AZ_printf_err( "this column number.\n");
      (void) AZ_printf_err( "Note: matrices must be square.\n");
      (void) AZ_printf_err( "The following columns were not found:\n");

      for (i = 0 ; i < num_sublist_elements; i++ ) {
        (void) AZ_printf_err( "%3d ", external[i * stride + first_extern]);
      }
      (void) AZ_printf_err("\n");
        exit(-1);
    }

    /* send a zero length message to processors in the neigh_list */

    type            = proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE)%AZ_NUM_MSGS + AZ_MSG_TYPE;


    /* post receives */

    for (i = 0; i < num_recv_neighbors; i++) {
      partner       = -1;
      (void) mdwrap_iread((void *) &k, 0, &partner, &type, request+i);
    }

    /* send messages */

    for (i = 0; i < num_send_neighbors; i++) {
      k = 1;
      mdwrap_write((void *) &k, 0, -extern_mess[i] - 1, type, &cflag);
    }

    /* receive the messages and create a new neighbor list */

    for (i = 0; i < num_recv_neighbors; i++) {
      partner  = -1;
      length   = mdwrap_wait((void *) &k, 0, &partner, &type, &cflag,request+i);
      neigh_list[i] = partner;
    }


    type            = proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE)%AZ_NUM_MSGS + AZ_MSG_TYPE;

    k               = N_external - first_extern;
    if (k > AZ_MAX_MESSAGE_SIZE) k = AZ_MAX_MESSAGE_SIZE;
    if (k > buf_size) k = buf_size;

    num_sent  = 0;
    exch_neighbor = -1;
    to_recv = num_recv_neighbors + num_send_neighbors;

    /* set variables to do a nonblocking sends */

    send_counter = 0; 
    exch_count   = 0;
    oldk         = k;
    ext_copy = (int *) AZ_allocate((k+1)*sizeof(int));
    exch_buffers = (int **) AZ_allocate((num_send_neighbors+1)*sizeof(int *));
    if ( (ext_copy == NULL) || (exch_buffers == NULL) ) {
       AZ_printf_err("Error: Out of memory in AZ_transform\n");
       exit(1);
    }
    for (iii = 0 ; iii < k ; iii++ ) ext_copy[iii] = external[first_extern+iii];


    while ( to_recv || (exch_neighbor != -1) ) {
      partner = -1;

      if (to_recv) mdwrap_iread((void *) extern_mess, buf_size*sizeof(int), 
                                 &partner, &type, request);
 
      if (exch_neighbor != -1) {
        mdwrap_iwrite((void *) exch_buffers[exch_count], exch_length, 
                       exch_neighbor, type, &cflag, 
                       &(send_request[send_counter]));
        send_counter++;
        exch_neighbor = -1;
        exch_count++;
      }
      else if (num_sent < num_recv_neighbors) {
        mdwrap_iwrite((void *) ext_copy, oldk * sizeof(int),
                      neigh_list[num_sent], type, &cflag, 
                      &(send_request[send_counter]));
        send_counter++;
        num_sent++;
      }
 
      if (to_recv) {
         length  = mdwrap_wait((void *) extern_mess, buf_size*sizeof(int), 
                                &partner, &type, &cflag, request);
         length = length/sizeof(int);
         to_recv--;
 
         if (extern_mess[0] >= 0 ) {
 
            /* sift through neighbors and reduce */
 
            kk = 0;
            for (j = 0 ; j < length ; j++ ) {
              if (AZ_find_index(extern_mess[j],update,N_update) != -1)
                 extern_mess[kk++] = extern_mess[j];
            }
            exch_buffers[exch_count] = (int *) AZ_allocate((kk+1)*sizeof(int));
            if (exch_buffers[exch_count] == NULL) {
               (void) AZ_printf_err( " Not enough memory to find the ");
               (void) AZ_printf_err( "external procs\n");
               exit(-1);
            }
            newbuffer = exch_buffers[exch_count];
            newbuffer[0] = -1;
            for (j = 0; j < kk ; j++ ) newbuffer[j+1] = extern_mess[j];


            exch_neighbor = partner;
            exch_length   = (kk+1)*sizeof(int); /* first element indicates */
 					        /* message type */
         }
         else {
 
            /* mark found external points */
 
            jj = 1;
            j  = first_extern;
            while( (j < N_external) && (jj < length) ) {
              if (external[j] == extern_mess[jj]) {
                (*extern_proc)[j] = partner;
                jj++;
              }
              j++;
            }
 
            /* Move things we did not find to the back of the list */
 
            ii = N_external - 1;
            for (j = N_external - 1; j >= first_extern ; j--) {
              if ((*extern_proc)[j] == -1) { /* swap */
                 kk = (*extern_proc)[ii];
                 (*extern_proc)[ii] = (*extern_proc)[j];
                 (*extern_proc)[j]  = kk;
                 kk = external[ii];
                 external[ii] = external[j];
                 external[j]  = kk;
                 ii--;
              }
            }
            first_extern += (jj-1);
            k   = N_external - first_extern;
            if (k > AZ_MAX_MESSAGE_SIZE) k = AZ_MAX_MESSAGE_SIZE;
 

         }
      }
    }

    /* did we finish */

    if (first_extern == N_external) flag = 1;
    else                            flag = 0;

    flag = AZ_gsum_int(flag, proc_config);

    if (flag == nprocs) flag = 1;
    else                flag = 0;

    for (iii = exch_count-1 ; iii >= 0; iii--) AZ_free(exch_buffers[iii]);
    AZ_free(exch_buffers);
    AZ_free(ext_copy);
    for (i = 0; i < send_counter; i++) 
       mdwrap_request_free( &(send_request[i]) );
  }

  AZ_free((char *) neigh_list);
  AZ_free((char *) tempneigh);
  AZ_free((char *) extern_mess);

  AZ_sort(external, N_external, *extern_proc, NULL);

} /* AZ_find_procs_for_externs */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_find_local_indices(int N_update, int bindx[], int update[],
                           int **external, int *N_external, int mat_type,
                           int bpntr[])

/*******************************************************************************

  Given the global column indices for the matrix and a list of elements updated
  on this processor, compute the external elements and store them in the list
  'external' and change the global column indices to local column indices. In
  particular,

  On input, the column index bindx[k] is converted to j on output where

          update[j] = bindx[k]
   or
          external[j - N_update] = bindx[k]

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        Number of elements updated on this processor.

  bindx:           MSR or VBR column indices. On input, they refer to global
                   column indices. On output, they refer to local column indices
                   as described above. See User's Guide for more information.

  update:          List (global indices) of elements updated on this node.

  external:        List (global indices) of external elements on this node.

  N_external:      Number of external elements on this processor.

  mat_type:        Indicates whether this is an MSR (= AZ_MSR_MATRIX) or a
                   VBR (= AZ_VBR_MATRIX).

  bnptr:           'bpntr[N_update]' indicates the location of the
                   last VBR nonzero.

*******************************************************************************/

{

  /* local variables */

  int  i, j, k, kk;
  int *temp_ele;
  int *bins,shift;
  int  start,end;

  /**************************** execution begins ******************************/

  /* set up some bins so that we will be able to use AZ_quick_find() */

  bins = (int *) AZ_allocate((unsigned) (N_update / 4 + 10)*sizeof(int));
  if  (bins == NULL) {
    (void) AZ_printf_err( "ERROR: Not enough temp space\n");
    exit(-1);
  }
  for (i = 0 ; i < N_update / 4 + 10 ; i++ ) bins[i] = 0;

  AZ_init_quick_find(update, N_update, &shift, bins);

  /*
   * Compute the location of the first and last column index that is stored in
   * the bindx[].
   */

  if (mat_type == AZ_MSR_MATRIX) {
    start = bindx[0]; end = bindx[N_update];
  }
  else {
    start = 0;        end = bpntr[N_update];
  }

  /*
   * Estimate the amount of space we will need by counting the number of
   * references to columns not found among 'update[]'.  At the same time replace
   * column indices found in update[] by the appropriate index into update[].
   * Add N_update to columns not found in 'update[]' (this effectively marks
   * them as external elements).
   *
   * Note the space estimate will be an over-estimate as we do not take into
   * account that there will be duplicates in the external list.
   */

  *N_external = 0;
  for (j = start; j < end; j++) {
    k = AZ_quick_find(bindx[j], update, N_update,shift,bins);

    if (k != -1) bindx[j] = k;
    else {
      if (bindx[j] < 0) {
        (void) AZ_printf_err( "ERROR: Negative column number found %d\n",
                       bindx[j]);
        exit(-1);
      }

      bindx[j] = N_update + bindx[j];
      (*N_external)++;
    }
  }

  AZ_free((char *) bins);

  /* temporarily record all external column indices in 'temp_ele' */

  temp_ele = (int *) AZ_allocate((unsigned) (*N_external+1)*sizeof(int));
  if (temp_ele == NULL) {
    (void) AZ_printf_err(
                   "Not enough temp space in AZ_find_local_indices()\n");
    exit(-1);
  }

  *N_external = 0;
  for (i = start; i < end; i++)
    if (bindx[i] >= N_update)
      temp_ele[(*N_external)++] = bindx[i] - N_update;

  /* sort the external elements */

  AZ_sort(temp_ele, *N_external, NULL, NULL);

  /* remove duplicates */

  kk = 0;
  for (k = 1; k < *N_external; k++) {
    if (temp_ele[kk] != temp_ele[k]) {
      kk++;
      temp_ele[kk] = temp_ele[k];
    }
  }

  if (*N_external != 0) kk++;
  *N_external = kk;

  /* Now store the external elements in 'external' */

  if (!AZ_using_fortran)
    *external = (int *) AZ_allocate((unsigned) (*N_external+1)*sizeof(int));
  if (*external == NULL) {
    (void) AZ_printf_err( "Not enough space for external in %s",
                   "AZ_find_local_indices()\n");
    exit(-1);
  }

  for (i = 0; i < *N_external; i++) (*external)[i] = temp_ele[i];
  AZ_free((char *) temp_ele);

  /*
   * Update 'bindx[]' so that external elements are replaced by an index into
   * 'external[]'.
   */

  for (i = start; i < end; i++) {
    if (bindx[i] >= N_update) {
      k = AZ_find_index(bindx[i]-N_update,*external,*N_external);
      bindx[i] = k + N_update;
    }
  }

} /* AZ_find_local_indices */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_read_update(int *N_update, int *update[], int proc_config[],
                    int N, int chunk, int input_option)

/*******************************************************************************

  This routine initializes update[] to the global indices updated by this
  processor and initializes N_update to the total number of elements to be
  updated.

  If input_option == AZ_linear Do a linear partitioning of the chunks.
     Specifically, proc 0 is assigned the first floor( (N+P-1)/P ) chunks,
     processor 1 is assigned the next floor( (N+P-2)/P ) chunks, etc. where P =
     proc_config[AZ_N_procs].
  If input_option == AZ_file Processor 0 reads the file '.update'.  This file
     should contain nprocs lists.  Each list consists of a number telling how
     many global indices are in the list followed by a list of global indices.
     The first list is then sent to processor 'nprocs-1', the second list is
     sent to processor 'nprocs-2', etc.
  If input_option == AZ_box
     we do a box partitioning of the unknowns (see comments below).

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N_update:        On Output, number of unknowns updated by this processor.

  update:          On Output, list of unknowns updated by this processor in
                   ascending order.

  proc_config:     proc_config[AZ_node] is node number.
                   proc_config[AZ_N_procs] is the number of processors.

  N:               Total number of chunks to be distributed.

  chunk:           Size of each chunk to be treated as a single unit.
                   The unknowns contained in the kth chunk are given
                   by {k*chunk, k*chunk + 1, ..... , (k+1)*chunk - 1}
                   and 'N*chunk' is the total number of unknowns to be
                   distributed.

  input_option:    AZ_linear   ==> perform linear partitioning
                   AZ_file     ==> read partioning from file '.update'
                   AZ_box      ==> perform a box partitioning.

*******************************************************************************/

{

  /* local variables */

  int   t1, t2, i;
  int   ii, j;
  int   allocated, length;
  int   cflag;
  int   partner;
  int   proc_x, proc_y, proc_z;
  int   pts_x, pts_y, pts_z;
  int   total_pts_x, total_pts_y;
  int   px, py, pz, k;
  int   start_x, start_y, start_z;
  int   end_x, end_y, end_z;
  int   pt_number;
  int   count, check;
  int   proc, nprocs;
  int   type, type2;
  FILE *fp = NULL;
  MPI_AZRequest request;  /* Message handle */


  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  /*
   * Figure out which chunks should be assigned to this processor using a box
   * decomposition.  That is, it is assumed that all the chunks are ordered
   * naturally corresponding to an m x m x m box where m = N^(1/3).  Boxes of
   * chunks are assigned to processors.
   *
   * NOTE: it is assumed that nprocs = 2^power and that the number of chunks in
   * each direction is divisible by the number of processors in each direction.
   */

  if (input_option == AZ_box) {

    /* determine the number of processors in each direction */

    if (proc == 0) {
       (void) printf("Input the dimensions of the processor cube\n\n");
       (void) printf("Enter the number of processors along x axis>");
       (void) fflush(stdout);
       scanf("%d",&proc_x);
       (void) printf("Enter the number of processors along y axis>");
       (void) fflush(stdout);
       scanf("%d",&proc_y);
       (void) printf("Enter the number of processors along z axis>");
       (void) fflush(stdout);
       scanf("%d",&proc_z);

       (void) printf("Input the grid dimensions\n\n");
       (void) printf("Enter the number of grid points along x axis>");
       (void) fflush(stdout);
       scanf("%d",&pts_x);
       (void) printf("Enter the number of grid points along y axis>");
       (void) fflush(stdout);
       scanf("%d",&pts_y);
       (void) printf("Enter the number of grid points along z axis>");
       (void) fflush(stdout);
       scanf("%d",&pts_z);
    }
    AZ_broadcast((char *) &proc_x, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &proc_y, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &proc_z, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &pts_x , sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &pts_y , sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &pts_z , sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) NULL   , 0          , proc_config, AZ_SEND);
 
    total_pts_x = pts_x;
    total_pts_y = pts_y;


    if ( proc_x*proc_y*proc_z != nprocs) {
        if (proc == 0) {
          (void) printf("Error: %d x %d x %d != %d ",
 			 proc_x, proc_y, proc_z, nprocs);
          (void) printf(" (total number of processors)\n");
        }
	exit(1);
    }

    if ( pts_x * pts_y * pts_z != N ) {
        if (proc == 0) {
          (void) printf("Error: %d x %d x %d != %d ",
 			 pts_x, pts_y, pts_z, N);
          (void) printf(" (total number of grid points)\n");
        }
	exit(1);
    }
    if ( pts_x%proc_x != 0 ) {
        if (proc == 0) {
          (void) printf("Error: grid points along x axis are not an ");
          (void) printf("even multiple of processors\n");
	  (void) printf("       along x axis.");
        }
	exit(1);
    }
    if ( pts_y%proc_y != 0 ) {
        if (proc == 0) {
          (void) printf("Error: grid points along y axis is not an ");
          (void) printf("even multiple of processors\n");
	  (void) printf("       along y axis.");
        }
	exit(1);
    }
    if ( pts_z%proc_z != 0 ) {
        if (proc == 0) {
          (void) printf("Error: grid points along z axis is not an ");
          (void) printf("even multiple of processors\n");
	  (void) printf("       along z axis.");
        }
	exit(1);
    }
    pts_x = pts_x/proc_x;
    pts_y = pts_y/proc_y;
    pts_z = pts_z/proc_z;

    /* compute the number of elements per processor in each direction */

    *N_update = pts_x * pts_y * pts_z * chunk;
    if (!AZ_using_fortran) 
       *update     = (int *) AZ_allocate((*N_update)*sizeof(int));

    /* compute the lower left corner and the upper right corner */

    px = proc % proc_x;
    pz = (proc-px) / proc_x;
    py = pz % proc_y;
    pz = (pz-py) / proc_y;

    start_x = px * pts_x;
    end_x   = start_x + pts_x;
    start_y = py * pts_y;
    end_y   = start_y + pts_y;
    start_z = pz * pts_z;
    end_z   = start_z + pts_z;

    /* set update[] */

    count = 0;
    for (k = start_z; k < end_z; k++ ) {
      for (j = start_y; j < end_y; j++ ) {
        for (i = start_x; i < end_x; i++ ) {
          for (ii = 0; ii < chunk; ii++ ) {
            pt_number = (i + j * total_pts_x + k * total_pts_x * total_pts_y) * 
                            chunk + ii;
            (*update)[count++] = pt_number;
          }
        }
      }
    }
  }

  else if (input_option == AZ_linear) {

    /*
     * Figure out which chunks should be assigned to this processor for linear
     * partitioning.  This means that processor 0 is assigned the chunks
     * approximately corresponding to 0, ... , N/nprocs and processor 1 is
     * approximately assigned the chunks 1+N/nprocs to 2*N/nprocs.
     */

    t1 = N/nprocs;
    t2 = N - t1 * nprocs;

    if ( proc >= t2) t2 += (proc * t1);
    else {
      t1++;
      t2    = proc*t1;
    }
    *N_update = t1*chunk;
    t2   *= chunk;

    if (!AZ_using_fortran) 
       *update = (int *) AZ_allocate((*N_update+1)*sizeof(int));

    if (*update == NULL) {
      (void) printf( "Not enough space to allocate 'update'\n");
      exit(-1);
    }

    for (i = 0; i < *N_update; i++) (*update)[i] = i + t2;
  }

  else if (input_option == AZ_file) {

    /* read the update elements from the file '.update' */

    type            = proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE)%AZ_NUM_MSGS +AZ_MSG_TYPE;
    type2           = proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] = (type2+1-AZ_MSG_TYPE)%AZ_NUM_MSGS +AZ_MSG_TYPE;

    /*
     * Processor 0 reads file '.update' and distributes the lists to the other
     * processors.
     */

    t1 = 0;            /* total number of points distributed */
    if (proc == 0) {
      (void) printf("Reading from file .update\n"); AZ_flush_out();

      if ( (fp = fopen(".update", "r")) == NULL) {
        (void) printf( "ERROR: file '.update' not found\n");
        exit(-1);
      }

      if (!AZ_using_fortran) *update = 0;
      allocated   = 0;
      for (i = nprocs - 1; i >= 0; i-- ) {

        /* read in list length and send to processor i  */

        fscanf(fp, "%d", &length);
        t1 += length;
        if (i != 0)
          mdwrap_write((void *) &length, sizeof(int), i, type, &cflag);

        /*
         * If this is the last list, we allocate the right amount of space and
         * keep the list instead of sending it off
         */

        if (i == 0) {
          *N_update = length;
          allocated       = 0;
        }

        /* allocate enough space for list */

        if (length > allocated ) {
          if ((*update != NULL) && (!AZ_using_fortran)) AZ_free(*update);
          allocated = length + 1;

          if (!AZ_using_fortran)
            *update = (int *) AZ_allocate(allocated*sizeof(int));
          if (*update == NULL) {
            (void) fprintf(stderr,
                           "Not enough space to allocate 'update'\n");
            exit(-1);
          }
        }

        /* read a list and send it off to proc i (if not last list) */

        for (j = 0; j < length; j++ ) fscanf(fp, "%d", *update + j);
        if (i != 0)
          mdwrap_write((void *) *update, length*sizeof(int), i, type2, &cflag);
      }
      fclose(fp);

      if (t1 != N*chunk) {
        (void) fprintf(stderr,"AZ_read_update() found %d points in file\n", t1);
        (void) fprintf(stderr,"'.update' instead of the requested %d\n",
                       N*chunk);
        exit(-1);
      }
    }

    else {

      /* read the update list from processor 0 */

      partner = 0;
      mdwrap_iread((void *) N_update, sizeof(int), &partner, &type, &request);
      mdwrap_wait((void *) N_update, sizeof(int), &partner, &type, &cflag, &request);

      if (!AZ_using_fortran)
        *update = (int *) AZ_allocate((*N_update+1)*sizeof(int));
      if (*update == NULL)  {
        (void) fprintf(stderr, "Not enough space to allocate 'update'\n");
        exit(-1);
      }

      partner = 0;
      mdwrap_iread((void *) *update, *N_update * sizeof(int), &partner, &type2,
            &request);
      mdwrap_wait((void *) *update, *N_update * sizeof(int), &partner, &type2,
            &cflag, &request);
    }

    AZ_sort(*update, *N_update, NULL, NULL);

    /* check that element '0' is contained on 1 processor. That is,  */
    /* make sure the user has numbered from 0 to n-1 instead of from */
    /* 1 to n                                                        */
    check = 0;
    if ( (*N_update > 0) && ((*update)[0] == 0) ) check = 1;
    check = AZ_gsum_int(check, proc_config);
    if (check != 1) {
       if (proc == 0) {
          (void) fprintf(stderr,"Error: In AZ_read_update(), the '.update'");
          (void) fprintf(stderr,"file does not contain\n       one ");
          (void) fprintf(stderr,"occurance of row 0. Make sure that rows are");
          (void) fprintf(stderr," numbered\n       from 0 to n-1.\n");
       }
       exit(1);
    }


  }
  else {
    (void) fprintf(stderr,"Unknown input option (%d) in AZ_read_update()\n",
                   input_option);
    exit(1);
  }


} /* AZ_read_update */

/*****************************************************************************/
/*****************************************************************************/
void AZ_input_update(char datafile[], int *N_update, int *update[], int proc_config[],
                    int N, int chunk, int input_option)

/*******************************************************************************

Exactly the same as AZ_read_update except it reads the update information from
a file speficied by the input argument datafile instead of .update

*******************************************************************************/

{

  /* local variables */

  int   t1, t2, i;
  int   ii, j;
  int   allocated, length;
  int   cflag;
  int   partner;
  int   proc_x, proc_y, proc_z;
  int   pts_x, pts_y, pts_z;
  int   total_pts_x, total_pts_y;
  int   px, py, pz, k;
  int   start_x, start_y, start_z;
  int   end_x, end_y, end_z;
  int   pt_number;
  int   count, check;
  int   proc, nprocs;
  int   type, type2;
  FILE *fp = NULL;
  MPI_AZRequest request;  /* Message handle */


  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  /*
   * Figure out which chunks should be assigned to this processor using a box
   * decomposition.  That is, it is assumed that all the chunks are ordered
   * naturally corresponding to an m x m x m box where m = N^(1/3).  Boxes of
   * chunks are assigned to processors.
   *
   * NOTE: it is assumed that nprocs = 2^power and that the number of chunks in
   * each direction is divisible by the number of processors in each direction.
   */

  if (input_option == AZ_box) {

    /* determine the number of processors in each direction */

    if (proc == 0) {
       (void) printf("Input the dimensions of the processor cube\n\n");
       (void) printf("Enter the number of processors along x axis>");
       (void) fflush(stdout);
       scanf("%d",&proc_x);
       (void) printf("Enter the number of processors along y axis>");
       (void) fflush(stdout);
       scanf("%d",&proc_y);
       (void) printf("Enter the number of processors along z axis>");
       (void) fflush(stdout);
       scanf("%d",&proc_z);

       (void) printf("Input the grid dimensions\n\n");
       (void) printf("Enter the number of grid points along x axis>");
       (void) fflush(stdout);
       scanf("%d",&pts_x);
       (void) printf("Enter the number of grid points along y axis>");
       (void) fflush(stdout);
       scanf("%d",&pts_y);
       (void) printf("Enter the number of grid points along z axis>");
       (void) fflush(stdout);
       scanf("%d",&pts_z);
    }
    AZ_broadcast((char *) &proc_x, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &proc_y, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &proc_z, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &pts_x , sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &pts_y , sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) &pts_z , sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast((char *) NULL   , 0          , proc_config, AZ_SEND);
 
    total_pts_x = pts_x;
    total_pts_y = pts_y;


    if ( proc_x*proc_y*proc_z != nprocs) {
        if (proc == 0) {
          (void) printf("Error: %d x %d x %d != %d ",
 			 proc_x, proc_y, proc_z, nprocs);
          (void) printf(" (total number of processors)\n");
        }
	exit(1);
    }

    if ( pts_x * pts_y * pts_z != N ) {
        if (proc == 0) {
          (void) printf("Error: %d x %d x %d != %d ",
 			 pts_x, pts_y, pts_z, N);
          (void) printf(" (total number of grid points)\n");
        }
	exit(1);
    }
    if ( pts_x%proc_x != 0 ) {
        if (proc == 0) {
          (void) printf("Error: grid points along x axis are not an ");
          (void) printf("even multiple of processors\n");
	  (void) printf("       along x axis.");
        }
	exit(1);
    }
    if ( pts_y%proc_y != 0 ) {
        if (proc == 0) {
          (void) printf("Error: grid points along y axis is not an ");
          (void) printf("even multiple of processors\n");
	  (void) printf("       along y axis.");
        }
	exit(1);
    }
    if ( pts_z%proc_z != 0 ) {
        if (proc == 0) {
          (void) printf("Error: grid points along z axis is not an ");
          (void) printf("even multiple of processors\n");
	  (void) printf("       along z axis.");
        }
	exit(1);
    }
    pts_x = pts_x/proc_x;
    pts_y = pts_y/proc_y;
    pts_z = pts_z/proc_z;

    /* compute the number of elements per processor in each direction */

    *N_update = pts_x * pts_y * pts_z * chunk;
    if (!AZ_using_fortran) 
       *update     = (int *) AZ_allocate((*N_update)*sizeof(int));

    /* compute the lower left corner and the upper right corner */

    px = proc % proc_x;
    pz = (proc-px) / proc_x;
    py = pz % proc_y;
    pz = (pz-py) / proc_y;

    start_x = px * pts_x;
    end_x   = start_x + pts_x;
    start_y = py * pts_y;
    end_y   = start_y + pts_y;
    start_z = pz * pts_z;
    end_z   = start_z + pts_z;

    /* set update[] */

    count = 0;
    for (k = start_z; k < end_z; k++ ) {
      for (j = start_y; j < end_y; j++ ) {
        for (i = start_x; i < end_x; i++ ) {
          for (ii = 0; ii < chunk; ii++ ) {
            pt_number = (i + j * total_pts_x + k * total_pts_x * total_pts_y) * 
                            chunk + ii;
            (*update)[count++] = pt_number;
          }
        }
      }
    }
  }

  else if (input_option == AZ_linear) {

    /*
     * Figure out which chunks should be assigned to this processor for linear
     * partitioning.  This means that processor 0 is assigned the chunks
     * approximately corresponding to 0, ... , N/nprocs and processor 1 is
     * approximately assigned the chunks 1+N/nprocs to 2*N/nprocs.
     */

    t1 = N/nprocs;
    t2 = N - t1 * nprocs;

    if ( proc >= t2) t2 += (proc * t1);
    else {
      t1++;
      t2    = proc*t1;
    }
    *N_update = t1*chunk;
    t2   *= chunk;

    if (!AZ_using_fortran) 
       *update = (int *) AZ_allocate((*N_update+1)*sizeof(int));

    if (*update == NULL) {
      (void) AZ_printf_err( "Not enough space to allocate 'update'\n");
      exit(-1);
    }

    for (i = 0; i < *N_update; i++) (*update)[i] = i + t2;
  }

  else if (input_option == AZ_file) {

    /* read the update elements from the data file */

    type            = proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE)%AZ_NUM_MSGS +AZ_MSG_TYPE;
    type2           = proc_config[AZ_MPI_Tag];
    proc_config[AZ_MPI_Tag] = (type2+1-AZ_MSG_TYPE)%AZ_NUM_MSGS +AZ_MSG_TYPE;

    /*
     * Processor 0 reads the data file and distributes the lists to the other
     * processors.
     */

    t1 = 0;            /* total number of points distributed */
    if (proc == 0) {
      (void) AZ_printf_out("reading from file %s\n", datafile); AZ_flush_out();

      if ( (fp = fopen(datafile, "r")) == NULL) {
        (void) AZ_printf_err( "ERROR: file %s not found\n", datafile);
        exit(-1);
      }

      if (!AZ_using_fortran) *update = 0;
      allocated   = 0;
      for (i = nprocs - 1; i >= 0; i-- ) {

        /* read in list length and send to processor i  */

        fscanf(fp, "%d", &length);
        t1 += length;
        if (i != 0)
          mdwrap_write((void *) &length, sizeof(int), i, type, &cflag);

        /*
         * If this is the last list, we allocate the right amount of space and
         * keep the list instead of sending it off
         */

        if (i == 0) {
          *N_update = length;
          allocated       = 0;
        }

        /* allocate enough space for list */

        if (length > allocated ) {
          if ((*update != NULL) && (!AZ_using_fortran)) AZ_free(*update);
          allocated = length + 1;

          if (!AZ_using_fortran)
            *update = (int *) AZ_allocate(allocated*sizeof(int));
          if (*update == NULL) {
            (void) AZ_printf_err(
                           "Not enough space to allocate 'update'\n");
            exit(-1);
          }
        }

        /* read a list and send it off to proc i (if not last list) */

        for (j = 0; j < length; j++ ) fscanf(fp, "%d", *update + j);
        if (i != 0)
          mdwrap_write((void *) *update, length*sizeof(int), i, type2, &cflag);
      }
      fclose(fp);

      if (t1 != N*chunk) {
        (void) AZ_printf_err("AZ_read_update() found %d points in file\n", t1);
        (void) AZ_printf_err("%s instead of the requested %d\n", datafile,
                       N*chunk);
        exit(-1);
      }
    }

    else {

      /* read the update list from processor 0 */

      partner = 0;
      mdwrap_iread((void *) N_update, sizeof(int), &partner, &type, &request);
      mdwrap_wait((void *) N_update, sizeof(int), &partner, &type, &cflag, &request);

      if (!AZ_using_fortran)
        *update = (int *) AZ_allocate((*N_update+1)*sizeof(int));
      if (*update == NULL)  {
        (void) AZ_printf_err( "Not enough space to allocate 'update'\n");
        exit(-1);
      }

      partner = 0;
      mdwrap_iread((void *) *update, *N_update * sizeof(int), &partner, &type2,
            &request);
      mdwrap_wait((void *) *update, *N_update * sizeof(int), &partner, &type2,
            &cflag, &request);
    }

    AZ_sort(*update, *N_update, NULL, NULL);

    /* check that element '0' is contained on 1 processor. That is,  */
    /* make sure the user has numbered from 0 to n-1 instead of from */
    /* 1 to n                                                        */
    check = 0;
    if ( (*N_update > 0) && ((*update)[0] == 0) ) check = 1;
    check = AZ_gsum_int(check, proc_config);
    if (check != 1) {
       if (proc == 0) {
          (void) AZ_printf_err("Error: In AZ_read_update(), the %s", datafile);
          (void) AZ_printf_err("file does not contain\n       one ");
          (void) AZ_printf_err("occurance of row 0. Make sure that rows are");
          (void) AZ_printf_err(" numbered\n       from 0 to n-1.\n");
       }
       exit(1);
    }


  }
  else {
    (void) AZ_printf_err("Unknown input option (%d) in AZ_read_update()\n",
                   input_option);
    exit(1);
  }


} /* AZ_input_update */


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
extern int AZ_read_external(int N_external, int external[],
                     int **extern_proc, FILE *fp, int proc_config[]);

int AZ_read_external(int N_external, int external[],
                     int **extern_proc, FILE *fp, int proc_config[])

/*******************************************************************************

  Processor 0 reads from the file pointer 'fp'.  This file should contain nprocs
  lists.  Each list consists of a number telling how many pairs (global index,
  proc) are in the list followed by the pairs.  The first list is sent to
  processor 'nprocs-1', the second list is sent to processor 'nprocs-2', etc.
  Each processor takes its list and assigns a processor to each external point.
  Routine returns 1.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     int, value 1.
  ============

  Parameter list:
  ===============

  N_external:      Number of external elements.

  external:        List of external elements needed by this processor.

  extern_proc:     On Output, the value of 'extern_proc[i]' is updating
                   processor for 'external[i]'.

  fp:              On input, openned file pointer where data will be read.

  proc_config:     proc_config[AZ_node] is node number.
                   proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int   i;
  int   loc;
  int  *temp_buffer, length;
  int   allocated = -1, cflag, j;
  int   partner;
  int   proc, nprocs;
  int   type, type2;
  MPI_AZRequest request;


  char *yo = "AZ_read_external: ";


  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;
  type2           = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type2+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  if (proc == 0) {
    temp_buffer = NULL;
    for (i = nprocs-1; i >= 0; i-- ) {

      /* read list length send to proc i */

      length = -1;
      fscanf (fp, "%d", &length);

      if (i == 0) {
        if (length != N_external) {
          (void) AZ_printf_err( "%sERROR: %d: number of extern elements in "
                         "file", yo, proc);
          (void) AZ_printf_err( " is not what I\n    found in matrix"
                         "(%d vs %d)\n", length, N_external);
          exit(-1);
        }
      }
      else {
        mdwrap_write((char *) &length, sizeof(int), i, type, &cflag);
      }

      /* allocate enough space to hold list */

      length *= 2;
      if (length > allocated ) {
        if (temp_buffer != NULL)  AZ_free(temp_buffer);

        allocated = length + 1;

        temp_buffer = (int *) AZ_allocate(allocated*sizeof(int));
        if (temp_buffer == NULL) {
          (void) AZ_printf_err( "%sERROR: not enough dynamic memory to "
                         "allocate 'temp_buffer'\n", yo);
          exit(-1);
        }
      }

      /* read list and send it to proc i */

      for (j = 0; j < length; j++ ) fscanf(fp, "%d", &(temp_buffer[j]));

      if (i != 0) {
        mdwrap_write((char *) temp_buffer, length*sizeof(int), i, type2, &cflag);
      }
    }
  }

  else {

    /* receive list from processor 0 */

    partner = 0;
    mdwrap_iread((void *) &length, sizeof(int), &partner, &type, &request);
    mdwrap_wait((void *) &length, sizeof(int), &partner, &type, &cflag,&request);

    if (length != N_external) {
      (void) AZ_printf_err( "%sERROR: %d:number of extern elements in file is",
                     yo, proc);
      (void) AZ_printf_err( " not what I\n   found in the matrix "
                     "(%d vs %d)\n", length, N_external);
      exit(-1);
    }

    length *= 2;

    temp_buffer = (int *) AZ_allocate((length+1)*sizeof(int));
    if (temp_buffer == NULL) {
      (void) AZ_printf_err( "%sERROR: not enough dynamic memory to allocate "
                      "'temp'\n", yo);
      exit(-1);
    }

    partner = 0;
    mdwrap_iread((void *) temp_buffer, length*sizeof(int), &partner, &type2, 
                  &request);
    mdwrap_wait((void *) temp_buffer, length*sizeof(int), &partner, &type2, 
                  &cflag, &request);
  }

  /*
   * Use list of pairs (global index, processor) to assign processors to
   * external elements.
   */

  if (!AZ_using_fortran) {
     *extern_proc = (int *) AZ_allocate((N_external+1)*sizeof(int));
  }
  if (*extern_proc == NULL) {
    (void) AZ_printf_err( "%sERROR: not enough dynamic memory for external "
                   "procs\n", yo);
    exit(-1);
  }
  if (!AZ_using_fortran) {
     for (i = 0 ; i < N_external ; i++) extern_proc[i] = 0;
  }

  for (i = 0 ; i < N_external ; i++ ) {
    loc = AZ_find_index(temp_buffer[2*i], external, N_external);
    if (loc == -1) {
      (void) AZ_printf_err( "%sERROR: external point (%d) in input \n", yo,
                     temp_buffer[2*i]);
      (void) AZ_printf_err( "       file was not found in the matrix \n");
                     exit(-1);
    }

    (*extern_proc)[loc] = temp_buffer[2*i+1];
  }

  AZ_free(temp_buffer);

  return 1;

} /* AZ_read_external */

/*******************************************************************************
 Program for reading in an input file and generating a distributed matrix
 across the processors.
*******************************************************************************/

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_read_msr_matrix(int update[], double **val, int **bindx, int N_update,
                        int proc_config[])

/*******************************************************************************

  Read an input file and create a distributed matrix.  Effectively, processor 0
  reads the input file.  If the new row to be added resides in processor 0's
  update, it is added to processor 0's matrix.  Otherwise, processor 0 checks
  which processor has requested this row and sends the row to this processor so
  that it can add it in its local matrix.

  The form of the input file is as follows:

         Number_of_matrix_rows
         col_num1 entry1 col_num2 entry2
         col_num3 entry3 -1
         col_num4 entry4 col_num5 entry5
         col_num6 entry6 -1

  This input corresponds to 2 rows (row 0 and row 1).  Row 0 contains entry1 in
  column col_num1, entry2 in column col_num2, and entry3 in column col_num3.
  Row 1 contains entry4 in column col_num4, entry5 in column col_num5, and
  entry6 in column col_num6.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  update:          'update[i]' is the global index of the ith point to be
                   updated on this processor.

  val, bindx:      MSR matrix arrays. See User's Guide.

  N_update:        Number of elements updated on this processor.

  nz_pointer:      On output, indicates the number of matrix nonzeros stored on
                   this processor.

  proc_config:     proc_config[AZ_node] is node number.
                   proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int    nz_ptr;
  char  *str;
  int    i,j,k, jjj;
  int    current;
  int    st, pnode;
  int    temp, *lil;
  double dtemp;
  int   *requests;
  int    total;
  FILE  *dfp = NULL;
  int    proc, nprocs;
  int    type, type2;
  unsigned int buf_len = 1000;
  int    msr_len;
  int    count = 0;
  int    kkk, m_one = -1, need_request = 0;
  int    column0 = 0;
  MPI_AZRequest request;  /* Message handle */
  int    totalN;

  char   *tchar;

  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;
  type2           = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type2+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  totalN = AZ_gsum_int(N_update, proc_config);
  str    = (char *) AZ_allocate((buf_len+1)*sizeof(char));
  if (str == NULL) {
    AZ_printf_err("ERROR: NOT enough dynamic memory in AZ_read_msr_matrix\n");
    exit(-1);
  }
  msr_len = 8*N_update+2;
  if (!AZ_using_fortran) {
    *bindx = (int *)    AZ_allocate(msr_len*sizeof(int));
    *val   = (double *) AZ_allocate(msr_len*sizeof(double));
  }

  if (*val == NULL) {
    (void) AZ_printf_err(
                   "ERROR: NOT enough dynamic memory in AZ_read_msr_matrix\n");
    exit(-1);
  }
  if (!AZ_using_fortran) {
     for (i = 0 ; i < msr_len; i++) (*bindx)[i] = 0;
     for (i = 0 ; i < msr_len; i++) (*val)[i] = 0;
  }

  nz_ptr      = N_update + 1;
  (*bindx)[0] = nz_ptr;
  current     = 0;

  if (proc != 0) {

    /*
     * Send requests to processor 0.  When the response is received add the
     * corresponding row to the matrix and make another request.  When all the
     * requests are done, send -1 as a request to signal processor 0 that we are
     * finished.
     */

    for (i = 0; i < N_update; i++ ) {
      mdwrap_write((void *) &(update[i]), sizeof(int), 0, type, &st);
      pnode = 0;
      mdwrap_iread(str, buf_len, &pnode, &type2, &request);
      j = mdwrap_wait(str, buf_len, &pnode, &type2, &st, &request);
      while (j == sizeof(int)) {
        lil = (int *) str;
        buf_len = (unsigned int) lil[0];
        str = (char *) AZ_realloc(str,buf_len*sizeof(char));
        if (str == 0) {
          (void) AZ_printf_err(
                         "ERROR: NoT enough dynamic memory in AZ_read_msr()\n");
          exit(-1);
        }
        mdwrap_iread(str, buf_len, &pnode, &type2, &request);
        j = mdwrap_wait(str, buf_len, &pnode, &type2, &st, &request);
      }

      AZ_add_new_row(update[i], &nz_ptr, &current, val, bindx, str, dfp,
                     &msr_len,&column0);
    }

    temp = -1;
    mdwrap_write((void *) &temp, sizeof(int), 0, type, &st);
  }

  else {
#ifdef binary
    dfp = fopen(".data", "rb");
#else
    dfp = fopen(".data", "r");
#endif
    if (dfp == NULL) {
      (void) AZ_printf_err( "ERROR: Matrix file '.data' not found\n");
      exit(-1);
    }
    (void) AZ_printf_out( " reading matrix (current version is very slow) .");
    (void) AZ_flush_out();

    /* read in number of blks */
    /*
      fscanf(dfp, "%d", &total);
      for (i = 0; i <= total; i++ ) fscanf(dfp, "%d", &temp);
      */

    /* read past cnptr info (not used) */

#ifdef binary
    kkk = fread(&total, sizeof(int), 1, dfp);
#else
    kkk = fscanf(dfp, "%d", &total);  /* read in number of elements */
#endif

    if (kkk <= 0) {
       (void) AZ_printf_err("\nfile '.data' is empty\n");
       exit(1);
    }

    if (total != totalN) {
       (void) AZ_printf_err("\nError:file '.data' contains %d rows ",total);
       (void) AZ_printf_err("while the user's input\n     requires %d rows.\n",
                      totalN);
       exit(1);
    }

    /* get the first requests from all the processors */

    requests    = (int *) AZ_allocate(nprocs*sizeof(int));
    requests[0] = -1;
    for (i = 1; i < nprocs; i++) {
      pnode = -1;
      mdwrap_iread((void *) &temp, sizeof(int), &pnode, &type, &request);
      mdwrap_wait((void *) &temp, sizeof(int), &pnode, &type, &st, &request);
      requests[pnode] = temp;
    }

    /*
     * Go through all the rows, for those rows that we own add them to our local
     * matrix.  Otherwise, read the row into a string, determine which processor
     * has requested the row, send the string to this processor, and receive
     * another request from this processor.
     */

    for (i = 0; i < total; i++) {
      count++;
      if (count%1000 == 0) {
        (void) AZ_printf_out( ".");
        (void) AZ_flush_out();
      }
      if ((current < N_update) && (i == update[current])) {
        AZ_add_new_row(i, &nz_ptr, &current, val, bindx, 0, dfp, &msr_len,
		       &column0);
      }
      else {
#ifdef binary
        kkk = fread(&temp, sizeof(int), 1, dfp);
#else
        kkk = fscanf(dfp, "%d", &temp);
#endif
        if (kkk <= 0) {
           (void) AZ_printf_err("\nError: AZ_read_msr(), end-of-file reached");
           (void) AZ_printf_err(" while reading row %d.\n",i);
           exit(1);
        }
        if (temp == 0) column0 = 1;

        j = 0;

        while (temp != -1) {
#ifdef binary
          kkk = fread(&dtemp, sizeof(double), 1, dfp);
#else
          kkk = fscanf(dfp, "%lf", &dtemp);
#endif
          if (kkk <= 0) {
           (void) AZ_printf_err("\nError: AZ_read_msr(), end-of-file reached");
           (void) AZ_printf_err(" while reading row %d.\n",i);
           exit(1);
          }

          if (j + 30 > (int) buf_len) {
            buf_len = 2*buf_len + 30;
            str = (char *) AZ_realloc(str,buf_len*sizeof(char));

            if (str == 0) {
              (void) AZ_printf_err("ERROR: Not Enough dynamic memory in "
                             "AZ_read_msr()\n");
              exit(-1);
            }
            if (need_request != 0)  {
               mdwrap_iread((void *) &(requests[need_request]), 
		        sizeof(int), &need_request,&type,&request);
               mdwrap_wait((void *) &(requests[need_request]), 
		        sizeof(int), &need_request,&type,&st,&request);
               need_request = 0;
            }
            for (jjj = 1; jjj < nprocs; jjj++) {
              if (requests[jjj] != -1) 
                 mdwrap_write((void *) &buf_len, sizeof(int), jjj, type2, &st);
	    }
          }

          /* put 'temp' and 'dtemp' into 'str' so that they can be sent */
          /* to another processor.                                      */

          tchar = (char *) &temp;
          for (kkk = 0 ; kkk < (int)sizeof(int); kkk++) str[j+kkk] = tchar[kkk];
          j += sizeof(int);
          tchar = (char *) &dtemp;
          for (kkk = 0 ; kkk < (int) sizeof(double); kkk++ ) 
             str[j+kkk] = tchar[kkk];
          j += sizeof(double);
#ifdef binary
          kkk = fread(&temp, sizeof(int), 1, dfp);
#else
          kkk = fscanf(dfp, "%d", &temp);
#endif
          if (kkk <= 0) {
           (void) AZ_printf_err("\nError: AZ_read_msr(), end-of-file reached");
           (void) AZ_printf_err(" while reading row %d.\n",i);
           exit(1);
          }
          if (temp == 0) column0 = 1;
        }
        tchar = (char *) &m_one;
        for (kkk = 0 ; kkk < (int)sizeof(int) ; kkk++ ) str[j+kkk] = tchar[kkk];
        j += sizeof(int);

        k = 0;
        if (need_request != 0)  {
           mdwrap_iread((void *) &(requests[need_request]), sizeof(int), 
		    &need_request,&type,&request);
           mdwrap_wait((void *) &(requests[need_request]), sizeof(int), 
		    &need_request,&type,&st, &request);
           need_request = 0;
        }

        while ((k < nprocs) && (requests[k] != i)) k++;
        if (k == nprocs) {
           (void) AZ_printf_err("\nError: AZ_read_msr(), input file contains");
           (void) AZ_printf_err(" a row (%d)\n       that is not ",i);
           (void) AZ_printf_err("assigned to any processor?\n");
           exit(1);
        }
        mdwrap_write((void *) str, (unsigned int) j, k, type2, &st);
        need_request = k;  /* read is deferred until we will need it */
      }
    }
    if (need_request != 0)  {
       mdwrap_iread((void *) &(requests[need_request]), sizeof(int), 
		&need_request,&type,&request);
       mdwrap_wait((void *) &(requests[need_request]), sizeof(int), 
		&need_request,&type,&st,&request);
       need_request = 0;
    }

    /* at this point all the requests should be -1 */

    for (k = 0 ; k < nprocs ; k++ ) {
       if (requests[k] != -1) {
          (void) AZ_printf_err("\nError: AZ_read_msr(), processor %d ",k);
          (void) AZ_printf_err("requested  but never received\n       ");
          (void) AZ_printf_err("matrix row %d. Check '.data' file.\n",
                         requests[k]);
          exit(1);
       }
    }

    if (column0 == 0) {
       (void) AZ_printf_err("\nWARNING: AZ_read_msr(), column 0 contains ");
       (void) AZ_printf_err("no off-diagonal elements.\n         Are you ");
       (void) AZ_printf_err("sure that you numbered the matrix rows/columns");
       (void) AZ_printf_err(" from\n         0 to n-1 (and not 1 to n).\n");
    }


    AZ_free(requests);
    fclose(dfp);
    (void) AZ_printf_out( "\n");
  }

  AZ_free(str);

} /* AZ_read_msr_matrix */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void AZ_input_msr_matrix(char datafile[], int update[], double **val, int **bindx, 
												 int N_update, int proc_config[])

/*******************************************************************************

Exactly the same as AZ_read_msr_matrix except it reads that data in from a 
file specified by the input argument datafile instead from a file called
.data

*******************************************************************************/

{

  /* local variables */

  int    nz_ptr;
  char  *str;
  int    i,j,k, jjj;
  int    current;
  int    st, pnode;
  int    temp, *lil;
  double dtemp;
  int   *requests;
  int    total;
  FILE  *dfp = NULL;
  int    proc, nprocs;
  int    type, type2;
  unsigned int buf_len = 1000;
  int    msr_len;
  int    count = 0;
  int    kkk, m_one = -1, need_request = 0;
  int    column0 = 0;
  MPI_AZRequest request;  /* Message handle */
  int    totalN;

  char   *tchar;

  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;
  type2           = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type2+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  totalN = AZ_gsum_int(N_update, proc_config);
  str    = (char *) AZ_allocate((buf_len+1)*sizeof(char));
  if (str == NULL) {
    AZ_printf_out("ERROR: NOT enough dynamic memory in AZ_input_msr_matrix\n");
    exit(-1);
  }
  msr_len = 8*N_update+2;
  if (!AZ_using_fortran) {
    *bindx = (int *)    AZ_allocate(msr_len*sizeof(int));
    *val   = (double *) AZ_allocate(msr_len*sizeof(double));
  }

  if (*val == NULL) {
    (void) AZ_printf_err(
                   "ERROR: NOT enough dynamic memory in AZ_input_msr_matrix\n");
    exit(-1);
  }
  if (!AZ_using_fortran) {
     for (i = 0 ; i < msr_len; i++) (*bindx)[i] = 0;
     for (i = 0 ; i < msr_len; i++) (*val)[i] = 0;
  }

  nz_ptr      = N_update + 1;
  (*bindx)[0] = nz_ptr;
  current     = 0;

  if (proc != 0) {

    /*
     * Send requests to processor 0.  When the response is received add the
     * corresponding row to the matrix and make another request.  When all the
     * requests are done, send -1 as a request to signal processor 0 that we are
     * finished.
     */

    for (i = 0; i < N_update; i++ ) {
      mdwrap_write((void *) &(update[i]), sizeof(int), 0, type, &st);
      pnode = 0;
      mdwrap_iread(str, buf_len, &pnode, &type2, &request);
      j = mdwrap_wait(str, buf_len, &pnode, &type2, &st, &request);
      while (j == sizeof(int)) {
        lil = (int *) str;
        buf_len = (unsigned int) lil[0];
        str = (char *) AZ_realloc(str,buf_len*sizeof(char));
        if (str == 0) {
          (void) AZ_printf_err(
                         "ERROR: Not enough dynamic memory in AZ_input_msr()\n");
          exit(-1);
        }
        mdwrap_iread(str, buf_len, &pnode, &type2, &request);
        j = mdwrap_wait(str, buf_len, &pnode, &type2, &st, &request);
      }

      AZ_add_new_row(update[i], &nz_ptr, &current, val, bindx, str, dfp,
                     &msr_len,&column0);
    }

    temp = -1;
    mdwrap_write((void *) &temp, sizeof(int), 0, type, &st);
  }

  else {
#ifdef binary
    dfp = fopen(datafile, "rb");
#else
    dfp = fopen(datafile, "r");
#endif
    if (dfp == NULL) {
      (void) AZ_printf_err( "ERROR: Matrix data file %s not found\n", datafile);
      exit(-1);
    }
    (void) AZ_printf_out( "Reading matrix from %s "
                           "(current version is very slow) .",datafile);
    (void) AZ_flush_out();

    /* read in number of blks */
    /*
      fscanf(dfp, "%d", &total);
      for (i = 0; i <= total; i++ ) fscanf(dfp, "%d", &temp);
      */

    /* read past cnptr info (not used) */

#ifdef binary
    kkk = fread(&total, sizeof(int), 1, dfp);
#else
    kkk = fscanf(dfp, "%d", &total);  /* read in number of elements */
#endif

    if (kkk <= 0) {
       (void) AZ_printf_err("data file %s is empty\n", datafile);
       exit(1);
    }

    if (total != totalN) {
       (void) AZ_printf_err("\nError: data file %s contains %d rows ",datafile, total);
       (void) AZ_printf_err("while the user's input\n     requires %d rows.\n",
                      totalN);
       exit(1);
    }

    /* get the first requests from all the processors */

    requests    = (int *) AZ_allocate(nprocs*sizeof(int));
    requests[0] = -1;
    for (i = 1; i < nprocs; i++) {
      pnode = -1;
      mdwrap_iread((void *) &temp, sizeof(int), &pnode, &type, &request);
      mdwrap_wait((void *) &temp, sizeof(int), &pnode, &type, &st, &request);
      requests[pnode] = temp;
    }

    /*
     * Go through all the rows, for those rows that we own add them to our local
     * matrix.  Otherwise, read the row into a string, determine which processor
     * has requested the row, send the string to this processor, and receive
     * another request from this processor.
     */

    for (i = 0; i < total; i++) {
      count++;
      if (count%1000 == 0) {
        (void) AZ_printf_out( ".");
        (void) AZ_flush_out();
      }
      if ((current < N_update) && (i == update[current])) {
        AZ_add_new_row(i, &nz_ptr, &current, val, bindx, 0, dfp, &msr_len,
		       &column0);
      }
      else {
#ifdef binary
        kkk = fread(&temp, sizeof(int), 1, dfp);
#else
        kkk = fscanf(dfp, "%d", &temp);
#endif
        if (kkk <= 0) {
           (void) AZ_printf_err("\nError: AZ_input_msr(), end-of-file reached");
           (void) AZ_printf_err(" while reading row %d.\n",i);
           exit(1);
        }
        if (temp == 0) column0 = 1;

        j = 0;

        while (temp != -1) {
#ifdef binary
          kkk = fread(&dtemp, sizeof(double), 1, dfp);
#else
          kkk = fscanf(dfp, "%lf", &dtemp);
#endif
          if (kkk <= 0) {
           (void) AZ_printf_err("\nError: AZ_input_msr(), end-of-file reached");
           (void) AZ_printf_err(" while reading row %d.\n",i);
           exit(1);
          }

          if (j + 30 > (int) buf_len) {
            buf_len = 2*buf_len + 30;
            str = (char *) AZ_realloc(str,buf_len*sizeof(char));

            if (str == 0) {
              (void) AZ_printf_err("ERROR: Not Enough dynamic memory in "
                             "AZ_input_msr()\n");
              exit(-1);
            }
            if (need_request != 0)  {
               mdwrap_iread((void *) &(requests[need_request]), 
		        sizeof(int), &need_request,&type,&request);
               mdwrap_wait((void *) &(requests[need_request]), 
		        sizeof(int), &need_request,&type,&st,&request);
               need_request = 0;
            }
            for (jjj = 1; jjj < nprocs; jjj++) {
              if (requests[jjj] != -1) 
                 mdwrap_write((void *) &buf_len, sizeof(int), jjj, type2, &st);
	    }
          }

          /* put 'temp' and 'dtemp' into 'str' so that they can be sent */
          /* to another processor.                                      */

          tchar = (char *) &temp;
          for (kkk = 0 ; kkk < (int)sizeof(int); kkk++) str[j+kkk] = tchar[kkk];
          j += sizeof(int);
          tchar = (char *) &dtemp;
          for (kkk = 0 ; kkk < (int) sizeof(double); kkk++ ) 
             str[j+kkk] = tchar[kkk];
          j += sizeof(double);
#ifdef binary
          kkk = fread(&temp, sizeof(int), 1, dfp);
#else
          kkk = fscanf(dfp, "%d", &temp);
#endif
          if (kkk <= 0) {
           (void) AZ_printf_err("\nError: AZ_input_msr(), end-of-file reached");
           (void) AZ_printf_err(" while reading row %d.\n",i);
           exit(1);
          }
          if (temp == 0) column0 = 1;
        }
        tchar = (char *) &m_one;
        for (kkk = 0 ; kkk < (int)sizeof(int) ; kkk++ ) str[j+kkk] = tchar[kkk];
        j += sizeof(int);

        k = 0;
        if (need_request != 0)  {
           mdwrap_iread((void *) &(requests[need_request]), sizeof(int), 
		    &need_request,&type,&request);
           mdwrap_wait((void *) &(requests[need_request]), sizeof(int), 
		    &need_request,&type,&st, &request);
           need_request = 0;
        }

        while ((k < nprocs) && (requests[k] != i)) k++;
        if (k == nprocs) {
           (void) AZ_printf_err("\nError: AZ_input_msr(), input file contains");
           (void) AZ_printf_err(" a row (%d)\n       that is not ",i);
           (void) AZ_printf_err("assigned to any processor?\n");
           exit(1);
        }
        mdwrap_write((void *) str, (unsigned int) j, k, type2, &st);
        need_request = k;  /* read is deferred until we will need it */
      }
    }
    if (need_request != 0)  {
       mdwrap_iread((void *) &(requests[need_request]), sizeof(int), 
		&need_request,&type,&request);
       mdwrap_wait((void *) &(requests[need_request]), sizeof(int), 
		&need_request,&type,&st,&request);
       need_request = 0;
    }

    /* at this point all the requests should be -1 */

    for (k = 0 ; k < nprocs ; k++ ) {
       if (requests[k] != -1) {
          (void) AZ_printf_err("\nError: AZ_input_msr(), processor %d ",k);
          (void) AZ_printf_err("requested  but never received\n       ");
          (void) AZ_printf_err("matrix row %d. Check data file.\n",
                         requests[k]);
          exit(1);
       }
    }

    if (column0 == 0) {
       (void) AZ_printf_err("\nWARNING: AZ_input_msr(), column 0 contains ");
       (void) AZ_printf_err("no off-diagonal elements.\n         Are you ");
       (void) AZ_printf_err("sure that you numbered the matrix rows/columns");
       (void) AZ_printf_err(" from\n         0 to n-1 (and not 1 to n).\n");
    }


    AZ_free(requests);
    fclose(dfp);
    (void) AZ_printf_out( "\n");
  }

  AZ_free(str);

} /* AZ_input_msr_matrix */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_add_new_row(int therow, int *nz_ptr, int *current, double **val,
                    int **bindx, char *input, FILE *dfp, int *msr_len,
		    int *column0)

/*******************************************************************************

  Add a new row to the matrix.  If input = 0, the new matrix is read from file
  pointer dfp.  Otherwise, it is read from the string 'input'.  The form of the
  input is as follows:

         col_num1 entry1 col_num2 entry2
         col_num3 entry3 -1

  On output, val[] and  bindx[] are updated to incorporate the new row.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  therow:          The global index of the row being added.

  nz_ptr:          The next available space in val[] and bindx[] for storing
                   nonzero offdiagonals.

  current:         The next available space in a[] to store the matrix diagonal.

  val, bindx:      MSR matrix arrays that will be updated to incorporate new
                   row. See User's Guide.

  input:           Contains the information describing the row to be added (if
                   input == 0, the row information is read from standard input).

*******************************************************************************/

{

  /* local variables */

  int    old_nz;
  double sum = 0.0;
  int    temp;
  double dtemp;
  int    k = 0, kk;
  char   *tchar;

  /**************************** execution begins ******************************/

  old_nz = *nz_ptr;

  if (input == 0) { 
#ifdef binary
    kk  = fread(&temp,sizeof(int),1,dfp);
#else
    kk  = fscanf(dfp, "%d", &temp);
#endif

    if (kk <= 0) {
         (void) AZ_printf_err("\nError: format error in '.data' file ");
         (void) AZ_printf_err("on row '%d'\n",*current);
         (void) AZ_printf_err("      This can be caused if exponents are\n");
         (void) AZ_printf_err("      given using 'D' instead of 'E'. \n");
       exit(1);
    }
    if (temp == 0) *column0 = 1;
  }
  else {
    tchar = (char *) &temp;
    for (kk = 0 ; kk < (int) sizeof(int) ; kk++ ) tchar[kk] = input[kk];
    k    += sizeof(int);
  }

  while (temp != -1) {
    if (input == 0) {
#ifdef binary
       kk = fread(&dtemp, sizeof(double), 1, dfp);
#else
       kk = fscanf(dfp, "%lf", &dtemp);
#endif
       if (kk <= 0) {
         (void) AZ_printf_err("\nError: format error in '.data' file ");
         (void) AZ_printf_err("on row '%d'\n",*current);
         (void) AZ_printf_err("       This can be caused if exponents are\n");
         (void) AZ_printf_err("       given using 'D' instead of 'E'. \n");
         exit(1);
       }
    }
    else {
      tchar = (char *) &dtemp;
      for (kk = 0 ; kk < (int) sizeof(double) ; kk++ ) tchar[kk] = input[k+kk];
      k += sizeof(double);
    }

    if (temp == therow) sum = dtemp;
    else {
      if (*nz_ptr >= *msr_len) {
        *msr_len = (int) ( 1.5 * (double) *msr_len);
        if (!AZ_using_fortran) {
          *bindx = (int *) AZ_realloc(*bindx,*msr_len*sizeof(int));
          *val   = (double *) AZ_realloc(*val,*msr_len*sizeof(double));
        }
        if (*val == 0) {
          (void) AZ_printf_err(
                         "ERROR: Not enough dynamic memory in AZ_read_msr()\n");
          exit(-1);
        }
      }
      (*bindx)[*nz_ptr] =  temp;
      (*val)[*nz_ptr]   = dtemp;
      (*nz_ptr)++;
    }

    if (input == 0) {
#ifdef binary
       kk  = fread(&temp,sizeof(int),1,dfp);
#else
       kk = fscanf(dfp, "%d", &temp);
#endif
       if (kk <= 0) {
         (void) AZ_printf_err("\nError: format error in '.data' file ");
         (void) AZ_printf_err("on row '%d'\n",*current);
         (void) AZ_printf_err("       This can be caused if exponents are\n");
         (void) AZ_printf_err("       given using 'D' instead of 'E'. \n");
         exit(1);
       }
       if (temp == 0) *column0 = 1;
    }
    else {
      tchar = (char *) &temp;
      for (kk = 0 ; kk < (int) sizeof(int) ; kk++ ) tchar[kk] = input[kk+k];
      k    += sizeof(int);
    }
  }

  (*val)[*current]     = sum;
  (*bindx)[*current+1] = (*bindx)[*current] + (*nz_ptr - old_nz);
  (*current)++;

} /* AZ_add_new_row */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_check_update(int update[], int N_update, int proc_config[])

/*******************************************************************************

  See that all the update elements are there by doing a cheap test.
  Specifically,
  1) calculate the total number of update elements
  2) Assuming that they range from 0 .... total-1 compute what the sum of the
     update elements should be mod a certain number (11576) to avoid overflow.
  3) calculate the sum of the update elements and see if they match.  If they
  match we assume that the update elements are okay.
  4) If they don't match then we do a long test where each processor will send
  all of its elements to processor 0 (sending just 1 point at a time) and
  processor 0 will determine which point is missing.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  update:          'update[i]' is the global index of the ith point to be
                   updated on this processor.

  N_update:        Number of elements updated on this processor.

  proc_config:     proc_config[AZ_node] is node number.
                   proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int current = 0, st, pnode;
  int i, j, k, temp;
  int total, *requests;
  int proc, nprocs;
  int type, type2;
  int sum, real_sum = 0;
  MPI_AZRequest request;  /* Message handle */



  /**************************** execution begins ******************************/

  AZ__MPI_comm_space_ok();
  proc   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;
  type2           = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type2+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* compute the total number of update elements */

  total = AZ_gsum_int(N_update, proc_config);

  /* compute theoretical sum (mod a number) */

  i    = total % (2 * 11576);
  j    = (total - 1) % (2 * 11576);
  sum  = i * j;
  sum /= 2;
  sum  = sum % (11576);

  /* compute actual sum */

  for (i = 0; i < N_update; i++) {
    real_sum += update[i];
    real_sum  = real_sum % 11576;
  }

  real_sum = AZ_gsum_int(real_sum, proc_config);
  real_sum = real_sum % (11576);

  if (real_sum != sum) {
    if (proc == 0) {
      (void) AZ_printf_err( "ERROR: update elements test failed\n");
      (void) AZ_printf_err( "       theor. sum of update = %d\n", sum);
      (void) AZ_printf_err( "       actual sum of update = %d\n",
                     real_sum);
    }
  }
  else return;

  if (proc == 0)
    (void) AZ_printf_err( "Performing a detailed (slow) check\n");

  /* Now do a slow and detailed check */

  if (proc != 0) {

    /*
     * Send requests to processor 0.  When the response is received make another
     * request.  When all the requests are done, send -1 as a request to signal
     * processor 0 that we are finished.
     */

    for (i = 0; i < N_update; i++ ) {
      mdwrap_write((void *) &(update[i]), sizeof(int), 0, type, &st);
      pnode = 0;
      mdwrap_iread((void *) &temp, sizeof(int), &pnode, &type2, &request);
      mdwrap_wait((void *) &temp, sizeof(int), &pnode, &type2, &st, &request);
    }
    temp = -1;
    mdwrap_write((void *) &temp, sizeof(int), 0, type, &st);
  }
  else {

    /* get the first requests from all the processors */

    requests    = (int *) AZ_allocate(nprocs*sizeof(int));
    requests[0] = -1;

    for (i = 1; i < nprocs; i++ ) {
      pnode = -1;
      mdwrap_iread((void *) &temp, sizeof(int), &pnode, &type, &request);
      mdwrap_wait((void*) &temp, sizeof(int), &pnode, &type, &st, &request);
      requests[pnode] = temp;
    }

    /*
     * Go through all the rows, for those rows that we own add them to our local
     * matrix.  Otherwise, read the row into a string, determine which processor
     * has requested the row, send the string to this processor, and receive
     * another request from this processor.
     */

    for (i = 0; i < total; i++ ) {
      if ((current < N_update) && (i == update[current])) {
        current++;
      }
      else {
        k = 0;
        while ((requests[k] != i) && (k < nprocs)) k++;

        if (k == nprocs) {
          (void) AZ_printf_err( "ERROR: A grid point (%d) was not found\n",
                         current);
          (void) AZ_printf_err( "       among the update elements\n");
          exit(-1);
        }

        mdwrap_write((void *) &temp, sizeof(int), k, type2, &st);
        mdwrap_iread((void *) &temp, sizeof(int), &k, &type, &request);
        mdwrap_wait((void *) &temp, sizeof(int), &k, &type, &st, &request);
        requests[k] = temp;
      }
    }

    AZ_free(requests);
  }

} /* AZ_check_update */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_vb2msr(int m, double val[], int indx[], int bindx[], int rpntr[],
               int cpntr[], int bpntr[], double msr_val[], int msr_bindx[])

/*******************************************************************************

  This routine converts the VBR matrix (stored in the arrays val[], indx[],
  bindx[], rpntr[], cpntr[], bpntr[]) into an MSR matrix (stored in the arrays:
  msr_val, msr_bindx)

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  m:               Number of (block) rows in the matrix.

  val:             Matrix A in sparse format (VBR).

  rnptr,
  bindx,
  indx,
  bnptr,
  cnptr:           Sparse matrix arrays. See User's Guide.
                   On input, the matrix corresponds to the initial ordering
                   (e.g. row i corresponds to global row update[i]).
                   On output, the matrix rows and columns are renumbered to
                   correspond to the ordering given by 'update_index' and
                   'extern_index'. (e.g. global row update[i] now appears
                   as row update_index[i] in the matrix).

  msr_val,
  msr_bindx:       On output, msr_val[] contains the nonzeros of the matrix
                   while msr_bindx[] contains pointer information of the
                   matrix in C-MSR format. That is,
                       msr_val[i] = A(i,i)               for i < N
                       msr_val[i] = A(j,msr_bindx[i])    for i > N
                   where A() refers to the matrix, N is the number of (point)
                   rows in the matrix, and j is given implicitly by
                       msr_bindx[j] <= i  < msr_bindx[j+1].

*******************************************************************************/

{

  /* local variables */

  register int indx_ptr, irow, iblk_row, jcol, jblk_col, icol_blk, iblk = 0;
  int          num_blk_rows, num_col_blks, num_blk_cols;
  int          nz_counter;
  int          realrow, realcol;

  /**************************** execution begins ******************************/

  /* loop over the block rows */

  nz_counter = rpntr[m] + 1;
  msr_bindx[0]    = nz_counter;

  for (iblk_row = 0; iblk_row < m; iblk_row++) {

    /* find out how many rows are in this block row */

    num_blk_rows = rpntr[iblk_row+1] - rpntr[iblk_row];

    /* find out how many block columns are in this block row */

    num_col_blks = bpntr[iblk_row+1] - bpntr[iblk_row];

    /* loop over all the rows in this block row */

    for (irow = 0; irow < num_blk_rows; irow++) {
      realrow = irow + rpntr[iblk_row];

      /* loop over all the column blocks in this block row */

      for (jblk_col = 0; jblk_col < num_col_blks; jblk_col++) {

        /* find out which column block this is */

        icol_blk = bindx[iblk];
        indx_ptr = indx[iblk++];

        /* find out how many columns are in this block */

        num_blk_cols = cpntr[icol_blk+1] - cpntr[icol_blk];

        /* loop over all the columns in this block */

        for (jcol = 0; jcol < num_blk_cols; jcol++) {
          realcol = jcol + cpntr[icol_blk];

          /* store val(realrow, realcol) into the MSR matrix */

          if (realcol == realrow) {
            msr_val[realrow] = val[indx_ptr + jcol*num_blk_rows + irow];
          }
          else {
            /* We need to stick in the zeros so that bilu */ 
            /* works from within 'az_subdomain_driver.c'. */
            /* if (val[indx_ptr + jcol*num_blk_rows + irow] != 0.0) { */

              msr_val[nz_counter] = val[indx_ptr + jcol*num_blk_rows + irow];
              msr_bindx[nz_counter++]  = realcol;
            /* } */
          }
        }
      }

      iblk -= num_col_blks;
      msr_bindx[realrow+1] = nz_counter;
    }

    iblk += num_col_blks;
  }

} /* vb2msr */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void AZ_transform(int proc_config[], int *external[], int bindx[], double val[],
                  int update[], int *update_index[], int *extern_index[],
                  int *data_org[], int N_update, int indx[], int bnptr[],
                  int rnptr[], int *cnptr[], int mat_type)

/*******************************************************************************

  Convert a global distributed matrix to a parallel local distributed matrix.
  This includes the following steps:
      1) reorder matrix rows so that all the rows corresponding to internal
         points preceed all the rows corresponding to border points.
      2) replace global indicies by local indicies.
      3) make a list of the external unknowns and store them in external[].
      4) make a list of processors which update each external unknown and store
         this list in extern_proc where extern_proc[i] is the processor that
         updates external[i].
      5) make 2 arrays (update_index[], extern_index[]) which define the mapping
         between global and local indicies. In particular, global index
         update[i] corresponds to the locally numbered variable update_index[i]
         and global index external[i] corresponds to the locally numbered
         variable extern_index[i].
      6) Initialize all the quanities in data_org[] to their appropriate values
         (so that communication can properly occur).

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  proc_config:     Processor information corresponding to:
                      proc_config[AZ_node] = name of this processor
                      proc_config[AZ_N_procs] = # of processors used

  external:        On output, list of external blocks

  update:          On input, list of pts to be updated on this node

  update_index,    On output, ordering of update and external locally on this
  extern_index:    processor. For example  'update_index[i]' gives the index
                   location of the block which has the global index 'update[i]'.

  data_org:        On output, indicates how the data is set out on this node.
                   For example, data_org[] contains information on how many
                   unknowns are internal, external and border unknowns as well
                   as which points need to be communicated. See Aztec User's
                   guide for more details.

  N_update         Number of points to be updated on this node.

  val,bindx        On input, global distributed matrix (MSR or VBR) arrays
  indx, bnptr,     holding matrix values. On output, local renumbered matrix
  rnptr, cnptr:    (DMSR or DMVBR).

  mat_type:        Type of matrix (AZ_MSR_MATRIX or AZ_VBR_MATRIX).

*******************************************************************************/

{
  int        i, ii, j;
  static int mat_name = 1;

  int         N_extern;   /* Number of pts needed by this processor for
                             matrix-vector multiply but not updated by this
                             processor.  */
  int         N_internal, /* Number of pts which can be updated without
                             communication */
              N_border;   /* Number of pts to be updated requiring communication
                           */
  int        *extern_proc;
  int        *tcnptr = NULL;

  AZ__MPI_comm_space_ok();
#ifdef AZTEC_MPI
  if ( proc_config[AZ_Comm_Set] != AZ_Done_by_User) {
      AZ_printf_err("Error: Communicator not set. Use AZ_set_comm()\n");
      AZ_printf_err("       (e.g. AZ_set_comm(proc_config,MPI_COMM_WORLD)).\n");
      exit(1);
  }
#endif

  /*
   * Compute the external points and change the global indices to
   * local indices. That is,
   *   On input:                        On output:
   *      bindx[k] = update[j]      ==>   bindx[k] = j
   *      bindx[k] = external[j]    ==>   bindx[k] = j + N_update
   */

  AZ_find_local_indices(N_update, bindx, update, external, &N_extern, mat_type,
                        bnptr);

  /* compute the cnptr array for VBR matrices */

  if (mat_type == AZ_VBR_MATRIX) {
    if (!AZ_using_fortran) {
      *cnptr = (int *) AZ_allocate((N_update + N_extern + 1)*sizeof(int));
      if (*cnptr == NULL) {
         AZ_printf_err("Out of memory in AZ_transform\n");
         exit(1);
      }
    }

    tcnptr = *cnptr;
    for (i = 0 ; i < N_update + N_extern + 1; i++) tcnptr[i] = 0;

    for (i = 0; i < N_update; i++) tcnptr[i] = rnptr[i+1] - rnptr[i];

    for (i = 0; i < N_update; i++) {
      for (j = bnptr[i]; j < bnptr[i+1]; j++) {
        ii = bindx[j];

        if ((ii >= N_update) && ( tcnptr[ii] == 0)) {
          tcnptr[ii] = (indx[j+1]-indx[j]) / (rnptr[i+1]-rnptr[i]);
        }
      }
    }

    AZ_convert_values_to_ptrs(tcnptr, N_update + N_extern, 0);
  }

  /*
   * Read or compute (and sort) the processor numbers of the processors which
   * update the external points.
   */

  i                = AZ_using_fortran;
  AZ_using_fortran = AZ_FALSE;

  AZ_find_procs_for_externs(N_update, update, *external, N_extern, proc_config,
                            &extern_proc);
  AZ_using_fortran = i;

  /*
   * Determine a new ordering for the points:
   *    a) lowest numbers for internal points,
   *    b) next lowest numbers for border points
   *    c) highest nubers for the external points
   *       NOTE: external points updated by the same processor are consecutively
   *             ordered.
   */

  if (!AZ_using_fortran) {
    *update_index = (int *) AZ_allocate((N_update + 1)*sizeof(int));
    *extern_index = (int *) AZ_allocate((N_extern + 1)*sizeof(int));
  }

  if (*extern_index == NULL)  {
    (void) AZ_printf_err(
                   "Error: Not enough space in main() for extern_index[]\n");
    exit(1);
  }

  AZ_order_ele(*update_index, *extern_index, &N_internal, &N_border, N_update,
#ifdef MB_MODIF
               bnptr, bindx, extern_proc, N_extern, AZ_EXTERNS, mat_type);
#else
               bnptr, bindx, extern_proc, N_extern, AZ_ALL, mat_type);
#endif

  /*
   * Permute the matrix using the new ordering.  IMPORTANT: This routine assumes
   * that update_index[] contains 2 sequencies that are ordered but
   * intertwined. See AZ_reorder_matrix().
   */

  AZ_reorder_matrix(N_update, bindx, val, *update_index, *extern_index,
                    indx, rnptr, bnptr, N_extern, tcnptr, AZ_ALL,mat_type);

  /*
   * Initialize 'data_org' so that local information can be exchanged to update
   * the external points.
   */

  AZ_set_message_info(N_extern, *extern_index, N_update, *external, extern_proc,
                      update, *update_index, proc_config, tcnptr, data_org,
                      mat_type);

  (*data_org)[AZ_name]       = mat_name;
  (*data_org)[AZ_N_int_blk]  = N_internal;
  (*data_org)[AZ_N_bord_blk] = N_border;
  (*data_org)[AZ_N_ext_blk]  = N_extern;

  if (mat_type == AZ_VBR_MATRIX) {
    (*data_org)[AZ_N_internal] = rnptr[N_internal];
    (*data_org)[AZ_N_external] = tcnptr[N_update + N_extern] - rnptr[N_update];
    (*data_org)[AZ_N_border]   = rnptr[N_update] - rnptr[N_internal];
  }

  else {
    (*data_org)[AZ_N_internal] = N_internal;
    (*data_org)[AZ_N_external] = N_extern;
    (*data_org)[AZ_N_border]   = N_border;
  }

  mat_name++;
  AZ_free(extern_proc);

} /* AZ_transform */

/**************************************************************/
/**************************************************************/
/**************************************************************/

void AZ_msr2vbr_mem_efficient(int N, int **ibindx,double **ival, 
	int **icpntr, int **ibpntr, int **iindx, int *N_blk_rows,
	int name, char *label, int special) {
/*
 * Convert a local msr matrix to a local vbr matrix.
 *
 *   1) This version is fairly memory efficient.
 *   2) The size of the blocks are deduced by comparing
 *      the sparsity pattern of adjacent rows. That is,
 *      two adjacent rows with the same sparsity pattern
 *      are grouped together.
 *
 * Author : Ray Tuminaro, SNL, 9222 (3/5/98)
  =======

  Parameter list:
  ===============

  N:              On input, the number of rows in the msr matrix.

  ibindx, ival    On input, the arrays contain the MSR matrix that will
                  be convert to vbr format. 
                  On output, ibindx and ival correspond to the VBR matrix.
                  Note: ibindx is reallocated (i.e. shortened) in this
                  routine.

  icpntr, ibpntr  On output, these correspond to the resulting VBR matrix 
  iindx           after conversion from MSR.
                  Note: These arrays are allocated inside this routine.

  N_blk_rows      On output, the number of block rows in the vbr matrix.

  name            On input, indicates the name of the matrix to be used with 
                  manage_memory calls associated with VBR matrix.

  label           On input, indicates tail of string to be used for manage_memory
                  call associated with reallocating bindx.
*/

int    *bindx, *bpntr, *indx, *cpntr;
double *val, *dtemp;
int    *itemp;
int    i,j, ii, jj, next, offset;
int    row_i,row_i_plus1,row_i_plus2;
int    sparsity_pattern_match;
int    row_size, col_size, block_size, nz_per_row;
char   string[80];
unsigned long int diff;

   if (N == 0) return;

   bindx = *ibindx;
   val   = *ival;

   /* allocate dtemp so that it can hold */
   /* either N+1 integers or N+1 doubles */

   i = sizeof(double);
   if (sizeof(int) > sizeof(double)) i = sizeof(int);
   dtemp = (double *) AZ_allocate(((unsigned int) i*((unsigned int) (N+1))));
   if (dtemp == NULL) {
      AZ_printf_err("Note enough space in msr2vbr_mem_eff()\n");
      exit(1);
   }
   itemp = (int    *) dtemp;


   /* move the diagonal into the matrix (CSR) */


   for (i = 0 ; i < N ; i++ )  dtemp[i] = val[i];

   next = 0;
   for (i = 0 ; i < N ; i++ ) {
      val[next++] = dtemp[i];
      for (j = bindx[i] ; j < bindx[i+1]; j++ ) val[next++] = val[j];
   }

   for (i = 0 ; i <= N ; i++ ) itemp[i] = bindx[i];

   next = 0; 
   for (i = 0 ; i < N ; i++ ) {
      bindx[next++] = i;
      for (j = itemp[i] ; j < itemp[i+1] ; j++ ) bindx[next++] = bindx[j];
   }
   for (i = 0 ; i <= N ; i++ ) itemp[i] += i - N - 1;

   /* sort the columns within each row */

   for (i = 0 ; i < N ; i++ ) {
      row_i = itemp[i];  row_i_plus1 = itemp[i+1];
      AZ_sort( &(bindx[row_i]), row_i_plus1 - row_i , NULL, &(val[row_i]));
   }

   /* To get rid of itemp, we encode the last element */
   /* element in each row using a negative number.    */

   for (i = 1 ; i <= N ; i++ )  bindx[itemp[i]-1] = -1-bindx[itemp[i]-1];
   AZ_free(itemp);

   sprintf(string,"cpntr %s",label);
   cpntr = (int *) AZ_manage_memory((unsigned int) (N+1)*sizeof(int),
                                 AZ_ALLOC, name,string,&i);
   cpntr[N] = -100;

   /* Calculate cpntr. We basically decide that two   */
   /* adjacent rows should be in the same block if    */
   /* they have the same sparsity pattern.            */
   /* Note: After this calculation cpntr[i] indicates */
   /*       which block row i is within.              */

   *N_blk_rows = 0;
   block_size = 1;
   row_i = 0;
   row_i_plus1 = row_i;
   while( bindx[row_i_plus1] >= 0) row_i_plus1++;
   row_i_plus1++;
   i = 0;
   cpntr[i] = 0;
   while (i < N-1) {
      row_i_plus2 = row_i_plus1;
      while(bindx[row_i_plus2] >= 0) row_i_plus2++;
      row_i_plus2++;

      sparsity_pattern_match = 1;
      if ( row_i_plus2-row_i_plus1!= row_i_plus1-row_i) sparsity_pattern_match = 0;
      else {
         for (jj = 0 ; jj < row_i_plus1-row_i ; jj++ ) {
            if (bindx[row_i +jj] != bindx[row_i_plus1+jj]) {
               sparsity_pattern_match = 0;
               break;
            }
         }
      }

      if (sparsity_pattern_match == 1) block_size++;
      else {
         (*N_blk_rows)++;
         block_size = 1;
      }
      i++;
      cpntr[i] = *N_blk_rows;
      row_i = row_i_plus1;
      row_i_plus1   = row_i_plus2;
   }
   (*N_blk_rows)++;

   AZ_check_block_sizes(bindx, cpntr, N, N_blk_rows);


   /* convert bindx */

   sprintf(string,"bpntr %s",label);
   bpntr = (int *) AZ_manage_memory((unsigned int) (*N_blk_rows + 1)*sizeof(int),
                                     AZ_ALLOC, name,string,&i);

   *N_blk_rows = 0;
   next = 0;
   i = 0 ; 
   bpntr[(*N_blk_rows)++] = 0;
   row_i = 0;
   while ( i < N ) {
      row_i_plus1 = row_i;
      while( bindx[row_i_plus1] >= 0) row_i_plus1++;
      bindx[row_i_plus1] = -1 - bindx[row_i_plus1];  /* decode value */
      row_i_plus1++;
      nz_per_row = row_i_plus1 - row_i;
      j = row_i;
      while( j < row_i_plus1) {
         jj = cpntr[bindx[j]];
         bindx[next++] = jj;
         j++;
         while ((j < row_i_plus1) && (cpntr[bindx[j]] == jj)) j++;
      }
      bpntr[(*N_blk_rows)++] = next;
      i++;
      row_i += nz_per_row;
      while ( cpntr[i] == cpntr[i-1] ) { i++; row_i += nz_per_row;}
   }
   (*N_blk_rows)--;
   if (special != 1) {
      sprintf(string,"bindx %s",label);
      bindx = (int    *) AZ_manage_memory(((unsigned int) next)*sizeof(int),
                                       AZ_REALLOC, name,string,&i);
   }
   else {
      sprintf(string,"val %s",label);
      diff = (unsigned long int) bindx - (unsigned long int) val;
      val = (double *) AZ_manage_memory((unsigned int) (diff + 
                                        (unsigned long int) next*sizeof(int)), 
                                        AZ_SPEC_REALLOC,name,string,&i);
      bindx = (int *) ( (unsigned long int) val + diff);
   }
   sprintf(string,"indx %s",label);
   indx  = (int *) AZ_manage_memory(((unsigned int)(next+1))*sizeof(int),
                                    AZ_ALLOC, name,string,&i);

   /* compress cpntr */

   next = 0;
   j = 0;
   for (i = 1; i < N ; i++ ) {
      j++;
      if ( cpntr[i-1] != cpntr[i] ) { cpntr[next++] = j; j = 0; }
   }
   cpntr[next] = j+1;
 
   sprintf(string,"cpntr %s",label);
   cpntr = (int *) AZ_manage_memory(((unsigned int) (next+2))*sizeof(int),
                                AZ_REALLOC, name,string,&i);


   /* compute the size of dtemp */

   jj = -1;
   for (i = 0 ; i < *N_blk_rows; i++ ) {
      row_size = cpntr[i];
      nz_per_row = 0;
      for (j = bpntr[i]; j < bpntr[i+1] ; j++ ) nz_per_row += cpntr[bindx[j]]; 
      ii = row_size*nz_per_row;
      if (ii > jj) jj = ii;
   }
   dtemp = (double *) AZ_allocate(sizeof(double)*((unsigned int) (jj+1)));
   if (dtemp == NULL) {
      AZ_printf_err("Note enough space in msr2vbr_mem_eff()\n");
      exit(1);
   }


   /* compute indx and shift val */

   indx[0] = 0;

   for (i = 0 ; i < *N_blk_rows; i++ ) {
      offset = indx[bpntr[i]];
      nz_per_row = 0;
      for (j = bpntr[i]; j < bpntr[i+1] ; j++ ) 
         nz_per_row += cpntr[bindx[j]]; 
      row_size = cpntr[i];
      next = 0;
      for (j = bpntr[i]; j < bpntr[i+1] ; j++ ) {
         col_size = cpntr[bindx[j]];
         indx[j+1] = indx[j] + row_size*col_size;
         for (ii = 0; ii < row_size ; ii++ ) {
            for (jj = 0 ; jj < col_size ; jj++ ) {
               dtemp[next++] = val[offset+ jj + ii*nz_per_row];
            }
         }
         offset += col_size;
      }
      next = 0;
      for (j = bpntr[i]; j < bpntr[i+1] ; j++ ) {
         offset = indx[j];
         col_size = cpntr[bindx[j]];
         for (ii = 0; ii < row_size ; ii++ ) {
            for (jj = 0 ; jj < col_size ; jj++ ) {
               val[offset + ii + jj*row_size] =  dtemp[next++];
            }
         }
      }

   }
   AZ_free(dtemp);

   /* convert cpntr[i] from the size of   */
   /* block i to the start row of block i */

   j = 0;
   for (i = 0 ; i < *N_blk_rows ; i++ ) {
      jj = cpntr[i];
      cpntr[i] = j;
      j += jj;
   }
   cpntr[*N_blk_rows] = j;

   *ibindx = bindx;
   *ival   = val;
   *icpntr = cpntr;
   *ibpntr = bpntr;
   *iindx  = indx;

}
void AZ_find_global_ordering(int proc_config[], AZ_MATRIX *Amat, 
                             int **global_bindx, int **update)

{
  int i, ii;
  int N_rows, N_cols, N_blk_rows = 0, N_external, N_ext_blks = 0;
  int n_nonzeros = 0, n_blk_nonzeros = 0;
  int *data_org;
  int *bindx, *rpntr, *externals = NULL;
  double *rownums;
  int is_VBR = 0, max_per_proc, offset;

  data_org = Amat->data_org;

  bindx = Amat->bindx;
  rpntr = Amat->rpntr;
  N_rows = data_org[AZ_N_internal] + data_org[AZ_N_border];
  N_external = data_org[AZ_N_external];
  N_cols = N_rows + N_external;

  if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
    {
      is_VBR = 1;
      N_blk_rows = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
      N_ext_blks = data_org[AZ_N_ext_blk];
      n_nonzeros = Amat->indx[Amat->bpntr[N_blk_rows]];
      n_blk_nonzeros = Amat->bpntr[N_blk_rows];
    }
  else if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX)
    {
      is_VBR = 0;
      N_blk_rows = N_rows;
      N_ext_blks = data_org[AZ_N_external];
      n_nonzeros = Amat->bindx[N_rows]-1;
      n_blk_nonzeros = n_nonzeros;
    }
  else
    {
      AZ_perror("Unsupported matrix type in AZ_find_global_ordering.");
    }

  
  *global_bindx = (int *) AZ_allocate((n_blk_nonzeros+1)*sizeof(int));
  if ((*global_bindx) == NULL) 
    AZ_perror("Error: Not enough space in AZ_find_global_ordering");
  
  if (N_ext_blks>0)
  externals      = (int    *) AZ_allocate(N_ext_blks*sizeof(int));
  rownums        = (double *) AZ_allocate(N_cols    *sizeof(double));
  if (rownums == NULL) 
    AZ_perror("Error: Not enough space in AZ_find_global_ordering");

  /* 
     Tranform the local matrices to a global matrix
     by using the exchange boundary routine to pass indices around.
   */

  max_per_proc = AZ_gmax_int(N_blk_rows,proc_config);
  max_per_proc++;
  offset       = max_per_proc*proc_config[AZ_node];
  
  if (is_VBR)
    {
      for (i = 0 ; i < N_cols; i++) rownums[i] = -1.0;
      for (i = 0 ; i < N_blk_rows; i++ ) 
	rownums[rpntr[i]] = (double) (offset + i);
    }
  else
    for (i = 0 ; i < N_blk_rows; i++ ) rownums[i] = (double) (offset + i);

  AZ_exchange_bdry(rownums, data_org, proc_config);

  if (is_VBR)
    {
      ii = 0;
      for (i = 0 ; i < N_external ; i++ )
	if (rownums[i + N_rows] >= 0.0)
	  externals[ii++] = (int) rownums[i + N_rows];
    }
  else
      for (i = 0 ; i < N_external ; i++ )
	externals[i] = (int) rownums[i + N_rows];

  /* change matrix columns to reflect global numbers */

  if (is_VBR)
    {
      for ( i = 0; i < n_blk_nonzeros; i++ ) 
	{
	  if ( bindx[i] < N_blk_rows) (*global_bindx)[i] = bindx[i] + offset;
	  else (*global_bindx)[i] = externals[bindx[i] - N_blk_rows];
	}
    }
  else
    {
      for ( i = 0; i < N_rows+1; i++ ) (*global_bindx)[i] = bindx[i];
      for ( i = N_rows+1; i < n_nonzeros+1; i++ ) 
	{
	  if ( bindx[i] < N_rows) (*global_bindx)[i] = bindx[i] + offset;
	  else (*global_bindx)[i] = externals[bindx[i] - N_rows];
	}
    }
  if (N_ext_blks>0)
  AZ_free((void *) externals);
  AZ_free((void *) rownums);

  /* Build update vector */

  *update = ( int *) AZ_allocate(N_rows*sizeof(int));

  for (i=0; i< N_blk_rows; i++) (*update)[i] = offset +i;

  /* end AZ_find_global_ordering */
}

void AZ_revert_to_global(int proc_config[], 
			 AZ_MATRIX *Amat, int **global_bindx, int **update) {
  int i, ii;
  int N_rows, N_cols, N_blk_rows = 0, N_external, N_ext_blks = 0;
  int n_nonzeros = 0, n_blk_nonzeros = 0;
  int *data_org;
  int *bindx, *rpntr, *externals = NULL;
  double *rownums;
  int is_VBR = 0, max_per_proc, offset;
  int * tmp_update;

  data_org = Amat->data_org;

  bindx = Amat->bindx;
  rpntr = Amat->rpntr;
  N_rows = data_org[AZ_N_internal] + data_org[AZ_N_border];
  N_external = data_org[AZ_N_external];
  N_cols = N_rows + N_external;

  if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
    {
      is_VBR = 1;
      N_blk_rows = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
      N_ext_blks = data_org[AZ_N_ext_blk];
      n_nonzeros = Amat->indx[Amat->bpntr[N_blk_rows]];
      n_blk_nonzeros = Amat->bpntr[N_blk_rows];
    }
  else if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX)
    {
      is_VBR = 0;
      N_blk_rows = N_rows;
      N_ext_blks = data_org[AZ_N_external];
      n_nonzeros = Amat->bindx[N_rows]-1;
      n_blk_nonzeros = n_nonzeros;
    }
  else
    {
      AZ_perror("Unsupported matrix type in AZ_find_global_ordering.");
    }

  
  *global_bindx = (int *) AZ_allocate((n_blk_nonzeros+1)*sizeof(int));
  if ((*global_bindx) == NULL) 
    AZ_perror("Error: Not enough space in AZ_find_global_ordering");
  
  if (N_ext_blks>0)
  externals      = (int    *) AZ_allocate(N_ext_blks*sizeof(int));
  rownums        = (double *) AZ_allocate(N_cols    *sizeof(double));
  if (rownums == NULL) 
    AZ_perror("Error: Not enough space in AZ_find_global_ordering");

  /* 
     Tranform the local matrices to a global matrix
     by using the exchange boundary routine to pass indices around.
   */

    tmp_update = ( int *) AZ_allocate(N_rows*sizeof(int));
  if (Amat->update==NULL) {
  
    /* Build update vector */
    
    max_per_proc = AZ_gmax_int(N_blk_rows,proc_config);
    max_per_proc++;
    offset       = max_per_proc*proc_config[AZ_node];
    
    for (i=0; i< N_blk_rows; i++) tmp_update[i] = offset +i;
    
  }
  else
    for (i=0; i< N_blk_rows; i++) tmp_update[i] = Amat->update[i];

  /* Now use boundary exchange to get global ids */
  if (is_VBR)
    {
      for (i = 0 ; i < N_cols; i++) rownums[i] = -1.0;
      for (i = 0 ; i < N_blk_rows; i++ ) 
	rownums[rpntr[i]] = (double) tmp_update[i];
    }
  else
    for (i = 0 ; i < N_blk_rows; i++ ) rownums[i] = (double)  tmp_update[i];

  AZ_exchange_bdry(rownums, data_org, proc_config);

  if (is_VBR)
    {
      ii = 0;
      for (i = 0 ; i < N_external ; i++ )
	if (rownums[i + N_rows] >= 0.0)
	  externals[ii++] = (int) rownums[i + N_rows];
    }
  else
      for (i = 0 ; i < N_external ; i++ )
	externals[i] = (int) rownums[i + N_rows];

  /* change matrix columns to reflect global numbers */

  if (is_VBR)
    {
      for ( i = 0; i < n_blk_nonzeros; i++ ) 
	{
	  if ( bindx[i] < N_blk_rows) (*global_bindx)[i] = tmp_update[bindx[i]];
	  else (*global_bindx)[i] = externals[bindx[i] - N_blk_rows];
	}
    }
  else
    {
      for ( i = 0; i < N_rows+1; i++ ) (*global_bindx)[i] = bindx[i];
      for ( i = bindx[0]; i < n_nonzeros+1; i++ ) 
	{
	  if ( bindx[i] < N_rows) (*global_bindx)[i] = tmp_update[bindx[i]];
	  else (*global_bindx)[i] = externals[bindx[i] - N_rows];
	}
    }
  AZ_free((void *) rownums);
  if (N_ext_blks>0) AZ_free((void *) externals);

  *update = tmp_update;


  /* end AZ_revert_to_global
   */
}
