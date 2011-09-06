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

/* System Include files */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "az_aztec.h"

#ifndef TRUE
#define TRUE  1
#define FALSE 0
#endif

extern int az_iterate_id;

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void AZ_exchange_bdry(double x[], int data_org[], int proc_config[])

/*******************************************************************************

  Routine to locally exchange the components of the vector "x". This routine
  gathers the necessary components of the vector and then sends the required
  "num_neighbors" messages. The messages which are received are placed
  contiguously in the external_nodes part of x.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:               Vector of unknowns defined on the current processor.
                   Indirect addressing will be used to gather those unknowns
                   that other processors need.

  data_org:        Array containing information on the distribution of the
                   matrix to this processor as well as communication parameters
                   (see Aztec User's Guide).

*******************************************************************************/

{

  /* local variables */

  char *message_send_add[AZ_MAX_NEIGHBORS];
  /* message_send_add[i] points to the beginning of the list of values to be
     sent to the ith neighbor (i.e. data_org[AZ_neighbors+i] or sometimes
     locally defined as proc_num_neighbor[i]). That is, *(message_send_add[i]+j)
     is the jth value to be sent to the ith neighbor. */
  unsigned int  message_send_length[AZ_MAX_NEIGHBORS];
  /* message_send_length[i] is the number of bytes to be sent to the ith
     neighbor (i.e. data_org[AZ_neighbors+i] or sometimes locally defined as
     proc_num_neighbor[i]). */
  char *message_recv_add[AZ_MAX_NEIGHBORS];
  /* message_recv_add[i] points to the beginning of the list of locations which
     are to receive values sent by the ith neighbor (i.e.
     data_org[AZ_neighbors+i] or sometimes locally defined as
     proc_num_neighbor[i]). That is, *(message_recv_add[i] + j) is the location
     where the jth value sent from the ith neighbor will be stored. */
  unsigned int  message_recv_length[AZ_MAX_NEIGHBORS];
  /* message_recv_length[i] is the number of bytes to be sent to the ith
     neighbor (i.e. data_org[AZ_neighbors+i] or sometimes locally defined as
     proc_num_neighbor[i]). */

  int              n;
  double          *ptr_send_list, *ptr_recv_list;
  register double *ptrd;
  register int    *ptr_int, i;
  int              size, num_send, num_recv;
  int              type;
  int              num_neighbors, *proc_num_neighbor, total_num_send_unknowns;
  int             *num_unknowns_send_neighbor, *list_send_unknowns;
  int              external_index, *num_unknowns_recv_neighbor;

  /**************************** execution begins ******************************/

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

#ifdef AZTEC_MPI 
  if ( proc_config[AZ_Comm_Set] != AZ_Done_by_User) {
      AZ_printf_err("Error: Communicator not set. Use AZ_set_comm()\n");
      AZ_printf_err("       (e.g. AZ_set_comm(proc_config,MPI_COMM_WORLD)).\n");
      exit(1);
  }
#endif

  num_neighbors              = data_org[AZ_N_neigh];
  if (num_neighbors == 0) return;

  proc_num_neighbor          = &data_org[AZ_neighbors];
  total_num_send_unknowns    = data_org[AZ_total_send];
  num_unknowns_send_neighbor = &data_org[AZ_send_length];
  list_send_unknowns         = &data_org[AZ_send_list];
  external_index             = data_org[AZ_N_internal] + data_org[AZ_N_border];
  num_unknowns_recv_neighbor = &data_org[AZ_rec_length];


  /* single processor case */

  if (num_neighbors == 0) return;

  /* Set up send messages: Gather send unknowns from "x" vector */

  ptrd = (double *) AZ_manage_memory(data_org[AZ_total_send]*sizeof(double),
                                     AZ_ALLOC, AZ_SYS+az_iterate_id, "comm buff", &n);
  ptr_send_list = ptrd;

  ptr_int = list_send_unknowns;
  for (i = total_num_send_unknowns; i > 0; i--) {
    *ptrd++ = x[*ptr_int++];
  }

  /* Define arrays for message passing */

  ptr_recv_list = &x[external_index];

  size = sizeof(double);
  for (n = 0; n < num_neighbors; n++) {
    num_send               = num_unknowns_send_neighbor[n];
    num_recv               = num_unknowns_recv_neighbor[n];
    message_send_add[n]    = (char *) ptr_send_list;
    message_recv_add[n]    = (char *) ptr_recv_list;
    message_send_length[n] = size * num_send;
    message_recv_length[n] = size * num_recv;
    ptr_send_list         += num_send;
    ptr_recv_list         += num_recv;
  }

  AZ_exchange_local_info(num_neighbors, proc_num_neighbor, message_send_add,
                         message_send_length, message_recv_add,
                         message_recv_length, type, proc_config);

} /* AZ_exchange_bdry */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_exchange_local_info(int num_neighbors, int proc_num_neighbor[],
                            char *message_send_add[], 
                            unsigned int message_send_length[],
                            char *message_recv_add[], 
                            unsigned int message_recv_length[],
                            int type, int proc_config[])

/*******************************************************************************

  Routine to communicate between a small number of processors by first writing
  all messages and then reading all messages. This generic routine can be called
  with the standard addresses, message lengths and message types for a wide
  variety of uses.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  num_neighbors:             Total number of neighboring processors.

  proc_num_neighbor:         Array containing the processor id for each of the
                             current processor's neighbors.

  message_send_add:          message_send_add[i] points to the beginning of
                             the list of values to be sent to the ith neighbor
                             (i.e. data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]). That is,
                             *(message_send_add[i] + j) is the jth value to be
                             sent to the ith neighbor.

  message_send_length:       message_send_length[i] is the number of bytes to
                             be sent to the ith neighbor (i.e.
                             data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]).

  message_recv_add:          message_recv_add[i] points to the beginning of the
                             list of locations which are to receive values sent
                             by the ith neighbor (i.e.  data_org[AZ_neighbors+i]
                             or sometimes locally defined as
                             proc_num_neighbor[i]). That is,
                             *(message_recv_add[i] + j) is the location where
                             the jth value sent from the ith neighbor will be
                             stored.

  message_recv_length:       message_recv_length[i] is the number of bytes to
                             be sent to the ith neighbor (i.e.
                             data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]).

  type:                      message type used when exchanging information.

*******************************************************************************/

{

  /* local declarations */

  register int n;
  int          rtype, st;

  MPI_AZRequest request[AZ_MAX_NEIGHBORS];  /* Message handle */

  /*********************** first executable statment *****************/

  /* post receives for all messages */

  for (n = 0; n < num_neighbors; n++) {
    rtype = type;
    (void) mdwrap_iread((void *) *(message_recv_add+n),
                         *(message_recv_length+n), proc_num_neighbor+n, &rtype,
                         request+n);
  }

  /* write out all messages */

  for (n = 0; n < num_neighbors; n++) {
    (void) mdwrap_write((void *) *(message_send_add+n),
                         *(message_send_length+n), *(proc_num_neighbor+n),
                         type, &st);
  }

  /* wait for all messages */

  for (n = 0; n < num_neighbors; n++) {
    rtype = type;
    (void) mdwrap_wait((void *) *(message_recv_add+n),
                        *(message_recv_length+n), proc_num_neighbor+n, &rtype,
                        &st, request+n);
  }

} /* AZ_exchange_local_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_gather_mesg_info(double x[] ,int data_org[], char *message_recv_add[],
                         char *message_send_add[], int message_recv_length[],
                         int message_send_length[])

/*******************************************************************************

  Routine to locally exchange the components of the vector "x". This routine
  gathers the necessary components of the vector and then sends the required
  "num_neighbors" messages. The messages which are received are placed
  contiguously in the external_nodes part of x.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  x:                         Vector of unknowns defined on the current
                             processor. Indirect addressing will be used to
                             gather those unknowns that other processors need.

  data_org:                  Array containing information on the distribution of
                             the matrix to this processor as well as
                             communication parameters (see Aztec User's Guide).

  message_send_add:          message_send_add[i] points to the beginning of
                             the list of values to be sent to the ith neighbor
                             (i.e. data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]). That is,
                             *(message_send_add[i] + j) is the jth value to be
                             sent to the ith neighbor.

  message_send_length:       message_send_length[i] is the number of bytes to
                             be sent to the ith neighbor (i.e.
                             data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]).

  message_recv_add:          message_recv_add[i] points to the beginning of the
                             list of locations which are to receive values sent
                             by the ith neighbor (i.e.  data_org[AZ_neighbors+i]
                             or sometimes locally defined as
                             proc_num_neighbor[i]). That is,
                             *(message_recv_add[i] + j) is the location where
                             the jth value sent from the ith neighbor will be
                             stored.

  message_recv_length:       message_recv_length[i] is the number of bytes to
                             be sent to the ith neighbor (i.e.
                             data_org[AZ_neighbors+i] or sometimes locally
                             defined as proc_num_neighbor[i]).

  num_neighbors:             Total number of neighboring processors.


  Important external definitions:
  ===============================

  Num_Neighbors:             Total number of neighboring processors.
                             (type - int value)

  Proc_Num_Neighbor[]:       Array containing the processor id for each of the
                             current processor's neighbors.
                             (type - int vector with fixed length
                             AZ_MAX_NEIGHBORS)

  Total_Num_Send_Unknowns:   Total number of unknowns which are sent to
                             neighboring processors (size of mesg buff).
                             (type - int value)

  Num_Unknowns_Send_Neighbor[]:
                             Vector containing the number of unknowns to be sent
                             to each processor.
                             (type - int vector with fixed length
                             AZ_MAX_NEIGHBORS)

  List_Send_Unknowns[]:      Vector of local node numbers for the unknowns to be
                             sent to each neighbor.
                             (type - pointer to variable length int vector)

  Num_Unknowns_Recv_Neighbor[]:
                             Array of number of unknowns to be received from
                             each of the neighboring processors.
                             (type - int vector with fixed length
                             AZ_MAX_NEIGHBORS)

*******************************************************************************/

{

  /* local variables */

  double          *ptr_send_list, *ptr_recv_list;
  register double *ptrd;
  register int    *ptr_int, i;
  int              size, num_send, num_recv;
  int              n;

  int              Num_Neighbors;
  int             *Num_Unknowns_Recv_Neighbor, *Num_Unknowns_Send_Neighbor;

  /*********************** first executable statement *****************/

  Num_Neighbors              =  data_org[AZ_N_neigh];
  Num_Unknowns_Recv_Neighbor = &data_org[AZ_rec_length];
  Num_Unknowns_Send_Neighbor = &data_org[AZ_send_length];

  /* single processor case */

  if (Num_Neighbors == 0) return;

  /* Set up send messages:  Gather send unknowns from "x" vector */

  ptrd = (double *) AZ_manage_memory(data_org[AZ_total_send]*sizeof(double),
                                     AZ_ALLOC, AZ_SYS+az_iterate_id, "ptrd", &n);
  ptr_send_list = ptrd;

  ptr_int = &data_org[AZ_send_list];
  for (i = data_org[AZ_total_send]; i > 0; i--) {
    *ptrd++ = x[*ptr_int++];
  }

  /* Define arrays for message passing */

  ptr_recv_list = &x[ data_org[AZ_N_internal] + data_org[AZ_N_border] ];

  size = sizeof(double);
  for (n = 0; n < Num_Neighbors; n++) {
    num_send               = Num_Unknowns_Send_Neighbor[n];
    num_recv               = Num_Unknowns_Recv_Neighbor[n];
    message_send_add[n]    = (char *) ptr_send_list;
    message_recv_add[n]    = (char *) ptr_recv_list;
    message_send_length[n] = size * num_send;
    message_recv_length[n] = size * num_recv;
    ptr_send_list         += num_send;
    ptr_recv_list         += num_recv;
  }

} /* AZ_gather_mesg_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_gsum_int(int val, int proc_config[])

/*******************************************************************************

  Global integer sum.

  Author:
  =======

  Return code:     int, result of global sum.
  ============

  Parameter list:
  ===============

  val:             Individual processor value to be summed.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int   type;             /* type of next message */
  int   partner;          /* processor I exchange with */
  int   mask;             /* bit pattern identifying partner */
  int   hbit;             /* largest nonzero bit in nprocs */
  int   nprocs_small;     /* largest power of 2 <= nprocs */
  int   val2;             /* arriving value to add */
  int   cflag;            /* dummy argument for compatability */
  int   node, nprocs;
  char *yo = "AZ_gsum_int: ";

  MPI_AZRequest request;  /* Message handle */

  /*********************** first executable statment *****************/

  node            = proc_config[AZ_node];
  nprocs          = proc_config[AZ_N_procs];
  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small * 2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (mdwrap_iread((void *) &val2, sizeof(int), &partner, &type, &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (mdwrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (mdwrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* sum values */

    val += val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;
      if (mdwrap_iread((void *) &val2, sizeof(int), &partner, &type,
                        &request)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                       &request) != sizeof(int)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      val += val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (mdwrap_iread((void *) &val, sizeof(int), &partner, &type, &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (mdwrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (mdwrap_wait((void *) &val, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gsum_int */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

#ifndef PUMA_GSUMD  /* Only compile the generic gsum_double function in */
                    /* if the "optimized" gsum_double for PUMA is not   */
                    /* desired.                                         */

double AZ_gsum_double(double val, int proc_config[])

/*******************************************************************************

  Global double sum.

  Author:
  =======

  Return code:     double, result of global sum.
  ============

  Parameter list:
  ===============

  val:             Individual processor value to be summed.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int    type;          /* type of next message */
  int    partner;       /* processor I exchange with */
  int    mask;          /* bit pattern identifying partner */
  int    hbit;          /* largest nonzero bit in nprocs */
  int    nprocs_small;  /* largest power of 2 <= nprocs */
  double val2;          /* arriving value to add */
  int    cflag;         /* dummy argument for compatability */
  int    node, nprocs;
  char  *yo = "AZ_gsum_double: ";

  MPI_AZRequest request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (mdwrap_iread((void *) &val2, sizeof(double), &partner, &type,
                      &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (mdwrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (mdwrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* sum values */

    val += val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (mdwrap_iread((void *) &val2, sizeof(double), &partner, &type,
                        &request)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                       &request) != sizeof(double)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      val += val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (mdwrap_iread((void *) &val, sizeof(double), &partner, &type,
                      &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (mdwrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (mdwrap_wait((void *) &val, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gsum_double */

#endif    /* ifndef PUMA_GSUMD */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gmax_double(double val, int proc_config[])

/*******************************************************************************

  Global max of type double.

  Author:
  =======

  Return code:     double, maximum value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int    type;          /* type of next message */
  int    partner;       /* processor I exchange with */
  int    mask;          /* bit pattern identifying partner */
  int    hbit;          /* largest nonzero bit in nprocs */
  int    nprocs_small;  /* largest power of 2 <= nprocs */
  double val2;          /* arriving value to add */
  int    cflag;         /* dummy argument for compatability */
  int    node, nprocs;
  char  *yo = "AZ_gmax_double: ";

  MPI_AZRequest request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (mdwrap_iread((void *) &val2, sizeof(double), &partner, &type,
                      &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (mdwrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (mdwrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* get max value */

    if (val2 > val) val = val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (mdwrap_iread((void *) &val2, sizeof(double), &partner, &type,
                        &request)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                       &request) != sizeof(double)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (val2 > val)
        val = val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (mdwrap_iread((void *) &val, sizeof(double), &partner, &type,
                      &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (mdwrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (mdwrap_wait((void *) &val, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gmax_double */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gmin_double(double val, int proc_config[])

/*******************************************************************************

  Global min of type double.

  Author:
  =======

  Return code:     double, minimum value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int    type;               /* type of next message */
  int    partner;            /* processor I exchange with */
  int    mask;               /* bit pattern identifying partner */
  int    hbit;               /* largest nonzero bit in nprocs */
  int    nprocs_small;       /* largest power of 2 <= nprocs */
  double val2;               /* arriving value to add */
  int    cflag;              /* dummy argument for compatability */
  int    node, nprocs;
  char  *yo = "AZ_gmin_double: ";

  MPI_AZRequest request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (mdwrap_iread((void *) &val2, sizeof(double), &partner, &type,
                      &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (mdwrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (mdwrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* get max value */

    if (val2 < val) val = val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (mdwrap_iread((void *) &val2, sizeof(double), &partner, &type,
                        &request)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_wait((void *) &val2, sizeof(double), &partner, &type, &cflag,
                       &request) != sizeof(double)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (val2 < val)
        val = val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (mdwrap_iread((void *) &val, sizeof(double), &partner, &type,
                      &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (mdwrap_write((void *) &val, sizeof(double), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (mdwrap_wait((void *) &val, sizeof(double), &partner, &type, &cflag,
                     &request) != sizeof(double)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gmin_double */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_gmax_int(int val, int proc_config[])

/*******************************************************************************

  Global max of type int.

  Author:
  =======

  Return code:     int, maximum value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int   type;                     /* type of next message */
  int   partner;                  /* processor I exchange with */
  int   mask;                     /* bit pattern identifying partner */
  int   hbit;                     /* largest nonzero bit in nprocs */
  int   nprocs_small;             /* largest power of 2 <= nprocs */
  int   val2;                     /* arriving value to add */
  int   cflag;                    /* dummy argument for compatability */
  int   node, nprocs;
  char *yo = "AZ_gmax_int: ";

  MPI_AZRequest request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (mdwrap_iread((void *) &val2, sizeof(int), &partner, &type, &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (mdwrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (mdwrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* get max value */

    if (val2 > val) val += val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;
      if (mdwrap_iread((void *) &val2, sizeof(int), &partner, &type,
                        &request)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                       &request) != sizeof(int)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (val2 > val) val = val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (mdwrap_iread((void *) &val, sizeof(int), &partner, &type, &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (mdwrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (mdwrap_wait((void *) &val, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gmax_int */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_gavg_double(double val, int proc_config[])

/*******************************************************************************

  Global average of type double.

  Author:
  =======

  Return code:     double, average value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  return (AZ_gsum_double(val, proc_config) / proc_config[AZ_N_procs]);

} /* AZ_gavg_double */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int AZ_gmin_int(int val, int proc_config[])

/*******************************************************************************

  Global min of type int.

  Author:
  =======

  Return code:     int, minimum value across all processors.
  ============

  Parameter list:
  ===============

  val:             Individual processor value.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int   type;             /* type of next message */
  int   partner;          /* processor I exchange with */
  int   mask;             /* bit pattern identifying partner */
  int   hbit;             /* largest nonzero bit in nprocs */
  int   nprocs_small;     /* largest power of 2 <= nprocs */
  int   val2;             /* arriving value to add */
  int   cflag;            /* dummy argument for compatability */
  int   node, nprocs;
  char *yo = "AZ_gmin_int: ";

  MPI_AZRequest request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (mdwrap_iread((void *) &val2, sizeof(int), &partner, &type, &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (mdwrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (mdwrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* get min value */

    if (val2 < val) val = val2;
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;
      if (mdwrap_iread((void *) &val2, sizeof(int), &partner, &type,
                        &request)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_wait((void *) &val2, sizeof(int), &partner, &type, &cflag,
                       &request) != sizeof(int)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (val2 < val) val = val2;
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (mdwrap_iread((void *) &val, sizeof(int), &partner, &type, &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (mdwrap_write((void *) &val, sizeof(int), partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (mdwrap_wait((void *) &val, sizeof(int), &partner, &type, &cflag,
                     &request) != sizeof(int)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  return val;

} /* AZ_gmin_int */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_gdot_vec(int N, double dots[], double dots2[], int proc_config[])

/*******************************************************************************

  Author:
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  N:               Length of vectors 'dots' and 'dots2'.

  dots:            On output, AZ_gdot_vec() produces N global sums. Specifically
                   dots[i] (on all processors) = dots[i] (on proc 0) +
                                                 dots[i] (on proc 1) +
                                                      .
                                                      .
                                                      .
                                                 dots[i] (on proc proc_config[
                                                          AZ_N_procs])
                   for 0 <= i < N .

  dots2:           double precision work vector of length N.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  /* local variables */

  int   type;             /* type of next message */
  unsigned int   msg_size;/* length of messages */
  int   partner;          /* processor I exchange with */
  int   mask;             /* bit pattern identifying partner */
  int   hbit;             /* largest nonzero bit in nprocs */
  int   nprocs_small;     /* largest power of 2 <= nprocs */
  int   cflag;            /* dummy argument for compatability */
  int   i;                /* loop counter */
  int   node, nprocs;
  char *yo = "AZ_gdot_vec: ";

  MPI_AZRequest request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small * 2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  msg_size = N * sizeof(double);
  partner  = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (mdwrap_iread((void *) dots2, msg_size, &partner, &type, &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (mdwrap_write((void *) dots, msg_size, partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (mdwrap_wait((void *) dots2, msg_size, &partner, &type, &cflag,
                     &request) < msg_size) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* sum values */

    for (i = 0; i < N; i++) dots[i] += dots2[i];
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (mdwrap_iread((void *) dots2, msg_size, &partner, &type, &request)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_write((void *) dots, msg_size, partner, type, &cflag)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_wait((void *) dots2, msg_size, &partner, &type, &cflag,
                       &request) < msg_size) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      for (i = 0; i < N; i++) dots[i] += dots2[i];
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (mdwrap_iread((void *) dots, msg_size, &partner, &type, &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (mdwrap_write((void *) dots, msg_size, partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (mdwrap_wait((void *) dots, msg_size, &partner, &type, &cflag,
                     &request) < msg_size) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

} /* AZ_gdot_vec */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_print_sync_start(int proc, int do_print_line, int proc_config[])

/*******************************************************************************

  Routine to allow IO between print_sync_start and print_sync_end to be printed
  by each processor entirely before the next processor begins its IO.  The
  printing sequence is from proc = 0 to the last processor,
  number_of_procs = nprocs - 1.

  NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  proc:            Current processor number.

  do_print_line:   Boolean variable.  If true, a line of # is printed to
                   indicate the start of a print_sync I/O block.

*******************************************************************************/

{

  /* local variables */

  int flag = 1, from, st, type;
  MPI_AZRequest request;

  /**************************** execution begins ******************************/

  type = proc_config[AZ_MPI_Tag];

  if (proc_config[AZ_node] != 0) {
    from = proc - 1;

    mdwrap_iread((void *) &flag, sizeof(int), &from, &type, &request);
    mdwrap_wait((void *) &flag, sizeof(int), &from, &type, &st, &request);
  }

  else {
    if (do_print_line) {
      (void) AZ_printf_out("\n");
      for (flag = 0; flag < 37; flag++) (void) AZ_printf_out("#");
      (void) AZ_printf_out(" PRINT_SYNC_START ");
      for (flag = 0; flag < 25; flag++) (void) AZ_printf_out("#");
      (void) AZ_printf_out("\n");
    }
  }

} /* AZ_print_sync_start */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_print_sync_end(int proc_config[], int do_print_line)

/*******************************************************************************

  Routine to allow IO between print_sync_start and print_sync_end to be printed
  by each processor entirely before the next processor begins its IO. The
  printing sequence is from proc = 0 to the last processor,
  number_of_procs = nprocs - 1.

  NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  proc:            Current processor number.

  nprocs:          Number of processors in the current machine configuration.

  do_print_line:   Boolean variable.  If true, a line of # is printed to
                   indicate the start of a print_sync I/O block.

*******************************************************************************/

{

  /* local variables */

  int st, flag = 1, from, type, to, proc, nprocs;
  MPI_AZRequest request, request2;

  /**************************** execution begins ******************************/

  proc = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];
  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  if (proc < nprocs -1) to = proc + 1;
  else {
    to = 0;
    if (do_print_line) {
      (void) AZ_printf_out("\n");
      for (flag = 0; flag < 37; flag++) (void) AZ_printf_out("#");
      (void) AZ_printf_out(" PRINT_SYNC_END__ ");
      for (flag = 0; flag < 25; flag++) (void) AZ_printf_out("#");
      (void) AZ_printf_out("\n\n");
    }
  }

  mdwrap_iwrite((void *) &flag, sizeof(int), to, type, &st, &request);

  if (proc == 0) {
    from = nprocs -1;
    mdwrap_iread((void *) &flag, sizeof(int), &from, &type, &request2);
    mdwrap_wait((void *) &flag, sizeof(int), &from, &type, &st, 
                 &request2);
  }

  mdwrap_request_free(&request); /* Free request object from iwrite call */
  /*
   * Do a final sync amongst all the processors, so that all of the other
   * processors must wait for Proc 0 to receive the final message from
   * Proc (Num_Proc-1).
   */

  AZ_sync(proc_config);

} /* AZ_print_sync_end */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_sync(int proc_config[])

/*******************************************************************************

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  node:            Current processor number.

  nprocs:          Number of processors in the current machine configuration.

*******************************************************************************/

{

  /* local variables */

  int   type;                     /* type of next message */
  int   partner;                  /* processor I exchange with */
  int   mask;                     /* bit pattern identifying partner */
  int   hbit;                     /* largest nonzero bit in nprocs */
  int   nprocs_small;             /* largest power of 2 <= nprocs */
  int   cflag;                    /* dummy argument for compatability */
  int   node, nprocs;
  char *yo = "sync: ";

  MPI_AZRequest request;  /* Message handle */

  /**************************** execution begins ******************************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];
  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /*  Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;
  if (nprocs_small*2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (mdwrap_iread((void *) NULL, 0, &partner, &type, &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (mdwrap_write((void *) NULL, 0, partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /*
     * Wait to receive the messages.  These messages will return length 1
     * because MPI will not necessarily send a zero-length message.
     */

    (void) mdwrap_wait((void *) NULL, 0, &partner, &type, &cflag, &request);
  }

  /*  Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small>>1; mask; mask >>= 1) {
      partner = node ^ mask;
      if (mdwrap_iread((void *) NULL, 0, &partner, &type, &request)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_write((void *) NULL, 0, partner, type, &cflag)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      /*
       * Wait to receive the messages.  These messages will return length 1
       * because MPI will not necessarily send a zero-length message.
       */

      (void) mdwrap_wait((void *) NULL, 0, &partner, &type, &cflag, &request);
    }
  }

  /*  Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (mdwrap_iread((void *) NULL, 0, &partner, &type, &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (mdwrap_write((void *) NULL, 0, partner, type, &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  /*
   * Wait to receive the messages.  These messages will return length 1
   * because MPI will not necessarily send a zero-length message.
   */

  if (node & nprocs_small) {
    (void) mdwrap_wait((void *) NULL, 0, &partner, &type, &cflag, &request);
  }

} /* AZ_sync */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_gsum_vec_int(int vals[], int vals2[], int length, int proc_config[])

/*******************************************************************************

  For each element in vals[], perform a global sum with the other processors.
  That is, on output vals[i] is equal to the sum of the input values in vals[i]
  on all the processors.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  vals:            On input, vals[i] on this processor is to be summed with
                   vals[i] on all the other processors.
                   On output, vals[i] is the sum of the input values in val[i]
                   defined on all processors.

  vals2:           Work space of size 'length'.

  node:            Current processor number.

  nprocs:          Number of processors in the current machine configuration.

  length:          Number of values in 'vals' (i.e. number of global sums).

*******************************************************************************/

{

  /* local variables */

  int   type;             /* type of next message */
  int   partner;          /* processor I exchange with */
  int   mask;             /* bit pattern identifying partner */
  int   hbit;             /* largest nonzero bit in nprocs */
  int   nprocs_small;     /* largest power of 2 <= nprocs */
  int   cflag;            /* dummy argument for compatability */
  int   k;
  int   node, nprocs;
  char *yo = "AZ_gsum_vec_int: ";

  MPI_AZRequest request;  /* Message handle */

  /*********************** first executable statment *****************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;

  if (nprocs_small * 2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;
  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if (mdwrap_iread((void *) vals2, length*sizeof(int), &partner, &type,
                      &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (mdwrap_write((void *) vals, length*sizeof(int), partner, type,
                      &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    if (mdwrap_wait((void *) vals2, length*sizeof(int), &partner, &type,
                     &cflag, &request) < length*sizeof(int)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }

    /* sum values */

    for (k = 0; k < length; k++) vals[k] += vals2[k];
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small >> 1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (mdwrap_iread((void *) vals2, length*sizeof(int), &partner, &type,
                        &request)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_write((void *) vals, length*sizeof(int), partner, type,
                        &cflag)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_wait((void *) vals2, length*sizeof(int), &partner, &type,
                       &cflag, &request) < length*sizeof(int)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      for (k = 0; k < length; k++) vals[k] += vals2[k];
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (mdwrap_iread((void *) vals, length*sizeof(int), &partner, &type,
                      &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (mdwrap_write((void *) vals, length*sizeof(int), partner, type,
                      &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    if (mdwrap_wait((void *) vals, length*sizeof(int), &partner, &type, &cflag,
                     &request) < length*sizeof(int)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }


} /* AZ_gsum_vec_int */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_gappend(int vals[], int *cur_length, int total_length,
                int proc_config[])

/*******************************************************************************

  Take the list contained in vals[] that is defined on this processor and append
  it to the other lists (vals[]) defined on all the other processors.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  vals:            On input, vals[i] on this processor contains list to be
                   appended with vals[i] on all the other processors.
                   On output, vals[i] is the concatenation of the input lists
                   in val[i] defined on all processors.

  cur_length:      On input, initial length of vals[].
                   On output, length of new vals[] (after concatenation).

  total_length:    Maximum allowable length for vals[].

  node:            Current processor number.

  nprocs:          Number of processors in the current machine configuration.

*******************************************************************************/

{

  /* local variables */

  int   type;         /* type of next message */
  int   partner;      /* processor I exchange with */
  int   mask;         /* bit pattern identifying partner */
  int   hbit;         /* largest nonzero bit in nprocs */
  int   nprocs_small; /* largest power of 2 <= nprocs */
  int   cflag;        /* dummy argument for compatability */
  int   length;
  int   node, nprocs;
  char *yo = "AZ_gappend: ";

  MPI_AZRequest request;  /* Message handle */

  /*********************** first executable statment *****************/

  node   = proc_config[AZ_node];
  nprocs = proc_config[AZ_N_procs];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find next lower power of 2. */

  for (hbit = 0; (nprocs >> hbit) != 1; hbit++);

  nprocs_small = 1 << hbit;

  if (nprocs_small * 2 == nprocs) {
    nprocs_small *= 2;
    hbit++;
  }

  partner = node ^ nprocs_small;

  if (node+nprocs_small < nprocs) {

    /* post receives on the hypercube portion of the machine partition */

    if ( mdwrap_iread((void *) &(vals[*cur_length]),
                       (total_length - *cur_length) * sizeof(int), &partner,
                       &type, &request) ) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }
  else if (node & nprocs_small) {

    /*
     * Send messages from the portion of the machine partition "above" the
     * largest hypercube to the hypercube portion.
     */

    if (mdwrap_write((void *) vals, (*cur_length)*sizeof(int), partner, type,
                      &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node+nprocs_small < nprocs) {

    /* wait to receive the messages */

    length = mdwrap_wait((void *) &(vals[*cur_length]),
                          (total_length - *cur_length)*sizeof(int), &partner,
                          &type, &cflag, &request);
    (*cur_length) += (length / sizeof(int));
  }

  /* Now do a binary exchange on nprocs_small nodes. */

  if (!(node & nprocs_small)) {
    for (mask = nprocs_small >> 1; mask; mask >>= 1) {
      partner = node ^ mask;

      if (mdwrap_iread((void *) &(vals[*cur_length]),
                        (total_length - *cur_length)*sizeof(int), &partner,
                        &type, &request)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      if (mdwrap_write((void *) vals, *cur_length*sizeof(int), partner, type,
                        &cflag)) {
        (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                       "type = %d\n", yo, node, type);
        exit(-1);
      }

      length = mdwrap_wait((void *) &(vals[*cur_length]),
                            (total_length - *cur_length)*sizeof(int), &partner,
                            &type, &cflag, &request);
      (*cur_length) += (length / sizeof(int));
    }
  }

  /* Finally, send message from lower half to upper half. */

  partner = node ^ nprocs_small;
  if (node & nprocs_small) {
    if (mdwrap_iread((void *) vals, total_length*sizeof(int), &partner, &type,
                      &request)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_iread failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  else if (node+nprocs_small < nprocs ) {
    if (mdwrap_write((void *) vals, *cur_length*sizeof(int), partner, type,
                      &cflag)) {
      (void) AZ_printf_err( "%sERROR on node %d\nmd_write failed, message "
                     "type = %d\n", yo, node, type);
      exit(-1);
    }
  }

  if (node & nprocs_small) {
    length = mdwrap_wait((void *) vals, total_length*sizeof(int), &partner,
                          &type, &cflag, &request);
    (*cur_length) = (length / sizeof(int));
  }

} /* AZ_gappend */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_broadcast(char *ptr, int length, int proc_config[], int action)

/*******************************************************************************

  Used to concatenate a buffer of information and to then broadcast this
  information from processor 0 to the other processors. The four possiblities
  are

    1) action == AZ_PACK  && proc_config[AZ_node] == 0
          Store ptr into internal buffer (brdcst_buffer).
    2) action == AZ_PACK  && proc_config[AZ_node] != 0
          Read from the internal buffer (brdcst_buffer) to ptr. If the
          internal buffer is empty, first receive broadcast information
    3) action == AZ_SEND  && proc_config[AZ_node] == 0
          Broadcast internal buffer (filled by previous AZ_broadcast
          calls) and then clear it.
    4) action == AZ_SEND  && proc_config[AZ_node] != 0
          Clear the internal buffer.


  Sample usage:

    if (proc_config[AZ_node] == 0) {
       a = 1; b = 2;
    }
    AZ_broadcast(&a, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast(&b, sizeof(int), proc_config, AZ_PACK);
    AZ_broadcast(NULL,         0, proc_config, AZ_SEND);

  Note: There can be no other communication calls between AZ_PACK
  and AZ_SEND calls to AZ_broadcast.

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  ptr:             If action == AZ_PACK,
                   If proc==0, ptr is string to be broadcast to other nodes.
                   If proc~=0, on output ptr is string received from node 0.

  length:          Length of string to be broadcast/received.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  action:          Determines what AZ_broadcast() does as described above.

*******************************************************************************/

{

  /* local variables */

  int i;
  static char *brdcst_buffer = 0;
  static int   buffer_start = 0;
  char   *temp;
  int    *tt;
  static unsigned int   buf_len = 0, buffer_end = 0;

  /**************************** execution begins ******************************/

  if (action == AZ_PACK) {

    /* allocate the buffer */

    if (brdcst_buffer == 0) {
      buf_len = 1000;    /* Note: this is coordinated with the     */
                         /* statement 'if (buf_len != 1000)' below */

      brdcst_buffer = (char *) AZ_allocate(buf_len * sizeof(char));
      if (brdcst_buffer == NULL) {

        (void) AZ_printf_err( "no space in AZ_broadcast: brdcst_buffer\n");
        exit(-1);
      }
    }

    /* If processor 0, pack buffer */

    if (proc_config[AZ_node] == 0) {

      if (buffer_end+length > buf_len)  {

        /* Buffer is not big enough. Allocate more space */

        buf_len += AZ_MAX(500, length);
        temp = (char *) AZ_allocate(buf_len * sizeof(char));
        if (temp == NULL) {
          (void) AZ_printf_err( "no space in AZ_broadcast: temp\n");
          exit(-1);
        }

        if (brdcst_buffer != 0) {
          for (i = 0; i < (int) buffer_end; i++) temp[i] = brdcst_buffer[i];
          AZ_free(brdcst_buffer);
        }
        brdcst_buffer = temp;
      }

      if (brdcst_buffer == 0) {
        (void) AZ_printf_err(
                       "Error: Not enough space in AZ_broadcast_pack\n");
        exit(-1);
      }

      for (i = 0; i < length; i++) brdcst_buffer[i + buffer_end] = ptr[i];
      buffer_end += length;
    }

    /* For processors other than 0 ... */

    else {

      /*
       * If the buffer is empty, do a broadcast (i.e. post a read) to obtain
       * broadcast information.
       */

      if (buffer_end == 0) {

        buffer_end = AZ_broadcast_info(brdcst_buffer, proc_config, buf_len);

        /*
         * A single broadcasted integer, indicates that processor 0 increased
         * its buffer size during the packing stage. The size of this buffer is
         * in fact the value received.
         */

        if (buffer_end == sizeof(int)) {

          /* allocate a new buffer */

          tt = (int *) brdcst_buffer;
          buf_len = tt[0];
          AZ_free(brdcst_buffer);
          brdcst_buffer = (char *) AZ_allocate(buf_len * sizeof(char));
          if ( brdcst_buffer == NULL) {
            (void) AZ_printf_err(
                           "no space in AZ_broadcast: brdcst_buffer \n");
            exit(-1);

          }

          /* Do the broadcast again with the larger buffer */

          buffer_end = AZ_broadcast_info(brdcst_buffer, proc_config, buf_len);
        }
      }

      /* take the top elements off of 'buffer' and store them into 'ptr'. */

      for (i = 0; i < length; i++) ptr[i] = brdcst_buffer[buffer_start+i];
      buffer_start += length;
    }
  }

  else {
    if (proc_config[AZ_node] == 0) {

      /*
       * If additional buffer space was needed, processor 0 tells the others by
       * sending just 1 integer (the length of the new buffer)
       */

      if (buf_len != 1000)
        (void) AZ_broadcast_info((char *) &buffer_end,proc_config,sizeof(int));

      /*
       * If only 1 integer needs to be sent, we increase the number of bytes
       * sent to distinguish this situation from the case above where 1 integer
       * which is the length of the new buffer is sent
       */

      if (buffer_end == sizeof(int)) buffer_end++;

      /* broadcast the data */

      (void) AZ_broadcast_info(brdcst_buffer, proc_config, buffer_end);
    }

    /* clear the internal buffer */

    if (brdcst_buffer != (char *) NULL) AZ_free(brdcst_buffer);
    brdcst_buffer = (char *) NULL;
    buf_len = buffer_end = buffer_start = 0;
  }

} /* AZ_broadcast */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

unsigned int AZ_broadcast_info(char buffer[], int proc_config[], 
                               unsigned int length)

/*******************************************************************************

  The information in 'buffer' on processor 0 is broadcast to all other
  processors.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     unsigned int, length of string broadcast.
  ============

  Parameter list:
  ===============

  buffer:          Buffer used to pack information to be broadcast.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

  length:          Length of string to be broadcast/received.

*******************************************************************************/

{

  /* local variables */

  int i;
  int type;         /* type of next message */
  int partner;      /* processor I exchange with */
  int hbit;         /* largest nonzero bit in nprocs */
  int my_lbit;      /* smallest nonzero bit in proc */
  int cflag;        /* dummy argument for compatability */
  int nprocs, proc;

  MPI_AZRequest request;  /* Message handle */

  /*********************** first executable statment *****************/

  nprocs = proc_config[AZ_N_procs];
  proc   = proc_config[AZ_node];

  type            = proc_config[AZ_MPI_Tag];
  proc_config[AZ_MPI_Tag] = (type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* Find location of largest bit. */

  for (hbit = 0; ((nprocs - 1) >> hbit) != 0; hbit++);

  /* Find location of smallest bit corresponding to my processor name */

  if (proc != 0)
    for (my_lbit = 1; (proc | (1 << (my_lbit - 1))) != proc; my_lbit++);
  else my_lbit = hbit + 1;

  /* Zero out lowest bit in proc ... and receive from that processor */

  if (proc != 0) {
    partner = proc ^ (1 << (my_lbit - 1));
    (void) mdwrap_iread((void *) buffer, length, &partner, &type, &request);

    /* wait for messages */

    length  = mdwrap_wait((void *) buffer, length, &partner, &type, &cflag,
                           &request);
  }

  /*
   * Send to neighbors. The neighbors are defined by putting a 1 in one of the
   * locations in 'proc' to the right of the lowest nonzero bit.
   */

  for (i = my_lbit - 1; i > 0; i--) {
    partner = proc | (1 << (i - 1));
    if (partner < nprocs)
      (void) mdwrap_write((void *) buffer, length, partner, type, &cflag);
  }

  return length;

} /* AZ_broadcast_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double AZ_sync_timer(int proc_config[])

/*******************************************************************************

  Use to synchronize the processors and take a timing.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     double, current time.
  ============

  Parameter list:
  ===============

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

*******************************************************************************/

{

  int i = 0;

  i = AZ_gsum_int(i, proc_config);

  return AZ_second();

} /* AZ_sync_timer */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_splitup_big_msg(int num_neighbors, char *ibuffer, char *obuffer,
			unsigned int element_size, int *start_send_proc,
                        int *actual_send_length, int *actual_recv_length, int
                        *proc_num_neighbor, int type, int *total_num_recv,
                        int *proc_config)

/*******************************************************************************

   Author:        John Shadid, SNL, 1421 Date: 8/1/94
  =======

  Return code:    (void).
  ============

  Parameter list:
  ===============

  num_neighbors:      total number of neighbors to communicate with
  buffer:             on input  - buffer that holds the information to be sent
                      on output - buffer contains the recieve information
  start_send_proc:    contains index of start of information that is to be
                      sent to each processor
  actual_send_length: contains number of double precision entries to be sent to
                      each processor
  actual_recv_length: number of entries to be recieve from each processor
  proc_num_neighbor:  contains the processor number for each neighbor
  type:               the message type to be used in the mesage
  total_num_recv:     on output - total number of actual recvieved entried
  proc_config:        Machine configuration.  proc_config[AZ_node] is the node
                      number.  proc_config[AZ_N_procs] is the number
                      of processors.

*******************************************************************************/

{

  /*
   * This function handshakes big messages between all the neighbors.  The
   * length of the messages are calculated conservatively to not allow overflow
   * of the message buffers.
   */

  int     m, n, st, rtype, j, dummy_int;
  int     max_neighbors, messg_size_doubles, doubles_sent;
  int     total_doubles_to_send, dest, flag, messg_from, messg_type;
  int     total_doubles_to_recv, total_send_size;
  int     finished_send_messg[AZ_MAX_NEIGHBORS];
  int     finished_recv_messg[AZ_MAX_NEIGHBORS];
  int     number_of_messages, start_recv_proc[AZ_MAX_NEIGHBORS];
  int     allowed_buff_size, num_recv;
  int     max_buffer_size = 0, max_messg_size;
  char    *send_buffer;
  char   *char_ptr;
  char   *yo = "AZ_splitup_big_msg ";
  int     split_up = FALSE;
  int     dummy_add;
  int     DEBUG = FALSE;
  unsigned int length, size;
  
  
  MPI_AZRequest request[AZ_MAX_NEIGHBORS];  /* Message handle */

  /**************************** execution begins ****************************/

  /* Compute the global maximum message buffer size needed */

  
  for (n = 0; n < num_neighbors; n++) {
    max_buffer_size += actual_recv_length[n];
  }
  max_buffer_size = AZ_gmax_int(max_buffer_size, proc_config);

  /* Determine if splitting of messages is necessary */

  if (max_buffer_size > (int) (AZ_MAX_MSG_BUFF_SIZE / (2 * element_size))) {

     /* Too big for message buffers */

     split_up = TRUE;
  }

  if (split_up == TRUE) {

    /*
     * Compute maximum total message size in bytes that any processor will
     * recieve and the maximum number of neighbors that any proc must
     * communicate with. Also initalize some logical arrays.
     */

    max_messg_size = 0;
    for (n = 0; n < num_neighbors; n++) {
      max_messg_size = AZ_MAX(max_messg_size, actual_recv_length[n]);
      finished_send_messg[n] = finished_recv_messg[n] = AZ_FALSE;
    }
    max_messg_size = AZ_gmax_int(max_messg_size, proc_config);
    max_neighbors  = AZ_gmax_int(num_neighbors, proc_config);

    /*
     * Total received nonzeros and starting location for each processors
     * message that will be received
     */

    num_recv = 0;
    for (n = 0; n < num_neighbors; n++) {
      start_recv_proc[n] = num_recv;
      num_recv          += actual_recv_length[n];
    }
    *total_num_recv = num_recv;

    /*
     * Compute the global maximum allowed message size and the maximum number of
     * messages to send all the required information.
     */

    allowed_buff_size  = (int) floor(((double) AZ_MAX_MSG_BUFF_SIZE /
                                      (double) (3*element_size)));

    messg_size_doubles = (int) floor((double) allowed_buff_size / 
                                      (double) max_neighbors);

    number_of_messages = (int) ceil((double) max_messg_size /
                                    (double) (messg_size_doubles));

    if (proc_config[AZ_node] == 0 && DEBUG == TRUE) {
      (void) AZ_printf_out("\n\t\tSplitting up messages in splitup_big_msg\n"
                    "\t\tmax_buffer_size required  (bytes): %d\n",
                    max_buffer_size*element_size);
      (void) AZ_printf_out("\t\tmax_buffer_size allocated (bytes): %d\n",
                    allowed_buff_size*element_size);
      (void) AZ_printf_out("\t\tindividual message size   (bytes): %d\n",
                    messg_size_doubles*element_size);
      (void) AZ_printf_out("\t\ttotal number of split messages to be sent: %d\n\n",
                    number_of_messages);
    }

    if (ibuffer == obuffer) {
       /*
        * The input and output buffers are the same. Allocate a temporary 
        * send buffer that can hold all out going messages.
        * Then copy all info to this buffer.
        */

        total_send_size = 0;
        for (n = 0; n < num_neighbors; n++) {
           total_send_size += actual_send_length[n];
        }

        send_buffer =(char *) AZ_allocate((total_send_size+1)*element_size);
        if (send_buffer == NULL) {
           (void) AZ_printf_err(
                          "no space in AZ_splitup_big_msg: send_buffer \n");
           exit(-1);
        }
        for (n = 0; n < (int) (total_send_size*element_size) ; n++)
          send_buffer[n] = ibuffer[n];
    }
    else send_buffer = ibuffer;

    /*
     * Send and receive messages in a series of communications. Each set of
     * exchanges is followed by a syncronization to not allow message buffers to
     * overflow.
     */

    doubles_sent = 0;

    for (m = 0; m < number_of_messages; m++) {

      /* post recieves for split messages */

      for (n = 0; n < num_neighbors; n++) {

        total_doubles_to_recv = actual_recv_length[n];
        messg_from            = proc_num_neighbor[n];
        dummy_int             = type;

        if (doubles_sent + messg_size_doubles < total_doubles_to_recv ) {

          /* read messg_size_doubles bytes */

          length = messg_size_doubles*element_size;

          char_ptr = (char *) (&obuffer[start_recv_proc[n]*element_size] + 
                                       doubles_sent*element_size);
          (void) mdwrap_iread((void *) char_ptr, length, &messg_from, 
                               &dummy_int,  request+n);
        }
        else if (doubles_sent+messg_size_doubles >= total_doubles_to_recv &&
                 finished_recv_messg[n] == AZ_FALSE) {

          /* read actual_recv_length[n] - doubles_sent bytes */

          length = (total_doubles_to_recv - doubles_sent)*element_size;

          char_ptr = (char *) (&obuffer[start_recv_proc[n]*element_size] + 
                               doubles_sent*element_size);
          (void) mdwrap_iread((void *) char_ptr, length, &messg_from, 
                               &dummy_int,  request+n);
        }
        else if (finished_recv_messg[n] == AZ_TRUE) {

          /* read integer dummy message */

          length = sizeof(int);
          (void) mdwrap_iread((void *) &dummy_add, length, &messg_from, 
                                &dummy_int,  request+n);
        }
      }

      /* write split messages */

      for (n = 0; n < num_neighbors; n++) {
        total_doubles_to_send = actual_send_length[n];
        dest                  = proc_num_neighbor[n];

        if (doubles_sent + messg_size_doubles < total_doubles_to_send) {

          /* send out messg_size_doubles bytes */

          length = messg_size_doubles*element_size;
          char_ptr = (char *) (&send_buffer[element_size*start_send_proc[n]] + 
                               doubles_sent*element_size);
          (void) mdwrap_write((void *) char_ptr, length, dest, type, &flag);
        }
        else if (doubles_sent + messg_size_doubles >= total_doubles_to_send &&
                 finished_send_messg[n] == AZ_FALSE) {

          /* send out actual_send_length[n] - doubles_sent bytes */

          length = (total_doubles_to_send - doubles_sent)*element_size;

          char_ptr = (char *) (&send_buffer[start_send_proc[n]*element_size] + 
                               doubles_sent*element_size);
          (void) mdwrap_write((void *) char_ptr, length, dest, type, &flag);

          finished_send_messg[n] = AZ_TRUE;
        }
        else if (finished_send_messg[n] == AZ_TRUE) {

          /* send out integer dummy message */

          length = sizeof(int);
          (void) mdwrap_write((void *) &dummy_add, length, dest, type, &flag);
        }
      }

      /* read split messages */

      for (n = 0; n < num_neighbors; n++) {
        total_doubles_to_recv = actual_recv_length[n];
        messg_from            = proc_num_neighbor[n];
        messg_type            = type;

        if (doubles_sent + messg_size_doubles < total_doubles_to_recv ) {

          /* read messg_size_doubles bytes */

          length = messg_size_doubles*element_size;
          char_ptr = (char *) (&obuffer[start_recv_proc[n]*element_size] + 
                               doubles_sent*element_size);
          size =  mdwrap_wait((void *) char_ptr, length, &messg_from,
                             &messg_type, &flag, request+n); 

          if (length > size) {
           (void) AZ_printf_err("%sERROR on node %d\nmd_wait failed, message "
                          "type = %d\n", yo, proc_config[AZ_node], messg_type);
           exit(-1);
          }
        }
        else if (doubles_sent+messg_size_doubles >= total_doubles_to_recv &&
                 finished_recv_messg[n] == AZ_FALSE) {

          /* read actual_recv_length[n] - doubles_sent bytes */

          length = (total_doubles_to_recv - doubles_sent)*element_size;
          char_ptr = (char *) (&obuffer[start_recv_proc[n]*element_size] + 
                               doubles_sent*element_size);
          size =  mdwrap_wait((void *) char_ptr, length, &messg_from,
                             &messg_type, &flag, request+n); 

          if (length > size) {
           (void) AZ_printf_err("%sERROR on node %d\nmd_wait failed, message "
                          "type = %d\n", yo, proc_config[AZ_node], messg_type);
           exit(-1);
          }

          finished_recv_messg[n] = AZ_TRUE;
        }
        else if (finished_recv_messg[n] == AZ_TRUE) {

          /* read integer dummy message */

          length = sizeof(int);
          size =  mdwrap_wait((void *) &dummy_add, length, &messg_from,
                             &messg_type, &flag, request+n); 

          if (length > size) {
           (void) AZ_printf_err("%sERROR on node %d\nmd_wait failed, message "
                          "type = %d\n", yo, proc_config[AZ_node], messg_type);
           exit(-1);
          }

        }
      }

      doubles_sent += messg_size_doubles;


      AZ_sync(proc_config);
    }

    if (ibuffer == obuffer) AZ_free(send_buffer);
    return;
  }

  else {
     
     if (ibuffer == obuffer ) {
        /* Allocate a send buffer, if the input */
        /* and output buffers are the same.     */
        total_send_size = 0;
        for (n = 0; n < num_neighbors; n++) {
           total_send_size += actual_send_length[n];
        }
        send_buffer = (char *) AZ_allocate((total_send_size+1)*element_size);
        if (send_buffer == NULL) {
           (void) AZ_printf_err("no space AZ_splitup_big_msg: send_buffer \n");
           exit(-1);
        }
   
        for (n = 0; n < (int) (total_send_size*element_size) ; n++) 
           send_buffer[n] = ibuffer[n];
     }
     else send_buffer = ibuffer;
     
     /* post receives for message */
     
     j = 0;
     for (n = 0; n < num_neighbors; n++) {
        messg_from = proc_num_neighbor[n];
        dummy_int = type;
        size      = actual_recv_length[n]*element_size;

        (void) mdwrap_iread((void *) &obuffer[j], size, &messg_from, 
                             &dummy_int, request+n);
        j += actual_recv_length[n]*element_size;
     }

     /* send messages to each neighbor */

     for (n = 0; n < num_neighbors; n++) {
        size = actual_send_length[n]*element_size;
        (void) mdwrap_write((void *) &send_buffer[start_send_proc[n]*
                             element_size], size, proc_num_neighbor[n], type, 
                             &st);
     }             

     /* wait for all messages */

     j = 0;
     for (n = 0; n < num_neighbors; n++) {
        messg_from = proc_num_neighbor[n];
        rtype     = type;
        size      = actual_recv_length[n]*element_size;
        length =  mdwrap_wait((void *) &obuffer[j], size, &messg_from,
                               &rtype, &st, request+n); 
        if ((length != size) && (size !=0) ) {
           (void) AZ_printf_err( "%sERROR on node %d\nmd_wait failed, message "
                          "type = %d\n", yo, proc_config[AZ_node] , rtype);
           exit(-1);
        }
        j += length;
     }
     *total_num_recv = j/element_size;
     if (ibuffer == obuffer) AZ_free(send_buffer);
  }

} /* AZ_splitup_big_msg */

int AZ_extract_comm_info(int **idata_org,int (*user_comm)(double *,AZ_MATRIX *),
	 AZ_MATRIX *Amat, int proc_config[], int N_cols, int Nghost) {

/*
 * Fill in the pre-communication struction of an ML_Operator's getrow by
 * using a communication routine supplied by the user.
 */

   double *data_vec;
   int    *procs, *tempo;
   int    i, j, index, N_rcv_procs, *rcv_neighbors, *send_neighbors, **rcv_list;
   int    *rcv_length, *send_length, N_send_procs;
   int    proc_id, nprocs, *data_org;
   int *send_list, *collect;
   MPI_AZRequest *request;
   int type, start_send;
   int status, *sort_ind, *start_index, k, tt, jj;

   if (user_comm == NULL) {
      data_org = (int  *) AZ_allocate(AZ_COMMLESS_DATA_ORG_SIZE*sizeof(int));
      if (data_org == NULL) {
         AZ_printf_err("Error: Not enough space in AZ_extract_comm_info().\n");
         exit(1);
      }
      data_org[AZ_N_neigh]    = 0;
      data_org[AZ_total_send] = 0;
      data_org[AZ_internal_use] = 0;
      *idata_org = data_org;
      return 0;
   }
   proc_id  = proc_config[AZ_node];
   nprocs   = proc_config[AZ_N_procs];

   data_vec = (double *) AZ_allocate((N_cols + Nghost + 1)*sizeof(double));
   procs    = (int    *) AZ_allocate(nprocs  * sizeof(int));
   tempo    = (int    *) AZ_allocate(nprocs * sizeof(int));
   if ( (data_vec == NULL) || (procs == NULL) || (tempo == NULL)) {
       AZ_printf_err("AZ_extract_comm_info: Out of space\n"); exit(1);
   }

   for (i = 0; i < N_cols+Nghost; i++) data_vec[i] = (double ) proc_id;
   user_comm(data_vec, Amat);

   /* Compute the number of elements recvd from the ith */
   /* processor in procs[i] and store the total number  */
   /* of processors from which we rcv in 'N_rcv_procs'. */

   for (i = 0; i < nprocs; i++) procs[i] = 0;
   for (i = 0; i < N_cols+Nghost ; i++) {
      procs[ (int) data_vec[i] ]++;
   }
   procs[ proc_id ] = 0;

   N_rcv_procs = 0;
   for (i = 0; i < nprocs; i++) {
      if ( procs[i] > 0) N_rcv_procs++;
   }
   /* Store the neighbors in rcv_neighbors[k] and the number of */
   /* elements rcvd in rcv_length[k]. Allocate an array to hold */
   /* the rcv list and finally store the index k in procs[i].   */

   rcv_neighbors  = (int * ) AZ_allocate( (N_rcv_procs+1)*sizeof(int));
   rcv_length     = (int * ) AZ_allocate( (N_rcv_procs+1)*sizeof(int));
   rcv_list       = (int **) AZ_allocate( (N_rcv_procs+1)*sizeof(int *));
   if ( rcv_list == NULL) {
       AZ_printf_err("AZ_extract_comm_info: Out of space\n"); exit(1);
   }
   N_rcv_procs    = 0;
   for (i = 0; i < nprocs; i++) {
      if ( procs[i] > 0) {
         rcv_neighbors[N_rcv_procs] = i;
         rcv_list[N_rcv_procs] = (int *) AZ_allocate(procs[i] * sizeof(int) );
         if ( rcv_list[N_rcv_procs] == NULL) {
            AZ_printf_err("AZ_extract_comm_info: Out of space\n"); exit(1);
         }
         procs[i] = N_rcv_procs++;
      }
   }

   /* store the rcv list */

   for (i = 0; i < N_rcv_procs; i++) rcv_length[i] = 0;

   for (i = 0; i < N_cols+Nghost; i++) {
      j = (int) data_vec[i];
      if ( j != proc_id ) {
          index = procs[j];
          rcv_list[   index    ][ rcv_length[index]++ ] = i;
          if (i < N_cols) {
             AZ_printf_err("AZ_extract_comm_info: Received elements must be stored ");
             AZ_printf_err("after\n                   all %d local elements\n",N_cols);
             exit(1);
          }
      }
   }

   /* figure out the number of neighbors that we send to */

   for (i = 0; i < N_rcv_procs; i++) procs[rcv_neighbors[i]] = 1;
   AZ_gsum_vec_int(procs, tempo, nprocs, proc_config);
   N_send_procs = procs[proc_id];
   AZ_free(tempo);
   AZ_free(procs);

   /* figure out to whom we send and how many we send to them */

   i = N_send_procs + N_rcv_procs + 1;
   send_neighbors  = (int *    ) AZ_allocate( (i)*sizeof(int));
   send_length     = (int *    ) AZ_allocate( (i)*sizeof(int));
   request         = (MPI_AZRequest *)  AZ_allocate((N_send_procs+N_rcv_procs+1)*
                                              sizeof(MPI_AZRequest));
   if ( (send_neighbors==NULL)||(send_length==NULL)||(request==NULL)) {
       AZ_printf_err("AZ_extract_comm_info: Out of space\n"); exit(1);
   }

   type = 4901;
   for (i = 0; i < N_send_procs ; i++) {
     send_neighbors[i] = -1; /* receive from anyone */
     mdwrap_iread((void *) &(send_length[i]), sizeof(int) ,
                &(send_neighbors[i]), &type, request+i);
   }
   for (i = 0; i < N_rcv_procs; i++) {
      mdwrap_write((void *) &(rcv_length[i]), sizeof(int),
                          rcv_neighbors[i], type, &status);
   }
   for (i = 0; i < N_send_procs ; i++) {
     mdwrap_wait((void *) &(send_length[i]), sizeof(int) ,
                &(send_neighbors[i]), &type, &status, request+i);
   }
   AZ_sort( send_neighbors, N_send_procs , send_length, NULL);
   /* Fill in the send list */

   j = 1;
   for (i = 0; i < N_send_procs; i++) j += send_length[i];
   send_list = (int *) AZ_allocate(sizeof(int)*(j+1));
   if ( send_list == NULL) {
       AZ_printf_err("AZ_extract_comm_info: Out of space\n"); exit(1);
   }
   j = 1;
   for (i = 0; i < N_rcv_procs; i++)
      if (j < rcv_length[i]) j = rcv_length[i];
   collect = (int *) AZ_allocate(sizeof(int)*j);
   if ( collect == NULL) {
       AZ_printf_err("AZ_extract_comm_info: Out of space\n"); exit(1);
   }


   for (i = 0; i < N_cols+Nghost ; i++) data_vec[i] = (double) i;
   user_comm(data_vec, Amat);

   type++;
   j = 0;
   for (i = 0; i < N_send_procs ; i++) {
      mdwrap_iread((void *) &(send_list[j]), sizeof(int)*
                               send_length[i], &(send_neighbors[i]), &type,
                               request+i);
      j += send_length[i];
   }
   for (i = 0; i < N_rcv_procs; i++) {
      for (j = 0; j < rcv_length[i]; j++)
         collect[j] = (int) data_vec[ rcv_list[i][j] ];
      mdwrap_write((void *) collect, rcv_length[i]*sizeof(int),
                          rcv_neighbors[i], type, &status);
   }
   j = 0;
   for (i = 0; i < N_send_procs ; i++) {
      mdwrap_wait((void *) &(send_list[j]), sizeof(int)*
                             send_length[i], &(send_neighbors[i]), &type,
                             &status, request+i);
      j += send_length[i];
   }
   AZ_free(collect);
   AZ_free(request);
   AZ_free(data_vec);


   data_org = (int *) AZ_allocate((AZ_send_list+j)*sizeof(int));
   *idata_org = data_org;
   if (data_org == NULL) {
      AZ_printf_err("No space for data_org\n");
      exit(1);
   }
   data_org[AZ_total_send] = j;


   sort_ind    = (int *) AZ_allocate( sizeof(int)*(1+N_rcv_procs));
   start_index = (int *) AZ_allocate( sizeof(int)*(1+N_rcv_procs));
   if ( start_index == NULL) {
       AZ_printf_err("AZ_extract_comm_info: Out of space\n"); exit(1);
   }
   start_index[0] = N_cols;

   for (i = 0; i < N_rcv_procs; i++) {
      AZ_sort(rcv_list[i],rcv_length[i],NULL,NULL);
      for (k = 1; k < rcv_length[i]; k++) {
         if ( rcv_list[i][k] != rcv_list[i][k-1]+1 ) {
            AZ_printf_err("AZ_extract_comm_info: elements received from a processor ");
            AZ_printf_err("are not stored contiguously.\n");
            exit(1);
         }
      }
      if (rcv_length[i] > 0) start_index[i] = rcv_list[i][0];
      else start_index[i] = -1;
      sort_ind[i] = i;
   }
   AZ_sort( start_index, N_rcv_procs, sort_ind, NULL);
   if (start_index[0] != N_cols ) {
      AZ_printf_err("AZ_extract_comm_info: elements received do not start ");
      AZ_printf_err("immediately after\n                   local components.\n");
      exit(1);
   }
   for (i = 1; i < N_rcv_procs; i++) {
      if ( start_index[i-1] + rcv_length[sort_ind[i-1]] != start_index[i] ) {
          AZ_printf_err("AZ_extract_comm_info: elements received from processors ");
          AZ_printf_err("are not stored contiguously.\n");
          exit(1);
      }
   }
   AZ_free(start_index);

   j = 0;
   data_org[AZ_N_neigh] = N_rcv_procs;
   for (i = 0; i < N_rcv_procs; i++) {
      tt = sort_ind[i];
      data_org[AZ_neighbors  + i] = rcv_neighbors[tt];
      data_org[AZ_rec_length + i] = rcv_length[tt];
      start_send = 0;
      for (k = 0; k < N_send_procs; k++) {
         if ( send_neighbors[k] == rcv_neighbors[tt]) break;
         start_send += send_length[k];
      }
      if (k < N_send_procs) {
         send_neighbors[k] = -1;
         data_org[AZ_send_length + i] = send_length[k];
         for (jj = 0; jj < send_length[k]; jj++) {
            data_org[AZ_send_list+j] = send_list[start_send+jj];
            j++;
         }
      }
      else data_org[AZ_send_length + i] = 0;
   }
   AZ_free(sort_ind);

   start_send = 0;
   for (i = 0; i < N_send_procs; i++) {
      if (send_neighbors[i] != -1) {
         k = data_org[AZ_N_neigh];
         data_org[AZ_N_neigh]++;
         data_org[AZ_neighbors  + k] = send_neighbors[i];
         data_org[AZ_rec_length + k] = 0;
         data_org[AZ_send_length + k] = send_length[i];
         for (jj = 0; jj < send_length[i]; jj++) {
            data_org[AZ_send_list+j] = send_list[start_send+jj];
            j++;
         }
      }
      start_send += send_length[i];
   }
   AZ_free(rcv_neighbors);
   AZ_free(send_list);
   AZ_free(send_length);
   AZ_free(send_neighbors);
   for (i = 0; i < N_rcv_procs; i++) AZ_free(rcv_list[i]);
   AZ_free(rcv_list);
   AZ_free(rcv_length);
   return 1;
}
