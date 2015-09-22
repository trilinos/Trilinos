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
#include <nx.h>
#include "az_aztec.h"

extern int az_iterate_id;

/*
 * Intel specific version of local communication routines.
 */

long int mid;
extern int AZ_sys_msg_type;

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_write_local_info(int data_org[], char *message_recv_add[],
                         char *message_send_add[], int message_recv_length[],
                         int message_send_length[])

/*******************************************************************************

  Author:
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  data_org:                  Array containing information on the distribution of
                             the matrix to this processor as well as
                             communication parameters (see file Aztec User's Guide).

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

*******************************************************************************/

{

  /* local variables */

  register int    n;
  long int        mid_new;
  static long int minfo[8];
  int             Num_Neighbors, *Proc_Num_Neighbor;
  int             type;

  /**************************** execution begins ******************************/

  Num_Neighbors     = data_org[AZ_N_neigh];
  Proc_Num_Neighbor = &data_org[AZ_neighbors];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* post irecvs for all messages and write out all messages */

  mid = -1;
  for (n = 0; n < Num_Neighbors; n++) {
    mid_new = irecvx(type, message_recv_add[n],
                     (long int) message_recv_length[n],
                     (long int) Proc_Num_Neighbor[n], -1, minfo);

    mid = msgmerge(mid, mid_new);
  }

  /* write messages */

  for (n = 0; n < Num_Neighbors; n++) {
    csend(type, message_send_add[n], (long int) message_send_length[n],
    (long int) Proc_Num_Neighbor[n], 0);
  }

} /* AZ_write_local_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_read_local_info(int data_org[], char *message_recv_add[],
                        int message_recv_length[])

/*******************************************************************************

  Author:
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  data_org:                  Array containing information on the distribution of
                             the matrix to this processor as well as
                             communication parameters (see file Aztec User's Guide).

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

*******************************************************************************/

{

  msgwait(mid);

} /* AZ_read_local_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/*
 * Intel specific version of local communication routines.
 */

long int mid;

void gather_and_send_mesgs(double x[], int data_org[], char *message_recv_add[],
                           char *message_send_add[], int message_recv_length[],
                           int message_send_length[])

/*******************************************************************************

  Routine to locally exchange the components of the vector "x". This routine
  gathers the necessary components of the vector and then sends the required
  "num_neighbors" messages. The messages which are received are placed
  contiguously in the external_nodes part of x.

  Author:          Bruce Hendrickson, SNL, 1422
  =======          John N. Shadid, SNL, 1421

  Return code:     void
  ============

  Parameter list:
  ===============

  x:                         Vector of unknowns defined on the current
                             processor. Indirect addressing will be used to
                             gather those unknowns that other processors need.

  data_org:                  Array containing information on the distribution of
                             the matrix to this processor as well as
                             communication parameters (see file Aztec User's Guide).

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

  Important external definitions:
  ===============================

  Num_Neighbors:             Total number of neighboring processors.
                             (type - int value)

  Proc_Num_Neighbor[]:       Array containing the processor id for each of the
                             current processor's neighbors.
                             (type - int vector with fixed length MAX_NEIGHBORS)

  Total_Num_Send_Unknowns:   Total number of unknowns which are sent to
                             neighboring processors (size of mesg buff).
                             (type - int value)

  Num_Unknowns_Send_Neighbor[]:
                             Vector containing the number of unknowns to be sent
                             to each processor.
                             (type - int vector with fixed length MAX_NEIGHBORS)

  List_Send_Unknowns[]:      Vector of local node numbers for the unknowns to be
                             sent to each neighbor.
                             (type - pointer to variable length int vector,
                             previously "AZ_allocate"ed)

  Num_Unknowns_Recv_Neighbor[]:
                             Array of number of unknowns to be received from
                             each of the neighboring processors.
                             (type - int vector with fixed length MAX_NEIGHBORS)

*******************************************************************************/

{

  /* local variables */

  double          *ptr_send_list,*ptrd, *ptr_recv_list;
  register double *ptrd2;
  register int    *ptr_int, i;
  int              size, num_send, num_recv;
  int              offset;
  int              n;
  long int         mid_new;
  static long int  minfo[8];
  int              type;
  int              Num_Neighbors, *Proc_Num_Neighbor,External_Index;
  int             *Num_Unknowns_Send_Neighbor, *Num_Unknowns_Recv_Neighbor;
  int             *List_Send_Unknowns;

  /*********************** first executable statement *****************/

  Num_Neighbors              = data_org[AZ_N_neigh];
  Proc_Num_Neighbor          = &data_org[AZ_neighbors];
  External_Index             = data_org[AZ_N_internal] + data_org[AZ_N_border];
  Num_Unknowns_Send_Neighbor = &data_org[AZ_send_length];
  Num_Unknowns_Recv_Neighbor = &data_org[AZ_rec_length];
  List_Send_Unknowns         = &data_org[AZ_send_list];

  type            = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* single processor case */

  if (Num_Neighbors == 0) return;

  /* Define arrays for message passing */

  ptrd = (double *) AZ_manage_memory(data_org[AZ_total_send]*sizeof(double),
                                     AZ_ALLOC, AZ_SYS+az_iterate_id, "ptrd", &n);
  ptr_send_list = ptrd;
  ptr_recv_list = &x[External_Index];

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

  /* Post irecvs for all messages */

  mid = -1;
  for (n = 0; n < Num_Neighbors; n++) {
    mid_new = irecvx(type, message_recv_add[n],
                     (long int) message_recv_length[n],
                     (long int) Proc_Num_Neighbor[n], -1, minfo);

    mid = msgmerge(mid, mid_new);
  }

  /* Set up send messages:   Gather send and send unknowns from "x" vector */

  for (n = 0; n < Num_Neighbors; n++) {

    /* First gather the data. */

    ptrd2   = (double *) message_send_add[n];
    offset  = ((int) ptrd2 - (int) ptrd) / sizeof(double);
    ptr_int = &(List_Send_Unknowns[offset]);
    for (i = message_send_length[n] / sizeof(double); i; i--) {
      *ptrd2++ = x[*ptr_int++];
    }

    /* Now ship it out. */

    csend(type, message_send_add[n], (long int) message_send_length[n],
    (long int) Proc_Num_Neighbor[n], 0);
  }

} /* gather_and_send_mesgs */
