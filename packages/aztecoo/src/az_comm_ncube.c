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
#include "az_aztec.h"

extern int AZ_sys_msg_type;
int AZ_little_type;

/*
 * Ncube specific version of local communication routines.
 */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_write_local_info(int data_org[], char *message_recv_add[],
                      char *message_send_add[], int message_recv_length[],
                      int message_send_length[])

/*******************************************************************************

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  data_org:                  Array containing information on the distribution of
                             the matrix to this processor as well as
                             communication parameters (see file Aztec User's Guide).

  message_send_add:          message_send_add[i] points to the beginning of the
                             list of values to be sent to the ith neighbor
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

  int Num_Neighbors, *Proc_Num_Neighbor;
  int n, st;

  /**************************** execution begins ******************************/

  Num_Neighbors     = data_org[AZ_N_neigh];
  Proc_Num_Neighbor = &data_org[AZ_neighbors];

  AZ_little_type  = AZ_sys_msg_type;
  AZ_sys_msg_type = (AZ_sys_msg_type+1-AZ_MSG_TYPE) % AZ_NUM_MSGS + AZ_MSG_TYPE;

  /* write out all messages */

  for (n = 0; n < Num_Neighbors; n++) {
    (void) nwrite((char *) message_send_add[n], message_send_length[n],
                  Proc_Num_Neighbor[n], AZ_little_type, &st);
  }

} /* AZ_write_local_info */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void AZ_read_local_info(int data_org[], char *message_recv_add[],
                     int message_recv_length[])

/*******************************************************************************

  Author:          John N. Shadid, SNL, 1421
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  data_org:                  Array containing information on the distribution of
                             the matrix to this processor as well as
                             communication parameters (see file Aztec User's Guide).

  message_send_add:          message_send_add[i] points to the beginning of the
                             list of values to be sent to the ith neighbor
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

  int n, mesg_from, st;
  int Num_Neighbors, *Proc_Num_Neighbor;

  /**************************** execution begins ******************************/

  Num_Neighbors     = data_org[AZ_N_neigh];
  Proc_Num_Neighbor = &data_org[AZ_neighbors];

  for (n = 0; n < Num_Neighbors; n++) {
    mesg_from = Proc_Num_Neighbor[n];
    (void) nread((char *) message_recv_add[n], message_recv_length[n],
                 &mesg_from, &AZ_little_type, &st);
  }

} /* AZ_read_local_info */

