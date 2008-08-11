/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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

