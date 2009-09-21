/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *****************  MPI_Intercomm_create.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Intercomm_create ( MPI_Comm local_comm, int local_leader, 
                         MPI_Comm peer_comm, int remote_leader, int tag, 
                         MPI_Comm *comm_out )
{
  _MPI_COVERAGE();
  return PMPI_Intercomm_create (local_comm, local_leader, peer_comm, remote_leader, tag, comm_out);
}

