/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  PMPI_Comm_split.c       ************************/
/****************************************************************************/
/* Author : Lisa Alano July 19 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Comm_split ( MPI_Comm comm, int color, int key, MPI_Comm *comm_out )
{
  if ( _MPI_Comm_check(comm) != MPI_SUCCESS ) return MPI_ERR_COMM;

  /* Returns a new communicator in groups of the same color --> Here, */
  /* it just returns a new communicator                               */
  return PMPI_Comm_create(comm,_MPI_COMM_WORLD_GROUP,comm_out);
}

