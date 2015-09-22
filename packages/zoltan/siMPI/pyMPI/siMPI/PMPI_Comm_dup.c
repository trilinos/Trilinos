/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:52 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************     PMPI_Comm_dup.c      ************************/
/****************************************************************************/
/* Author : Lisa Alano July 18 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int PMPI_Comm_dup ( 
        MPI_Comm comm, 
        MPI_Comm *comm_out )
{
  _MPI_COVERAGE();
  if ( comm == MPI_COMM_NULL ) {
  _MPI_COVERAGE();
    if ( comm_out ) *comm_out = MPI_COMM_NULL;
    return MPI_SUCCESS;
  }
  return PMPI_Comm_create(comm,_MPI_COMM_WORLD_GROUP,comm_out);
}

