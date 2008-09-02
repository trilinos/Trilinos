/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  PMPI_Reduce_scatter.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int PMPI_Reduce_scatter ( void *sendbuf, void *recvbuf, int *recvcnts, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm ) {
  if ( !recvcnts ) {
    _MPI_ERR_ROUTINE (MPI_ERR_COUNT, "MPI_ERR_COUNT : Invalid recv count argument");
    MPI_Abort(comm, MPI_ERR_OTHER);
    return _MPI_NOT_OK;
  }

  return PMPI_Reduce(sendbuf,recvbuf,*recvcnts,datatype,op,0,comm);
}

