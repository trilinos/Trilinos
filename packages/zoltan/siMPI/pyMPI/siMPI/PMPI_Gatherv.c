/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  PMPI_Gatherv.c   ****************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int PMPI_Gatherv ( void *sendbuf, int sendcnt, MPI_Datatype sendtype, 
                  void *recvbuf, int *recvcnts, int *displs, 
                 MPI_Datatype recvtype, 
                  int root, MPI_Comm comm )
{
  _MPI_COVERAGE();
  return PMPI_Gather(sendbuf,sendcnt,sendtype,recvbuf,*recvcnts,recvtype,root,comm);
}

