/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************     MPI_Allgatherv.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 17 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_Allgatherv ( void *sendbuf, int sendcount, MPI_Datatype sendtype, 
                     void *recvbuf, int *recvcounts, int *displs, 
                    MPI_Datatype recvtype, MPI_Comm comm )
{
  _MPI_COVERAGE();
  return PMPI_Allgatherv ( sendbuf, sendcount, sendtype, 
                     recvbuf, recvcounts, displs, 
                     recvtype, comm );
}

