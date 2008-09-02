/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  **********************  MPI_Scatterv.c   ***************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Scatterv ( 
        void *sendbuf, 
        int *sendcnts, 
        int *displs, 
        MPI_Datatype sendtype, 
        void *recvbuf, 
        int recvcnt,  
        MPI_Datatype recvtype, 
        int root, 
        MPI_Comm comm )
{
  _MPI_COVERAGE();
  return PMPI_Scatterv (sendbuf, sendcnts, displs, sendtype, recvbuf, recvcnt, recvtype, root, comm);
}

