/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: knogard $
 *    Date: 2007/07/02 16:48:06 $
 *    Revision: 1.3 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************     PMPI_Alltoallv.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 17 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Alltoallv ( 
        void *sendbuf, 
        int *sendcnts, 
        int *sdispls, 
        MPI_Datatype sendtype, 
        void *recvbuf, 
        int *recvcnts, 
        int *rdispls, 
        MPI_Datatype recvtype, 
        MPI_Comm comm )
{
  _MPI_COVERAGE();
  return PMPI_Alltoall(sendbuf,*sendcnts,sendtype,
                       recvbuf,*recvcnts,recvtype,comm);
}
