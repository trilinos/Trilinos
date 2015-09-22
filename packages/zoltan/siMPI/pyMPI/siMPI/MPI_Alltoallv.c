/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:49 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************      MPI_Alltoallv.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 17 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_Alltoallv ( 
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
  return PMPI_Alltoallv( sendbuf, sendcnts, sdispls, sendtype, recvbuf, recvcnts, rdispls, recvtype, comm); 
}

