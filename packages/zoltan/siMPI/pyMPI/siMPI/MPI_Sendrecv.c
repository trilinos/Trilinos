/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************       MPI_Sendrecv.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 9 2002                                          */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int MPI_Sendrecv( void *sendbuf, int sendcount, MPI_Datatype sendtype, 
    int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, 
    int source, int recvtag, MPI_Comm comm, MPI_Status *status )
{
  _MPI_COVERAGE();
  return PMPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, 
                       recvcount, recvtype, source, recvtag, comm, status);
}
/*==========================================================================*/
