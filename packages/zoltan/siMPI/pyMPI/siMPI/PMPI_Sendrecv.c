/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************      PMPI_Sendrecv.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 9 2002                                          */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int PMPI_Sendrecv( void *sendbuf, int sendcount, MPI_Datatype sendtype, 
    int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, 
    int source, int recvtag, MPI_Comm comm, MPI_Status *status )
{
  int retval;
 _MPI_COVERAGE();

  _MPI_CHECK_STATUS(&comm);
  retval = PMPI_Send(sendbuf, sendcount, sendtype, dest, sendtag, comm);
  if (retval == MPI_SUCCESS)
  {
    return PMPI_Recv(recvbuf, recvcount, recvtype, source, recvtag, comm, status);
  }
  return retval;
}
/*==========================================================================*/
