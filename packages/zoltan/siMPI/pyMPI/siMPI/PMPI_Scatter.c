/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************      PMPI_Scatter.c      ************************/
/****************************************************************************/
/* Author : Lisa Alano July 16 2002                                          */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int PMPI_Scatter (void *sendbuf, int sendcnt, MPI_Datatype sendtype, void *recvbuf, 
    int recvcnt, MPI_Datatype recvtype, int root, MPI_Comm comm) 
{
  int sendsize, recvsize, retval;
 _MPI_COVERAGE();

  _MPI_CHECK_STATUS(&comm);
  retval = _MPI_checks(sendbuf, sendcnt, sendtype, _MPI_RANK, MPI_ANY_TAG, comm);
  if (retval != MPI_SUCCESS) return retval;

  retval = _MPI_checks(recvbuf,recvcnt,recvtype, _MPI_RANK, MPI_ANY_TAG, comm);
  if (retval == MPI_SUCCESS) {
    recvsize = _MPI_calculateSize(recvcnt, recvtype);
    sendsize = _MPI_calculateSize(sendcnt, sendtype);
    if (recvsize < sendsize) {/*MESSAGE IS TRUNCATED*/
      recvbuf = memcpy(recvbuf, sendbuf, recvsize);
      printf("MPI_RECV : Message truncated.\n");
      MPI_Abort(comm, MPI_ERR_COUNT);
      return MPI_ERR_COUNT;
    } else {
      recvbuf = memcpy(recvbuf, sendbuf, sendsize);
    }
    return MPI_SUCCESS;
  }

 _MPI_COVERAGE();
  return _MPI_NOT_OK;
}
/*==========================================================================*/

