/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************     PMPI_Alltoall.c      ************************/
/****************************************************************************/
/* Author : Lisa Alano July 17 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"
#include <string.h>

/* STUB */
int PMPI_Alltoall( void *sendbuf, int sendcount, MPI_Datatype sendtype, 
                  void *recvbuf, int recvcount, MPI_Datatype recvtype, 
                 MPI_Comm comm )
{
  int send_size;
  int recv_size;

 _MPI_COVERAGE();

  if ( sendbuf == 0 ) return MPI_ERR_ARG;
  if ( sendcount <= 0 ) return MPI_ERR_ARG;
  if ( recvbuf == 0 ) return MPI_ERR_ARG;
  if ( recvcount <= 0 ) return MPI_ERR_ARG;

  _MPI_CHECK_STATUS(&comm);


  /* TODO: This is not really correct */
  switch( _MPI_checkSendType(sendtype) ) {
  case _MPI_DEFAULT:
    _MPI_COVERAGE();
    switch( _MPI_checkSendType(recvtype) ) {
    case _MPI_DEFAULT:
      _MPI_COVERAGE();
      send_size = _MPI_calculateSize(sendcount, sendtype);  
      recv_size = _MPI_calculateSize(recvcount, recvtype);  
      if ( send_size < recv_size ) {
        _MPI_COVERAGE();
	memcpy(recvbuf,sendbuf,send_size);
      } else {
        _MPI_COVERAGE();
	memcpy(recvbuf,sendbuf,recv_size);
      }
      return MPI_SUCCESS;
    default:
      return MPI_Abort(comm, MPI_UNDEFINED); 
    }
  default:
    return MPI_Abort(comm, MPI_UNDEFINED); 
  }
}

