/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ********************  PMPI_Request_free.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Request_free( MPI_Request *request ) {
  if (request == 0 ||
      _MPI_checkRequest(*request) !=MPI_SUCCESS) {
    _MPI_ERR_ROUTINE (MPI_ERR_REQUEST, "MPI_REQUEST_FREE: argument error");
    MPI_Abort (MPI_COMM_NULL, MPI_ERR_REQUEST);
    return MPI_ERR_REQUEST;
  }
  _MPI_Req_Invalid(*request);
  *request = MPI_REQUEST_NULL;
  return MPI_SUCCESS;
}

