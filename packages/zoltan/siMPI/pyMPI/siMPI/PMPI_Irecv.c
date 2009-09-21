/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.3 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************       PMPI_Irecv.c       ************************/
/****************************************************************************/
/* Author : Lisa Alano June 27 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int PMPI_Irecv (void* message, int count, MPI_Datatype datatype, int source,
        int tag, MPI_Comm comm, MPI_Request* request)
{
  int retval, req;
  retval=_MPI_checks(message, count, datatype, source, tag, comm);

  if (retval != MPI_SUCCESS) {
    _MPI_ERR_ROUTINE (retval, "MPI_IRECV: argument error");
    MPI_Abort (comm, retval);
    return retval;
  }

  if (retval == MPI_SUCCESS && request != 0) {
    *request = _MPI_New_Request(message, count, datatype, tag, comm,_MPI_FALSE);
    if (*request == MPI_REQUEST_NULL) return _MPI_NOT_OK;
  }
  return retval;
}
/*==========================================================================*/


