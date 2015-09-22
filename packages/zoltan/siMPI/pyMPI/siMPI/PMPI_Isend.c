/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.2 $
 ****************************************************************************/
/***********************************************************************************************/
/* FILE  **************************      PMPI_Isend.c          *********************************/
/***********************************************************************************************/
/* Author : Lisa Alano June 27 2002                                                            */
/* Copyright (c) 2002 University of California Regents                                         */
/***********************************************************************************************/
#include "mpi.h"

/*=============================================================================================*/
int PMPI_Isend (void* message, int count, MPI_Datatype datatype, int dest,
        int tag, MPI_Comm comm, MPI_Request* request)
{
  int retval;
  retval = PMPI_Send(message, count, datatype, dest, tag, comm); 

  /* Fill in the request object so that we can query it later */
  if (retval == MPI_SUCCESS && request != 0) {
    *request = _MPI_New_Request(message, count, datatype, tag, comm,_MPI_TRUE);
    if (*request == MPI_REQUEST_NULL) return _MPI_NOT_OK;
  }
  return retval;
}
/*=============================================================================================*/

