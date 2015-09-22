/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: rrdrake $
 *    Date: 2009/07/17 15:14:49 $
 *    Revision: 1.4 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************        PMPI_Wait.c       ************************/
/****************************************************************************/
/* Author : Lisa Alano July 9 2002                                          */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int PMPI_Wait(MPI_Request* request, MPI_Status* status) {
  int retval;
  MPI_Status recv_status;

  if ( !request ) {
    _MPI_ERR_ROUTINE(MPI_ERR_REQUEST,"Request pointer is null");
    return _MPI_NOT_OK;
  }

  /* ----------------------------------------------- */
  /* A null request requires no wait                 */
  /* ----------------------------------------------- */
  if ((*request) == MPI_REQUEST_NULL) return MPI_SUCCESS;

  /* ----------------------------------------------- */
  /* Send requests are always ready (eager)          */
  /* We may have to actually do a recv() here if it  */
  /* is not a send request.                          */
  /* ----------------------------------------------- */
  if (!(*request)->send) {
    retval = PMPI_Recv((*request)->buffer,
                       (*request)->count,
                       (*request)->type,
                       _MPI_RANK, 
                       (*request)->tag,
                       (*request)->comm,
                       &recv_status); 
    if ( retval == MPI_ERR_TAG && (*request)->cancel )
    {
      /* no matching send and the recv request has been cancelled */
      _MPI_Req_Invalid((*request));
      *request = MPI_REQUEST_NULL;
      return MPI_SUCCESS;
    }
    else if (retval != MPI_SUCCESS) {
      return retval;
    }
  }
  
  /* Copy in the status */
  if ( status && status != MPI_STATUS_IGNORE) {
    status->MPI_SOURCE = _MPI_RANK; 
    status->MPI_TAG = (*request)->tag;
    status->MPI_ERROR = MPI_SUCCESS;
    if ((*request)->send) {
      status->__count = _MPI_calculateSize((*request)->count, (*request)->type);
    } else {
      status->__count = recv_status.__count;
    }
  }
  
  /* ----------------------------------------------- */
  /* Mark the request available in the pool and then */
  /* write REQUEST_NULL back into the original req   */
  /* so that subsequent requests will immediately    */
  /* succeed.                                        */
  /* ----------------------------------------------- */
  _MPI_Req_Invalid((*request));
  *request = MPI_REQUEST_NULL;
  return MPI_SUCCESS;
}
/*==========================================================================*/

