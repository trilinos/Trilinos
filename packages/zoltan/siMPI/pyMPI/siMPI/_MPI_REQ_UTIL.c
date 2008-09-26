/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/******************************************************************************/
/* FILE  ********************     _MPI_REQ_UTIL.c      ************************/
/******************************************************************************/
/* Author : Lisa Alano July 12 2002                                           */
/* Copyright (c) 2002 University of California Regents                        */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

MPI_Request _MPI_New_Request(void* buffer,
                             int count,
                             MPI_Datatype datatype,
                             int tag,
                             MPI_Comm comm,
                             int send) {
  MPI_Request request = 0;
  int i;
  int j;

  /* ----------------------------------------------- */
  /* We search through the array of lists looking    */
  /* for a request object that is not in use.        */
  /* ----------------------------------------------- */
  for(i=0;i<_MPI_REQ_ARRAY_SIZE;++i) {
    for(j=0;j<_MPI_PREALLOCATION_SIZE;++j) {
      if (_MPI_REQ_LIST_OF_LISTS[i][j].valid == _MPI_NOT_VALID) {
        request = _MPI_REQ_LIST_OF_LISTS[i]+j;
        break;
      }
    }
    if (request) break;
  }

  /* ----------------------------------------------- */
  /* Allocate another strip of the request list      */
  /* ----------------------------------------------- */
  if (!request) {
    /* Extend the list of lists (keeps old requests in same place */
    i = _MPI_REQ_ARRAY_SIZE++;
    _MPI_REQ_LIST_OF_LISTS = (_MPI_REQUEST_OBJECT **)
      _MPI_safeRealloc(_MPI_REQ_LIST_OF_LISTS,_MPI_REQ_ARRAY_SIZE*sizeof(_MPI_REQUEST_OBJECT*),
                       "Error with malloc of REQ_LIST");
    if (!_MPI_REQ_LIST_OF_LISTS) return MPI_REQUEST_NULL; /* malloc failure only */

    /* Allocate a new strip */
    _MPI_REQ_LIST_OF_LISTS[i] =
      _MPI_safeMalloc(_MPI_PREALLOCATION_SIZE*sizeof(_MPI_REQUEST_OBJECT),
                      "Error malloc a new REQ strip");

    /* Zero out the new strip and select a new "good" req object */
    for(j=0;j<_MPI_PREALLOCATION_SIZE;++j) {
      _MPI_Req_Invalid(_MPI_REQ_LIST_OF_LISTS[i]+j);
    }
    request = _MPI_REQ_LIST_OF_LISTS[i];
  }

  /* ----------------------------------------------- */
  /* Fill in the request object                      */
  /* ----------------------------------------------- */
  request->buffer = buffer;
  request->count = count;
  request->type = datatype;
  request->tag = tag;
  request->comm = comm;
  request->send = send;
  request->valid = _MPI_VALID;

  return request;
}

int _MPI_Req_Invalid (MPI_Request request) {
  _MPI_COVERAGE();

  request->buffer = (void *)0;
  request->count = _MPI_NOT_VALID;
  request->tag = _MPI_NOT_VALID;
  request->comm = MPI_COMM_NULL;
  request->send = _MPI_NOT_VALID;
  request->valid = _MPI_NOT_VALID;

  return MPI_SUCCESS;
}

int _MPI_Check_Request_Array(int count, MPI_Request array[]) {
  int index;
  for(index = 0; index<count; index++) {
    if (array[index] != MPI_REQUEST_NULL) return MPI_SUCCESS;
  }
  return _MPI_NOT_OK;
}

