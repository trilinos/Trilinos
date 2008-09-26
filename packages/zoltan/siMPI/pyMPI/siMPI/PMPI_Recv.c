/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************       PMPI_Recv.c        ************************/
/****************************************************************************/
/* Author : Lisa Alano June 27 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int PMPI_Recv (void* message, int count, MPI_Datatype datatype, int source,
        int tag, MPI_Comm comm, MPI_Status* status)
{
  int retval, receiver_size, message_size, index;

  retval=_MPI_checks(message, count, datatype, source, tag, comm);
  if (retval == MPI_SUCCESS) {
    receiver_size = _MPI_calculateSize(count, datatype);

    /* ----------------------------------------------- */
    /* Look up the send that matches this receive      */
    /* ----------------------------------------------- */
    index = _MPI_Buff_Find(tag, comm);
    if (index == _MPI_NOT_OK) return MPI_ERR_TAG;
        
    message_size = _MPI_calculateSize(_MPI_DATA_BUFF[index].count,
                                      _MPI_DATA_BUFF[index].type);

    if (status != MPI_STATUS_IGNORE) {
      status->MPI_SOURCE = _MPI_RANK;
      status->MPI_TAG = _MPI_DATA_BUFF[index].tag;
      status->__count = message_size;
    }
    if (message_size > receiver_size) {
      _MPI_ERR_ROUTINE(MPI_ERR_COUNT, "MPI_RECV : Message buffer too small for message");
      if (status != MPI_STATUS_IGNORE) status->MPI_ERROR = MPI_ERR_COUNT;
      _MPI_Data_Invalid(index);
      return MPI_ERR_COUNT;
    }

    memcpy(message, _MPI_DATA_BUFF[index].buffer, message_size);

    if (status != MPI_STATUS_IGNORE) status->MPI_ERROR = MPI_SUCCESS;
    _MPI_Data_Invalid(index);
    return MPI_SUCCESS;
  }

  return retval;  
}

/*==========================================================================*/
