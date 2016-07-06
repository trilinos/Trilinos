/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.3 $
 ****************************************************************************/
/***********************************************************************************************/
/* FILE  **************************       PMPI_Send.c        ***********************************/
/***********************************************************************************************/
/* Author : Lisa Alano June 25 2002                                                            */
/* Copyright (c) 2002 University of California Regents                                         */
/***********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*=============================================================================================*/
int PMPI_Send (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
  int size, retval, sendType, index;
  MPI_Aint position, copiedPointer;
  _MPI_TYPE_INFO *info;
  char* p;

  _MPI_COVERAGE();

  _MPI_CHECK_STATUS(&comm);
  retval = _MPI_checks(message, count, datatype, dest, tag, comm);
  if (retval == MPI_SUCCESS) {
    sendType = _MPI_checkSendType(datatype);
    switch (sendType) {
    case _MPI_DEFAULT:
      {
        size = _MPI_calculateSize(count, datatype);  
        p = (char *)_MPI_safeMalloc(size, "Error with malloc for send buffer."); 
        p = memcpy(p, message, size);
        retval =_MPI_Buff_Insert(p, count, datatype, tag, comm);
        return retval;
      }
    case _MPI_CONTIG:
      {
        sendType = _MPI_FindType(datatype);
        size = _MPI_TYPE_LIST[sendType].extent;
        p = (char *)_MPI_safeMalloc(size, "Error with malloc for send buffer."); 
        p = memcpy(p, message, size);
        retval =_MPI_Buff_Insert(p, count, datatype, tag, comm);
        return retval;
      }
    case _MPI_INDEXED:
      {
        sendType = _MPI_FindType(datatype);
        size = _MPI_TYPE_LIST[sendType].extent;
        p = (char *)_MPI_safeMalloc(size, "Error with malloc for send buffer."); 

        /* ================================================== */
        /* Determine the correct parts to save to the buffers */
        info = _MPI_TYPE_LIST[sendType].info;
        position = (MPI_Aint) 0;
        copiedPointer = (MPI_Aint) 0;
        for (index = 0; index < info->count; index++)
          {
            position = info->stride[index]*sizeof(info->types[0]);
            p = memcpy(p+copiedPointer, ((char*)message)+position, info->blocklen[index]*sizeof(info->types[0])); 
            copiedPointer += info->blocklen[index]*sizeof(info->types[0]);
          }
        retval =_MPI_Buff_Insert(p, count, datatype, tag, comm);
        return retval;
      }
    case _MPI_VECTOR:
      {
        sendType = _MPI_FindType(datatype);
        size = _MPI_TYPE_LIST[sendType].extent;
        p = (char *)_MPI_safeMalloc(size, "Error with malloc for send buffer."); 
        /* =================================== */
        /* Figure out the correct ones to pass */
        retval =_MPI_Buff_Insert(p, count, datatype, tag, comm);
        return retval;
      }
    case _MPI_STRUCT:
      {
        sendType = _MPI_FindType(datatype);
        size = _MPI_TYPE_LIST[sendType].extent;
        p = (char *)_MPI_safeMalloc(size, "Error with malloc for send buffer."); 
        /* =================================== */
        /* Figure out the correct ones to pass */
        retval =_MPI_Buff_Insert(p, count, datatype, tag, comm);
        return retval;
      }
    default:
      {
        fprintf(stderr,"SEND: mpi_Hindexed or mpi_Hvector not implemented\n");
        MPI_Abort (comm, _MPI_NOT_OK);
      }
    }
  } else { 
    _MPI_ERR_ROUTINE (retval, "MPI_SEND / MPI_ISEND: argument error");
    MPI_Abort (comm, retval);
  }

  _MPI_COVERAGE();
  return _MPI_NOT_OK;
}
/*=============================================================================================*/
