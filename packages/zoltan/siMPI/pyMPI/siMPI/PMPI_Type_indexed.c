/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ********************  PMPI_Type_indexed.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano August 7 2002                                        */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"
#include <string.h>

int PMPI_Type_indexed( 
        int count, 
        int blocklens[], 
        int indices[], 
        MPI_Datatype old_type, 
        MPI_Datatype *newtype )
{
  int index, i;
  MPI_Aint size;
  _MPI_TYPE_DES* currType;
  _MPI_TYPE_DES* prevType;
  _MPI_TYPE_INFO* info;
  index = _MPI_Find_free();
  /* ===================================== */
  /* Calculate total size of parts to copy */
  size = 0;
  for (i=0; i<count; i++)
  {
    size += (MPI_Aint)_MPI_calculateSize(blocklens[i], old_type);
  }

  /* ============================== */
  /* Give new id/unique of datatype */
  *newtype = _MPI_TYPE_LIST[index].id = _MPI_TYPE_COUNT+_MPI_TYPE_OFFSET;
  _MPI_TYPE_COUNT++;

  /* ====================== */
  /* Save Query information */
  _MPI_TYPE_LIST[index].size = size;
  _MPI_TYPE_LIST[index].extent = size;
  _MPI_TYPE_LIST[index].ub = size;
  _MPI_TYPE_LIST[index].lb = (MPI_Aint) 0;
  _MPI_TYPE_LIST[index].sendType = _MPI_INDEXED;
  _MPI_TYPE_LIST[index].next = 0;
  _MPI_TYPE_LIST[index].info = (_MPI_TYPE_INFO *) _MPI_safeMalloc(sizeof(_MPI_TYPE_INFO), "MPI_TYPE_INDEXED: Error with malloc");

  /* ========================================================== */
  /* Save information for packing and unpacking to/from buffers */
  info = _MPI_TYPE_LIST[index].info;
  info->count = count;
  info->blocklen = (int *) _MPI_safeMalloc(sizeof(int)*count, "MPI_TYPE_INDEXED: Error with malloc");;
  info->blocklen = memcpy(info->blocklen, blocklens, sizeof(int)*count); 
  info->stride = (int *) _MPI_safeMalloc(sizeof(int)*count, "MPI_TYPE_INDEXED: Error with malloc");;
  info->stride = memcpy(info->stride, indices, sizeof(int)*count); 
  info->types = (int *) _MPI_safeMalloc(sizeof(MPI_Datatype), "MPI_TYPE_INDEXED: Error with malloc");;
  info->types[0] = old_type;  

  /* ================================ */
  /* Create linked list of structures */
  prevType = &_MPI_TYPE_LIST[index];
  size = indices[0]*_MPI_getSize(old_type);
  for (index=0; index<count; index++)
  {
    currType = (_MPI_TYPE_DES *) _MPI_safeMalloc(sizeof(_MPI_TYPE_DES), "MPI_TYPE_INDEXED: Error with malloc");
    prevType->next = currType;
    currType->id = old_type;
    currType->size = blocklens[index]*_MPI_getSize(old_type);
    size += currType->size;
    currType->extent = size;
    currType->next = 0;
    prevType = currType;
  }

  return MPI_SUCCESS;
}

