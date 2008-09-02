/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  PMPI_Type_vector.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano August 7 2002                                        */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Type_vector( 
        int count, 
        int blocklen, 
        int stride, 
        MPI_Datatype old_type, 
        MPI_Datatype *newtype )
{
  int index;
  MPI_Aint size;
  _MPI_TYPE_DES* currType;
  _MPI_TYPE_DES* prevType;
  _MPI_TYPE_INFO* info;
  index = _MPI_Find_free();

  size = (MPI_Aint)count*_MPI_calculateSize(blocklen, old_type);
  *newtype = _MPI_TYPE_LIST[index].id = _MPI_TYPE_COUNT+_MPI_TYPE_OFFSET;
  _MPI_TYPE_COUNT++;
  _MPI_TYPE_LIST[index].size = size;
  _MPI_TYPE_LIST[index].extent = size;
  _MPI_TYPE_LIST[index].ub = size;
  _MPI_TYPE_LIST[index].lb = (MPI_Aint) 0;
  _MPI_TYPE_LIST[index].sendType = _MPI_VECTOR;
  _MPI_TYPE_LIST[index].next = 0;
  _MPI_TYPE_LIST[index].info = (_MPI_TYPE_INFO *) _MPI_safeMalloc(sizeof(_MPI_TYPE_INFO), "MPI_TYPE_INDEXED: Error with malloc");

  /* ========================================================== */
  /* Save information for packing and unpacking to/from buffers */
  info = _MPI_TYPE_LIST[index].info;
  info->count = count;
  info->blocklen = (int *) _MPI_safeMalloc(sizeof(int), "MPI_TYPE_INDEXED: Error with malloc");;
  *(info->blocklen) = blocklen;
  info->stride = (int *) _MPI_safeMalloc(sizeof(int), "MPI_TYPE_INDEXED: Error with malloc");;
  *(info->stride) = stride;

  prevType = &_MPI_TYPE_LIST[index];
  size = 0;
  for (index=0; index<count; index++)
  {
    currType = (_MPI_TYPE_DES *) _MPI_safeMalloc(sizeof(_MPI_TYPE_DES), "MPI_TYPE_VECTOR: Error with malloc.");
    prevType->next = currType;
    currType->id = old_type;
    currType->size = blocklen*_MPI_getSize(old_type);
    size += blocklen*_MPI_getSize(old_type);
    currType->extent = size;
    currType->next = 0;
    prevType = currType;
  }

  return MPI_SUCCESS;
}

