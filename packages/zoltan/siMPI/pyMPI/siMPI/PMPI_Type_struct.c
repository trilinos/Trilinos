/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  PMPI_Type_struct.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano August 8 2002                                        */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"
#include <string.h>

int PMPI_Type_struct( 
        int count, 
        int blocklens[], 
        MPI_Aint indices[], 
        MPI_Datatype old_types[], 
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
    size += (MPI_Aint)_MPI_calculateSize(blocklens[i], old_types[i]);
  }

  /* ============================== */
  /* Give new id/unique of datatype */
  *newtype = _MPI_TYPE_LIST[index].id = _MPI_TYPE_COUNT+_MPI_TYPE_OFFSET;
  _MPI_TYPE_COUNT++;

  /* ====================== */
  /* Save Query information */
  _MPI_TYPE_LIST[index].size = size;
  _MPI_TYPE_LIST[index].extent = indices[count-1]+(blocklens[count-1]*_MPI_getSize(old_types[count-1]));
  _MPI_TYPE_LIST[index].ub = indices[count-1]+(blocklens[count-1]*_MPI_getSize(old_types[count-1]));
  _MPI_TYPE_LIST[index].lb = indices[0];
  _MPI_TYPE_LIST[index].sendType = _MPI_STRUCT;
  _MPI_TYPE_LIST[index].next = 0;
  _MPI_TYPE_LIST[index].info = (_MPI_TYPE_INFO *) _MPI_safeMalloc(sizeof(_MPI_TYPE_INFO), "MPI_TYPE_INDEXED: Error with malloc");

  /* ========================================================== */
  /* Save information for packing and unpacking to/from buffers */
  info = _MPI_TYPE_LIST[index].info;
  info->count = count;
  info->blocklen = (int *) _MPI_safeMalloc(sizeof(int)*count, "MPI_TYPE_STRUCT: Error with malloc");;
  info->blocklen = memcpy(info->blocklen, blocklens, sizeof(int)*count);
  info->stride = (int *) _MPI_safeMalloc(sizeof(int)*count, "MPI_TYPE_STRUCT: Error with malloc");;
  info->stride = memcpy(info->stride, indices, sizeof(int)*count);
  info->types = (MPI_Datatype *) _MPI_safeMalloc(sizeof(MPI_Datatype)*count, "MPI_TYPE_STRUCT: Error with malloc");;
  info->types = memcpy(info->types, old_types, sizeof(int)*count);

  /* ================================ */
  /* Create linked list of structures */
  prevType = &_MPI_TYPE_LIST[index];
  for (index=0; index<count; index++)
  {
    currType = (_MPI_TYPE_DES *) _MPI_safeMalloc(sizeof(_MPI_TYPE_DES), "MPI_TYPE_INDEXED: Error with malloc");
    prevType->next = currType;
    currType->id = old_types[index];
    currType->size = blocklens[index]*_MPI_getSize(old_types[index]);
    currType->extent = indices[index]+currType->size;
    currType->next = 0;
    prevType = currType;
  }
  /* =============================================== */
  /* Add the MPI_UB at the end of the structure list */
  currType = (_MPI_TYPE_DES *) _MPI_safeMalloc(sizeof(_MPI_TYPE_DES), "MPI_TYPE_CONTIGUOUS: Error with malloc.");
  prevType->next = currType;
  currType->id = MPI_UB;
  currType->size = _MPI_TYPE_LIST[index].size;
  currType->extent = _MPI_TYPE_LIST[index].extent;
  currType->next = 0;

  return MPI_SUCCESS;  
}

