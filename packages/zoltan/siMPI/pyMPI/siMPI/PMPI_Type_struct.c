/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
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
  int nfields;
  MPI_Aint size;
  _MPI_TYPE_DES* currType;
  _MPI_TYPE_DES* prevType;
  _MPI_TYPE_INFO* info;

  /* KDD 6/2/16 Last entry of old_types must be MPI_UB to 
     get the padding correct. */
  nfields = count - 1;
  if (old_types[nfields] != MPI_UB) {
    _MPI_ERR_ROUTINE (MPI_ERR_TYPE, 
                      "MPI_Type_struct:  Must terminate type list with MPI_UB");
    MPI_Abort(MPI_COMM_WORLD, MPI_ERR_TYPE);
    return _MPI_NOT_OK;
  }

  index = _MPI_Find_free();
  /* ===================================== */
  /* Calculate total size of parts to copy */
  /* KDD 6/2/16  This size calc doesn't account for padding in the struct. 
     KDD 6/2/16  Instead, we'll compute size based on the indices.
     KDD 6/2/16  We assume that the indices are terminated with MPI_UB,
     KDD 6/2/16  as recommended in MPI-1 documentation, so that indices
     KDD 6/2/16  is an offsets array in CRS sense.
  for (i=0; i<nfields; i++)
  {
     size += (MPI_Aint)_MPI_calculateSize(blocklens[i], old_types[i]);
  }
  */
  size = indices[nfields] - indices[0];

  /* ============================== */
  /* Give new id/unique of datatype */
  *newtype = _MPI_TYPE_LIST[index].id = _MPI_TYPE_COUNT+_MPI_TYPE_OFFSET;
  _MPI_TYPE_COUNT++;

  /* ====================== */
  /* Save Query information */
  /* KDD 6/2/16  To account for padding in the structure, the 
   * KDD 6/2/16  extent and ub should be related to the calculated size 
  _MPI_TYPE_LIST[index].extent = indices[count-1]+(blocklens[count-1]*_MPI_getSize(old_types[count-1]));
  _MPI_TYPE_LIST[index].ub = indices[count-1]+(blocklens[count-1]*_MPI_getSize(old_types[count-1]));
  _MPI_TYPE_LIST[index].lb = indices[0];
  */
  _MPI_TYPE_LIST[index].extent = size;
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
  for (i=0; i<nfields; i++)
  {
    currType = (_MPI_TYPE_DES *) _MPI_safeMalloc(sizeof(_MPI_TYPE_DES), "MPI_TYPE_STRUCT: Error with malloc");
    prevType->next = currType;
    currType->id = old_types[i];
    /* KDD 6/2/16 use the actual extent provided by the indices 
    currType->extent = indices[i]+currType->size;
    */
    currType->extent = indices[i+1]-indices[i];
    currType->next = 0;
    prevType = currType;
  }
  /* =============================================== */
  /* Add the MPI_UB at the end of the structure list */
  currType = (_MPI_TYPE_DES *) _MPI_safeMalloc(sizeof(_MPI_TYPE_DES), "MPI_TYPE_STRUCT: Error with malloc.");
  prevType->next = currType;
  currType->id = MPI_UB;
  /* KDD 6/2/16  Not sure why MPI_UB should have a size or extent.
  currType->size = _MPI_TYPE_LIST[index].size;
  currType->extent = _MPI_TYPE_LIST[index].extent;
  */
  currType->extent = 0;
  currType->next = 0;

  return MPI_SUCCESS;  
}

