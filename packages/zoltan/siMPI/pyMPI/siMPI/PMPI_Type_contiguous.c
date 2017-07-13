/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  PMPI_Type_contiguous.c  ************************/
/****************************************************************************/
/* Author : Lisa Alano August 2 2002                                        */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Type_contiguous( 
        int count,
        MPI_Datatype old_type,
        MPI_Datatype *newtype)
{
  int index;
  MPI_Aint size;
  _MPI_TYPE_DES* currType;
  _MPI_TYPE_DES* prevType;
  index = _MPI_Find_free();
  size = (MPI_Aint)_MPI_calculateSize(count, old_type);

  *newtype = _MPI_TYPE_LIST[index].id = _MPI_TYPE_COUNT+_MPI_TYPE_OFFSET;
  _MPI_TYPE_COUNT++;
  _MPI_TYPE_LIST[index].extent = size;
  _MPI_TYPE_LIST[index].sendType = _MPI_CONTIG;
  _MPI_TYPE_LIST[index].next = 0;
  _MPI_TYPE_LIST[index].info = 0;
  prevType = &_MPI_TYPE_LIST[index];

  size = 0;
  /* ================================================================= */
  /* DONT THINK WE ACUTALLY NEED THIS ================================ */
  /* WE CAN JUST USE THE INFO STRUCT - same with vector and contiguous */
  for (index=0; index<count; index++)
  {
    currType = (_MPI_TYPE_DES *) _MPI_safeMalloc(sizeof(_MPI_TYPE_DES), "MPI_TYPE_CONTIGUOUS: Error with malloc.");
    prevType->next = currType;
    currType->id = old_type;
    size += _MPI_getSize(old_type);
    currType->extent = size;
    currType->next = 0;
    prevType = currType;
  }
  
  return MPI_SUCCESS;
}

