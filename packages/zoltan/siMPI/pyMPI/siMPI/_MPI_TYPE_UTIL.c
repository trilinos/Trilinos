/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************   _MPI_TYPE_UTIL.c    ************************/
/****************************************************************************/
/* Author : Lisa Alano August 2 2002                                        */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/*==========================================================================*/
int _MPI_BasicType (MPI_Datatype datatype)
{
  _MPI_COVERAGE();
  if ( (datatype >= _MPI_TYPE_START)&&(datatype <= _MPI_TYPE_END) )
  {
  _MPI_COVERAGE();
    return MPI_SUCCESS;
  }
  return _MPI_FALSE;
}
/*==========================================================================*/
int _MPI_FindType (MPI_Datatype datatype)
{
  int index;
  _MPI_COVERAGE();
  for (index = 0; index < _MPI_TYPE_COUNT; index++)
  {
  _MPI_COVERAGE();
    if (_MPI_TYPE_LIST[index].id == datatype)
    {
  _MPI_COVERAGE();
      return index;
    }
  } 
  return _MPI_NOT_OK;
}
/*==========================================================================*/
int _MPI_Free_datatype (MPI_Datatype datatype)
{
  int index;
  _MPI_TYPE_DES* child;
  _MPI_COVERAGE();
  index = _MPI_FindType(datatype);
  if (index == _MPI_NOT_OK) {
  _MPI_COVERAGE();
    return _MPI_NOT_OK;
  }
#ifdef KDD_REMOVED_DEBUG
  printf("%s:%d: BOGUS value??? %d %d 0x%x\n",__FILE__,__LINE__,datatype, index, (int)_MPI_TYPE_LIST[index].next);
#endif
  child = _MPI_TYPE_LIST[index].next;
  if ( child ) _MPI_deleteAll(child);
  _MPI_safeFree(_MPI_TYPE_LIST[index].info,"type info");
  /* 
  _MPI_safeFree(_MPI_TYPE_LIST[index].next,"type next");
  */

  _MPI_TYPE_LIST[index].id = _MPI_NOT_VALID;
  _MPI_TYPE_LIST[index].info = 0;
  _MPI_TYPE_LIST[index].next = 0;

  return MPI_SUCCESS;
}

/*==========================================================================*/
void _MPI_deleteAll (_MPI_TYPE_DES* parent)
{
  _MPI_COVERAGE();
  if (parent->next != NULL) {
  _MPI_COVERAGE();
    _MPI_deleteAll (parent->next);
    parent->next = 0;
  }
  _MPI_safeFree (parent,"Parent");
}
/*==========================================================================*/
int _MPI_Find_free (void) {
  int i;
  _MPI_COVERAGE();
  for (i=0; i<_MPI_TYPE_COUNT; i++) {
    _MPI_COVERAGE();
    if (_MPI_TYPE_LIST[i].id == _MPI_NOT_VALID) return i;
  }
  _MPI_TYPE_LIST = (_MPI_TYPE_DES*) _MPI_safeRealloc
    (_MPI_TYPE_LIST,
     (_MPI_TYPE_ARRAY_SIZE+_MPI_PREALLOCATION_SIZE)*sizeof(_MPI_TYPE_DES),
     "Error in _MPI_Find_free for reallocation");
  _MPI_TYPE_ARRAY_SIZE += _MPI_PREALLOCATION_SIZE;

  return _MPI_TYPE_COUNT;
}
/*==========================================================================*/
/* int _MPI_calculatePadding (MPI_Datatype firstType, MPI_Aint currentSize) */
/*==========================================================================*/
int _MPI_checkSendType (MPI_Datatype type) {
  int index;
  _MPI_COVERAGE();
  index = _MPI_FindType(type);
  if (index == _MPI_NOT_OK) {
    _MPI_COVERAGE();
    return _MPI_DEFAULT;
  }
  return _MPI_TYPE_LIST[index].sendType;
}
/*==========================================================================*/



