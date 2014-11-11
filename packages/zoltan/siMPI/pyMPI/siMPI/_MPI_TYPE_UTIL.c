/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
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
  /* KDDKDD 11/4/14:  changed loop max from _MPI_TYPE_COUNT to
   * _MPI_TYPE_ARRAY_SIZE; _MPI_TYPE_LIST entries are freed and 
   * reused, so _MPI_TYPE_ARRAY_SIZE accurately gives the max 
   * entries to search */
  for (index = 0; index < _MPI_TYPE_ARRAY_SIZE; index++)
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
  int ret;
  _MPI_COVERAGE();
  /* KDDKDD 11/4/14:  changed loop max from _MPI_TYPE_COUNT to
   * _MPI_TYPE_ARRAY_SIZE; _MPI_TYPE_LIST entries are freed and 
   * reused, so _MPI_TYPE_ARRAY_SIZE accurately gives the max 
   * entries to search */
  for (i=0; i<_MPI_TYPE_ARRAY_SIZE; i++) {
    _MPI_COVERAGE();
    if (_MPI_TYPE_LIST[i].id == _MPI_NOT_VALID) return i;
  }
  _MPI_TYPE_LIST = (_MPI_TYPE_DES*) _MPI_safeRealloc
    (_MPI_TYPE_LIST,
     (_MPI_TYPE_ARRAY_SIZE+_MPI_PREALLOCATION_SIZE)*sizeof(_MPI_TYPE_DES),
     "Error in _MPI_Find_free for reallocation");
  ret = _MPI_TYPE_ARRAY_SIZE;
  _MPI_TYPE_ARRAY_SIZE += _MPI_PREALLOCATION_SIZE;

  /* KDDKDD 11/4/14:  changed return value from_MPI_TYPE_COUNT to
   * the previous _MPI_TYPE_ARRAY_SIZE */
  return ret;
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



