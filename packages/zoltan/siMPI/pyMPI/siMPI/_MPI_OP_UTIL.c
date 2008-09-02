/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/* FILE  ******************     _MPI_OP_UTIL.c       ************************/
/****************************************************************************/
/* Author : Lisa Alano July 1 2002                                          */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
void _MPI_Default_Op ( void* a, void* b, int* len, MPI_Datatype* type )
{
  int size;
  _MPI_COVERAGE();
  size = _MPI_calculateSize(*len, *type); 
  memcpy(b, a, size);
}

/*==========================================================================*/
void _MPI_Ops_init(int i)
{
  _MPI_COVERAGE();
  _MPI_OP_LIST[i].valid = _MPI_NOT_VALID;
  _MPI_OP_LIST[i].function = _MPI_Default_Op;
  _MPI_OP_LIST[i].commute = 1;
  _MPI_OP_COUNT++;
}

/*==========================================================================*/
int _MPI_Op_insert (MPI_User_function* function, int commute, MPI_Op* id) {

  int i,j;

  _MPI_COVERAGE();

  if ( function == 0 ) {
    _MPI_COVERAGE();
    return MPI_ERR_ARG;
  }
  if ( id == 0 ) {
    _MPI_COVERAGE();
    return MPI_ERR_ARG;
  }

  for (i=0; i<_MPI_OP_ARRAY_SIZE; i++) {
    _MPI_COVERAGE();
    if (_MPI_OP_LIST[i].valid != _MPI_VALID) break;
  }

  if (i>=_MPI_COMM_ARRAY_SIZE) {
    _MPI_COVERAGE();
    /* ----------------------------------------------- */
    /* This realloc doesn't orphan ops because MPI_Op  */
    /* is just an int (index into list).               */
    /* ----------------------------------------------- */
    _MPI_OP_LIST = (_MPI_OP_TYPE*) _MPI_safeRealloc
      (_MPI_OP_LIST,
       (_MPI_OP_ARRAY_SIZE+_MPI_PREALLOCATION_SIZE)*sizeof(_MPI_OP_TYPE),
       "Error in _MPI_Op_Insert for reallocation"); 

    for(j=0;j<_MPI_PREALLOCATION_SIZE;j++) {
      _MPI_COVERAGE();
      _MPI_Ops_init(_MPI_OP_ARRAY_SIZE+j);
    }
    _MPI_OP_ARRAY_SIZE+=_MPI_PREALLOCATION_SIZE;
  }

  _MPI_OP_LIST[i].valid = _MPI_VALID;
  _MPI_OP_LIST[i].function = function;
  _MPI_OP_LIST[i].commute = commute;

  *id = i+_MPI_OP_OFFSET; 
  _MPI_OP_COUNT++;
  _MPI_COVERAGE();
  return MPI_SUCCESS;
} 
/*==========================================================================*/
int _MPI_checkOp(MPI_Op op)
{
  if ( (op >= MPI_MAX) && (op<_MPI_OP_OFFSET) ) {
    _MPI_COVERAGE();
    return MPI_SUCCESS;
  }

  else if( (op>=_MPI_OP_OFFSET)&&((op-_MPI_OP_OFFSET) < _MPI_OP_ARRAY_SIZE) ) {
    _MPI_COVERAGE();
    if (_MPI_OP_LIST[op-_MPI_OP_OFFSET].valid == _MPI_VALID) {
      _MPI_COVERAGE();
      return MPI_SUCCESS;
    }
  }

  printf("%s:%d: FAILS\n",__FILE__,__LINE__); /* HACK: */
  /* Else */
  _MPI_COVERAGE();
  return MPI_ERR_OP;
}
/*==========================================================================*/
int _MPI_Function_check(MPI_User_function* function)
{
  _MPI_COVERAGE();
  if (function == 0) return _MPI_NOT_OK;
  _MPI_COVERAGE();
  return MPI_SUCCESS;
}
/*==========================================================================*/
int _MPI_Commute_check(int commute)
{
  _MPI_COVERAGE();
  if (commute == 0) return MPI_ERR_OTHER;
  _MPI_COVERAGE();
  return MPI_SUCCESS;
}
/*==========================================================================*/
int _MPI_Op_invalid(int index)
{
  _MPI_COVERAGE();
  if ((index >=0 ) && (index < _MPI_OP_ARRAY_SIZE) ) {
    _MPI_COVERAGE();
    _MPI_OP_LIST[index].valid = _MPI_NOT_VALID;
    _MPI_OP_LIST[index].function = 0;
    _MPI_OP_LIST[index].commute = 0;
    return MPI_SUCCESS;
  }
  _MPI_COVERAGE();
  return MPI_ERR_OP;
}
/*==========================================================================*/
int _MPI_check_overlap(void* sendbuf, void* recvbuf, int size)
{
  _MPI_COVERAGE();
  if(sendbuf == recvbuf) { 
    _MPI_COVERAGE();
    return _MPI_NOT_OK;
  }
  if( (((char*)sendbuf+size)>(char*)recvbuf)&&(sendbuf<recvbuf) ) {
    _MPI_COVERAGE();
    return _MPI_NOT_OK;
  }
  if( (((char*)recvbuf+size)>(char*)sendbuf)&&(recvbuf<sendbuf) ) {
    _MPI_COVERAGE();
    return _MPI_NOT_OK;
  }
  _MPI_COVERAGE();
  return MPI_SUCCESS;
}
/*==========================================================================*/

