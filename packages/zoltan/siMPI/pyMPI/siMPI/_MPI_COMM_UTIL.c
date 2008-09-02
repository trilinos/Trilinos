/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********     _MPI_COMM_UTIL.c      ********************/
/******************************************************************/
/* Author : Lisa Alano July 12 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

/*============================================================================*/
int _MPI_Comm_Insert (int index)
{
  _MPI_COVERAGE();
  if (_MPI_Comm_check(index)==MPI_SUCCESS) return MPI_ERR_COMM;
  return _MPI_Comm_Insert0 (index, index, MPI_ERRORS_RETURN);
}
/*============================================================================*/
int _MPI_Comm_Insert0 (int index, MPI_Comm set_comm, MPI_Handler_function_pointer function)
{
  _MPI_COVERAGE();
  _MPI_COMM_LIST[index].valid = _MPI_VALID;
  _MPI_COMM_LIST[index].comm = set_comm;
  _MPI_COMM_COUNT ++;

  _MPI_COMM_LIST[index].handler = function;
  _MPI_COMM_LIST[index].group = _MPI_COMM_WORLD_GROUP;

#ifdef _MPI_DEBUG
  ;/*printf("Inserted a communicator at index %d.\n", index);*/
#endif

  return MPI_SUCCESS;
}
/*============================================================================*/
int _MPI_Comm_Invalid (int index)
{
  _MPI_COVERAGE();
  if (index == (int) MPI_COMM_WORLD) {
    _MPI_COVERAGE();
    #ifdef _MPI_DEBUG
    printf("MPI_COMM_WORLD was asked to be freed. Not executed\n");
    #endif
    return MPI_SUCCESS;
  }
  _MPI_COVERAGE();
  _MPI_COMM_LIST[index].valid = _MPI_NOT_VALID;
  _MPI_COMM_LIST[index].comm = _MPI_NOT_VALID;
  _MPI_COMM_LIST[index].handler = 0;
  _MPI_COMM_LIST[index].group = MPI_GROUP_NULL;
  return MPI_SUCCESS;
}
/*============================================================================*/
int _MPI_Comm_check (int index)
{
  _MPI_COVERAGE();
/* KDDKDDKDD  WHY DO WE CHANGE THE VALIDITY HERE?  
  _MPI_COMM_LIST[index].valid = _MPI_NOT_VALID;
*/
  if ( (index == MPI_COMM_NULL)||(index==0) ) {
    _MPI_COVERAGE();
    return MPI_ERR_COMM;
  }
  if (index == MPI_COMM_WORLD) {
    _MPI_COVERAGE();
    index = 0;
  }
  if ( (index >= _MPI_COMM_ARRAY_SIZE)||(index<0 ) ) {
    _MPI_COVERAGE();
    #ifdef _MPI_DEBUG
    _MPI_ERR_ROUTINE(MPI_ERR_COMM, "MPI_Comm is not valid.");
    #endif
    return MPI_ERR_COMM;
  }
  return (_MPI_COMM_LIST[index].valid == _MPI_VALID)?MPI_SUCCESS:MPI_ERR_COMM;
}
/*============================================================================*/
int _MPI_Group_check (MPI_Group group)
{
  if (group == MPI_GROUP_NULL) {
     _MPI_COVERAGE();
     return MPI_ERR_GROUP;
  }
  _MPI_COVERAGE();
  return MPI_SUCCESS;
}
/*============================================================================*/
int _MPI_Comm_check_legal (MPI_Comm comm, int *index)
{
  _MPI_COVERAGE();
  if ( (comm == MPI_COMM_NULL)||(comm==0) ) {
    _MPI_COVERAGE();
    *index = -1;
    return MPI_ERR_COMM;
  }
  if ( comm == MPI_COMM_WORLD ) {
    _MPI_COVERAGE();
    *index = 0;
    return MPI_SUCCESS;
  }
  if ( (comm >= _MPI_COMM_ARRAY_SIZE)||(comm<0 ) ) {
    _MPI_COVERAGE();
#ifdef _MPI_DEBUG
    _MPI_ERR_ROUTINE(MPI_ERR_COMM, "MPI_Comm is not valid.");
#endif
    *index = -1;
    return MPI_ERR_COMM;
  }

  else if (_MPI_COMM_LIST[comm].valid == _MPI_VALID) {
    _MPI_COVERAGE();
    *index = comm;
    return MPI_SUCCESS;
  }
  #ifdef _MPI_DEBUG
  _MPI_ERR_ROUTINE(MPI_ERR_COMM, "MPI_Comm is not valid.");
  #endif
  return MPI_ERR_COMM;
}
/*============================================================================*/

