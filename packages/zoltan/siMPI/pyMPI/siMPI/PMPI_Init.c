/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.2 $
 ****************************************************************************/
/***********************************************************************************************/
/* FILE  ***********************        PMPI_Init.c         ************************************/
/***********************************************************************************************/
/* Author : Lisa Alano June 18 2002                                                            */
/* Copyright (c) 2002 University of California Regents                                         */
/***********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


int PMPI_Init ( int *argc, char **argv[]) {
  int i;
  int retval = MPI_ERR_IN_STATUS;
  _MPI_COVERAGE();
  _MPI_INITIALIZED_FLAG = _MPI_TRUE;
  
  /* ----------------------------------------------- */
  /*  Check for the current status of MPI            */
  /* ----------------------------------------------- */
  if ( (_MPI_INIT_STATUS == _MPI_ENDED) || (_MPI_INIT_STATUS == _MPI_STARTED) ) {   
    if (_MPI_INIT_STATUS == _MPI_ENDED) {
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_IN_STATUS, "MPI was already finalized");
    } else {
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_IN_STATUS, "MPI was already initialized");	
    }
    MPI_ERRORS_ARE_FATAL((MPI_Comm)0, &retval);
    return MPI_ERR_OTHER;
  } 

  /* ----------------------------------------------- */
  /* If the status is fine, initialize the internal  */
  /* data structures                                 */ 
  /* ----------------------------------------------- */
  _MPI_COVERAGE();
  _MPI_INIT_STATUS = _MPI_STARTED; 

  /* ---------- */ 
  /* Allocation */
  _MPI_COMM_ARRAY_SIZE = _MPI_PREALLOCATION_SIZE;
  _MPI_DATA_ARRAY_SIZE = _MPI_PREALLOCATION_SIZE;
  _MPI_TYPE_ARRAY_SIZE = _MPI_PREALLOCATION_SIZE;
  _MPI_OP_ARRAY_SIZE = _MPI_PREALLOCATION_SIZE;
  _MPI_REQ_ARRAY_SIZE = _MPI_PREALLOCATION_SIZE;


  _MPI_COMM_LIST = (_MPI_COMM_IMPL *) _MPI_safeMalloc (_MPI_COMM_ARRAY_SIZE*sizeof(_MPI_COMM_IMPL), "Error with malloc of COMM_LIST");
  _MPI_DATA_BUFF = (_MPI_DATA_ENTRY *) _MPI_safeMalloc (_MPI_DATA_ARRAY_SIZE*sizeof(_MPI_DATA_ENTRY), "Error with malloc of DATA_BUFF");
  _MPI_TYPE_LIST = (_MPI_TYPE_DES *) _MPI_safeMalloc (_MPI_TYPE_ARRAY_SIZE*sizeof(_MPI_TYPE_DES), "Error with malloc of TYPE_LIST");
  _MPI_OP_LIST = (_MPI_OP_TYPE *) _MPI_safeMalloc (_MPI_OP_ARRAY_SIZE*sizeof(_MPI_OP_TYPE), "Error with malloc of OP_LIST");
  _MPI_REQ_ARRAY_SIZE = 0;
  _MPI_REQ_LIST_OF_LISTS = 0;

  /* ----------------------------------------------- */
  /* Communicators are not set up                    */
  /* ----------------------------------------------- */
  for (i=1; i<_MPI_COMM_ARRAY_SIZE; i++) {
    _MPI_COMM_LIST[i].valid = _MPI_NOT_VALID;
  }
  for (i=1; i<_MPI_DATA_ARRAY_SIZE; i++) {
    _MPI_DATA_BUFF[i].valid = _MPI_NOT_VALID;
  }
  for (i=1; i<_MPI_OP_ARRAY_SIZE; i++) {
    _MPI_OP_LIST[i].valid = _MPI_NOT_VALID;
  }
  for (i=1; i<_MPI_TYPE_ARRAY_SIZE; i++) {
    _MPI_TYPE_LIST[i].id = _MPI_NOT_VALID;
    _MPI_TYPE_LIST[i].info = 0;
    _MPI_TYPE_LIST[i].next = 0;
  }

  _MPI_COMM_COUNT = 0;
  _MPI_DATA_BUFF_COUNT = 0;
  _MPI_TYPE_COUNT = 0;
  _MPI_OP_COUNT = 0;
  _MPI_REQ_COUNT = 0;
    
  /* ------------------------- */
  /* Set entries all to "null" */ 
  for (i=0; i<_MPI_PREALLOCATION_SIZE; i++) {
    _MPI_COVERAGE();
    _MPI_Data_Invalid(i);
    _MPI_Comm_Invalid(i);
    _MPI_Type_Invalid(i);
  }                                                           /* --------------------------------------------------- */
  _MPI_Comm_Insert0(0, MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL); /* This inserts MPI_COMM_WORLD as the 1st communicator */
  return MPI_SUCCESS;

}
/*=============================================================================================*/
