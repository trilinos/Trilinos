/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:52 $
 *    Revision: 1.2 $
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********        PMPI_Finalize.c     ********************/
/******************************************************************/
/* Author : Lisa Alano June 18 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int PMPI_Finalize (void) {
  int i;  
  _MPI_FINALIZED_FLAG = _MPI_TRUE; 
  
  if ( (_MPI_INIT_STATUS)!=(_MPI_STARTED) ) {
    _MPI_ERR_ROUTINE(MPI_ERR_IN_STATUS, "MPI_FINALIZE: MPI has not been initialized.");    
    MPI_Abort (MPI_COMM_NULL, MPI_ERR_IN_STATUS);
    return _MPI_NOT_OK; 
  } else {
    _MPI_INIT_STATUS = _MPI_ENDED; 
    for (i = 1; i < _MPI_PREALLOCATION_SIZE; i++) {
      _MPI_Comm_Invalid(i);
      _MPI_Data_Invalid(i);
      _MPI_Type_Invalid(i);
    }
    free(_MPI_COMM_LIST); _MPI_COMM_LIST = 0;
    free(_MPI_DATA_BUFF);_MPI_DATA_BUFF = 0;
    free(_MPI_TYPE_LIST);_MPI_TYPE_LIST = 0;
    free(_MPI_OP_LIST);_MPI_OP_LIST = 0;
    for(i=0;i<_MPI_REQ_ARRAY_SIZE;++i) {
      free(_MPI_REQ_LIST_OF_LISTS[i]);
      _MPI_REQ_LIST_OF_LISTS[i] = 0;
    }
    free(_MPI_REQ_LIST_OF_LISTS);_MPI_REQ_LIST_OF_LISTS = 0;
    free(_MPI_REQNULL);_MPI_REQNULL = 0;
    return MPI_SUCCESS;
  } 
}
