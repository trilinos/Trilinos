/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/******************************************************************************/
/* FILE  ****************    PMPI_Comm_create.c     ***************************/
/******************************************************************************/
/* Author : Lisa Alano June 19 2002                                           */
/* Copyright (c) 2002 University of California Regents                        */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


/*============================================================================*/
int PMPI_Comm_create(MPI_Comm comm, MPI_Group new_group, MPI_Comm* new_comm) {
  int i,g,c,j;
  _MPI_COVERAGE();
  if (_MPI_CHECK_STATUS(new_comm) == _MPI_OK) {
    _MPI_COVERAGE();
    g = _MPI_Group_check(new_group);
    c = _MPI_Comm_check(comm);

    if ( (g==MPI_SUCCESS)&&(c==MPI_SUCCESS) ) {
      _MPI_COVERAGE();
      for (i=1; i<_MPI_COMM_ARRAY_SIZE; i++) {
        _MPI_COVERAGE();
        if (_MPI_COMM_LIST[i].valid != _MPI_VALID)
          break;
      }
      if (i>=_MPI_COMM_ARRAY_SIZE) {
        _MPI_COVERAGE();
        _MPI_COMM_LIST = (_MPI_COMM_IMPL*) _MPI_safeRealloc(_MPI_COMM_LIST, (_MPI_COMM_ARRAY_SIZE+_MPI_PREALLOCATION_SIZE)*sizeof(_MPI_COMM_IMPL), "Error in MPI_Comm_create reallocation");
        for(j=0;j<_MPI_PREALLOCATION_SIZE;++j) {
          _MPI_COVERAGE();
          _MPI_COMM_LIST[_MPI_COMM_ARRAY_SIZE+j].valid = _MPI_NOT_VALID;
        }
        _MPI_COMM_ARRAY_SIZE+=_MPI_PREALLOCATION_SIZE;
      }
      _MPI_Comm_Insert(i);  
      *new_comm = _MPI_COMM_LIST[i].comm;
      return MPI_SUCCESS;
    } else {
      _MPI_COVERAGE();
      if (g ==MPI_SUCCESS)
        return MPI_ERR_COMM;
      return MPI_ERR_GROUP;
    }
  }
  return MPI_ERR_INTERN;
}
/*============================================================================*/

