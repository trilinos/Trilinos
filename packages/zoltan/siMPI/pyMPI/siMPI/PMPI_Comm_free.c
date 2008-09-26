/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********     PMPI_Comm_free.c      ********************/
/******************************************************************/
/* Author : Lisa Alano June 20 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


/*============================================================================*/


int PMPI_Comm_free(MPI_Comm* comm)
{
int index;
  if (_MPI_CHECK_STATUS(comm) == _MPI_OK)
  {
    if (_MPI_Comm_check_legal(*comm, &index)==MPI_SUCCESS)
    {
      _MPI_Comm_Invalid( (int) *comm);
      *comm = _MPI_NOT_VALID;
      _MPI_COMM_COUNT --;
      return MPI_SUCCESS;
    } else {
       return MPI_ERR_COMM;
    }
  } else {
    _MPI_ERR_ROUTINE(MPI_ERR_IN_STATUS, "MPI Status Invalid.");
    return MPI_ERR_IN_STATUS;
  }
}


/*============================================================================*/
