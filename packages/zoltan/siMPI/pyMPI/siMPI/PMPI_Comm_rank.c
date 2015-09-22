/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:52 $
 *    Revision: 1.2 $
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********    PMPI_Comm_rank.c        ********************/
/******************************************************************/
/* Author : Lisa ALano June 19 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


int PMPI_Comm_rank(MPI_Comm comm, int* rank)
{
int index;
  if (_MPI_CHECK_STATUS(&comm) == _MPI_OK)
  {
    if (_MPI_Comm_check_legal(comm, &index) == MPI_SUCCESS)
    {
      *rank = _MPI_RANK;
      return MPI_SUCCESS;
    }
    _MPI_ERR_ROUTINE(MPI_ERR_COMM, "MPI_COMM_RANK: Null communicator.");
    MPI_Abort (comm, MPI_ERR_COMM); 
    return MPI_ERR_COMM;
  }
  else
  {
    _MPI_ERR_ROUTINE(MPI_ERR_IN_STATUS, "MPI_COMM_RANK: MPI initialization error.");
    MPI_Abort (comm, MPI_ERR_IN_STATUS); 
    return MPI_ERR_ARG;
  }
}

