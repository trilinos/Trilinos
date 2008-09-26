/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************      PMPI_Barrier.c      ************************/
/****************************************************************************/
/* Author : Lisa Alano July 18 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Barrier ( MPI_Comm comm)
{
  int retval;
  _MPI_COVERAGE();
  retval = _MPI_checkCommunicator (comm);
  if (retval != MPI_SUCCESS)
  {
  _MPI_COVERAGE();
    _MPI_ERR_ROUTINE (retval, "MPI_Barrier: Invalid communicator.");
     MPI_Abort (comm, retval);
  }
  return retval;
}

