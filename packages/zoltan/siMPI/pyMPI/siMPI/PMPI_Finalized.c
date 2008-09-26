/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************     PMPI_Finalized.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Finalized( int *flag )
{
  int retval;
  retval = _MPI_checkIntP (flag);
  if (retval!=MPI_SUCCESS)
  {
    _MPI_ERR_ROUTINE(MPI_ERR_IN_STATUS, "MPI_FINALIZED: Invalid pointer."); 
    MPI_Abort((MPI_Comm)NULL, MPI_ERR_OTHER);
    return retval;
  }
  if ( (_MPI_INIT_STATUS == _MPI_ENDED) && (_MPI_FINALIZED_FLAG) )
    *flag = _MPI_TRUE;
  else 
    *flag = _MPI_FALSE;
  return MPI_SUCCESS; 
}

