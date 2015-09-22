/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  PMPI_Initialized.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Initialized( int *flag )
{
  int retval;
  retval = _MPI_checkIntP (flag);
  if (retval!=MPI_SUCCESS)
  {
    MPI_Abort((MPI_Comm)0, MPI_ERR_OTHER);
    return retval;
  }
  *flag = _MPI_INITIALIZED_FLAG;
  return MPI_SUCCESS; 
}

