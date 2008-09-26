/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  PMPI_Get_count.c   *****************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Get_count( 
        MPI_Status *status, 
        MPI_Datatype datatype, 
        int *count )
{
  int size_of_one;
  if (!count) return MPI_SUCCESS;

  /* TODO: Verify datatype and return MPI_ERR_TYPE */

  size_of_one = _MPI_calculateSize(1,datatype);
  if (size_of_one == 0) {
    *count = 0;
    return MPI_SUCCESS;
  }
  *count = status->__count / size_of_one ;
  if ( status->__count % size_of_one != 0) {
    *count = MPI_UNDEFINED;
  }
  return MPI_SUCCESS;
}

