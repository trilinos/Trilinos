/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  PMPI_Group_compare.c    ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Group_compare ( MPI_Group group1, MPI_Group group2, int *result )
{
  _MPI_COVERAGE();
  *result = MPI_UNEQUAL;
  if ( _MPI_Group_check(group1) == MPI_SUCCESS &&
       _MPI_Group_check(group2) == MPI_SUCCESS )
  {
    _MPI_COVERAGE();
    if ( group1 == group2 )
      *result = MPI_IDENT;
    return MPI_SUCCESS;
  }
  return MPI_ERR_GROUP;
}

