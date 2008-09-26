/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************      MPI_Waitany.c       ************************/
/****************************************************************************/
/* Author : Lisa Alano July 12 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int MPI_Waitany( int count, MPI_Request array_of_requests[], int* index, MPI_Status* status)
{
  _MPI_COVERAGE();
  return PMPI_Waitany(count, array_of_requests, index, status);
}
/*==========================================================================*/

