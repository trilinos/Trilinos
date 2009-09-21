/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
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

