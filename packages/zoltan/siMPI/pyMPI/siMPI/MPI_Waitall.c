/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************      MPI_Waitall.c       ************************/
/****************************************************************************/
/* Author : Lisa Alano July 12 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int MPI_Waitall(
        int count,
        MPI_Request array_of_requests[],
        MPI_Status array_of_statuses[] )
{
  _MPI_COVERAGE();
  return PMPI_Waitall(count, array_of_requests, array_of_statuses);
}
/*==========================================================================*/

