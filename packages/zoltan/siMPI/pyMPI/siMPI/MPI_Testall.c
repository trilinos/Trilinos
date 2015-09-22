/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *************************   MPI_Testall.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Testall( 
        int count, 
        MPI_Request array_of_requests[], 
        int *flag, 
        MPI_Status array_of_statuses[] )
{
  _MPI_COVERAGE();
  return PMPI_Testall (count, array_of_requests, flag, array_of_statuses);
}

