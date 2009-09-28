/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  MPI_Testany.c   *****************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Testany( 
        int count, 
        MPI_Request array_of_requests[], 
        int *index, int *flag, 
        MPI_Status *status )
{
  _MPI_COVERAGE();
  return PMPI_Testany (count, array_of_requests, index, flag, status);
}

