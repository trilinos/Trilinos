/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ************************  PMPI_Waitsome.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int PMPI_Waitsome( 
        int incount, 
        MPI_Request array_of_requests[], 
        int *outcount, 
        int array_of_indices[], 
        MPI_Status array_of_statuses[] )
{
  return MPI_Abort((MPI_Comm)0, MPI_UNDEFINED); 
}

