/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ***********************   MPI_Iprobe.c   ***************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Iprobe( int source, int tag, MPI_Comm comm, int *flag, 
               MPI_Status *status )
{
  _MPI_COVERAGE();
  return PMPI_Iprobe (source, tag, comm, flag, status);
}

