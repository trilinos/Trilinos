/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
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

