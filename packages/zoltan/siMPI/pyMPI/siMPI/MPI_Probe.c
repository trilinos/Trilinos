/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *************************  MPI_Probe.c   ***************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Probe( int source, int tag, MPI_Comm comm, MPI_Status *status )
{
  _MPI_COVERAGE();
  return PMPI_Probe (source, tag, comm, status);
}

