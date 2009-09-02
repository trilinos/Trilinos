/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  **************  MPI_Status_set_cancelled.c  ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Status_set_cancelled( MPI_Status *status, int flag )
{
  _MPI_COVERAGE();
  return PMPI_Status_set_cancelled (status, flag);
}

