/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************   MPI_Errhandler_set.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Errhandler_set( MPI_Comm comm, MPI_Errhandler errhandler )
{
  _MPI_COVERAGE();
  return PMPI_Errhandler_set (comm, errhandler);
}

