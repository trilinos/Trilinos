/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:50 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *************** MPI_Graph_neighbors_count.c ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Graph_neighbors_count ( MPI_Comm comm, int rank, int *nneighbors )
{
  _MPI_COVERAGE();
  return PMPI_Graph_neighbors_count (comm, rank, nneighbors);
}

