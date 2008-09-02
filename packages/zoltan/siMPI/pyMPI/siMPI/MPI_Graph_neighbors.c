/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_Graph_neighbors.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Graph_neighbors ( MPI_Comm comm, int rank, int maxneighbors, 
                        int *neighbors )
{
  _MPI_COVERAGE();
  return PMPI_Graph_neighbors (comm, rank, maxneighbors, neighbors);
}

