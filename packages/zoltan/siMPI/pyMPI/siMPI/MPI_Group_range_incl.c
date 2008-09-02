/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_Group_range_incl.c  ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Group_range_incl ( MPI_Group group, int n, int ranges[][3], 
                         MPI_Group *newgroup )
{
  _MPI_COVERAGE();
  return PMPI_Group_range_incl (group, n, ranges, newgroup);
}

