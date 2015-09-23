/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:50 $
 *    Revision: 1.2 $
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

