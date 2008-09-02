/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  MPI_Group_union.c    ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Group_union ( MPI_Group group1, MPI_Group group2, 
                     MPI_Group *group_out )
{
  _MPI_COVERAGE();
  return PMPI_Group_union (group1, group2, group_out);
}

