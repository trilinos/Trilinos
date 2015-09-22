/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:50 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_Group_difference.c  ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Group_difference ( MPI_Group group1, MPI_Group group2, 
                         MPI_Group *group_out )
{
  _MPI_COVERAGE();
  return PMPI_Group_difference (group1, group2, group_out);
}

