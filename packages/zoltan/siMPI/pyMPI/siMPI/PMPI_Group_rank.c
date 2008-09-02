/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ********************   PMPI_Group_rank.c    ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int PMPI_Group_rank ( MPI_Group group, int *rank )
{
  if (group == MPI_GROUP_EMPTY) {
    *rank = MPI_ERR_RANK;
    return MPI_ERR_RANK;
  }

  *rank = _MPI_RANK;
  return MPI_SUCCESS;
}

