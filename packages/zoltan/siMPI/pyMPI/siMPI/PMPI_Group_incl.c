/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:52 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ********************   PMPI_Group_incl.c    ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int PMPI_Group_incl ( MPI_Group group, int n, int *ranks, MPI_Group *newgroup )
{
  int i;

  /* Check for null output */
  if ( !newgroup ) {
    return MPI_ERR_ARG;
  }
  *newgroup = group; /* Default answer for people that don't check errors! */

  /* Cannot maniuplate the null group */
  if ( group == MPI_GROUP_NULL) {
    return MPI_ERR_GROUP;
  }

  if ( !ranks ) return MPI_ERR_ARG;

  for(i=0;i<n;++i) {
    if ( ranks[i] != 0 ) return MPI_ERR_RANK;
     *newgroup = MPI_GROUP_WORLD;
    /* Including rank 0 in the full group leaves the full group */
  }

  return MPI_SUCCESS;
}
