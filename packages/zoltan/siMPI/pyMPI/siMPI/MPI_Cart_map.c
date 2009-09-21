/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:49 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************      MPI_Cart_map.c      ************************/
/****************************************************************************/
/* Author : Lisa Alano July 18 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_Cart_map (
        MPI_Comm comm_old,
        int ndims,
        int *dims,
        int *periods,
        int *newrank)
{
  _MPI_COVERAGE();
  return PMPI_Cart_map(comm_old, ndims, dims, periods, newrank); 
}

