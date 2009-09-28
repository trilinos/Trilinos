/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:49 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_Dims_create.c       ************************/
/****************************************************************************/
/* Author : Lisa Alano July 19 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_Dims_create(
        int nnodes, 
        int ndims, 
        int *dims)
{
  _MPI_COVERAGE();
  return PMPI_Dims_create ( nnodes, ndims, dims); 
}

