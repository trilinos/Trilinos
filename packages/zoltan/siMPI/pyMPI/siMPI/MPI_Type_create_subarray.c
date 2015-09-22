/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  **************** MPI_Type_create_subarray.c ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Type_create_subarray(
        int ndims, 
        int *array_of_sizes, 
        int *array_of_subsizes, 
        int *array_of_starts, 
        int order, 
        MPI_Datatype oldtype, 
        MPI_Datatype *newtype)
{
  _MPI_COVERAGE();
  return PMPI_Type_create_subarray (ndims, array_of_sizes, array_of_subsizes, array_of_starts, 
                                    order, oldtype, newtype);
}

