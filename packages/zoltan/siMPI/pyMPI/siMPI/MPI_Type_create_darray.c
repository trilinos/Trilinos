/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ****************  MPI_Type_create_darray.c  ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Type_create_darray(int size, int rank, int ndims, 
                           int *array_of_gsizes, int *array_of_distribs, 
                           int *array_of_dargs, int *array_of_psizes, 
                           int order, MPI_Datatype oldtype, 
                           MPI_Datatype *newtype) 
{
  _MPI_COVERAGE();
  return PMPI_Type_create_darray (size, rank, ndims, array_of_gsizes, array_of_distribs, 
                           array_of_dargs, array_of_psizes, order, oldtype, newtype);
}

