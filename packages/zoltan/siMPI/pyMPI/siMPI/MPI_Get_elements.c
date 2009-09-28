/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:50 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  MPI_Get_elements.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Get_elements ( MPI_Status *status, MPI_Datatype datatype, 
                      int *elements )
{
  _MPI_COVERAGE();
  return PMPI_Get_elements (status, datatype, elements);
}

