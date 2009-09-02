/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  **************************  MPI_Type_lb.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Type_lb ( MPI_Datatype datatype, MPI_Aint *displacement )
{
  _MPI_COVERAGE();
  return PMPI_Type_lb (datatype, displacement);
}

