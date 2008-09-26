/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_Type_contiguous.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Type_contiguous( 
        int count,
        MPI_Datatype old_type,
        MPI_Datatype *newtype)
{
  _MPI_COVERAGE();
  return PMPI_Type_contiguous (count, old_type, newtype);
}

