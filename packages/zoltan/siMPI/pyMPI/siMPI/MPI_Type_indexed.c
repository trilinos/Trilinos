/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  MPI_Type_indexed.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Type_indexed( 
        int count, 
        int blocklens[], 
        int indices[], 
        MPI_Datatype old_type, 
        MPI_Datatype *newtype )
{
  _MPI_COVERAGE();
  return PMPI_Type_indexed (count, blocklens, indices, old_type, newtype);
}

