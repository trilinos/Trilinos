/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  **********************  MPI_Type_struct.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Type_struct( 
        int count, 
        int blocklens[], 
        MPI_Aint indices[], 
        MPI_Datatype old_types[], 
        MPI_Datatype *newtype )
{
  _MPI_COVERAGE();
  return PMPI_Type_struct (count, blocklens, indices, old_types, newtype);
}

