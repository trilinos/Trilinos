/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
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

