/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  MPI_Type_hvector.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Type_hvector( 
        int count, 
        int blocklen, 
        MPI_Aint stride, 
        MPI_Datatype old_type, 
        MPI_Datatype *newtype )
{
  _MPI_COVERAGE();
  return PMPI_Type_hvector (count, blocklen, stride, old_type, newtype);
}

