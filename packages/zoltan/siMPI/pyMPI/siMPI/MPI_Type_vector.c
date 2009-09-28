/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  **********************  MPI_Type_vector.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Type_vector( 
        int count, 
        int blocklen, 
        int stride, 
        MPI_Datatype old_type, 
        MPI_Datatype *newtype )
{
  _MPI_COVERAGE();
  return PMPI_Type_vector (count, blocklen, stride, old_type, newtype);
}

