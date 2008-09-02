/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *****************  MPI_Type_get_envelope.c  ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Type_get_envelope(
        MPI_Datatype datatype, 
        int *num_integers, 
        int *num_addresses, 
        int *num_datatypes, 
        int *combiner)
{
  _MPI_COVERAGE();
  return PMPI_Type_get_envelope (datatype, num_integers, num_addresses, num_datatypes, combiner);
}

