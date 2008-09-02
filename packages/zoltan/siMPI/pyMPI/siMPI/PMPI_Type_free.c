/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ***********************  PMPI_Type_free.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Type_free ( MPI_Datatype *datatype ) {
  int retval;

  retval = _MPI_Free_datatype (*datatype);
  if ( datatype ) *datatype = _MPI_NOT_VALID;
  return retval;
}

