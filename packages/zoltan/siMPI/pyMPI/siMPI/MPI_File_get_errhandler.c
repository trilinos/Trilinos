/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ****************** MPI_File_get_errhandler.c  **********************/
/****************************************************************************/
/* Author : Lisa Alano July 22 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_File_get_errhandler(MPI_File fh, MPI_Errhandler *errhandler)
{
  _MPI_COVERAGE();
  return PMPI_File_get_errhandler(fh, errhandler); 
}

