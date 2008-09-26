/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ****************  MPI_File_set_errhandler.c ************************/
/****************************************************************************/
/* Author : Lisa Alano July 22 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_File_set_errhandler(MPI_File fh, MPI_Errhandler errhandler)
{
  _MPI_COVERAGE();
  return PMPI_File_set_errhandler(fh, errhandler); 
}

