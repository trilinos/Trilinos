/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_File_set_atomicity.c  **********************/
/****************************************************************************/
/* Author : Lisa Alano July 22 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_File_set_atomicity(MPI_File fh, int flag)
{
  _MPI_COVERAGE();
  return PMPI_File_set_atomicity (fh, flag); 
}

