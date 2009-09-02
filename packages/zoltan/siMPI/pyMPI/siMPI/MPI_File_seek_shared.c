/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************   MPI_File_seek_shared.c   **********************/
/****************************************************************************/
/* Author : Lisa Alano July 22 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence)
{
  _MPI_COVERAGE();
  return PMPI_File_seek_shared (fh, offset, whence); 
}

