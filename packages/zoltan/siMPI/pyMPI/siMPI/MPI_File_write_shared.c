/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_File_write_shared.c  ***********************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_File_write_shared(MPI_File fh, void *buf, int count, 
                          MPI_Datatype datatype, MPI_Status *status)
{
  _MPI_COVERAGE();
  return PMPI_File_write_shared (fh, buf, count, datatype, status);
}

