/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********  MPI_File_write_ordered_begin.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_File_write_ordered_begin(MPI_File fh, void *buf, int count, 
                              MPI_Datatype datatype)
{
  _MPI_COVERAGE();
  return PMPI_File_write_ordered_begin (fh, buf, count, datatype);
}

