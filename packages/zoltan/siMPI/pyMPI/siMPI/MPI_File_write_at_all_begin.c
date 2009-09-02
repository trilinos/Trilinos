/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ************** MPI_File_write_at_all_begin.c ***********************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf,
                         int count, MPI_Datatype datatype)
{
  _MPI_COVERAGE();
  return PMPI_File_write_at_all_begin (fh, offset, buf, count, datatype); 
}

