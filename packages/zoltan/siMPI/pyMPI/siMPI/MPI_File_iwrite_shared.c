/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ****************   MPI_File_iwrite_shared.c   **********************/
/****************************************************************************/
/* Author : Lisa Alano July 22 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_File_iwrite_shared(MPI_File fh, void *buf, int count, 
                       MPI_Datatype datatype, MPIO_Request *request)
{
  _MPI_COVERAGE();
  return PMPI_File_iwrite_shared(fh, buf, count, datatype, request); 
}

