/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *************** MPI_File_get_type_extent.c  ************************/
/****************************************************************************/
/* Author : Lisa Alano July 22 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype, 
                             MPI_Aint *extent)
{
  _MPI_COVERAGE();
  return PMPI_File_get_type_extent (fh, datatype, extent); 
}

