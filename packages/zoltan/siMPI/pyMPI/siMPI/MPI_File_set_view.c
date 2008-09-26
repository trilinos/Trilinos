/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_File_set_view.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 22 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype,
                MPI_Datatype filetype, char *datarep, MPI_Info info)
{
  _MPI_COVERAGE();
  return PMPI_File_set_view (fh, disp, etype, filetype, datarep, info); 
}

