/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:50 $
 *    Revision: 1.2 $
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

