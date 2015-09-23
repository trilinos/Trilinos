/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:50 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************   MPI_File_get_byte_offset.c   ******************/
/****************************************************************************/
/* Author : Lisa Alano July 22 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset, MPI_Offset *disp)
{
  _MPI_COVERAGE();
  return PMPI_File_get_byte_offset(fh, offset, disp); 
}

