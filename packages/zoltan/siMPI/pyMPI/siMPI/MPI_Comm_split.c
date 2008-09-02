/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************    MPI_Comm_split.c      ************************/
/****************************************************************************/
/* Author : Lisa Alano July 19 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_Comm_split ( MPI_Comm comm, int color, int key, MPI_Comm *comm_out )
{
  _MPI_COVERAGE();
  return PMPI_Comm_split (comm, color, key, comm_out); 
}

