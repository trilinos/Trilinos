/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:49 $
 *    Revision: 1.2 $
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

