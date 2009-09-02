/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_Comm_remote_group.c ************************/
/****************************************************************************/
/* Author : Lisa Alano July 18 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_Comm_remote_group ( MPI_Comm comm, MPI_Group *group )
{
  _MPI_COVERAGE();
  return PMPI_Comm_remote_group (comm, group); 
}

