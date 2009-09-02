/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/***********************************************************************************************/
/* FILE  *************************         MPI_Send.c        ***********************************/
/***********************************************************************************************/
/* Author : Lisa Alano June 26 2002                                                            */
/* Copyright (c) 2002 University of California Regents                                         */
/***********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "mpi.h"

int MPI_Send (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
  _MPI_COVERAGE();
  return PMPI_Send (message, count, datatype, dest, tag, comm);
}
