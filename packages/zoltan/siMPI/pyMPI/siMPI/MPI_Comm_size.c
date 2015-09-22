/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:49 $
 *    Revision: 1.2 $
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********    MPI_Comm_size.c        ********************/
/******************************************************************/
/* Author : Lisa Alano June 19 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"


int MPI_Comm_size(MPI_Comm comm, int* size)
{
  _MPI_COVERAGE();
  return PMPI_Comm_size(comm, size);
}

