/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:49 $
 *    Revision: 1.2 $
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********      MPI_Comm_create.c      ********************/
/******************************************************************/
/* Author : Lisa Alano June 21 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#include <stdio.h>
#include "mpi.h"

int MPI_Comm_create(MPI_Comm comm, MPI_Group new_group, MPI_Comm* new_comm)
{
  _MPI_COVERAGE();
  return  PMPI_Comm_create(comm, new_group, new_comm);
}

