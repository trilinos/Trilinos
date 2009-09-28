/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/* FILE  ******************      MPI_BCAST.c        ************************/
/****************************************************************************/
/* Author : Lisa Alano July 9 2002                                          */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int MPI_Bcast ( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm ) {
  _MPI_COVERAGE();
  return PMPI_Bcast(buffer, count, datatype, root, comm);
}
/*==========================================================================*/

