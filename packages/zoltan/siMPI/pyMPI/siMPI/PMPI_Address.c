/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************       PMPI_Address.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 17 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int PMPI_Address( void *location, MPI_Aint *address)
{
  _MPI_COVERAGE();
  if ( address ) (*address) = (MPI_Aint)(location);
  return MPI_SUCCESS;
}

