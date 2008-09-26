/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ************************   MPI_Info_get.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Info_get(MPI_Info info, char *key, int valuelen, char *value, int *flag)
{
  _MPI_COVERAGE();
  return PMPI_Info_get (info, key, valuelen, value, flag);
}

