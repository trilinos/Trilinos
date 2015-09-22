/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ****************************   MPI_Pack.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Pack ( void *inbuf, int incount, MPI_Datatype datatype, 
               void *outbuf, int outcount, int *position, MPI_Comm comm )
{
  _MPI_COVERAGE();
  return PMPI_Pack (inbuf, incount, datatype, outbuf, outcount, position, comm);
}

