/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ***************************   MPI_Unpack.c  ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Unpack ( void *inbuf, int insize, int *position, 
                void *outbuf, int outcount, MPI_Datatype datatype, 
                MPI_Comm comm )
{
  _MPI_COVERAGE();
  return PMPI_Unpack (inbuf, insize, position, outbuf, outcount, datatype, comm);
}

