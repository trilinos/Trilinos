/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  **********************  MPI_Recv_init.c   **************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Recv_init( void *buf, int count, MPI_Datatype datatype, int source, 
                  int tag, MPI_Comm comm, MPI_Request *request )
{
  _MPI_COVERAGE();
  return PMPI_Recv_init (buf, count, datatype, source, tag, comm, request);
}

