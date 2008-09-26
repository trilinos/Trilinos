/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************       MPI_Irecv.c        ************************/
/****************************************************************************/
/* Author : Lisa Alano June 27 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Irecv (void* message, int count, MPI_Datatype datatype, int source, 
        int tag, MPI_Comm comm, MPI_Request* request)
{
  _MPI_COVERAGE();
  return PMPI_Irecv (message,count, datatype, source, tag, comm, request);
}
