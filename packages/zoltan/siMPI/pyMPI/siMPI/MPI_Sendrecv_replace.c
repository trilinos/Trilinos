/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************  MPI_Sendrecv_replace.c  ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Sendrecv_replace( void *buf, int count, MPI_Datatype datatype, 
                        int dest, int sendtag, int source, int recvtag, 
                        MPI_Comm comm, MPI_Status *status )
{
  _MPI_COVERAGE();
  return PMPI_Sendrecv_replace (buf, count, datatype, dest, sendtag, source, recvtag, comm, status);
}

