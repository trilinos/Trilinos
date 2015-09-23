/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
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

