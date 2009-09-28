/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************     PMPI_Allreduce.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 17 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int PMPI_Allreduce ( void *sendbuf, void *recvbuf, int count, 
                   MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )
{
  _MPI_COVERAGE();
  return PMPI_Reduce( sendbuf,recvbuf,count,datatype,op,0,comm);
}

