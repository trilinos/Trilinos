/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/***********************************************************************************************/
/* FILE  **************************      PMPI_Ibsend.c         *********************************/
/***********************************************************************************************/
/* Author : Lisa Alano July 15 2002                                                            */
/* Copyright (c) 2002 University of California Regents                                         */
/***********************************************************************************************/

#include "mpi.h"

int PMPI_Ibsend (void* message, int count, MPI_Datatype datatype, int dest,
                 int tag, MPI_Comm comm, MPI_Request* request) {
  return PMPI_Isend(message,count,datatype,dest,tag,comm,request);
}

