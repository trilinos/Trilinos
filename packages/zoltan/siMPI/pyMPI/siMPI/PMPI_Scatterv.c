/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  PMPI_Scatterv.c   ***************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int PMPI_Scatterv ( 
        void *sendbuf, 
        int *sendcnts, 
        int *displs, 
        MPI_Datatype sendtype, 
        void *recvbuf, 
        int recvcnt,  
        MPI_Datatype recvtype, 
        int root, 
        MPI_Comm comm )
{
  fprintf(stderr,"%s:%d: NOT IMPLEMENTED\n",__FILE__,__LINE__);
  return MPI_Abort((MPI_Comm)NULL, MPI_UNDEFINED); 
}

