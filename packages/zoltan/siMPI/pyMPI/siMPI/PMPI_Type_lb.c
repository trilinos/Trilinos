/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *************************  PMPI_Type_lb.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano August 7 2002                                        */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Type_lb ( MPI_Datatype datatype, MPI_Aint *displacement )
{
  fprintf(stderr,"%s:%d: NOT IMPLEMENTED\n",__FILE__,__LINE__);
  return MPI_Abort((MPI_Comm)0, MPI_UNDEFINED);
}

