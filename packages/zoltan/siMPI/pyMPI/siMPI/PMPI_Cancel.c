/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************     PMPI_Cancel.c        ************************/
/****************************************************************************/
/* Author : Lisa Alano July 18 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Cancel( MPI_Request *request )
{
  _MPI_COVERAGE();
  if ( !request )
    return MPI_ERR_REQUEST;
  
  (*request)->cancel = 1;  /* true */
  return MPI_SUCCESS;
}
