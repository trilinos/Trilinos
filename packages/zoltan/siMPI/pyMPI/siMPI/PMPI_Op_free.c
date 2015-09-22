/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************    PMPI_Op_free.c        ************************/
/****************************************************************************/
/* Author : Lisa Alano July 8 2002                                          */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int PMPI_Op_free( MPI_Op *op )
{
  int index;
  if ( !op ) return MPI_ERR_ARG;

  index = *op-_MPI_OP_OFFSET;
  if ((index >=0 ) && (index < _MPI_OP_ARRAY_SIZE) ) {
    if (_MPI_OP_LIST[index].valid != _MPI_VALID) return MPI_ERR_OP;
    return MPI_SUCCESS;
  }

  if (_MPI_Op_invalid(index)==MPI_SUCCESS) {
    *op = MPI_OP_NULL;
    _MPI_OP_COUNT--;
    return MPI_SUCCESS;
  }
  return MPI_ERR_OP;
}
/*==========================================================================*/
