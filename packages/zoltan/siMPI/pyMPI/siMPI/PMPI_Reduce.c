/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************      PMPI_Reduce.c       ************************/
/****************************************************************************/
/* Author : Lisa Alano July 1 2002                                          */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int PMPI_Reduce ( void *sendbuf, void *recvbuf, int count, 
   MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm )
{
  int size;
  int retval;

 _MPI_COVERAGE();

  /* ----------------------------------------------- */
  /* We shortcut on length 0 sends since there is no */
  /* work to do...                                   */
  /* ----------------------------------------------- */
  if ( count == 0 ) return MPI_SUCCESS;

  /* ----------------------------------------------- */
  /* First we verify that the sendbuf, recvbuf, and  */
  /* operator are OK.                                */
  /* ----------------------------------------------- */
  retval = _MPI_checks(sendbuf, count, datatype, root, 1, comm);
  if ( retval != MPI_SUCCESS ) {
    _MPI_ERR_ROUTINE (MPI_ERR_OTHER, "MPI_REDUCE : Invalid argument");
    MPI_Abort(comm, retval);
    return _MPI_NOT_OK;
  }
      
  if ( _MPI_checkBuffer(recvbuf) != MPI_SUCCESS ) {
    _MPI_ERR_ROUTINE (MPI_ERR_BUFFER, "MPI_REDUCE : Invalid buffer pointer");
    MPI_Abort(comm, MPI_ERR_BUFFER);
    return _MPI_NOT_OK;
  }
   
  if ( _MPI_checkOp(op) != MPI_SUCCESS ) {
    _MPI_ERR_ROUTINE(MPI_ERR_OP, "MPI_REDUCE : Invalid MPI_Op");
    MPI_Abort(comm, MPI_ERR_OP);
    return _MPI_NOT_OK;
  }

  /* ----------------------------------------------- */
  /* Guard against buffer overlap...                 */
  /* ----------------------------------------------- */
  size = _MPI_calculateSize(count, datatype);
  if (  _MPI_check_overlap(sendbuf, recvbuf, size) != MPI_SUCCESS ) {
    _MPI_ERR_ROUTINE (MPI_ERR_BUFFER, "MPI_REDUCE : Invalid buffer pointer: Arguments must specify different buffers (no aliasing)");
    MPI_Abort(comm, MPI_ERR_BUFFER);
    return _MPI_NOT_OK;
  }

  /* ----------------------------------------------- */
  /* All functions (at 1 proc), just copy over the   */
  /* data.                                           */
  /* ----------------------------------------------- */
  _MPI_Default_Op(sendbuf, recvbuf, &count, &datatype);
  return MPI_SUCCESS;
}
/*==========================================================================*/

