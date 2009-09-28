/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********    _MPI_CHECK_STATUS.c      ********************/
/******************************************************************/
/* Author : Lisa ALano June 19 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int _MPI_CHECK_STATUS(MPI_Comm* comm)
{
  int errorcode = MPI_ERR_IN_STATUS;
  _MPI_COVERAGE();
  if (_MPI_INIT_STATUS == _MPI_STARTED) {
    return MPI_SUCCESS;
  } else {
    _MPI_ERR_ROUTINE(MPI_ERR_IN_STATUS, "MPI is not in good status.");
    MPI_Abort(0, errorcode);
    /*(_MPI_COMM_LIST[*comm].handler)(comm, &errorcode);*/
    return _MPI_NOT_OK;
  }
}
