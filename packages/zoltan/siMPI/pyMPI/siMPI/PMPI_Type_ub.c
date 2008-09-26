/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *************************  PMPI_Type_ub.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano August 6 2002                                        */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Type_ub ( MPI_Datatype datatype, MPI_Aint *displacement )
{
  int index;
  index = _MPI_FindType(datatype);
  if (datatype == _MPI_NOT_OK)
  {
    _MPI_ERR_ROUTINE (MPI_ERR_TYPE, "MPI_TYPE_UB: invalid datatype.");
    MPI_Abort (MPI_COMM_NULL, MPI_ERR_TYPE);
  }
  *displacement = _MPI_TYPE_LIST[index].ub;
  return MPI_SUCCESS;
}

