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
  int index;                                      
  index = _MPI_FindType(datatype);
  if (datatype == _MPI_NOT_OK)
  {
    _MPI_ERR_ROUTINE (MPI_ERR_TYPE, "MPI_TYPE_LB: invalid datatype.");
    MPI_Abort (MPI_COMM_NULL, MPI_ERR_TYPE);
  }
  *displacement = _MPI_TYPE_LIST[index].lb;
  return MPI_SUCCESS;
}

