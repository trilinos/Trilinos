/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ***********************  PMPI_Type_size.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Type_size ( MPI_Datatype datatype, int *size )
{
  int index;
  if ( _MPI_BasicType(datatype) == MPI_SUCCESS )
    return _MPI_getSize(datatype);
  else
  {
    index = _MPI_FindType (datatype);
    if (index == _MPI_NOT_OK)
    {
      _MPI_ERR_ROUTINE (MPI_ERR_TYPE, "MPI_TYPE_SIZE: datatype error");
      MPI_Abort (MPI_COMM_NULL, MPI_ERR_TYPE); 
    }
    return _MPI_TYPE_LIST[index].size;
  }
}

