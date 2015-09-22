/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ***************  PMPI_Type_get_contents.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Type_get_contents(
        MPI_Datatype datatype, 
        int max_integers, 
        int max_addresses, 
        int max_datatypes, 
        int *array_of_integers, 
        MPI_Aint *array_of_addresses, 
        MPI_Datatype *array_of_datatypes)
{
  int index, position;
  position = _MPI_FindType(datatype);
  if (position == _MPI_NOT_OK)
  {
    position = _MPI_BasicType (datatype);
    if (position == MPI_SUCCESS)
    {
      array_of_integers[0] = 1;
      array_of_addresses[0] = (MPI_Aint) 0;
      array_of_datatypes[0] = datatype;
      return MPI_SUCCESS;
    } else {
      _MPI_ERR_ROUTINE (MPI_ERR_TYPE, "MPI_TYPE_GET_CONTENTS: datatype error");
      MPI_Abort (MPI_COMM_NULL, MPI_ERR_TYPE);
    }
  }
  if ( (_MPI_TYPE_LIST[position].info->count < max_addresses) || (_MPI_TYPE_LIST[position].info->count < max_integers) ||
       (_MPI_TYPE_LIST[position].info->count < max_datatypes) )
  {
    _MPI_ERR_ROUTINE (MPI_ERR_TYPE, "MPI_TYPE_GET_CONTENTS: invalid max_* error");
    MPI_Abort (MPI_COMM_NULL, MPI_ERR_TYPE);    
  }
  if (_MPI_TYPE_LIST[position].sendType == _MPI_STRUCT)
  {
    for (index=0; index<max_integers; index++)
    {
      array_of_integers[index] = _MPI_TYPE_LIST[position].info->blocklen[index];
    }
    for (index=0; index<max_addresses; index++)
    {
      array_of_addresses[index] = _MPI_TYPE_LIST[position].info->stride[index];
    }
    for (index=0; index<max_datatypes; index++)
    {
      array_of_addresses[index] = _MPI_TYPE_LIST[position].info->types[index];
    }
  } else {
    _MPI_ERR_ROUTINE (MPI_ERR_TYPE, "MPI_TYPE_GET_CONTENTS: Not Struct datatype error");
    MPI_Abort (MPI_COMM_NULL, MPI_ERR_TYPE);
  }
  return MPI_SUCCESS;
}

