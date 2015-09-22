/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:53 $
 *    Revision: 1.2 $
 ****************************************************************************/
/***********************************************************************************************/
/* FILE  **************************       PMPI_Msend.c        ***********************************/
/***********************************************************************************************/
/* Author : Lisa Alano June 25 2002                                                            */
/* Copyright (c) 2002 University of California Regents                                         */
/***********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/* This is the old PMPI_Send.c */
/*=============================================================================================*/
int PMPI_Msend (void* message, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
{
  int size, retval;
  char* p;

  /* option 1: Assert( _MPI_IS_INITIALIZED() ); */
 _MPI_COVERAGE();
  
  _MPI_CHECK_STATUS(&comm);
  retval = _MPI_checks(message, count, datatype, dest, tag, comm);
  /* option 2: STANDARD_MPI_CHECK(comm); */
  if (retval == MPI_SUCCESS)
  {
    size = _MPI_calculateSize(count, datatype);  
    p = (char *)_MPI_safeMalloc(size, "Error with malloc for send buffer."); 
    p = memcpy(p, message, size);
    retval =_MPI_Buff_Insert(p, count, datatype, tag, comm);
    return retval;
  } else { 
     _MPI_ERR_ROUTINE (retval, "MPI_SEND / MPI_ISEND: argument error");
     MPI_Abort (comm, retval);
  }

 _MPI_COVERAGE();
  return _MPI_NOT_OK;
}
/*=============================================================================================*/
