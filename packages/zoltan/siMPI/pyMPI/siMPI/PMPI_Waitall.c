/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************       PMPI_Waitall.c     ************************/
/****************************************************************************/
/* Author : Lisa Alano July 11 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int PMPI_Waitall(
        int count,
        MPI_Request array_of_requests[],
        MPI_Status  array_of_statuses[] )
{
  int retval, i;

  retval = _MPI_Check_Request_Array(count, array_of_requests); 
  if ( retval != MPI_SUCCESS ) return MPI_UNDEFINED;

  /* POLL THROUGH REQUESTS */
  for(i=0; i<count; i++) {
    if (array_of_requests[i] != MPI_REQUEST_NULL) {
      retval = PMPI_Wait(array_of_requests+i,array_of_statuses+i);
      if (array_of_requests[i] != MPI_REQUEST_NULL) return retval;
    }
  }  
  return MPI_SUCCESS;
  
}
/*==========================================================================*/

