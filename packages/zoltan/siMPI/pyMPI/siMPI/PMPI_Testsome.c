/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ************************  PMPI_Testsome.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Testsome( 
        int incount, 
        MPI_Request array_of_requests[], 
        int *outcount, 
        int array_of_indices[], 
        MPI_Status array_of_statuses[] )
{
  int i, j, num_active, flag, t;
  *outcount = 0;
  j = 0;
  num_active = 0;
  for ( i = 0; i < incount; ++i )
  {
    if ( array_of_requests[i] != MPI_REQUEST_NULL )
    {
      ++num_active;
      if ( array_of_statuses == MPI_STATUSES_IGNORE )
        t = PMPI_Test( array_of_requests+i, &flag, MPI_STATUS_IGNORE );
      else
        t = PMPI_Test( array_of_requests+i, &flag, array_of_statuses+j );
      if ( t != MPI_SUCCESS )
        return t;
      
      if ( flag )
      {
        array_of_indices[j] = i;
        ++(*outcount);
        ++j;
      }
    }
  }
  
  if ( num_active == 0 )
    *outcount = MPI_UNDEFINED;
  
  return MPI_SUCCESS;
}

