/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  *********************  PMPI_Test.c   *******************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Test ( 
        MPI_Request  *request,
        int          *flag,
        MPI_Status   *status)
{
  int err;
  
  /* defer to the code in Wait() */
  err = PMPI_Wait(request,status);
  
  if ( err == MPI_ERR_TAG )
  {
    /* could not match the request up, so assume not done */
    *flag = 0;
    return MPI_SUCCESS;
  }
  
  /* Flag indicates message was sent */
  if ( flag ) *flag = 1;
  
  return err;
}
