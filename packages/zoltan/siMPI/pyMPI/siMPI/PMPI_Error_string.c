/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:52 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *******************   PMPI_Error_string.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"
#include <string.h>

/* STUB */
int PMPI_Error_string( int errorcode, char *string, int *resultlen )
{
  char* p;
  int n;

  switch(errorcode) {
  case MPI_SUCCESS		: p = "MPI_SUCCESS "; break;
  case MPI_ERR_BUFFER		: p = "MPI_ERR_BUFFER"; break;
  case MPI_ERR_COUNT		: p = "MPI_ERR_COUNT"; break;
  case MPI_ERR_TYPE		: p = "MPI_ERR_TYPE"; break;
  case MPI_ERR_TAG		: p = "MPI_ERR_TAG"; break;
  case MPI_ERR_COMM		: p = "MPI_ERR_COMM"; break;
  case MPI_ERR_RANK		: p = "MPI_ERR_RANK"; break;
  case MPI_ERR_ROOT		: p = "MPI_ERR_ROOT"; break;
  case MPI_ERR_GROUP		: p = "MPI_ERR_GROUP"; break;
  case MPI_ERR_OP		: p = "MPI_ERR_OP"; break;
  case MPI_ERR_TOPOLOGY		: p = "MPI_ERR_TOPOLOGY"; break;
  case MPI_ERR_DIMS		: p = "MPI_ERR_DIMS"; break;
  case MPI_ERR_ARG		: p = "MPI_ERR_ARG"; break;
  case MPI_ERR_UNKNOWN		: p = "MPI_ERR_UNKNOWN"; break;
  case MPI_ERR_TRUNCATE		: p = "MPI_ERR_TRUNCATE"; break;
  case MPI_ERR_OTHER		: p = "MPI_ERR_OTHER"; break;
  case MPI_ERR_IN_STATUS	: p = "MPI_ERR_IN_STATUS"; break;
  case MPI_ERR_PENDING		: p = "MPI_ERR_PENDING"; break;
  case MPI_ERR_REQUEST		: p = "MPI_ERR_REQUEST"; break;
  case MPI_ERR_LASTCODE		: p = "MPI_ERR_LASTCODE"; break;
  case MPI_ERR_INTERN		: p = "MPI_ERR_INTERN"; break;
  default			: p = 0;
  }

  if ( string == 0 ) return MPI_ERR_ARG;
  if ( resultlen == 0 ) return MPI_ERR_ARG;
  n = *resultlen;
  if ( n <= 0 ) return MPI_ERR_ARG;

  if ( p ) {
    strncpy(string,p,n);
    *resultlen = strlen(p);
    if ( *resultlen < n ) *resultlen = n;
    return MPI_SUCCESS;
  }
  string[0] = 0;
  *resultlen = 0;
  return MPI_ERR_ARG;
}

