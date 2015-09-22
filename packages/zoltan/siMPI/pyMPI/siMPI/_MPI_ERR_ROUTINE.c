/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:54 $
 *    Revision: 1.2 $
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********    _MPI_ERR_ROUTINE.c     ********************/
/******************************************************************/
/* Author : Lisa Alano June 19 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#ifndef _MPI_ERR_ROUT
#define _MPI_ERR_ROUT

int _MPI_ERR_ROUTINE(MPI_Error_Class error, char* message)
{
  _MPI_COVERAGE();
#ifdef _MPI_DEBUG
  printf("ERROR: %d %s\n",error, message);
#endif 
  return 0;
}

#endif
