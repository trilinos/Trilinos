/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:49 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  *****************  MPI_Errhandler_create.c  ************************/
/****************************************************************************/
/* Author : Lisa Alano July 24 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int MPI_Errhandler_create( 
        MPI_Handler_function *function,
        MPI_Errhandler       *errhandler)
{
  _MPI_COVERAGE();
  return PMPI_Errhandler_create (function, errhandler);
}

