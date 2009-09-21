/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:51 $
 *    Revision: 1.2 $
 ****************************************************************************/
/****************************************************************************/
/* FILE  ***********************  MPI_Int2handle.c   ************************/
/****************************************************************************/
/* Author : Lisa Alano July 23 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

MPI_Handle_type MPI_Int2handle( MPI_Fint f_handle, MPI_Handle_enum handle_kind )
{
  _MPI_COVERAGE();
  return PMPI_Int2handle (f_handle, handle_kind);
}

