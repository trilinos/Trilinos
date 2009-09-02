/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************      MPI_Attr_get.c      ************************/
/****************************************************************************/
/* Author : Lisa Alano July 18 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

/* STUB */
int MPI_Attr_get ( 
        MPI_Comm comm, 
        int keyval, 
        void *attr_value, 
        int *flag )
{
  _MPI_COVERAGE();
  return PMPI_Attr_get(comm, keyval, attr_value, flag); 
}

