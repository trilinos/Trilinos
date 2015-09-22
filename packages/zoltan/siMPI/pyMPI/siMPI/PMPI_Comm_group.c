/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:52 $
 *    Revision: 1.2 $
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********    PMPI_Comm_group.c      ********************/
/******************************************************************/
/* Author : Lisa Alano June 20 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

/*============================================================================*/
int PMPI_Comm_group(MPI_Comm comm, MPI_Group* group)
{
int index;
  if (_MPI_CHECK_STATUS(&comm) == _MPI_OK)
  {
    if  (_MPI_Comm_check_legal(comm, &index)==MPI_SUCCESS)
    {
      if (comm == MPI_COMM_NULL)
         return MPI_ERR_COMM;
      *group = _MPI_COMM_LIST[index].group;
      return MPI_SUCCESS;
    }
  }
  return MPI_ERR_COMM;
}

/*============================================================================*/

