/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************   PMPI_Comm_compare.c    ************************/
/****************************************************************************/
/* Author : Lisa Alano July 18 2002                                         */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include "mpi.h"

int PMPI_Comm_compare ( 
        MPI_Comm  comm1,
        MPI_Comm  comm2,
        int *result)
{
  int index1, index2;  /* treat an index as the MPI communicator "context" */
  MPI_Group group1, group2;
  int gcmp;
  _MPI_COVERAGE();
  *result = MPI_UNEQUAL;
  if (_MPI_CHECK_STATUS(&comm1) == _MPI_OK &&
      _MPI_CHECK_STATUS(&comm2) == _MPI_OK)
  {
    _MPI_COVERAGE();
    if  (_MPI_Comm_check_legal(comm1, &index1)==MPI_SUCCESS &&
         _MPI_Comm_check_legal(comm2, &index2)==MPI_SUCCESS)
    {
      _MPI_COVERAGE();
      if (comm1 == MPI_COMM_NULL || comm2 == MPI_COMM_NULL)
         return MPI_ERR_COMM;
      
      if ( index1 == index2 ) {
        _MPI_COVERAGE();
        *result = MPI_IDENT;
      }
      else
      {
        _MPI_COVERAGE();
        group1 = _MPI_COMM_LIST[index1].group;
        group2 = _MPI_COMM_LIST[index2].group;
        
        if ( PMPI_Group_compare( group1, group2, &gcmp ) != MPI_SUCCESS )
          return MPI_ERR_GROUP;
        
        if ( gcmp == MPI_IDENT )
          *result = MPI_CONGRUENT;
        
        if ( gcmp == MPI_SIMILAR )
          *result = MPI_SIMILAR;
      }
      
      return MPI_SUCCESS;
    }
  }
  return MPI_ERR_COMM;
}

