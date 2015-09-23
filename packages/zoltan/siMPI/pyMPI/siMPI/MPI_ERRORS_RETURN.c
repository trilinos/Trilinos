/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    Author: patmiller $
 *    Date: 2007/06/11 14:12:49 $
 *    Revision: 1.2 $
 ****************************************************************************/
/******************************************************************/
/* FILE  ***********    MPI_ERRORS_RETURN.c    ********************/
/******************************************************************/
/* Author : Lisa Alano June 24 2002                               */
/* Copyright (c) 2002 University of California Regents            */
/******************************************************************/

#include "mpi.h"

void MPI_ERRORS_RETURN (MPI_Comm* comm, int* error_code, ...)
{
  _MPI_COVERAGE();
  switch(*error_code)
  {
    case MPI_ERR_BUFFER:          
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_BUFFER, "Error with Buffer.");
      break;
    case MPI_ERR_COUNT:          
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_COUNT, "Error with count value.");
      break;
    case MPI_ERR_TYPE:            
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_TYPE, "Error with datatype.");
      break;
    case MPI_ERR_TAG:       
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_TAG, "Error with tag.");
      break;
    case MPI_ERR_COMM:     
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_COMM, "Error with communicator.");
      break;
    case MPI_ERR_RANK:   
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_RANK, "Error with rank.");
      break;
    case MPI_ERR_ROOT:           
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_ROOT, "Error with root.");
      break;
    case MPI_ERR_GROUP:         
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_GROUP, "Error with group.");
      break;
    case MPI_ERR_OP:         
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_OP, "Error  with Op");
      break;
    case MPI_ERR_TOPOLOGY:      
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_TOPOLOGY, "Error with topology.");
      break;
    case MPI_ERR_DIMS:         
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_DIMS, "Error with Dims.");
      break;
    case MPI_ERR_ARG:         
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_ARG, "Error with argument.");
      break;
    case MPI_ERR_UNKNOWN:    
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_UNKNOWN, "Error unknown.");
      break;
    case MPI_ERR_TRUNCATE:      
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_TRUNCATE, "Error with truncate.");
      break;
    case MPI_ERR_OTHER:        
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_OTHER, "Error with other."); 
      break;
    case MPI_ERR_IN_STATUS:   
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_IN_STATUS, "Error with Init status.");
      break;
    case MPI_ERR_PENDING:    
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_PENDING, "Error pending.");
      break;
    case MPI_ERR_REQUEST:       
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_REQUEST, "Error with request.");
      break;
    case MPI_ERR_LASTCODE:   
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_LASTCODE, "Error with Last code.");
      break;
    case MPI_ERR_INTERN:   
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(MPI_ERR_INTERN, "Error with Internal.");
      break;
    default:
      _MPI_COVERAGE();
      _MPI_ERR_ROUTINE(0, "Unknown Error.");
      break;
  }
  return;
}

