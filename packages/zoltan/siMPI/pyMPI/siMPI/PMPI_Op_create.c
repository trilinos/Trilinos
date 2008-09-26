/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/
/****************************************************************************/
/* FILE  ******************    PMPI_Op_create.c      ************************/
/****************************************************************************/
/* Author : Lisa Alano July 8 2002                                          */
/* Copyright (c) 2002 University of California Regents                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

/*==========================================================================*/
int PMPI_Op_create( MPI_User_function* function, int commute, MPI_Op* op) {
  int f,c;
  int status;

  f = _MPI_Function_check(function);
  c = _MPI_Commute_check(commute);

  if ( (f==MPI_SUCCESS)&&(c==MPI_SUCCESS) ) {
#ifdef KDD_REMOVED_DEBUG
    printf("%s:%d: CREATE %d %d\n",__FILE__,__LINE__,f,c); /* HACK: */
#endif
    status = _MPI_Op_insert (function, commute, op);
#ifdef KDD_REMOVED_DEBUG
    printf("%s:%d: created %d\n",__FILE__,__LINE__,*op); /* HACK: */
#endif
  } else {
#ifdef KDD_REMOVED_DEBUG
    printf("%s:%d: cannot CREATE %d %d\n",__FILE__,__LINE__,f,c); /* HACK: */
#endif
    if ( op ) *op = MPI_OP_NULL;
    status = MPI_ERR_OTHER;
  }

  return status;
}
/*==========================================================================*/
