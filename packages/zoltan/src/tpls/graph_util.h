/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifndef __COMMON_H
#define __COMMON_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zoltan_util.h"
#include "third_library_const.h"

extern int Zoltan_Verify_Graph(MPI_Comm, indextype *, indextype *,
       indextype *, weighttype *, weighttype *,
       int, int, int, int, int);

extern int Zoltan_Scatter_Graph(indextype **, indextype **,
       indextype **, weighttype **, indextype **, weighttype **,
       realtype **, int, int, ZZ *, ZOLTAN_COMM_OBJ **);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
