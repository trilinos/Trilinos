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
 *    Revision: 1.3 $
 ****************************************************************************/


#ifndef __COMM_CONST_H
#define __COMM_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


#include <mpi.h>

/* Interface to the Zoltan Communication Package. */
/* This file should be included by user's of the  */
/* Communication package.                         */

struct Zoltan_Comm_Obj;
typedef struct Zoltan_Comm_Obj ZOLTAN_COMM_OBJ;

/* function prototypes */

extern int Zoltan_Comm_Create(ZOLTAN_COMM_OBJ **, int, int *, MPI_Comm, 
    int, int *);

extern int Zoltan_Comm_Destroy(ZOLTAN_COMM_OBJ **);

extern int Zoltan_Comm_Invert_Map(int *, int *, int, int, int **, int **, int *,
    int, int, int, int, MPI_Comm);

extern int Zoltan_Comm_Sort_Ints(int *, int *, int);

extern int Zoltan_Comm_Exchange_Sizes(int *, int *, int, int, int *, int *,
    int, int *, int, int, MPI_Comm);

extern int Zoltan_Comm_Resize(ZOLTAN_COMM_OBJ *, int *, int, int *);

extern int Zoltan_Comm_Do(ZOLTAN_COMM_OBJ *, int, char *, int, char *);

extern int Zoltan_Comm_Do_Reverse(ZOLTAN_COMM_OBJ *, int, char *, int, int *, 
    char *);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
