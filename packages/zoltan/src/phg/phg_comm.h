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
#ifndef __ZOLTAN_PHG_COMM_H
#define __ZOLTAN_PHG_COMM_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
    
/********************************************/
/* Communication and Distribution variables */
/********************************************/
struct PHGCommStruct {
    /* the following three (Communicator, Proc and Num_Proc) are copies from ZZ
       just for convenience */
  MPI_Comm Communicator;          /*  The MPI Communicator.                  */
  int Proc;                       /*  The processor's ID within the MPI
                                      Communicator.                          */
  int Num_Proc;                   /*  The number of processors in the MPI
                                      Communicator.                          */
  int nProc_x;    /* number of processors in x-direction of 2D data distrib.  */
  int nProc_y;    /* number of processors in y-direction of 2D data distrib.  */
                  /* nProc_x * nProc_y should equal number of processors!     */
  int myProc_x;   /* my processor's row block number in [0,nProc_x-1] */
  int myProc_y;   /* my processor's column block number in [0,nProc_y-1] */
  MPI_Comm row_comm; /* my processor's row communicator */
  MPI_Comm col_comm; /* my processor's column communicator */
};

typedef struct PHGCommStruct PHGComm;

    
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_PHG_COMM_H */
    
