// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef __ZOLTAN_PHG_COMM_H
#define __ZOLTAN_PHG_COMM_H

#include "zz_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/********************************************/
/* Communication and Distribution variables */
/********************************************/

/* Smallest prime number of processors allowed when using 2D decompositions.
   Larger prime numbers will be adjusted. */
#define SMALL_PRIME 7

struct PHGCommStruct {
  MPI_Comm Communicator;  /* MPI Communicator for all procs partitioning
			     this HG.  May not equal zz->Communicator when
			     splitting HG among procs. */
  int myProc;     /* my processor's rank within Communicator. */
  int nProc;      /* number of proc in Communicator.
		     nProc = nProc_x * nProc_y */
  unsigned int RNGState;  /* State for random-number generator
			      w.r.t. Communicator */
  int nProc_x;    /* number of processors in x-direction of 2D data distrib.  */
  int nProc_y;    /* number of processors in y-direction of 2D data distrib.  */
		  /* nProc_x * nProc_y should equal number of processors!     */
  int myProc_x;   /* my processor's row block number in [0,nProc_x-1] */
  int myProc_y;   /* my processor's column block number in [0,nProc_y-1] */
  MPI_Comm row_comm; /* my processor's row communicator */
  MPI_Comm col_comm; /* my processor's column communicator */
  unsigned int RNGState_row;  /* State for random-number generator w.r.t.
			       row_comm */
  unsigned int RNGState_col;  /* State for random-number generator w.r.t.
			       row_comm */
  ZZ  *zz;        /* for debugging purpose */
};

typedef struct PHGCommStruct PHGComm;

void
Zoltan_PHGComm_Destroy(PHGComm* comm);

void
Zoltan_PHGComm_Init(PHGComm* comm);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_PHG_COMM_H */
