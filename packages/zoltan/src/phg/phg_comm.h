/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */
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
