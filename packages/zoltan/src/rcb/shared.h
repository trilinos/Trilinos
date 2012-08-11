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


#ifndef __SHARED_CONST_H
#define __SHARED_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


/* Definitions shared by parallel RCB and RIB */
#define DEFAULT_CHECK_GEOM 1
#define TINY 1.0e-6
#define RB_MAX_WGTS 4 

/* Data structures shared by parallel RCB and RIB */
/* dot to balance on for RCB and RIB */ 
#if 0
struct Dot_Struct {	        /* dot = point in 3-space */
  double  X[3];			/* location of dot */
  double  Weight[RB_MAX_WGTS];  /* weight of dot - if used must be >= 0 */
  int Proc;                     /* Processor ID for processor owning a dot.
                                   For now, we'll keep it with a dot, even 
                                   though the global and local ids are now
                                   stored separately.                       */
  int Input_Part;               /* Partition to which the dot is assigned
                                   initially (input).  */
  int Part;                     /* New partition to which the dot is 
                                   assigned.  */
  int Size;                     /* Migration size */
};
#else
struct Dot_Struct {	        /* dot = point in 3-space */
  double  *X, *Y, *Z;
  double *Weight;                /* Weight-tuples of dots, in dot-major order */
  int nWeights;                  /* number of weights per dot in Weight array */
  double uniformWeight;      /* if non-zero, then Weight is NULL and all weights are the same */
  int *Proc;                     /* Processor ID for processor owning a dot.
                                   For now, we'll keep it with a dot, even 
                                   though the global and local ids are now
                                   stored separately.                       */
  int *Input_Part;               /* Partition to which the dot is assigned
                                   initially (input).  */
  int *Part;                     /* New partition to which the dot is 
                                   assigned.  */
  int *Size;                     /* Migration size, NULL means all are size 1 */
};
#endif


extern int Zoltan_RB_Build_Structure(ZZ *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  struct Dot_Struct *, int *, int *, int *, int, double, int, int);

extern void Zoltan_RB_Print_All(ZZ *, ZOLTAN_ID_PTR , struct Dot_Struct *,
  int , int , ZOLTAN_ID_PTR , int *);

extern int Zoltan_RB_Send_Outgoing(ZZ *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  struct Dot_Struct *, int **, int *, int *, int *, int, int *, double, int,
  ZOLTAN_GNO_TYPE *, int, MPI_Comm, int, int, int, int);

extern int Zoltan_RB_Send_To_Part(ZZ *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  struct Dot_Struct *, int **, int *, int *, int *, int *, double, int,
  ZOLTAN_GNO_TYPE *, int);

extern int Zoltan_RB_Send_Dots(ZZ *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  struct Dot_Struct *, int **, int *, int, int *, int *, int, int *, double,
  int, ZOLTAN_GNO_TYPE *, int, MPI_Comm);

extern int Zoltan_RB_Send_Dots_less_memory(ZZ *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  struct Dot_Struct *, int **, int *, int, int *, int *, int, int *, double,
  int, ZOLTAN_GNO_TYPE *, int, MPI_Comm);

extern void Zoltan_Free_And_Reset_Dot_Structure(struct Dot_Struct *);

extern int Zoltan_RB_Remap(ZZ *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, 
  struct Dot_Struct *, int *, int *, int *, double, int , ZOLTAN_GNO_TYPE *, 
  int);

extern int Zoltan_RB_Return_Arguments(ZZ *, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, 
  struct Dot_Struct *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *, int **, int **, 
  int);

extern int Zoltan_RB_check_geom_input(ZZ *, struct Dot_Struct *, int);

extern int Zoltan_RB_check_geom_output(ZZ *, struct Dot_Struct *, float *,
  int, int, int, int, void *);

extern void Zoltan_RB_stats(ZZ *, double, struct Dot_Struct *, int , 
                            float *, double *, ZOLTAN_GNO_TYPE *, int, 
                            ZOLTAN_GNO_TYPE *, void *, int);

extern int Zoltan_RB_Use_IDs(ZZ *);

extern int Zoltan_RB_Tree_Gatherv(ZZ *, int, int *, int *, int *);

extern int Zoltan_RB_Candidates_Copy_Input(ZZ *, int, ZOLTAN_ID_PTR,
  ZOLTAN_ID_PTR, struct Dot_Struct *, int *, ZOLTAN_ID_PTR *, ZOLTAN_ID_PTR *,
  int **, int **);

extern int Zoltan_RB_Candidates_Output(ZZ *, int, int *, 
  ZOLTAN_ID_PTR, ZOLTAN_ID_PTR,
  struct Dot_Struct *, int, ZOLTAN_ID_PTR, int *, ZOLTAN_ID_PTR *);



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
