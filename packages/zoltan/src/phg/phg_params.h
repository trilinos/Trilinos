// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __PHG_PARAMS_H
#define __PHG_PARAMS_H
 
#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "params_const.h"
#include "zz_const.h"

/*  Parameters structure for parallel HG method.  */
static PARAM_VARS PHG_params[] = {
  /* Add parameters here. */
  {"HYPERGRAPH_PACKAGE",              NULL,  "STRING", 0},
    /* Software package: PHG (Zoltan) or Patoh */
  {"PHG_MULTILEVEL",                  NULL,  "INT", 0},
    /* Indicate whether or not to use multilevel method (1/0) */
  {"PHG_CUT_OBJECTIVE",               NULL,  "STRING", 0},
    /* connectivity or hyperedges */
  {"PHG_OUTPUT_LEVEL",                NULL,  "INT",    0},
    /* Higher value -> more output */
  {"FINAL_OUTPUT",                    NULL,  "INT",    0},
    /* Output related to final partitioning only */ 
  {"CHECK_GRAPH",                     NULL,  "INT",    0},
    /* Same as CHECK_HYPERGRAPH */
  {"CHECK_HYPERGRAPH",                NULL,  "INT",    0},
    /* Same as CHECK_GRAPH */
  {"PHG_NPROC_VERTEX",                NULL,  "INT",    0},
    /* No. of processors along vertex direction initially in 2D layout */
  {"PHG_NPROC_EDGE",                  NULL,  "INT",    0},
    /* No. of processors along hyperedge direction initially in 2D layout */
  {"PHG_COARSENING_LIMIT",            NULL,  "INT",    0},
    /* When to stop coarsening (global no. of vertices) */
  {"PHG_COARSENING_NCANDIDATE",       NULL,  "INT",    0},  
    /* Max no. of candidate vertices in a round of matching */
  {"PHG_COARSENING_METHOD",           NULL,  "STRING", 0},
    /* Coarsening method (ipm, agg, etc. ) */
  {"PHG_COARSENING_METHOD_FAST",      NULL,  "STRING", 0},
    /* Used by a-ipm to alternate between full ipm and a faster method */
  {"PHG_VERTEX_VISIT_ORDER",          NULL,  "INT",    0},
    /* Vertex ordering for greedy matching (linear, random,  ) */
  {"PHG_EDGE_SCALING",                NULL,  "INT",    0},
    /* Edge scaling schemes to tweak inner product similarity in matching */
  {"PHG_VERTEX_SCALING",              NULL,  "INT",    0},
    /* Vertex scaling schemes to tweak inner product similarity in matching */
  {"PHG_COARSEPARTITION_METHOD",      NULL,  "STRING", 0},
    /* Coarse partitioning method: linear, random, greedy, auto */
  {"PHG_REFINEMENT_METHOD",           NULL,  "STRING", 0},
    /* Only 2-way FM (fm2) for now */
  {"PHG_DIRECT_KWAY",                 NULL,  "INT",    0},
    /* Direct k-way partitioning not yet implemented! */
  {"PHG_REFINEMENT_LOOP_LIMIT",       NULL,  "INT",    0},
    /* Max no. of loops in KL/FM. */
  {"PHG_REFINEMENT_MAX_NEG_MOVE",     NULL,  "INT",    0},    
    /* Max. no. of negative moves allowed before exiting refinement. */
  {"PHG_REFINEMENT_QUALITY",          NULL,  "FLOAT",  0},
    /* 1.0 is default; higher (lower) value gives more (less) refinement. */
  {"PHG_USE_TIMERS",                  NULL,  "INT",    0},    
    /* Same as USE_TIMERS. */
  {"USE_TIMERS",                      NULL,  "INT",    0},    
    /* Same as PHG_USE_TIMERS. */
  {"PHG_EDGE_SIZE_THRESHOLD",         NULL,  "FLOAT",  0},
    /* Ignore hyperedges larger than this threshold times nvertex */
    /* If PHG_EDGE_SIZE_THRESHOLD>1, interpret it as absolute value. */
  {"PHG_MATCH_EDGE_SIZE_THRESHOLD",   NULL,  "INT",    0},
    /* Ignore hyperedges larger than this threshold, in local processor, during matching */
  {"PHG_BAL_TOL_ADJUSTMENT",          NULL,  "FLOAT",  0},  
    /* Adjustment factor for balance in recursive bisection. */
  {"PHG_EDGE_WEIGHT_OPERATION",       NULL,  "STRING",  0},
    /* How to handle inconsistent edge weights across processors */
  {"PARKWAY_SERPART",                 NULL,  "STRING", 0},
    /* Serial partitioner for ParKway (PaToH or hMetis) */
  {"ADD_OBJ_WEIGHT",                  NULL,  "STRING", 0},
    /* Add implicit vertex weight, like no. of pins (nonzeros)? */
  {"PHG_RANDOMIZE_INPUT",             NULL,  "INT",    0},    
    /* Randomizing input often improves load balance within PHG but destroys 
       locality, so may produce lower quality partitions  */
  {"PHG_PROCESSOR_REDUCTION_LIMIT",   NULL,  "FLOAT",  0},
    /* When to move data to fewer processors. */
  {"PHG_REPART_MULTIPLIER",           NULL,  "FLOAT",  0},
    /* Multiplier for communication to migration trade-off in repartitioning. */
  {"HYBRID_REDUCTION_FACTOR",        NULL,  "FLOAT",    0},
    /* Factor by which to reduce the # of vtx when using geometric matching. */
  {"HYBRID_REDUCTION_LEVELS",        NULL,  "INT",    0},
    /* # of levels on which to use geometric matching; INT_MAX-->all levels. */
  {"PATOH_ALLOC_POOL0",               NULL,  "INT",    0},
    /* Memory allocation parameter for Patoh. */
  {"PATOH_ALLOC_POOL1",               NULL,  "INT",    0},   
  /* Memory allocation parameter for Patoh. */
  {"PHG_KEEP_TREE",                   NULL,  "INT",    0},
  /* Keep bisection tree */
  {NULL,                              NULL,  NULL,     0}     
};


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
