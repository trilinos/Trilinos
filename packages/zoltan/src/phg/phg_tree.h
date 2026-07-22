// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef __ZOLTAN_PHG_TREE_H
#define __ZOLTAN_PHG_TREE_H

#include "zz_const.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*********************************************************/
/* Structure for storing bisection tree for phg          */
/*********************************************************/
typedef struct Zoltan_PHG_Tree_ {
  int size;
  int * array;
} Zoltan_PHG_Tree;

/*********************************************************/
/* Structure for storing timer indices from call to call */
/*********************************************************/
struct phg_timer_indices {
  int all;  /* All of Zoltan_PHG; this timer includes other timers and their
               synchronization time, so it will be high. */
  int build;  /* hypergraph build time */
  int setupvmap;
  int parkway;  /* ParKway time */
  int patoh;  /* PaToH time */
  int retlist;  /* Building return lists time */
  int finaloutput;  /* printing final output time */
  int match;   /* Matching time */
  int coarse;  /* Coarsening time */
  int refine; /* Refinement time */
  int coarsepart;  /* Coarse partitioning time */
  int project;    /* Project coarse-to-fine */
  int procred;    /* Processor reduction */
  int vcycle;    /* Vcycle time */
  int comerge;   /* Part of coarsening */
  int coshuffle;   /* Part of coarsening */
  int coremove;   /* Part of coarsening */
  int cotheend;   /* Part of coarsening */
  int matchstage[7];  /* Matching stages */
  int rdrdivide;  /* Rdivide time. */
  int rdbefore;   /* Part of Rdivide */
  int rdafter;    /* Part of Rdivide */
  int rdsplit;    /* Part of Rdivide */
  int rdredist;    /* Part of Rdivide */
  int rdsend;    /* Part of Rdivide */
  int rdwait;    /* Part of Rdivide */
  int rfrefine;    /* Refinement time */
  int rfpins;   /* Part of Refinement */
  int rfiso;   /* Part of Refinement */
  int rfgain;   /* Part of Refinement */
  int rfheap;   /* Part of Refinement */
  int rfpass;   /* Part of Refinement */
  int rfroll;   /* Part of Refinement */
  int rfnonroot;   /* Part of Refinement */
  int cpart;   /* Coarse partitioning time */
  int cpgather;   /* Part of Coarse Partitioning */
  int cprefine;   /* Part of Coarse Partitioning */
};

typedef struct Zoltan_PHG_LB_Data_ {
  struct phg_timer_indices * timers;
  Zoltan_PHG_Tree * tree;
#ifdef CEDRIC_2D_PARTITIONS
  struct Zoltan_DD_Struct *ddHedge;
  int * partTree; /* Not used yet */
  int numParts;
  int *sizeParts;
#endif /* CEDRIC_2D_PARTITIONS */
} Zoltan_PHG_LB_Data;

struct phg_timer_indices *
Zoltan_PHG_LB_Data_timers(ZZ const * zz);

Zoltan_PHG_Tree *
Zoltan_PHG_LB_Data_tree(ZZ const * zz);

void
Zoltan_PHG_LB_Data_free_timers(ZZ* zz);

void
Zoltan_PHG_LB_Data_free_tree(ZZ* zz);

/* Build a centralized tree */
int
Zoltan_PHG_Tree_centralize(ZZ *zz);

int
Zoltan_PHG_Tree_create(int part_number, ZZ* zz);

void
Zoltan_PHG_Tree_Set(ZZ* zz, int father, int lo, int hi);


Zoltan_PHG_Tree *
get_tree(ZZ* zz);

int
get_tree_size(Zoltan_PHG_Tree * tree);

#define get_tree_size(tree) ((tree)->size)

/* Find interval (int [2]) in tree. Return the position of the smallest node
   that contains this interval, or -1 if none is found */
int
find_interval_in_tree(Zoltan_PHG_Tree *tree, int *interval);

int*
compute_part_number(Zoltan_PHG_Tree *tree);

int
Zoltan_PHG_Timers_init(ZZ* zz);


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __ZOLTAN_PHG_TREE_H */
