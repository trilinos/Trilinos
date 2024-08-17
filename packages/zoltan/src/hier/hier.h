// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __HIER_H
#define __HIER_H

#include "zoltan_dd.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* header file for hierarchical balancing */

#define ZOLTAN_PLATFORM_MAX_LEVELS 8

typedef struct _spec{
  /*
   * name of predefined topologies, or null if topology given by parameter
   */
  char *platform_name;

  /*
   * size of num_siblings and my_part arrays
   */
  int numLevels;

  /*
   * number of objects (cores, caches, sockets, etc), or number of
   *  children of the parent, of this level
   */
  int num_siblings[ZOLTAN_PLATFORM_MAX_LEVELS];

  /*
   * the part computed by this process at this level
   */
  int my_part[ZOLTAN_PLATFORM_MAX_LEVELS];
} zoltan_platform_specification;


/* Parameters to hierarchical balancing */
struct HierPartParamsStruct {
  int output_level;                  /* amount of debugging info */
  int checks;                        /* should we do sanity checks? */
  int gen_files;                      /* call Zoltan_Generate_Files */

  int num_levels;                    /* number of levels I do */
  int level;                         /* level currently being processed */
  MPI_Comm hier_comm;                /* MPI communicator for each level */
 
  ZZ *origzz;                        /* Zoltan struct passed into top level */
  ZZ *hierzz;                        /* internal zoltan struct for balancing 
					within the hierarchy */

  int part_to_compute;               /* part to compute at each level */
  int num_parts;                     /* number of parts to compute */

  int use_geom, use_graph;           /* flags for whether methods to be
					used will require geometric
					and/or graph information */
  int num_obj;                       /* number of local objects at start */
  int obj_wgt_dim, edge_wgt_dim;     /* object and edge weight dimensions */
  ZOLTAN_GNO_TYPE invalid_gno;  /* a value guaranteed not to be one of the gnos */
  ZOLTAN_GNO_TYPE *gno;              /* global vertex ids */
  float *vwgt;                       /* vector of vertex weights */
  int *xadj;                         /* intermediate graph structure storage */
  ZOLTAN_GNO_TYPE *adjncy;              /*    see Zoltan_Build_Graph */
  float *ewgts;                      /* edge weights for intermediate struct */
  int *adjproc;                      /* adjacent proc (in current MPI group) */

  int ndims;                         /* number of dimensions for geom data */
  double *geom_vec;                  /* geometry of objects in intermediate */
  int use_timers;                    /* control degree of timing done with hier*/

  zoltan_platform_specification *spec;   /* levels based on network topology */
};
typedef struct HierPartParamsStruct HierPartParams;

/*
 * Hierarchical balancing output levels.
 */
#define HIER_DEBUG_NONE 0
#define HIER_DEBUG_LIST 1
#define HIER_DEBUG_ALL  2
#define HIER_DEBUG_PRINT 3

/* Macro for error handling */
#define ZOLTAN_HIER_ERROR(error,str) {ierr = error ; \
 ZOLTAN_PRINT_ERROR(zz->Proc, yo, str) ; goto End ;}

/* prototype for set_param function needed by params/set_param.c */
extern int Zoltan_Hier_Set_Param(char *name, char *val);
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
