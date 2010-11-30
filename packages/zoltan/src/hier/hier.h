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

#ifndef __HIER_H
#define __HIER_H

#include "third_library_const.h"
#include "zoltan_dd.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/* header file for hierarchical balancing */

/* Parameters to hierarchical balancing */
struct HierPartParamsStruct {
  int output_level;                  /* flag indicating amount of debugging
					output from hierarchical balancing.
					See levels in HIER_DEBUG_* below */
  int num_levels;                    /* number of levels of hierarchy to be
					dealt with at local processor */
  int global_num_levels;             /* max num_levels across all procs */
  int level;                         /* level currently being processed */
  MPI_Comm hier_comm;                /* MPI communicator for each level */
  int *hier_ranks_of_orig;           /* orig ranks of the procs in hier_comm */
 
  ZZ *origzz;                        /* Zoltan struct passed into top level */
  ZZ *hierzz;                        /* internal zoltan struct for balancing 
					within the hierarchy */
  int part_to_compute;               /* partition to compute at each level */
  int num_parts;                     /* number of partitions to compute */
  int use_geom, use_graph;           /* flags for whether methods to be
					used will require geometric
					and/or graph information */
  int checks;                        /* should we do sanity checks? */
  int init_num_obj;                  /* number of local objects at start */
  int hier_num_obj;                  /* number of local objects during
					hierarchical balancing */
  ZOLTAN_ID_PTR local_ids;           /* / lists of local ids and global ids */
  ZOLTAN_ID_PTR global_ids;          /* \ for initial partitioning */
  char *migrated;                    /* array of flags indicating whether
					each gid in global_ids has been
					migrated somewhere else */
  struct Zoltan_DD_Struct *dd;       /* distributed data to track migrated 
					gnos during hierarchical balancing */
  int allocsize_gnos_of_interest;    /* size of gnos_of_interest array */
  int num_gnos_of_interest;          /* num gnos in gnos_of_interest */
  ZOLTAN_GNO_TYPE *gnos_of_interest;    /* list of gnos of interest, used
					when looking up remote proc locations
					for graph edge callbacks */
  int *gnos_of_interest_procs;       /* list of procs where gnos of interest
					are located */
  int obj_wgt_dim, edge_wgt_dim;     /* object and edge weight dimensions */
  float *vwgt;                       /* vector of vertex weights */
  int *input_parts;                  /* Initial partitions for objects. */
  ZOLTAN_GNO_TYPE *vtxdist;          /* intermediate graph structure storage */
  int *xadj;                         /* intermediate graph structure storage */
  ZOLTAN_GNO_TYPE *adjncy;              /*    see Zoltan_Build_Graph */
  float *ewgts;                      /* edge weights for intermediate struct */
  int *adjproc;                      /* adjacent proc for each edge */
  int ndims;                         /* number of dimensions for geom data */
  int num_edges;                     /* number of edges in graph rep */
  double *geom_vec;                  /* geometry of objects in intermediate */
  int num_migrated_in_gnos;          /* number of gnos migrated to this proc */
  int alloc_migrated_in_gnos;        /* size of allocated array of migrated
					in gnos */
  ZOLTAN_GNO_TYPE *migrated_in_gnos;    /* ordered array of gnos migrated in */
  void **migrated_in_data;           /* data migrated in, parallel array to
					migrated_in_gnos */

                                         /* If non-null, Zoltan is computing the */
  zoltan_platform_specification *spec;   /* levels based on network topology */
  int gen_files;                      /* call Zoltan_Generate_Files */
                                         
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
