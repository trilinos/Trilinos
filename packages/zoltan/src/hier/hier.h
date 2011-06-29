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
  int output_level;                  /* amount of debugging info */
  int checks;                        /* should we do sanity checks? */
  int gen_files;                      /* call Zoltan_Generate_Files */

  int num_levels;                    /* number of levels I do */
  int global_num_levels;             /* max number of levels */
  int level;                         /* level currently being processed */
  MPI_Comm hier_comm;                /* MPI communicator for each level */
 
  ZZ *origzz;                        /* Zoltan struct passed into top level */
  ZZ *hierzz;                        /* internal zoltan struct for balancing 
					within the hierarchy */

  int part_to_compute;               /* partition to compute at each level */
  int num_parts;                     /* number of partitions to compute */

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
