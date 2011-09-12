/*
 * Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.  
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile: ne_nemesisI.h,v $
 *
 * $Author: gdsjaar $
 *
 * $Date: 2008/02/07 13:40:50 $
 *
 * $Revision: 1.17 $
 *
 * $Name:  $
 *====================================================================*/

/****************************************************************************
 * This file contains prototypes for the functions found in the NEMESIS
 * library.
 ****************************************************************************/

#ifndef _NE_NEMESIS_H
#define _NE_NEMESIS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define NEMESIS_API_VERSION		3.12
#define NEMESIS_API_VERSION_NODOT	312

/*=============================================================================
 *     Initial Information Routines
 *===========================================================================*/
extern int
ne_get_init_info(int   neid,		/* NemesisI file ID */
                 int  *num_proc,	/* Number of processors */
                 int  *num_proc_in_f,	/* Number of procs in this file */
                 char *ftype
                 );

extern int
ne_put_init_info(int   neid,		/* NemesisI file ID */
                 int   num_proc,	/* Number of processors */
                 int   num_proc_in_f,	/* Number of procs in this file */
                 char *ftype
                 );

extern int
ne_get_init_global(int   neid, 		  /* NemesisI file ID */
                   int  *num_nodes_g,	  /* Number of global FEM nodes */
                   int  *num_elems_g,	  /* Number of global FEM elements */
                   int  *num_elem_blks_g, /* Number of global elem blocks */
                   int  *num_node_sets_g, /* Number of global node sets */
                   int  *num_side_sets_g  /* Number of global side sets */
                   );
extern int
ne_put_init_global(int neid, 		/* NemesisI file ID */
                   int num_nodes_g,	/* Number of global FEM nodes */
                   int num_elems_g,	/* Number of global FEM elements */
                   int num_elem_blks_g,	/* Number of global elem blocks */
                   int num_node_sets_g,	/* Number of global node sets */
                   int num_side_sets_g	/* Number of global side sets */
                   );

/*=============================================================================
 *     Loadbalance Parameter Routines
 *===========================================================================*/
extern int
ne_get_loadbal_param(int   neid, 	/* NetCDF/Exodus file ID */
                     int  *num_int_nodes,  /* Number of internal FEM nodes */
                     int  *num_bor_nodes,  /* Number of border FEM nodes */
                     int  *num_ext_nodes,  /* Number of external FEM nodes */
                     int  *num_int_elems,  /* Number of internal FEM elems */
                     int  *num_bor_elems,  /* Number of border FEM elems */
                     int  *num_node_cmaps, /* Number of nodal comm maps */
                     int  *num_elem_cmaps, /* Number of elemental comm maps */
                     int   processor         /* Processor ID */
                     );

extern int
ne_put_loadbal_param(int   neid, 	  /* NemesisI file ID  */
                     int   num_int_nodes, /* Number of internal FEM nodes */
                     int   num_bor_nodes, /* Number of border FEM nodes */
                     int   num_ext_nodes, /* Number of external FEM nodes */
                     int   num_int_elems, /* Number of internal FEM elems */
                     int   num_bor_elems, /* Number of border FEM elems */
                     int   num_node_cmaps,/* Number of nodal comm maps */
                     int   num_elem_cmaps,/* Number of elemental comm maps */
                     int   processor	  /* Processor ID */
                     );

extern int
ne_put_loadbal_param_cc(int   neid,		/* NetCDF/Exodus file ID */
                        int  *num_int_nodes,  /* Number of internal node IDs */
                        int  *num_bor_nodes,  /* Number of border node IDs */
                        int  *num_ext_nodes,  /* Number of external node IDs */
                        int  *num_int_elems,  /* Number of internal elem IDs */
                        int  *num_bor_elems,  /* Number of border elem IDs */
                        int  *num_node_cmaps, /* Number of nodal comm maps */
                        int  *num_elem_cmaps  /* Number of elem comm maps */
                        );

/*=============================================================================
 *     NS, SS & EB Global Parameter Routines
 *===========================================================================*/
extern int
ne_get_ns_param_global(int neid,	     /* NetCDF/Exodus file ID */
                       int *ns_ids_glob,     /* Global IDs of node sets */
                       int *ns_n_cnt_glob,   /* Count of nodes in node sets */
                       int *ns_df_cnt_glob   /* Count of dist. factors in ns */
                       );

extern int
ne_put_ns_param_global(int neid, 	    /* NemesisI file ID */
                       int *global_ids,	    /* Vector of global node-set IDs */
                       int *global_n_cnts,  /* Vector of node counts in */
                                            /* node-sets */
                       int *global_df_cnts  /* Vector of dist factor */
                                            /* counts in node-sets */
                       );

extern int
ne_get_ss_param_global(int neid,	    /* NetCDF/Exodus file ID */
                       int *ss_ids_glob,    /* Global side-set IDs */
                       int *ss_s_cnt_glob,  /* Global side count */
                       int *ss_df_cnt_glob  /* Global dist. factor count */
                       );

extern int
ne_put_ss_param_global(int neid, 	    /* NemesisI file ID */
                       int *global_ids,	    /* Vector of global side-set IDs */
                       int *global_el_cnts, /* Vector of element/side */
					    /* counts in each side set */
                       int *global_df_cnts  /* Vector of dist. factor */
					    /* counts in each side set */
                       );

extern int
ne_get_eb_info_global(int neid,		/* NemesisI file ID                 */
                      int *el_blk_ids,	/* Vector of global element IDs     */
                      int *el_blk_cnts	/* Vector of global element counts  */
                      );

extern int
ne_put_eb_info_global(int neid,		/* NemesisI file ID */
                      int *el_blk_ids,	/* Vector of global element IDs     */
                      int *el_blk_cnts	/* Vector of global element counts  */
                      );

/*=============================================================================
 *     NS, SS & EB Subset Routines
 *===========================================================================*/
extern int
ne_get_n_side_set(int  neid,		    /* NetCDF/Exodus file ID */
                  int  side_set_id,	    /* Side-set ID to read */
                  int  start_side_num,      /* Starting element number */
                  int  num_sides,	    /* Number of sides to read */
                  int *side_set_elem_list,  /* List of element IDs */
                  int *side_set_side_list   /* List of side IDs */
                  );

extern int
ne_put_n_side_set(int  neid,                /* NetCDF/Exodus file ID */
                  int  side_set_id,         /* Side-set ID to write */
                  int  start_side_num,      /* Starting element number */
                  int  num_sides,           /* Number of sides to write */
                  int *side_set_elem_list,  /* List of element IDs */
                  int *side_set_side_list   /* List of side IDs */
                  );

extern int
ne_get_n_side_set_df(int   neid,		/* NetCDF/Exodus file ID */
                     int   side_set_id,		/* Side-set ID */
                     int   start_num,		/* Starting df number */
                     int   num_df_to_get,	/* Number of df's to read */
                     void *side_set_df 		/* Distribution factors */
                     );

extern int
ne_put_n_side_set_df(int   neid,                /* NetCDF/Exodus file ID */
                     int   side_set_id,         /* Side-set ID */
                     int   start_num,           /* Starting df number */
                     int   num_df_to_get,       /* Number of df's to write */
                     void *side_set_df          /* Distribution factors */
                     );

extern int
ne_get_n_node_set(int  neid,		   /* NetCDF/Exodus file ID */
                  int  node_set_id,	   /* Node set ID */
                  int  start_node_num,	   /* Node index to start reading at */
                  int  num_node,	   /* Number of nodes to read */
                  int *node_set_node_list  /* List of nodes in node set */
                  );

extern int
ne_put_n_node_set(int  neid,		   /* NetCDF/Exodus file ID */
                  int  node_set_id,	   /* Node set ID */
                  int  start_node_num,	   /* Node index to start writing at */
                  int  num_node,	   /* Number of nodes to write */
                  int *node_set_node_list  /* List of nodes in node set */
                  );

extern int
ne_get_n_node_set_df(int   neid,		/* NetCDF/Exodus file ID */
                     int   node_set_id,		/* Node-set ID */
                     int   start_num,		/* Starting df number */
                     int   num_df_to_get,	/* Number of df's to read */
                     void *node_set_df 		/* Distribution factors */
                     );

extern int
ne_put_n_node_set_df(int   neid,		/* NetCDF/Exodus file ID */
                     int   node_set_id,		/* Node-set ID */
                     int   start_num,		/* Starting df number */
                     int   num_df_to_get,	/* Number of df's to write */
                     void *node_set_df 		/* Distribution factors */
                     );

extern int
ne_get_n_coord(int   neid,		/* NetCDF/Exodus file ID */
               int   start_node_num,	/* Starting position to read from */
               int   num_nodes,		/* Number of coords to read */
               void *x_coor,		/* Vector of X coordinates */
               void *y_coor,		/* Vector of Y coordinates */
               void *z_coor		/* Vector of Z coordinates */
               );

extern int
ne_put_n_coord(int   neid,              /* NetCDF/Exodus file ID */
               int   start_node_num,    /* Starting position to write to */
               int   num_nodes,         /* Number of coords to write */
               void *x_coor,            /* Vector of X coordinates */
               void *y_coor,            /* Vector of Y coordinates */
               void *z_coor             /* Vector of Z coordinates */
               );

extern int
ne_get_n_elem_conn (int   neid,		  /* NetCDF/Exodus file ID */
                    int   elem_blk_id,	  /* Element block ID */
                    int   start_elem_num, /* Starting position to read from */
                    int   num_elems,	  /* Number of elements to read */
                    int  *connect	  /* Connectivity vector */
                    );

extern int
ne_put_n_elem_conn (int   neid,           /* NetCDF/Exodus file ID */
                    int   elem_blk_id,    /* Element block ID */
                    int   start_elem_num, /* Starting position to write to */
                    int   num_elems,      /* Number of elements to write */
                    int  *connect         /* Connectivity vector */
);

extern int
ne_get_n_elem_attr (int   neid,		   /* NetCDF/Exodus file ID */
                    int   elem_blk_id,	   /* Element block ID */
                    int   start_elem_num,  /* Starting position to read from */
                    int   num_elems,	   /* Number of elements to read */
                    void *attrib	   /* Attribute */
                    );

extern int
ne_put_n_elem_attr (int   neid,            /* NetCDF/Exodus file ID */
                    int   elem_blk_id,     /* Element block ID */
                    int   start_elem_num,  /* Starting position to write to */
                    int   num_elems,       /* Number of elements to write */
                    void *attrib           /* Attribute */
                    );

extern int
ne_get_elem_type(int   neid,            /* NetCDF/Exodus file ID */
                 int   elem_blk_id,     /* Element block ID */
                 char *elem_type        /* The name of the element type */
                 );

/*=============================================================================
 *     Variable Routines
 *===========================================================================*/
extern int
ne_get_n_elem_var (int   neid,              /* NetCDF/Exodus file ID */
                   int   time_step,         /* time index */
                   int   elem_var_index,    /* elemental variable index */
                   int   elem_blk_id,       /* elemental block id */
                   int   num_elem_this_blk, /* number of elements in block */
                   int   start_elem_num,    /* Starting position for read */
                   int   num_elem,          /* Number of elements to read */
                   void *elem_var_vals      /* variable values */
                   );

extern int
ne_put_elem_var_slab (int   neid,           /* NetCDF/Exodus file ID */
                      int   time_step,      /* time index */
                      int   elem_var_index, /* elemental variable index */
                      int   elem_blk_id,    /* elemental block id */
                      int   start_pos,      /* Starting position to write to */
                      int   num_vals,       /* Number of elements to write */
                      void *elem_var_vals   /* variable values */
                      );

extern int
ne_get_n_nodal_var(int   neid,               /* NetCDF/Exodus file ID */
                   int   time_step,          /* whole time step number */
                   int   nodal_var_index,    /* index of desired nodal var */
                   int   start_node_num,     /* starting node number */
                   int   num_nodes,          /* number of nodes to read */
                   void *nodal_vars          /* array of nodal var values */
                   );

extern int
ne_put_nodal_var_slab(int   neid,            /* NetCDF/Exodus file ID */
                      int   time_step,       /* The time step index */
                      int   nodal_var_index, /* Nodal variable index */
                      int   start_pos,       /* Start position for write */
                      int   num_vals,        /* Number of nodal variables */
                      void *nodal_var_vals   /* Nodal variable values */
                      );

/*=============================================================================
 *     Number Map Routines
 *===========================================================================*/
extern int
ne_get_n_elem_num_map (int  neid,           /* NetCDF/Exodus file ID */
                       int  start_ent,      /* Starting position to read from */
                       int  num_ents,       /* Number of elements to read */
                       int *elem_map        /* element map numbers */
                       );

extern int
ne_put_n_elem_num_map (int  neid,           /* NetCDF/Exodus file ID */
                       int  start_ent,      /* Starting position to read from */
                       int  num_ents,       /* Number of elements to read */
                       int *elem_map        /* element map numbers */
                       );

extern int
ne_get_n_node_num_map(int   neid,	     /* NetCDF/Exodus file ID */
                      int   start_ent,       /* starting node number */
                      int   num_ents,        /* number of nodes to read */
                      int  *node_map         /* vector for node map */
                      );

extern int
ne_put_n_node_num_map(int   neid,	     /* NetCDF/Exodus file ID */
                      int   start_ent,       /* starting node number */
                      int   num_ents,        /* number of nodes to read */
                      int  *node_map         /* vector for node map */
                      );

extern int
ne_get_node_map(int   neid,		/* NetCDF/Exodus file ID */
                int  *node_mapi,	/* Internal FEM node IDs */
                int  *node_mapb,	/* Border FEM node IDs */
                int  *node_mape,	/* External FEM node IDs */
                int   processor		/* Processor IDs */
                );

extern int
ne_put_node_map(int   neid,		/* NetCDF/Exodus file ID */
                int  *node_mapi,	/* Internal FEM node IDs */
                int  *node_mapb,	/* Border FEM node IDs */
                int  *node_mape,	/* External FEM node IDs */
                int   processor		/* This processor ID */
                );

extern int
ne_get_elem_map(int   neid,		/* NetCDF/Exodus file ID */
                int  *elem_mapi,	/* Internal element IDs */
                int  *elem_mapb,	/* Border element IDs */
                int   processor		/* Processor ID */
                );

extern int
ne_put_elem_map(int   neid,		/* NetCDF/Exodus file ID */
                int  *elem_mapi,	/* Internal FEM element IDs */
                int  *elem_mapb,	/* Border FEM element IDs */
                int   processor		/* This processor ID */
                );


/*=============================================================================
 *     Communications Maps Routines
 *===========================================================================*/

extern int
ne_get_cmap_params(int neid,                  /* NetCDF/Exodus file ID */
                   int *node_cmap_ids,        /* Nodal comm. map IDs */
                   int *node_cmap_node_cnts,  /* Number of nodes in each map */
                   int *elem_cmap_ids,        /* Elemental comm. map IDs */
                   int *elem_cmap_elem_cnts,  /* Number of elems in each map */
                   int  processor             /* This processor ID */
                   );

extern int
ne_put_cmap_params(int  neid,			/* NetCDF/Exodus file ID */
                   int *node_map_ids,		/* Node map IDs */
                   int *node_map_node_cnts,	/* Nodes in nodal comm */
                   int *elem_map_ids,		/* Elem map IDs */
                   int *elem_map_elem_cnts,	/* Elems in elemental comm */
                   int  processor		/* This processor ID */
                   );

extern int
ne_put_cmap_params_cc(int  neid,		/* NetCDF/Exodus file ID */
                      int *node_map_ids,	/* Node map IDs */
                      int *node_map_node_cnts,	/* Nodes in nodal comm */
                      int *node_proc_ptrs,      /* Pointer into array for */
						/* node maps		  */
                      int *elem_map_ids,	/* Elem map IDs */
                      int *elem_map_elem_cnts,	/* Elems in elemental comm */
                      int *elem_proc_ptrs	/* Pointer into array for */
						/* elem maps		  */
                      );

extern int
ne_get_node_cmap(int  neid,             /* NetCDF/Exodus file ID */
                 int  map_id,           /* Map ID */
                 int *node_ids,         /* FEM node IDs */
                 int *proc_ids,         /* Processor IDs */
                 int  processor         /* This processor ID */
                 );

extern int
ne_put_node_cmap(int  neid,	/* NetCDF/Exodus file ID */
                 int  map_id,	/* Nodal comm map ID */
                 int *node_ids,	/* FEM node IDs */
                 int *proc_ids, /* Processor IDs */
                 int  processor	/* This processor ID */
                 );

extern int
ne_get_elem_cmap(int  neid,     /* NetCDF/Exodus file ID */
                 int  map_id,   /* Elemental comm map ID */
                 int *elem_ids, /* Element IDs */
                 int *side_ids, /* Element side IDs */
                 int *proc_ids, /* Processor IDs */
                 int  processor /* This processor ID */
                 );

extern int
ne_put_elem_cmap(int  neid,	/* NetCDF/Exodus file ID */
                 int  map_id,	/* Elemental comm map ID */
                 int *elem_ids,	/* Vector of element IDs */
                 int *side_ids, /* Vector of side IDs */
                 int *proc_ids, /* Vector of processor IDs */
                 int  processor	/* This processor ID */
                 );

/* Only used internally by NEMESIS */

extern int
ne_leavedef(int neid, 		/* NemesisI file ID         */
            char *func_name	/* Name of calling function */
            );

extern int
ne_get_file_type(int neid,	/* NetCDF/Exodus file ID */
                 char *ftype	/* Nemesis file type */
                 );

extern char *
ne_catstr2(char *name,	/* The name to attach num1 and num2 to */
           int   num1,	/* First number to tack to name */
           int   num2	/* Second number to tack to name */
           );

extern int
ne_id_lkup(int   neid,		/* NetCDF/Exodus file ID */
           char    *var_name,	/* Nemesis variable name */
           int64_t *idx,         /* index variable for variable, length 2 */
           int   ne_var_id	/* NetCDF variable ID */
           );

extern int
ne_get_map_status(int neid,
		  char *stat_var,
		  int proc_id,
		  int *stat);

extern int
ne_put_version(int neid		/* NetCDF/Exodus file ID */
               );

extern int
ne_check_file_version(int neid	/* NetCDF/Exodus file ID */
                      );

extern int
ne_get_idx(int      neid,	   /* NetCDF/Exodus file ID */
           char    *ne_var_name,   /* Nemesis index variable name */
           int64_t *index,	   /* array of length 2 to hold results */
           int      pos		   /* position of this proc/cmap in index */
           );

#ifdef __cplusplus
}
#endif

#endif /* _NE_NEMESIS_H */
