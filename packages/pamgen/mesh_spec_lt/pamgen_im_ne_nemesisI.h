// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

/****************************************************************************
 * This file contains prototypes for the functions found in the NEMESIS
 * library.
 ****************************************************************************/

#ifndef _IM_NE_NEMESIS_H
#define _IM_NE_NEMESIS_H

#ifdef __cplusplus
extern "C" {
#endif

/*=============================================================================
 *     Initial Information Routines
 *===========================================================================*/
extern int
im_ne_get_init_info(int   neid,		/* NemesisI file ID */
                 int  *num_proc,	/* Number of processors */
                 int  *num_proc_in_f,	/* Number of procs in this file */
                 char *ftype
                 );

extern int
im_ne_get_init_global(int   neid, 		  /* NemesisI file ID */
                   int  *num_nodes_g,	  /* Number of global FEM nodes */
                   int  *num_elems_g,	  /* Number of global FEM elements */
                   int  *num_elem_blks_g, /* Number of global elem blocks */
                   int  *num_node_sets_g, /* Number of global node sets */
                   int  *num_side_sets_g  /* Number of global side sets */
                   );

/*=============================================================================
 *     Loadbalance Parameter Routines
 *===========================================================================*/
extern int
im_ne_get_loadbal_param(int   neid, 	/* NetCDF/Exodus file ID */
                     int  *num_int_nodes,  /* Number of internal FEM nodes */
                     int  *num_bor_nodes,  /* Number of border FEM nodes */
                     int  *num_ext_nodes,  /* Number of external FEM nodes */
                     int  *num_int_elems,  /* Number of internal FEM elems */
                     int  *num_bor_elems,  /* Number of border FEM elems */
                     int  *num_node_cmaps, /* Number of nodal comm maps */
                     int  *num_elem_cmaps, /* Number of elemental comm maps */
                     int   processor         /* Processor ID */
                     );


/*=============================================================================
 *     NS, SS & EB Global Parameter Routines
 *===========================================================================*/
extern int
im_ne_get_ns_param_global(int neid,	     /* NetCDF/Exodus file ID */
                       int *ns_ids_glob,     /* Global IDs of node sets */
                       int *ns_n_cnt_glob,   /* Count of nodes in node sets */
                       int *ns_df_cnt_glob   /* Count of dist. factors in ns */
                       );

extern int
im_ne_get_ss_param_global(int neid,	    /* NetCDF/Exodus file ID */
                       int *ss_ids_glob,    /* Global side-set IDs */
                       int *ss_s_cnt_glob,  /* Global side count */
                       int *ss_df_cnt_glob  /* Global dist. factor count */
                       );

extern int
im_ne_get_eb_info_global(int neid,		/* NemesisI file ID                 */
                      int *el_blk_ids,	/* Vector of global element IDs     */
                      int *el_blk_cnts	/* Vector of global element counts  */
                      );


/*=============================================================================
 *     Number Map Routines
 *===========================================================================*/

extern int
im_ne_get_node_map(int   neid,		/* NetCDF/Exodus file ID */
                int  *node_mapi,	/* Internal FEM node IDs */
                int  *node_mapb,	/* Border FEM node IDs */
                int  *node_mape,	/* External FEM node IDs */
                int   processor		/* Processor IDs */
                );

extern int
im_ne_get_elem_map(int   neid,		/* NetCDF/Exodus file ID */
                int  *elem_mapi,	/* Internal element IDs */
                int  *elem_mapb,	/* Border element IDs */
                int   processor		/* Processor ID */
                );


/*=============================================================================
 *     Communications Maps Routines
 *===========================================================================*/

extern int
im_ne_get_cmap_params(int neid,                  /* NetCDF/Exodus file ID */
                   int *node_cmap_ids,        /* Nodal comm. map IDs */
                   int *node_cmap_node_cnts,  /* Number of nodes in each map */
                   int *elem_cmap_ids,        /* Elemental comm. map IDs */
                   int *elem_cmap_elem_cnts,  /* Number of elems in each map */
                   int  processor             /* This processor ID */
                   );


extern int
im_ne_get_node_cmap(int  neid,             /* NetCDF/Exodus file ID */
                 int  map_id,           /* Map ID */
                 int *node_ids,         /* FEM node IDs */
                 int *proc_ids,         /* Processor IDs */
                 int  processor         /* This processor ID */
                 );

extern int
im_ne_get_elem_cmap(int  neid,     /* NetCDF/Exodus file ID */
                 int  map_id,   /* Elemental comm map ID */
                 int *elem_ids, /* Element IDs */
                 int *side_ids, /* Element side IDs */
                 int *proc_ids, /* Processor IDs */
                 int  processor /* This processor ID */
                 );



#ifdef __cplusplus
}
#endif

#endif /* _IM_NE_NEMESIS_H */
