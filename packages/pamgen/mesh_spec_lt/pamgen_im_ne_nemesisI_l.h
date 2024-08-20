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
im_ne_get_init_info_l(int   neid,		/* NemesisI file ID */
                 long long  *num_proc,	/* Number of processors */
                 long long  *num_proc_in_f,	/* Number of procs in this file */
                 char *ftype
                 );

extern int
im_ne_get_init_global_l(int   neid, 		  /* NemesisI file ID */
                   long long  *num_nodes_g,	  /* Number of global FEM nodes */
                   long long  *num_elems_g,	  /* Number of global FEM elements */
                   long long  *num_elem_blks_g, /* Number of global elem blocks */
                   long long  *num_node_sets_g, /* Number of global node sets */
                   long long  *num_side_sets_g  /* Number of global side sets */
                   );

/*=============================================================================
 *     Loadbalance Parameter Routines
 *===========================================================================*/
extern int
im_ne_get_loadbal_param_l(int   neid, 	/* NetCDF/Exodus file ID */
                     long long  *num_int_nodes,  /* Number of internal FEM nodes */
                     long long  *num_bor_nodes,  /* Number of border FEM nodes */
                     long long  *num_ext_nodes,  /* Number of external FEM nodes */
                     long long  *num_int_elems,  /* Number of internal FEM elems */
                     long long  *num_bor_elems,  /* Number of border FEM elems */
                     long long  *num_node_cmaps, /* Number of nodal comm maps */
                     long long  *num_elem_cmaps, /* Number of elemental comm maps */
                     long long   processor         /* Processor ID */
                     );


/*=============================================================================
 *     NS, SS & EB Global Parameter Routines
 *===========================================================================*/
extern int
im_ne_get_ns_param_global_l(int neid,	     /* NetCDF/Exodus file ID */
                       long long *ns_ids_glob,     /* Global IDs of node sets */
                       long long *ns_n_cnt_glob,   /* Count of nodes in node sets */
                       long long *ns_df_cnt_glob   /* Count of dist. factors in ns */
                       );

extern int
im_ne_get_ss_param_global_l(int neid,	    /* NetCDF/Exodus file ID */
                       long long *ss_ids_glob,    /* Global side-set IDs */
                       long long *ss_s_cnt_glob,  /* Global side count */
                       long long *ss_df_cnt_glob  /* Global dist. factor count */
                       );

extern int
im_ne_get_eb_info_global_l(int neid,		/* NemesisI file ID                 */
                      long long *el_blk_ids,	/* Vector of global element IDs     */
                      long long *el_blk_cnts	/* Vector of global element counts  */
                      );


/*=============================================================================
 *     Number Map Routines
 *===========================================================================*/

extern int
im_ne_get_node_map_l(int   neid,		/* NetCDF/Exodus file ID */
                long long  *node_mapi,	/* Internal FEM node IDs */
                long long  *node_mapb,	/* Border FEM node IDs */
                long long  *node_mape,	/* External FEM node IDs */
                long long   processor		/* Processor IDs */
                );

extern int
im_ne_get_elem_map_l(int   neid,		/* NetCDF/Exodus file ID */
                long long  *elem_mapi,	/* Internal element IDs */
                long long  *elem_mapb,	/* Border element IDs */
                long long   processor		/* Processor ID */
                );


/*=============================================================================
 *     Communications Maps Routines
 *===========================================================================*/

extern int
im_ne_get_cmap_params_l(int neid,                  /* NetCDF/Exodus file ID */
                   long long *node_cmap_ids,        /* Nodal comm. map IDs */
                   long long *node_cmap_node_cnts,  /* Number of nodes in each map */
                   long long *elem_cmap_ids,        /* Elemental comm. map IDs */
                   long long *elem_cmap_elem_cnts,  /* Number of elems in each map */
                   long long  processor             /* This processor ID */
                   );


extern int
im_ne_get_node_cmap_l(int  neid,             /* NetCDF/Exodus file ID */
                 long long  map_id,           /* Map ID */
                 long long *node_ids,         /* FEM node IDs */
                 long long *proc_ids,         /* Processor IDs */
                 long long  processor         /* This processor ID */
                 );

extern int
im_ne_get_elem_cmap_l(int  neid,     /* NetCDF/Exodus file ID */
                 long long  map_id,   /* Elemental comm map ID */
                 long long *elem_ids, /* Element IDs */
                 long long *side_ids, /* Element side IDs */
                 long long *proc_ids, /* Processor IDs */
                 long long  processor /* This processor ID */
                 );


/*=============================================================================
 *     Get IJK
 *===========================================================================*/
extern int
im_ne_get_global_ijk_l(int neid,
                long long *global_ijk
                );
extern int
im_ne_get_num_ijk_l(int neid,
                long long *num_ijk
                );
extern int
im_ne_get_local_num_ijk_l(int neid,
                long long *local_num_ijk
                );

#ifdef __cplusplus
}
#endif

#endif /* _IM_NE_NEMESIS_H */
