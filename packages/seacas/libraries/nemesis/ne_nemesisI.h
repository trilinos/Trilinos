/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/****************************************************************************
 * This file contains prototypes for the functions found in the NEMESIS
 * library.
 ****************************************************************************/

#ifndef _NE_NEMESIS_H
#define _NE_NEMESIS_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef EXODUS_II_HDR
#include <stdint.h>
/** Specifies that this argument is the id of an entity: element block, nodeset, sideset, ... */
typedef int64_t ex_entity_id;

/** The mechanism for passing double/float and int/int64_t both use a
   void*; to avoid some confusion as to whether a function takes an
   integer or a float/double, the following typedef is used for the
   integer argument
*/
typedef void void_int;
#endif

/*=============================================================================
 *     Initial Information Routines
 *===========================================================================*/
extern int ne_get_init_info(int   neid,          /* NemesisI file ID */
                            int * num_proc,      /* Number of processors */
                            int * num_proc_in_f, /* Number of procs in this file */
                            char *ftype);

extern int ne_put_init_info(int   neid,          /* NemesisI file ID */
                            int   num_proc,      /* Number of processors */
                            int   num_proc_in_f, /* Number of procs in this file */
                            char *ftype);

extern int ne_get_init_global(int       neid,            /* NemesisI file ID */
                              void_int *num_nodes_g,     /* Number of global FEM nodes */
                              void_int *num_elems_g,     /* Number of global FEM elements */
                              void_int *num_elem_blks_g, /* Number of global elem blocks */
                              void_int *num_node_sets_g, /* Number of global node sets */
                              void_int *num_side_sets_g  /* Number of global side sets */
);
extern int ne_put_init_global(int     neid,            /* NemesisI file ID */
                              int64_t num_nodes_g,     /* Number of global FEM nodes */
                              int64_t num_elems_g,     /* Number of global FEM elements */
                              int64_t num_elem_blks_g, /* Number of global elem blocks */
                              int64_t num_node_sets_g, /* Number of global node sets */
                              int64_t num_side_sets_g  /* Number of global side sets */
);

extern int ne_put_version(int neid);

/*=============================================================================
 *     Loadbalance Parameter Routines
 *===========================================================================*/
extern int ne_get_loadbal_param(int       neid,           /* NetCDF/Exodus file ID */
                                void_int *num_int_nodes,  /* Number of internal FEM nodes */
                                void_int *num_bor_nodes,  /* Number of border FEM nodes */
                                void_int *num_ext_nodes,  /* Number of external FEM nodes */
                                void_int *num_int_elems,  /* Number of internal FEM elems */
                                void_int *num_bor_elems,  /* Number of border FEM elems */
                                void_int *num_node_cmaps, /* Number of nodal comm maps */
                                void_int *num_elem_cmaps, /* Number of elemental comm maps */
                                int       processor       /* Processor ID */
);

extern int ne_put_loadbal_param(int     neid,           /* NemesisI file ID  */
                                int64_t num_int_nodes,  /* Number of internal FEM nodes */
                                int64_t num_bor_nodes,  /* Number of border FEM nodes */
                                int64_t num_ext_nodes,  /* Number of external FEM nodes */
                                int64_t num_int_elems,  /* Number of internal FEM elems */
                                int64_t num_bor_elems,  /* Number of border FEM elems */
                                int64_t num_node_cmaps, /* Number of nodal comm maps */
                                int64_t num_elem_cmaps, /* Number of elemental comm maps */
                                int     processor       /* Processor ID */
);

extern int ne_put_loadbal_param_cc(int       neid,           /* NetCDF/Exodus file ID */
                                   void_int *num_int_nodes,  /* Number of internal node IDs */
                                   void_int *num_bor_nodes,  /* Number of border node IDs */
                                   void_int *num_ext_nodes,  /* Number of external node IDs */
                                   void_int *num_int_elems,  /* Number of internal elem IDs */
                                   void_int *num_bor_elems,  /* Number of border elem IDs */
                                   void_int *num_node_cmaps, /* Number of nodal comm maps */
                                   void_int *num_elem_cmaps  /* Number of elem comm maps */
);

/*=============================================================================
 *     NS, SS & EB Global Parameter Routines
 *===========================================================================*/
extern int ne_get_ns_param_global(int       neid,          /* NetCDF/Exodus file ID */
                                  void_int *ns_ids_glob,   /* Global IDs of node sets */
                                  void_int *ns_n_cnt_glob, /* Count of nodes in node sets */
                                  void_int *ns_df_cnt_glob /* Count of dist. factors in ns */
);

extern int
ne_put_ns_param_global(int       neid,          /* NemesisI file ID */
                       void_int *global_ids,    /* Vector of global node-set IDs */
                       void_int *global_n_cnts, /* Vector of node counts in node-sets */
                       void_int *global_df_cnts /* Vector of dist factor counts in node-sets */
);

extern int ne_get_ss_param_global(int       neid,          /* NetCDF/Exodus file ID */
                                  void_int *ss_ids_glob,   /* Global side-set IDs */
                                  void_int *ss_s_cnt_glob, /* Global side count */
                                  void_int *ss_df_cnt_glob /* Global dist. factor count */
);

extern int ne_put_ss_param_global(int       neid,           /* NemesisI file ID */
                                  void_int *global_ids,     /* Vector of global side-set IDs */
                                  void_int *global_el_cnts, /* Vector of element/side */
                                                            /* counts in each side set */
                                  void_int *global_df_cnts  /* Vector of dist. factor */
                                                            /* counts in each side set */
);

extern int ne_get_eb_info_global(int       neid,       /* NemesisI file ID                 */
                                 void_int *el_blk_ids, /* Vector of global element IDs     */
                                 void_int *el_blk_cnts /* Vector of global element counts  */
);

extern int ne_put_eb_info_global(int       neid,       /* NemesisI file ID */
                                 void_int *el_blk_ids, /* Vector of global element IDs     */
                                 void_int *el_blk_cnts /* Vector of global element counts  */
);

/*=============================================================================
 *     NS, SS & EB Subset Routines
 *===========================================================================*/
extern int ne_get_n_side_set(int          neid,               /* NetCDF/Exodus file ID */
                             ex_entity_id side_set_id,        /* Side-set ID to read */
                             int64_t      start_side_num,     /* Starting element number */
                             int64_t      num_sides,          /* Number of sides to read */
                             void_int *   side_set_elem_list, /* List of element IDs */
                             void_int *   side_set_side_list  /* List of side IDs */
);

extern int ne_put_n_side_set(int             neid,               /* NetCDF/Exodus file ID */
                             ex_entity_id    side_set_id,        /* Side-set ID to write */
                             int64_t         start_side_num,     /* Starting element number */
                             int64_t         num_sides,          /* Number of sides to write */
                             const void_int *side_set_elem_list, /* List of element IDs */
                             const void_int *side_set_side_list  /* List of side IDs */
);

extern int ne_get_n_side_set_df(int          neid,          /* NetCDF/Exodus file ID */
                                ex_entity_id side_set_id,   /* Side-set ID */
                                int64_t      start_num,     /* Starting df number */
                                int64_t      num_df_to_get, /* Number of df's to read */
                                void *       side_set_df    /* Distribution factors */
);

extern int ne_put_n_side_set_df(int          neid,          /* NetCDF/Exodus file ID */
                                ex_entity_id side_set_id,   /* Side-set ID */
                                int64_t      start_num,     /* Starting df number */
                                int64_t      num_df_to_get, /* Number of df's to write */
                                void *       side_set_df    /* Distribution factors */
);

extern int ne_get_n_node_set(int          neid,              /* NetCDF/Exodus file ID */
                             ex_entity_id node_set_id,       /* Node set ID */
                             int64_t      start_node_num,    /* Node index to start reading at */
                             int64_t      num_node,          /* Number of nodes to read */
                             void_int *   node_set_node_list /* List of nodes in node set */
);

extern int ne_put_n_node_set(int             neid,              /* NetCDF/Exodus file ID */
                             ex_entity_id    node_set_id,       /* Node set ID */
                             int64_t         start_node_num,    /* Node index to start writing at */
                             int64_t         num_node,          /* Number of nodes to write */
                             const void_int *node_set_node_list /* List of nodes in node set */
);

extern int ne_get_n_node_set_df(int          neid,          /* NetCDF/Exodus file ID */
                                ex_entity_id node_set_id,   /* Node-set ID */
                                int64_t      start_num,     /* Starting df number */
                                int64_t      num_df_to_get, /* Number of df's to read */
                                void *       node_set_df    /* Distribution factors */
);

extern int ne_put_n_node_set_df(int          neid,          /* NetCDF/Exodus file ID */
                                ex_entity_id node_set_id,   /* Node-set ID */
                                int64_t      start_num,     /* Starting df number */
                                int64_t      num_df_to_get, /* Number of df's to write */
                                void *       node_set_df    /* Distribution factors */
);

extern int ne_get_n_coord(int     neid,           /* NetCDF/Exodus file ID */
                          int64_t start_node_num, /* Starting position to read from */
                          int64_t num_nodes,      /* Number of coords to read */
                          void *  x_coor,         /* Vector of X coordinates */
                          void *  y_coor,         /* Vector of Y coordinates */
                          void *  z_coor          /* Vector of Z coordinates */
);

extern int ne_put_n_coord(int     neid,           /* NetCDF/Exodus file ID */
                          int64_t start_node_num, /* Starting position to write to */
                          int64_t num_nodes,      /* Number of coords to write */
                          void *  x_coor,         /* Vector of X coordinates */
                          void *  y_coor,         /* Vector of Y coordinates */
                          void *  z_coor          /* Vector of Z coordinates */
);

extern int ne_get_n_elem_conn(int          neid,           /* NetCDF/Exodus file ID */
                              ex_entity_id elem_blk_id,    /* Element block ID */
                              int64_t      start_elem_num, /* Starting position to read from */
                              int64_t      num_elems,      /* Number of elements to read */
                              void_int *   connect         /* Connectivity vector */
);

extern int ne_put_n_elem_conn(int             neid,           /* NetCDF/Exodus file ID */
                              ex_entity_id    elem_blk_id,    /* Element block ID */
                              int64_t         start_elem_num, /* Starting position to write to */
                              int64_t         num_elems,      /* Number of elements to write */
                              const void_int *connect         /* Connectivity vector */
);

extern int ne_get_n_elem_attr(int          neid,           /* NetCDF/Exodus file ID */
                              ex_entity_id elem_blk_id,    /* Element block ID */
                              int64_t      start_elem_num, /* Starting position to read from */
                              int64_t      num_elems,      /* Number of elements to read */
                              void *       attrib          /* Attribute */
);

extern int ne_put_n_elem_attr(int          neid,           /* NetCDF/Exodus file ID */
                              ex_entity_id elem_blk_id,    /* Element block ID */
                              int64_t      start_elem_num, /* Starting position to write to */
                              int64_t      num_elems,      /* Number of elements to write */
                              void *       attrib          /* Attribute */
);

extern int ne_get_elem_type(int          neid,        /* NetCDF/Exodus file ID */
                            ex_entity_id elem_blk_id, /* Element block ID */
                            char *       elem_type    /* The name of the element type */
);

/*=============================================================================
 *     Variable Routines
 *===========================================================================*/
extern int ne_get_n_elem_var(int          neid,              /* NetCDF/Exodus file ID */
                             int          time_step,         /* time index */
                             int          elem_var_index,    /* elemental variable index */
                             ex_entity_id elem_blk_id,       /* elemental block id */
                             int64_t      num_elem_this_blk, /* number of elements in block */
                             int64_t      start_elem_num,    /* Starting position for read */
                             int64_t      num_elem,          /* Number of elements to read */
                             void *       elem_var_vals      /* variable values */
);

extern int ne_put_elem_var_slab(int          neid,           /* NetCDF/Exodus file ID */
                                int          time_step,      /* time index */
                                int          elem_var_index, /* elemental variable index */
                                ex_entity_id elem_blk_id,    /* elemental block id */
                                int64_t      start_pos,      /* Starting position to write to */
                                int64_t      num_vals,       /* Number of elements to write */
                                void *       elem_var_vals   /* variable values */
);

extern int ne_get_n_nodal_var(int     neid,            /* NetCDF/Exodus file ID */
                              int     time_step,       /* whole time step number */
                              int     nodal_var_index, /* index of desired nodal var */
                              int64_t start_node_num,  /* starting node number */
                              int64_t num_nodes,       /* number of nodes to read */
                              void *  nodal_vars       /* array of nodal var values */
);

extern int ne_put_nodal_var_slab(int     neid,            /* NetCDF/Exodus file ID */
                                 int     time_step,       /* The time step index */
                                 int     nodal_var_index, /* Nodal variable index */
                                 int64_t start_pos,       /* Start position for write */
                                 int64_t num_vals,        /* Number of nodal variables */
                                 void *  nodal_var_vals   /* Nodal variable values */
);

/*=============================================================================
 *     Number Map Routines
 *===========================================================================*/
extern int ne_get_n_elem_num_map(int       neid,      /* NetCDF/Exodus file ID */
                                 int64_t   start_ent, /* Starting position to read from */
                                 int64_t   num_ents,  /* Number of elements to read */
                                 void_int *elem_map   /* element map numbers */
);

extern int ne_put_n_elem_num_map(int             neid,      /* NetCDF/Exodus file ID */
                                 int64_t         start_ent, /* Starting position to read from */
                                 int64_t         num_ents,  /* Number of elements to read */
                                 const void_int *elem_map   /* element map numbers */
);

extern int ne_get_n_node_num_map(int       neid,      /* NetCDF/Exodus file ID */
                                 int64_t   start_ent, /* starting node number */
                                 int64_t   num_ents,  /* number of nodes to read */
                                 void_int *node_map   /* vector for node map */
);

extern int ne_put_n_node_num_map(int             neid,      /* NetCDF/Exodus file ID */
                                 int64_t         start_ent, /* starting node number */
                                 int64_t         num_ents,  /* number of nodes to read */
                                 const void_int *node_map   /* vector for node map */
);

extern int ne_get_node_map(int       neid,      /* NetCDF/Exodus file ID */
                           void_int *node_mapi, /* Internal FEM node IDs */
                           void_int *node_mapb, /* Border FEM node IDs */
                           void_int *node_mape, /* External FEM node IDs */
                           int       processor  /* Processor IDs */
);

extern int ne_put_node_map(int       neid,      /* NetCDF/Exodus file ID */
                           void_int *node_mapi, /* Internal FEM node IDs */
                           void_int *node_mapb, /* Border FEM node IDs */
                           void_int *node_mape, /* External FEM node IDs */
                           int       processor  /* This processor ID */
);

extern int ne_get_elem_map(int       neid,      /* NetCDF/Exodus file ID */
                           void_int *elem_mapi, /* Internal element IDs */
                           void_int *elem_mapb, /* Border element IDs */
                           int       processor  /* Processor ID */
);

extern int ne_put_elem_map(int       neid,      /* NetCDF/Exodus file ID */
                           void_int *elem_mapi, /* Internal FEM element IDs */
                           void_int *elem_mapb, /* Border FEM element IDs */
                           int       processor  /* This processor ID */
);

/*=============================================================================
 *     Communications Maps Routines
 *===========================================================================*/

extern int ne_get_cmap_params(int       neid,                /* NetCDF/Exodus file ID */
                              void_int *node_cmap_ids,       /* Nodal comm. map IDs */
                              void_int *node_cmap_node_cnts, /* Number of nodes in each map */
                              void_int *elem_cmap_ids,       /* Elemental comm. map IDs */
                              void_int *elem_cmap_elem_cnts, /* Number of elems in each map */
                              int       processor            /* This processor ID */
);

extern int ne_put_cmap_params(int       neid,               /* NetCDF/Exodus file ID */
                              void_int *node_map_ids,       /* Node map IDs */
                              void_int *node_map_node_cnts, /* Nodes in nodal comm */
                              void_int *elem_map_ids,       /* Elem map IDs */
                              void_int *elem_map_elem_cnts, /* Elems in elemental comm */
                              int64_t   processor           /* This processor ID */
);

extern int ne_put_cmap_params_cc(int       neid,               /* NetCDF/Exodus file ID */
                                 void_int *node_map_ids,       /* Node map IDs */
                                 void_int *node_map_node_cnts, /* Nodes in nodal comm */
                                 void_int *node_proc_ptrs,     /* Pointer into array for */
                                                               /* node maps               */
                                 void_int *elem_map_ids,       /* Elem map IDs */
                                 void_int *elem_map_elem_cnts, /* Elems in elemental comm */
                                 void_int *elem_proc_ptrs      /* Pointer into array for */
                                                               /* elem maps               */
);

extern int ne_get_node_cmap(int          neid,     /* NetCDF/Exodus file ID */
                            ex_entity_id map_id,   /* Map ID */
                            void_int *   node_ids, /* FEM node IDs */
                            void_int *   proc_ids, /* Processor IDs */
                            int          processor /* This processor ID */
);

extern int ne_put_node_cmap(int          neid,     /* NetCDF/Exodus file ID */
                            ex_entity_id map_id,   /* Nodal comm map ID */
                            void_int *   node_ids, /* FEM node IDs */
                            void_int *   proc_ids, /* Processor IDs */
                            int          processor /* This processor ID */
);

extern int ne_get_elem_cmap(int          neid,     /* NetCDF/Exodus file ID */
                            ex_entity_id map_id,   /* Elemental comm map ID */
                            void_int *   elem_ids, /* Element IDs */
                            void_int *   side_ids, /* Element side IDs */
                            void_int *   proc_ids, /* Processor IDs */
                            int          processor /* This processor ID */
);

extern int ne_put_elem_cmap(int          neid,     /* NetCDF/Exodus file ID */
                            ex_entity_id map_id,   /* Elemental comm map ID */
                            void_int *   elem_ids, /* Vector of element IDs */
                            void_int *   side_ids, /* Vector of side IDs */
                            void_int *   proc_ids, /* Vector of processor IDs */
                            int          processor /* This processor ID */
);

/* Should be internal use only, but was in external include file for
   nemesis and some codes are using the function
*/
int ne_get_idx(int         neid,        /* NetCDF/Exodus file ID */
               const char *ne_var_name, /* Nemesis index variable name */
               int64_t *   index,       /* array of length 2 to hold results */
               int         pos          /* position of this proc/cmap in index */
);

#ifdef __cplusplus
}
#endif

#endif /* _NE_NEMESIS_H */
