// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <cstring>
#include "pamgen_im_ne_nemesisI.h"
#include "pamgen_mesh_specification.h"
#include <strings.h>


/*****************************************************************************/
int im_ne_get_init_global(int   /* neid */, 		  /* NemesisI file ID */
			  int  *num_nodes_g,	  /* Number of global FEM nodes */
			  int  *num_elems_g,	  /* Number of global FEM elements */
			  int  *num_elem_blks_g, /* Number of global elem blocks */
			  int  *num_node_sets_g, /* Number of global node sets */
			  int  *num_side_sets_g  /* Number of global side sets */
			  )
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  *num_nodes_g     = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODES_GLOBAL);
  *num_elems_g     = ms->getMSI(ms_lt::Mesh_Specification::NUM_ELEMS_GLOBAL);
  *num_elem_blks_g = ms->getMSI(ms_lt::Mesh_Specification::NUM_ELM_BLKS_GLOBAL);
  *num_node_sets_g = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_SETS_GLOBAL);
  *num_side_sets_g = ms->getMSI(ms_lt::Mesh_Specification::NUM_SIDE_SETS_GLOBAL);
  
  return 0;
}

/*****************************************************************************/
int im_ne_get_init_info(int   /* neid */,		/* NemesisI file ID */
			int  *num_proc,	/* Number of processors */
			int  *num_proc_in_f,	/* Number of procs in this file */
			char *ftype)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;    
  *num_proc = ms->getMSI(ms_lt::Mesh_Specification::NUM_TOTAL_PROC);
  *num_proc_in_f = ms->getMSI(ms_lt::Mesh_Specification::NUM_PROC_IN_FILE);
  const std::string ft = ms->File_Type();
  strcpy(ftype,ft.c_str());
  return 0;
}

/*****************************************************************************/
int im_ne_get_eb_info_global(int /* neid */,		/* NemesisI file ID                 */
			     int *el_blk_ids,	/* Vector of global element IDs     */
			     int *el_blk_cnts	/* Vector of global element counts  */
			     )
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  int num_eb = ms->getMSI(ms_lt::Mesh_Specification::NUM_ELM_BLKS_GLOBAL);
  long long * ebids = ms->getMSP(ms_lt::Mesh_Specification::ELEM_BLK_IDS_GLOBAL);
  long long * ebcts = ms->getMSP(ms_lt::Mesh_Specification::ELEM_BLK_CNTS_GLOBAL);
  for(int i = 0; i < num_eb; i ++){
    el_blk_ids[i] = ebids[i];
    el_blk_cnts[i] = ebcts[i];
  }
  return 0;
}

/*****************************************************************************/
int im_ne_get_loadbal_param(int   /* neid */, 	/* NetCDF/Exodus file ID */
			    int  *num_int_nodes,  /* Number of internal FEM nodes */
			    int  *num_bor_nodes,  /* Number of border FEM nodes */
			    int  *num_ext_nodes,  /* Number of external FEM nodes */
			    int  *num_int_elems,  /* Number of internal FEM elems */
			    int  *num_bor_elems,  /* Number of border FEM elems */
			    int  *num_node_cmaps, /* Number of nodal comm maps */
			    int  *num_elem_cmaps, /* Number of elemental comm maps */
			    int   /* processor */         /* Processor ID */
			    )
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  *num_int_nodes  = ms->getMSI(ms_lt::Mesh_Specification::NUM_INTERNAL_NODES);
  *num_bor_nodes  = ms->getMSI(ms_lt::Mesh_Specification::NUM_BORDER_NODES);
  *num_ext_nodes  = ms->getMSI(ms_lt::Mesh_Specification::NUM_EXTERNAL_NODES);
  *num_int_elems  = ms->getMSI(ms_lt::Mesh_Specification::NUM_INTERNAL_ELEMS);
  *num_bor_elems  = ms->getMSI(ms_lt::Mesh_Specification::NUM_BORDER_ELEMS);
  *num_node_cmaps = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_COMM_MAPS);
  *num_elem_cmaps = ms->getMSI(ms_lt::Mesh_Specification::NUM_ELEM_COMM_MAPS);
  
  
  return 0;
}

/*****************************************************************************/
int im_ne_get_elem_map(int   /* neid */,		/* NetCDF/Exodus file ID */
		       int  *elem_mapi,	/* Internal element IDs */
		       int  *elem_mapb,	/* Border element IDs */
		       int   /* processor */		/* Processor ID */
		       )
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  
  int int_elems = ms->getMSI(ms_lt::Mesh_Specification::NUM_INTERNAL_ELEMS);
  int bor_elems = ms->getMSI(ms_lt::Mesh_Specification::NUM_BORDER_ELEMS);
  
  long long * int_el_arr = ms->getMSP(ms_lt::Mesh_Specification::INTERNAL_ELEMENTS);
  long long * bor_el_arr = ms->getMSP(ms_lt::Mesh_Specification::BORDER_ELEMENTS);
  
  for(int i = 0; i < int_elems; i ++)elem_mapi[i] = int_el_arr[i];
  for(int i = 0; i < bor_elems; i ++)elem_mapb[i] = bor_el_arr[i];
  
  return 0;
}

/*****************************************************************************/
int im_ne_get_node_map(int   /* neid */,		/* NetCDF/Exodus file ID */
		       int  *node_mapi,	/* Internal FEM node IDs */
		       int  *node_mapb,	/* Border FEM node IDs */
		       int  *node_mape,	/* External FEM node IDs */
		       int   /* processor */		/* Processor IDs */
		       )
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1; 
  
  int int_nodes = ms->getMSI(ms_lt::Mesh_Specification::NUM_INTERNAL_NODES);
  int bor_nodes = ms->getMSI(ms_lt::Mesh_Specification::NUM_BORDER_NODES);
  int ext_nodes = ms->getMSI(ms_lt::Mesh_Specification::NUM_EXTERNAL_NODES);

  long long * int_nd_arr = ms->getMSP(ms_lt::Mesh_Specification::INTERNAL_NODES);
  long long * bor_nd_arr = ms->getMSP(ms_lt::Mesh_Specification::BORDER_NODES);
  long long * ext_nd_arr = ms->getMSP(ms_lt::Mesh_Specification::EXTERNAL_NODES);

  for(int i = 0; i < int_nodes; i ++)node_mapi[i] = int_nd_arr[i];
  for(int i = 0; i < bor_nodes; i ++)node_mapb[i] = bor_nd_arr[i];
  for(int i = 0; i < ext_nodes; i ++)node_mape[i] = ext_nd_arr[i];

  return 0;
}

/*****************************************************************************/
int im_ne_get_cmap_params(int /* neid */,                  /* NetCDF/Exodus file ID */
			  int *node_cmap_ids,        /* Nodal comm. map IDs */
			  int *node_cmap_node_cnts,  /* Number of nodes in each map */
			  int *elem_cmap_ids,        /* Elemental comm. map IDs */
			  int *elem_cmap_elem_cnts,  /* Number of elems in each map */
			  int  /* processor */             /* This processor ID */
			  )
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1; 
  
  int num_node_cmaps = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_COMM_MAPS);
  int num_elem_cmaps = ms->getMSI(ms_lt::Mesh_Specification::NUM_ELEM_COMM_MAPS);

  long long * nd_cmap_id_arr = ms->getMSP(ms_lt::Mesh_Specification::NODE_CMAP_IDS);
  const long long * nd_cmap_ct_arr = ms->getMSP(ms_lt::Mesh_Specification::NODE_CMAP_NODE_CNTS);
  long long * el_cmap_id_arr = ms->getMSP(ms_lt::Mesh_Specification::ELEM_CMAP_IDS);
  const long long * el_cmap_ct_arr = ms->getMSP(ms_lt::Mesh_Specification::ELEM_CMAP_ELEM_CNTS);

  for(int i = 0; i < num_node_cmaps; i ++){
    node_cmap_ids[i]      = nd_cmap_id_arr[i];
    node_cmap_node_cnts[i] = nd_cmap_ct_arr[i];
  }

  for(int i = 0; i < num_elem_cmaps; i ++){
    elem_cmap_ids[i]       = el_cmap_id_arr[i];
    elem_cmap_elem_cnts[i] = el_cmap_ct_arr[i];
  }

  return 0;
}

/*****************************************************************************/
int im_ne_get_node_cmap(int  /* neid */,             /* NetCDF/Exodus file ID */
			int  map_id,           /* Map ID */
			int *node_ids,         /* FEM node IDs */
			int *proc_ids,         /* Processor IDs */
			int  /* processor */         /* This processor ID */
			)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int map_index = -1;  
  // get the map index
  int num_node_cmaps = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_COMM_MAPS);
  long long * nd_cmap_id_arr = ms->getMSP(ms_lt::Mesh_Specification::NODE_CMAP_IDS);
  for(int i = 0; i < num_node_cmaps; i ++)if(map_id == nd_cmap_id_arr[i])map_index = i;
  if(map_index == -1)return -1;


  const long long * nd_cmap_ct_arr = ms->getMSP(ms_lt::Mesh_Specification::NODE_CMAP_NODE_CNTS);
  int num_comm_nodes = nd_cmap_ct_arr[map_index];
  
  long long * const * nd_id_arr = ms->getMSPP(ms_lt::Mesh_Specification::COMM_NODE_IDS);//Comm_Node_Ids();
  long long * const * nd_proc_arr = ms->getMSPP(ms_lt::Mesh_Specification::COMM_NODE_PROC_IDS);//Comm_Node_Proc_Ids();

  for(int i = 0; i < num_comm_nodes; i ++){
    node_ids[i] = nd_id_arr[map_index][i];
    proc_ids[i] = nd_proc_arr[map_index][i];
  }
  
 return 0;
}

/*****************************************************************************/
int im_ne_get_elem_cmap(int  /* neid */,     /* NetCDF/Exodus file ID */
			int  map_id,   /* Elemental comm map ID */
			int *elem_ids, /* Element IDs */
			int *side_ids, /* Element side IDs */
			int *proc_ids, /* Processor IDs */
			int  /* processor */ /* This processor ID */
			)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int map_index = -1;
  // get the map index
  int num_elem_cmaps = ms->getMSI(ms_lt::Mesh_Specification::NUM_ELEM_COMM_MAPS);
  long long * elem_cmap_id_arr = ms->getMSP(ms_lt::Mesh_Specification::ELEM_CMAP_IDS);
  for(int i = 0; i < num_elem_cmaps; i ++)if(map_id == elem_cmap_id_arr[i])map_index = i;
  if(map_index == -1)return -1;


  const long long * el_cmap_ct_arr = ms->getMSP(ms_lt::Mesh_Specification::ELEM_CMAP_ELEM_CNTS);
  int num_comm_elements = el_cmap_ct_arr[map_index];

  long long * const * el_id_arr = ms->getMSPP(ms_lt::Mesh_Specification::COMM_ELEM_IDS);
  long long * const * el_side_arr = ms->getMSPP(ms_lt::Mesh_Specification::COMM_SIDE_IDS);
  long long * const * el_proc_arr = ms->getMSPP(ms_lt::Mesh_Specification::COMM_ELEM_PROC_IDS);

  for(int i = 0; i < num_comm_elements; i ++){
    elem_ids[i] = el_id_arr[map_index][i];
    side_ids[i] = el_side_arr[map_index][i];
    proc_ids[i] = el_proc_arr[map_index][i];
  }

  return 0;
}

/*****************************************************************************/
int im_ne_get_ns_param_global(int /* neid */,	     /* NetCDF/Exodus file ID */
			      int *ns_ids_glob,     /* Global IDs of node sets */
			      int *ns_n_cnt_glob,   /* Count of nodes in node sets */
			      int *ns_df_cnt_glob   /* Count of dist. factors in ns */
			      )
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  
  int nns = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_SETS_GLOBAL);
  
  const long long * nsids = ms->getMSP(ms_lt::Mesh_Specification::NS_IDS_GLOBAL);
  const long long * ncts = ms->getMSP(ms_lt::Mesh_Specification::NS_CNTS_GLOBAL);
  const long long * dfs = ms->getMSP(ms_lt::Mesh_Specification::NS_DF_CNTS_GLOBAL);
  for(int i = 0; i < nns; i ++){
    ns_ids_glob[i] = nsids[i];
    ns_n_cnt_glob[i] = ncts[i];
    ns_df_cnt_glob[i] = dfs[i];
  }
  return 0;
}

/*****************************************************************************/
int im_ne_get_ss_param_global(int /* neid */,	    /* NetCDF/Exodus file ID */
			      int *ss_ids_glob,    /* Global side-set IDs */
			      int *ss_s_cnt_glob,  /* Global side count */
			      int *ss_df_cnt_glob  /* Global dist. factor count */
			      )
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  
  int nss = ms->getMSI(ms_lt::Mesh_Specification::NUM_SIDE_SETS_GLOBAL);
  
  const long long * nsids = ms->getMSP(ms_lt::Mesh_Specification::SS_IDS_GLOBAL);
  const long long * ncts = ms->getMSP(ms_lt::Mesh_Specification::SS_CNTS_GLOBAL);
  const long long * dfs = ms->getMSP(ms_lt::Mesh_Specification::SS_DF_CNTS_GLOBAL);
  for(int i = 0; i < nss; i ++){
    ss_ids_glob[i] = nsids[i];
    ss_s_cnt_glob[i] = ncts[i];
    ss_df_cnt_glob[i] = dfs[i];
  }
  return 0;
}
