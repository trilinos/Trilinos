// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <cstring>
#include <stdlib.h>
#include "pamgen_im_exodusII.h"
#include "pamgen_mesh_specification.h"
#include <strings.h>
#include <stdio.h>


/*****************************************************************************/
int im_ex_get_init (int   /* im_exoid */,
		    char *title,
		    int  *num_dim,
		    int  *num_nodes,
		    int  *num_elem, 
		    int  *num_elem_blk,
		    int  *num_node_sets,
		    int  *num_side_sets)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  const std::string tit = ms->Title();
  strcpy(title,tit.c_str());
  *num_dim = ms->getMSI(ms_lt::Mesh_Specification::DIM);
  *num_nodes = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODES);
  *num_elem = ms->getMSI(ms_lt::Mesh_Specification::NUM_ELEMENTS);
  *num_elem_blk = ms->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);
  *num_node_sets = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_SETS);
  *num_side_sets = ms->getMSI(ms_lt::Mesh_Specification::NUM_SIDE_SETS);

  return 0;
}

/*****************************************************************************/
int im_ex_inquire (int   /* exoid */,
		   int   req_info,
		   int  *ret_int,
		   float *ret_float,
		   char *ret_char)
/*****************************************************************************/
{
  char  errmsg[MAX_ERR_LENGTH];
  
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  
  
  switch (req_info)
    {
      
    case IM_EX_INQ_API_VERS:
      *ret_float = 0.1;
      break;
      
    case IM_EX_INQ_TITLE:
      {
      /*     returns the title of the database */
      const std::string tit = ms->Title();
      strcpy(ret_char,tit.c_str());
      
      break;
      }
    case IM_EX_INQ_DIM:
      *ret_int = ms->getMSI(ms_lt::Mesh_Specification::DIM);
      
      break;
      
    case IM_EX_INQ_NODES:
      *ret_int = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODES);
      
      break;
      
    case IM_EX_INQ_ELEM:
      *ret_int = ms->getMSI(ms_lt::Mesh_Specification::NUM_ELEMENTS);
      
      break;
      
    case IM_EX_INQ_ELEM_BLK:
      *ret_int = ms->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);
      
      break;
      
    case IM_EX_INQ_NODE_SETS:
      *ret_int = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_SETS);
      break;
      
    case IM_EX_INQ_NS_NODE_LEN:
      *ret_int = 0;
      /*     returns the length of the concatenated node sets node list */      
      break;
      
    case IM_EX_INQ_NS_DF_LEN:
      /*     returns the length of the concatenated node sets dist factor list */
      *ret_int = 0;
      break;
      
    case IM_EX_INQ_SIDE_SETS:
      
      /*     returns the number of side sets */
      *ret_int = ms->getMSI(ms_lt::Mesh_Specification::NUM_SIDE_SETS);
      
      
      break;
      
    case IM_EX_INQ_SS_ELEM_LEN:
      *ret_int = 0;
      break;

    case IM_EX_INQ_SS_NODE_LEN:
      *ret_int = ms->getMSI(ms_lt::Mesh_Specification::NUM_SIDE_SET_NODES);
      /*     returns the length of the concatenated side sets node list */
      break;

    case IM_EX_INQ_SS_DF_LEN:
      *ret_int = 0;
      /*     returns the length of the concatenated side sets df factor */
      break;
      
      
    case IM_EX_INQ_QA:
      /*     returns the number of QA records */
      *ret_int = ms->getMSI(ms_lt::Mesh_Specification::NUM_QA_RECORDS);
      
      
      break;
      
    case IM_EX_INQ_INFO:
      /*     returns the number of information records */
      *ret_int = ms->getMSI(ms_lt::Mesh_Specification::NUM_INFO_RECORDS);
      
      break;
      
    case IM_EX_INQ_TIME:
      
      return -1;
      break;
    case IM_EX_INQ_EB_PROP:
      /*     returns the number of element block properties */
      *ret_int = 0;
      break;
      
    case IM_EX_INQ_NS_PROP:
      *ret_int = 0;
      /*     returns the number of node set properties */
      break;
      
    case IM_EX_INQ_SS_PROP:
      *ret_int = 0;
      /*     returns the number of side set properties */
      break;
      
    case IM_EX_INQ_ELEM_MAP:
      
      /*     returns the number of element maps */
      return -1;
      
      break;
      
    case IM_EX_INQ_EM_PROP:
      
      /*     returns the number of element map properties */
      return -1;
      break;
      
    case IM_EX_INQ_NODE_MAP:
      
      /*     returns the number of node maps */
      return -1;
      break;
      
    case IM_EX_INQ_NM_PROP:
      /*     returns the number of element map properties */
      return -1;
      break;
      
    default:
      *ret_int = 0;
      sprintf(errmsg, "Error: invalid inquiry %d", req_info);
      return(-1);
    }
  return (0);
}

/*****************************************************************************/
int im_ex_get_coord (int /* exoid */,
		     void *x_coor,
		     void *y_coor,
		     void *z_coor)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int num_nds = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODES);
  int dim = ms->getMSI(ms_lt::Mesh_Specification::DIM);

  double * x = (double *)x_coor;
  double * y = (double *)y_coor;
  double * z = (double *)z_coor;

  const double  * the_coords = ms->Coord();
  for(int i = 0; i < num_nds; i ++){
    x[i] = the_coords[i];
    y[i] = the_coords[i+num_nds];
    if(dim == 3)z[i] = the_coords[i+2*num_nds];
  }

  return 0;
}

/*****************************************************************************/
int im_ex_get_coord_names (int    /* exoid */,
			   char **coord_names)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  std::string * cnames =  ms->getMSPSA(ms_lt::Mesh_Specification::COORDINATE_NAMES);//ms->Coordinate_Names();
  int dim = ms->getMSI(ms_lt::Mesh_Specification::DIM);
  if(dim >=1)strcpy(coord_names[0],cnames[0].c_str());
  if(dim >=2)strcpy(coord_names[1],cnames[1].c_str());
  if(dim >=3)strcpy(coord_names[2],cnames[2].c_str());

  return 0;
}

/*****************************************************************************/
int im_ex_get_map (int  /* exoid */, int *elem_map)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  const long long * the_map = ms->getMSP(ms_lt::Mesh_Specification::ELEM_ORDER_MAP);
  int num_elem = ms->getMSI(ms_lt::Mesh_Specification::NUM_ELEMENTS);
  for(int i = 0; i < num_elem; i ++)elem_map[i] = the_map[i];

  return 0;
}

/*****************************************************************************/
int im_ex_get_elem_num_map (int  /* exoid */,
			    int *elem_map)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int num_elem = ms->getMSI(ms_lt::Mesh_Specification::NUM_ELEMENTS);
  const long long * gen = ms->getMSP(ms_lt::Mesh_Specification::GLOBAL_ELEMENT_NUMBERS);

  for(int i = 0; i < num_elem; i ++)elem_map[i] = gen[i];
  
  return 0;
}

/*****************************************************************************/
int im_ex_get_node_num_map (int  /* exoid */,
			    int *node_map)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int num_nds = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODES);
  const long long * gnn = ms->getMSP(ms_lt::Mesh_Specification::GLOBAL_NODE_NUMBERS);//Global_Node_Numbers();

  for(int i = 0; i < num_nds; i ++)node_map[i] = gnn[i];

  return 0;
}

/*****************************************************************************/
int im_ex_get_elem_blk_ids (int  /* exoid */, int *ids)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int nb = ms->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);
  const long long * bid_arr = ms->getMSP(ms_lt::Mesh_Specification::BLOCK_ID);

  for(int i = 0; i < nb; i ++)ids[i] = bid_arr[i];
  
  return 0;
}

/*****************************************************************************/
int im_ex_get_elem_block (int   /* exoid */,
			  int   elem_blk_id,
			  char *elem_type,
			  int  *num_elem_this_blk, 
			  int  *num_nodes_per_elem,
			  int  *num_attr)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int nb = ms->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);
  const long long * bid_arr = ms->getMSP(ms_lt::Mesh_Specification::BLOCK_ID);
  
  // get the index
  int el_blk_index = -1;
  for(int i = 0; i < nb; i ++)if(elem_blk_id == bid_arr[i])el_blk_index = i;
  if(el_blk_index == -1)return -1;

  const std::string * el_typ = ms->getMSPSA(ms_lt::Mesh_Specification::ELEMENT_TYPES);
  strcpy(elem_type,el_typ[el_blk_index].c_str());

  const long long * nbe = ms->getMSP(ms_lt::Mesh_Specification::ELEMENTS_IN_BLOCK);
  *num_elem_this_blk = nbe[el_blk_index];
  
  const long long * nnpe = ms->getMSP(ms_lt::Mesh_Specification::NODES_PER_ELEMENT);
  *num_nodes_per_elem = nnpe[el_blk_index];

  const long long * natr = ms->getMSP(ms_lt::Mesh_Specification::ELEMENT_ATTRIBUTES);
  *num_attr = natr[el_blk_index];
  
  return 0;
}



// int im_ex_get_elem_attr (int   exoid,
// 			 int   elem_blk_id,
// 			 void *attrib)
// {
//   ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
//   if(!ms)return -1;
//   int nb = ms->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);
 
  
//   return 0;
// }

/*****************************************************************************/
int im_ex_get_elem_conn (int   /* exoid */,
			 int   elem_blk_id,
			 int  *connect)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int nb = ms->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);
  const long long * bid_arr = ms->getMSP(ms_lt::Mesh_Specification::BLOCK_ID);
  
  // get the index
  int el_blk_index = -1;
  for(int i = 0; i < nb; i ++)if(elem_blk_id == bid_arr[i])el_blk_index = i;
  if(el_blk_index == -1)return -1;

  const long long * const * enl = ms->getMSPP(ms_lt::Mesh_Specification::ELMT_NODE_LINKAGE);//Element_Node_Linkage();
  const long long * nbe = ms->getMSP(ms_lt::Mesh_Specification::ELEMENTS_IN_BLOCK);
  const long long * nnpe_arr = ms->getMSP(ms_lt::Mesh_Specification::NODES_PER_ELEMENT);
  int nnpe = nnpe_arr[el_blk_index];
  int nel = nbe[el_blk_index];
  
  for(int i = 0; i < nel*nnpe; i ++)connect[i] = enl[el_blk_index][i];

  return 0;
}

/*****************************************************************************/
int im_ex_get_node_set_ids (int  /* exoid */,
			    int *ids)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int nns = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_SETS);
  const long long * ids_arr = ms->getMSP(ms_lt::Mesh_Specification::NODE_SET_ID);
  for(int i = 0; i < nns; i ++)ids[i] = ids_arr[i];
  
  return 0;
}

/*****************************************************************************/
int get_ss_index(int ssid)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int nss = ms->getMSI(ms_lt::Mesh_Specification::NUM_SIDE_SETS);
  const long long * ids_arr = ms->getMSP(ms_lt::Mesh_Specification::SIDE_SET_ID);
  for(int i = 0; i < nss; i ++)
    if(ssid == ids_arr[i])
      return i;

  return -1;  
}

/*****************************************************************************/
int get_ns_index(int nsid)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int nns = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_SETS);
  const long long * ids_arr = ms->getMSP(ms_lt::Mesh_Specification::NODE_SET_ID);
  for(int i = 0; i < nns; i ++)
    if(nsid == ids_arr[i])
      return i;
  
  return -1;  
}

/*****************************************************************************/
int im_ex_get_side_set_ids (int  /* exoid */,
			    int *ids)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  int nss = ms->getMSI(ms_lt::Mesh_Specification::NUM_SIDE_SETS);
  
  const long long * ids_arr = ms->getMSP(ms_lt::Mesh_Specification::SIDE_SET_ID);
  for(int i = 0; i < nss; i ++)ids[i] = ids_arr[i];
  return 0;
}

/*****************************************************************************/
int im_ex_get_node_set_param (int /* exoid */,
			      int node_set_id,
			      int * num_nodes_in_set,
			      int * num_dist_in_set)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  
  int ns_index = get_ns_index(node_set_id);
  if(ns_index == -1)return -1;

  const long long * nn_in_set = ms->getMSP(ms_lt::Mesh_Specification::NUM_NODES_IN_NODE_SET);
  *num_nodes_in_set = nn_in_set[ns_index];
  
  const long long * nn_df_in_set = ms->getMSP(ms_lt::Mesh_Specification::NUM_DF_IN_NODE_SET);
  *num_dist_in_set = nn_df_in_set[ns_index];

  return 0;
}

/*****************************************************************************/
int im_ex_get_side_set_param (int /* exoid */,
			      int side_set_id,
			      int * num_side_in_set,
			      int * num_dist_in_set)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  
  int ss_index = get_ss_index(side_set_id);
  if(ss_index == -1)return -1;

  const long long * ns_in_set = ms->getMSP(ms_lt::Mesh_Specification::NUM_ELEMENTS_IN_SIDE_SET);
  *num_side_in_set = ns_in_set[ss_index];
  
  const long long * ss_df_in_set = ms->getMSP(ms_lt::Mesh_Specification::NUM_DF_IN_SIDE_SET);
  *num_dist_in_set = ss_df_in_set[ss_index];

  return 0;
}

/*****************************************************************************/
int im_ex_get_node_set (int   /* exoid */,
			int   node_set_id,
			int  *node_set_node_list)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  
  int ns_index = get_ns_index(node_set_id);
  if(ns_index == -1)return -1;
  
  const long long * nn_in_set = ms->getMSP(ms_lt::Mesh_Specification::NUM_NODES_IN_NODE_SET);
  int num_nds_in_set = nn_in_set[ns_index];
  
  const long long * const * ns_arr = ms->getMSPP(ms_lt::Mesh_Specification::NODE_SET_NODES);//Node_Set_Nodes();
  for(int i = 0; i < num_nds_in_set; i ++)node_set_node_list[i] = ns_arr[ns_index][i];

  return 0;
}

/*****************************************************************************/
int im_ex_get_side_set (int   /* exoid */,
			int   side_set_id,
			int  *side_set_elem_list,
			int  *side_set_side_list)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  
  int ss_index = get_ss_index(side_set_id);
  if(ss_index == -1)return -1;
  
  const long long * ne_in_set = ms->getMSP(ms_lt::Mesh_Specification::NUM_ELEMENTS_IN_SIDE_SET);
  int num_fcs_in_set = ne_in_set[ss_index];
  
  const long long * const * ss_arr = ms->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_ELEMENTS);
  for(int i = 0; i < num_fcs_in_set; i ++)side_set_elem_list[i] = ss_arr[ss_index][i];

  const long long * const * ssf_arr = ms->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_FACES);
  for(int i = 0; i < num_fcs_in_set; i ++)side_set_side_list[i] = ssf_arr[ss_index][i];

  return 0;
}

/*****************************************************************************/
int im_ex_get_side_set_node_list(int /* exoid */,
				 int side_set_id,
				 int *side_set_node_cnt_list,
				 int *side_set_node_list)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  
  int ss_index = get_ss_index(side_set_id);
  if(ss_index == -1)return -1;

  const long long * ne_in_set = ms->getMSP(ms_lt::Mesh_Specification::NUM_ELEMENTS_IN_SIDE_SET);
  int num_fcs_in_set = ne_in_set[ss_index];

  const long long * const * ssnc = ms->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_NODE_COUNTER);
  for(int i = 0; i < num_fcs_in_set; i ++)side_set_node_cnt_list[i] = ssnc[ss_index][i];
  
  const long long * nn_in_sset = ms->getMSP(ms_lt::Mesh_Specification::NUM_NODES_IN_SIDE_SET);
  int nn = nn_in_sset[ss_index];
  
  const long long * const * ssn = ms->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_NODES);
  for(int i = 0; i < nn; i ++)side_set_node_list[i] = ssn[ss_index][i]; 

  return 0;
}

/*****************************************************************************/
int im_ex_get_qa (int /* exoid */,
		  char *qa_record[][4])
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;
  
  int num_qa = ms->getMSI(ms_lt::Mesh_Specification::NUM_QA_RECORDS);
  
  typedef std::string QA_Record[4];

  const QA_Record * qa_recs = ms->QA_Records();
  
  for(int i = 0; i < num_qa; i++){
    for(int j = 0; j < 4; j++){
      strcpy(qa_record[i][j],qa_recs[i][j].c_str());
    }
  }

  return 0;
}

/*****************************************************************************/
int im_ex_get_info (int /* exoid */, char **info)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int num_info = ms->getMSI(ms_lt::Mesh_Specification::NUM_INFO_RECORDS);
  const std::string * info_strings =  ms->getMSPSA(ms_lt::Mesh_Specification::INFO_STRINGS);
  
  for(int i = 0; i < num_info; i ++){
    strcpy(info[i],info_strings[i].c_str());
  }
  
  
  return 0;
}
/*****************************************************************************/
int im_ex_get_elem_blk_parent_mesh (int  /* exoid */, int *ids)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int nb = ms->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);
  const long long * bid_arr = ms->getMSP(ms_lt::Mesh_Specification::BLOCK_PARENT_MESHES);

  for(int i = 0; i < nb; i ++)ids[i] = bid_arr[i];
  
  return 0;
}

/*****************************************************************************/
int im_ex_get_ns_parent_mesh (int  /* exoid */, int *ids)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int nb = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_SETS);
  const long long * bid_arr = ms->getMSP(ms_lt::Mesh_Specification::NODESET_PARENT_MESHES);

  for(int i = 0; i < nb; i ++)ids[i] = bid_arr[i];
  
  return 0;
}
/*****************************************************************************/
int im_ex_get_ss_parent_mesh (int  /* exoid */, int *ids)
/*****************************************************************************/
{
  ms_lt::Mesh_Specification * ms = ms_lt::Mesh_Specification::first_ms_static_storage;
  if(!ms)return -1;

  int nb = ms->getMSI(ms_lt::Mesh_Specification::NUM_SIDE_SETS);
  const long long * bid_arr = ms->getMSP(ms_lt::Mesh_Specification::SIDESET_PARENT_MESHES);

  for(int i = 0; i < nb; i ++)ids[i] = bid_arr[i];
  
  return 0;
}
