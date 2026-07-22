// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "../mesh_spec_lt/pamgen_mesh_specification.h"
#include "inline_mesh_desc.h"
#include "uns_inline_decomp.h"
#include <iostream>
#include <strings.h>
#include <time.h>
#include <cstring>

/****************************************************************************/
ms_lt::Mesh_Specification * buildMeshSpecification_LT(PAMGEN_NEVADA::Inline_Mesh_Desc* imd,
    long long rank,
    long long num_procs)
  /****************************************************************************/
{
  imd->my_rank = rank;
  imd->num_processors = num_procs;

  long long num_nodes_per_element = 8;
  Element_Type the_element_type = HEX8;
  long long num_nodes_per_face = 4;
  long long dim = imd->dimension;
  if(dim == 3){}
  else if(dim == 2){
    num_nodes_per_element = 4;
    the_element_type = QUAD4;
    num_nodes_per_face = 2;
  }
  else{
  }

  imd->setStrides();


  //Pre-Process BC's
  // set up nnx/nny/nnz for calculating bc loop limits
  long long nnx = imd->nel_tot[0]+1;
  long long nny = imd->nel_tot[1]+1;
  long long nnz = 1;
  if(dim == 3){
    nnz = imd->nel_tot[2]+1;
  }

  imd->Size_BC_Sets(nnx,nny,nnz);

  if(imd->inline_geometry_type == RADIAL && imd->periodic_j){
    nny = imd->nel_tot[1];
  }

  ms_lt::Mesh_Specification * nemesis_db = new ms_lt::Mesh_Specification();
  nemesis_db->setMSI(ms_lt::Mesh_Specification::PROC_ID, imd->my_rank);


  // The strategy is to implement serial with a trivial decomposition.
  // The trivial decomposition is to disperse the elements based on their
  // 'global' ordering in the entire i,j,k domain.
  // The strategy is to implement serial with a trivial decomposition.
  // The trivial decomposition is to disperse the elements based on their
  // 'global' ordering in the entire i,j,k domain.

  //make up list of global elements on processor
  std::set <long long> global_el_ids;
  long long local_ijk[3];
  long long global_ijk[6];
  for(int i=0; i<3; i++) local_ijk[i] = 0;
  for(int i=0; i<6; i++) global_ijk[i] = 0;
  {
    long long error_code;
    if(imd->inline_geometry_type == RADIAL_TRISECTION){
      error_code = imd->Decompose(global_el_ids);
    } else {
      error_code = imd->Decompose(global_el_ids,local_ijk,global_ijk);
    }
    if(error_code){
      delete nemesis_db;
      nemesis_db = NULL;
      return NULL;
    }
    nemesis_db->Set_Local_Num_IJK(local_ijk[0],local_ijk[1],local_ijk[2]);
    nemesis_db->Set_Global_IJK(global_ijk[0],global_ijk[1],global_ijk[2],
                               global_ijk[3],global_ijk[4],global_ijk[5]);
  }
  // Set IJK
  nemesis_db->Set_Total_Num_IJK(nnx,nny,nnz);

  // reads in all the serial component of the mesh specification
  // including the maps
  std::vector <long long> element_vector;
  std::vector <long long> global_node_vector;
  std::list <long long> global_node_list;
  std::map <long long, long long> global_node_map;//maps global node id to local ordinal
  std::map <long long, long long> global_element_map;//maps global node id to local ordinal

  imd->Build_Global_Lists(global_el_ids,
      element_vector,
      global_node_list,
      global_node_vector,
      global_node_map,
      global_element_map);

  nemesis_db->Specify_Global_Information(std::string ("PAMGEN Inline Mesh"),//title_string
      dim,//dimension
      global_node_list.size(),//total num nodes that are local
      element_vector.size(),//total num elements that are local
      imd->numBlocks(),//total num blocks in entire problem
      imd->nodeset_list.size(),//num nodesets in entire problem
      imd->sideset_list.size(),//num_sidesets in entire problem
      1,//num_qa records
      0// num_info_records
      );


  time_t tim = time(NULL);
  char * s = ctime(&tim);
  s[strlen(s)-1]=0;

  typedef std::string QA_Record[4];
  QA_Record * qa_recs = nemesis_db->QA_Records();
  qa_recs[0][0] = "PAMGEN";
  qa_recs[0][1] = "PArallel Mesh GENerator";
  qa_recs[0][2] = s;
  qa_recs[0][3] = s;


  std::string *coord_names = nemesis_db->getMSPSA(ms_lt::Mesh_Specification::COORDINATE_NAMES);
  coord_names[0] = "X";
  coord_names[1] = "Y";
  if(dim ==3){
    coord_names[2] = "Z";
  }

  imd->Calc_Serial_Component(global_el_ids,
      global_node_vector);

  long long err_code = imd->Calc_Coord_Vectors();
  if(err_code){
    delete nemesis_db;
    nemesis_db = NULL;
    return NULL;
  }

  imd->Populate_Coords(nemesis_db->Coord(),
      global_node_vector,
      global_node_map,
      global_node_list.size());

  imd->Offset_Coords(nemesis_db->Coord(),
      global_node_list.size(),
      dim);

  if(!imd->getErrorString().empty()){
    delete nemesis_db;
    nemesis_db = NULL;
    return NULL;
  }
  
  imd->Customize_Coords(nemesis_db->Coord(),
      global_node_list.size(),
      dim);

  imd->Display_Class(imd->info_stream,"");


  std::list < PG_BC_Specification *> ::iterator setit;
  long long nsct = 0;
  long long the_num_side_set_nodes = 0;
  for(setit = imd->sideset_list.begin(); setit != imd->sideset_list.end();setit++,nsct ++){
    nemesis_db->Specify_Side_Set_Information(nsct,//the index
        (*setit)->id,//the id
        imd->sideset_vectors[nsct].size(),// number of faces
        imd->sideset_vectors[nsct].size()*num_nodes_per_face,// number of ss nodes
        0 );
    the_num_side_set_nodes += imd->sideset_vectors[nsct].size()*num_nodes_per_face;
  }
  nemesis_db->setMSI(ms_lt::Mesh_Specification::NUM_SIDE_SET_NODES,the_num_side_set_nodes);

  err_code = imd->Populate_Sideset_Info( global_element_map,
      global_node_map,
      nemesis_db->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_ELEMENTS),
      nemesis_db->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_FACES),
      nemesis_db->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_NODES),
      nemesis_db->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_NODE_COUNTER));

  if(!imd->getErrorString().empty()){
    delete nemesis_db;
    nemesis_db = NULL;
    return NULL;
  }
  if(err_code){
    delete nemesis_db;
    nemesis_db = NULL;
    return NULL;
  }

  nsct = 0;

  for(setit = imd->nodeset_list.begin(); setit != imd->nodeset_list.end();setit++,nsct ++){
    nemesis_db->Specify_Node_Set_Information(nsct,//the index
        (*setit)->id,//the id
        imd->nodeset_vectors[nsct].size(),//
        0/*number_of_df*/);
  }

  imd->Populate_Nodeset_Info(nemesis_db->getMSPP(ms_lt::Mesh_Specification::NODE_SET_NODES),//Node_Set_Nodes(),
      global_node_map);
  if(!imd->getErrorString().empty()){
    delete nemesis_db;
    nemesis_db = NULL;
    return NULL;
  }
  //CONNECTIVITY
  for(long long i = 0;i <nemesis_db->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);i++ ){
    nemesis_db->Specify_Block_Information(i,//index
        i+imd->inline_block_start,//id
        imd->element_block_lists[i].size(),//number of elements
        num_nodes_per_element,//num nodes per element
        0,//number of element_attributes
        the_element_type);
  }

  imd->Populate_Connectivity(nemesis_db->getMSPP(ms_lt::Mesh_Specification::ELMT_NODE_LINKAGE),
      global_node_map);
  if(!imd->getErrorString().empty()){
    delete nemesis_db;
    nemesis_db = NULL;
    return NULL;
  }



  std::string *el_types = nemesis_db->getMSPSA(ms_lt::Mesh_Specification::ELEMENT_TYPES);
  for(long long bct = 0; bct < nemesis_db->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);bct ++ ){
    el_types[bct] = "QUAD";
    if(dim==3){
      el_types[bct] = "HEX";
    }
  }

  long long * the_map = nemesis_db->getMSP(ms_lt::Mesh_Specification::ELEM_ORDER_MAP);
  long long * global_element_numbers = nemesis_db->getMSP(ms_lt::Mesh_Specification::GLOBAL_ELEMENT_NUMBERS);
  imd->Populate_Map_and_Global_Element_List(the_map,
      global_element_numbers);

  //   Read_Global_Numbers();
  long long * global_node_numbers = nemesis_db->getMSP(ms_lt::Mesh_Specification::GLOBAL_NODE_NUMBERS);
  for(unsigned gnv = 0;gnv < global_node_vector.size();gnv ++){
    global_node_numbers[gnv] = global_node_vector[gnv]+1;
  }
  //   Read_Global_Info();
  {
    long long telc;
    long long tnoc;
    long long tect;
    imd->calculateSize(telc,
        tnoc,
        tect);

    nemesis_db->Global_Data_Size(tnoc,
        imd->GlobalNumElements(),
        imd->numBlocks(),
        imd->nodeset_list.size(),//num_nodesets
        imd->sideset_list.size(),
        imd->num_processors,
        imd->my_rank);

  }

  long long* elem_blk_ids_global =    nemesis_db->getMSP(ms_lt::Mesh_Specification::ELEM_BLK_IDS_GLOBAL);//Elem_Blk_Ids_Global();
  for(long long bct = 0; bct <  imd->numBlocks();bct++)elem_blk_ids_global[bct] = bct + imd->inline_block_start;// add 1 for block index

  long long* elem_blk_cnts_global =   nemesis_db->getMSP(ms_lt::Mesh_Specification::ELEM_BLK_CNTS_GLOBAL);
  imd->getGlobal_Element_Block_Totals(elem_blk_cnts_global);

  long long* ns_ids_global =          nemesis_db->getMSP(ms_lt::Mesh_Specification::NS_IDS_GLOBAL);
  long long* ns_cnts_global =         nemesis_db->getMSP(ms_lt::Mesh_Specification::NS_CNTS_GLOBAL);
  long long* ns_df_cnts_global =      nemesis_db->getMSP(ms_lt::Mesh_Specification::NS_DF_CNTS_GLOBAL);
  nsct = 0;

  for(setit = imd->nodeset_list.begin(); setit != imd->nodeset_list.end();setit++,nsct ++){
    ns_ids_global[nsct]=(*setit)->id;
    ns_cnts_global[nsct] = 0;
    for(unsigned ict = 0;ict < (*setit)->the_locs.size();ict ++){
      PAMGEN_NEVADA::LoopLimits ll = (*setit)->the_locs[ict].limits;
      ns_cnts_global[nsct] += ll.total;
    }
    ns_df_cnts_global[nsct] = 0;
  }


  long long* ss_ids_global =          nemesis_db->getMSP(ms_lt::Mesh_Specification::SS_IDS_GLOBAL);
  long long* ss_cnts_global =         nemesis_db->getMSP(ms_lt::Mesh_Specification::SS_CNTS_GLOBAL);
  long long* ss_df_cnts_global =      nemesis_db->getMSP(ms_lt::Mesh_Specification::SS_DF_CNTS_GLOBAL);
  nsct = 0;
  for(setit = imd->sideset_list.begin(); setit != imd->sideset_list.end();setit++,nsct ++){
    ss_ids_global[nsct]=(*setit)->id;
    ss_cnts_global[nsct] = imd->sideset_global_count[nsct];
    ss_df_cnts_global[nsct] = 0;
  }
  //Declare containers for par info Calculation
  std::list <long long> internal_node_list;
  std::list < long long > border_nodes_list;
  std::list <long long> internal_element_list;
  std::list < long long > border_elements_list;
  std::list <long long> node_proc_id_list;
  std::list <long long> element_proc_id_list;
  std::vector <long long> node_neighbor_vector;
  std::list <long long> * boundary_node_list = NULL;
  std::vector <long long> element_neighbor_vector;
  std::list <std::pair <long long ,Topo_Loc > > *boundary_element_list = NULL;

  imd->Calc_Parallel_Info(element_vector,
      global_node_vector,
      global_node_map,
      internal_node_list,
      border_nodes_list,
      internal_element_list,
      border_elements_list,
      node_proc_id_list,
      element_proc_id_list,
      node_neighbor_vector,
      boundary_node_list,
      element_neighbor_vector,
      boundary_element_list);

  nemesis_db->Parallel_Data_Size(internal_node_list.size(),
      border_nodes_list.size(),
      0,//num_external_nodes,
      internal_element_list.size(),
      border_elements_list.size(),
      node_proc_id_list.size(),
      element_proc_id_list.size()
      );


  imd->Populate_Border_Nodes_Elements( nemesis_db->getMSP(ms_lt::Mesh_Specification::INTERNAL_ELEMENTS),
      nemesis_db->getMSP(ms_lt::Mesh_Specification::INTERNAL_NODES),
      nemesis_db->getMSP(ms_lt::Mesh_Specification::BORDER_ELEMENTS),
      nemesis_db->getMSP(ms_lt::Mesh_Specification::BORDER_NODES),
      internal_node_list,
      border_nodes_list,
      internal_element_list,
      border_elements_list,
      global_node_map,
      global_element_map);

  imd->Populate_Cmap( nemesis_db->getMSP(ms_lt::Mesh_Specification::NODE_CMAP_NODE_CNTS),
      nemesis_db->getMSP(ms_lt::Mesh_Specification::NODE_CMAP_IDS),
      nemesis_db->getMSP(ms_lt::Mesh_Specification::ELEM_CMAP_ELEM_CNTS),
      nemesis_db->getMSP(ms_lt::Mesh_Specification::ELEM_CMAP_IDS),
      node_neighbor_vector,
      element_neighbor_vector,
      boundary_node_list,
      boundary_element_list);

  nemesis_db->Allocate_Parallel_Data();

  imd->Populate_Parallel_Info( nemesis_db->getMSPP(ms_lt::Mesh_Specification::COMM_NODE_IDS),
      nemesis_db->getMSPP(ms_lt::Mesh_Specification::COMM_NODE_PROC_IDS),
      nemesis_db->getMSPP(ms_lt::Mesh_Specification::COMM_ELEM_IDS),
      nemesis_db->getMSPP(ms_lt::Mesh_Specification::COMM_SIDE_IDS),
      nemesis_db->getMSPP(ms_lt::Mesh_Specification::COMM_ELEM_PROC_IDS),
      node_neighbor_vector,
      element_neighbor_vector,
      boundary_node_list,
      global_node_map,
      boundary_element_list,
      global_element_map);


  if(node_proc_id_list.size())
    delete [] boundary_node_list;
  if(element_proc_id_list.size())
    delete [] boundary_element_list;

  return nemesis_db;
}


/****************************************************************************/
ms_lt::Mesh_Specification * consolidateMeshSpecification_LT(ms_lt::Mesh_Specification *)
  /****************************************************************************/
{
  return NULL;
}
