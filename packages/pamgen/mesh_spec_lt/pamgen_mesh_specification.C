// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <ctype.h>
#include <assert.h>

#include "pamgen_element_dictionary.h"
#include "pamgen_mesh_specification.h"
#include <string.h>
#include <time.h>
#include <map>
#include <vector>

namespace ms_lt{

  Mesh_Specification * ms_lt::Mesh_Specification::first_ms_static_storage = NULL;


/*****************************************************************************/
Mesh_Specification::Mesh_Specification()
/*****************************************************************************/
{
  Zero_Set();
  suppress_warnings = false;  //Must be here not in Zero_Set()
}

/*****************************************************************************/
Mesh_Specification::~Mesh_Specification()
/*****************************************************************************/
{
  if(next)delete next;
  Free();
}

/*****************************************************************************/
bool Mesh_Specification::Are_Warnings_Suppressed() const
/*****************************************************************************/
{
  return suppress_warnings;
}

/*****************************************************************************/
void Mesh_Specification::Suppress_Warnings(long long logical)
/*****************************************************************************/
{
  suppress_warnings=logical;
}

/*****************************************************************************/
void Mesh_Specification::Resize_Info_Store(long long value)
/*****************************************************************************/
//  Resize the space for info records so more information
//    can be appended to the database
//  value is the new total size of the information section
{
  int j;

  int old_num_info_records = msia[NUM_INFO_RECORDS];
  msia[NUM_INFO_RECORDS] = value;

  std::string* old_info_strings = mspsa[INFO_STRINGS];
  mspsa[INFO_STRINGS] = new std::string[msia[NUM_INFO_RECORDS]];

  for(j = 0; j < old_num_info_records; j++) {
    mspsa[INFO_STRINGS][j] = old_info_strings[j];
  }

  delete [] old_info_strings;
}


/*****************************************************************************/
void Mesh_Specification::Specify_Global_Information(const std::string &titl,
                                                    long long dimen,
                                                    long long n_nodes,
                                                    long long n_elements,
                                                    long long n_blocks,
                                                    long long num_ns,
                                                    long long num_ss,
                                                    long long num_qa,
                                                    long long num_info)
/*****************************************************************************/
{
  int i;

  Free();

  title                    = titl;
  msia[DIM]                = dimen;
  msia[NUM_NODES]                = n_nodes;
  msia[NUM_ELEMENTS]             = n_elements;
  msia[NUM_BLOCKS]               = n_blocks;
  msia[NUM_NODE_SETS]            = num_ns;
  msia[NUM_SIDE_SETS]            = num_ss;


  msia[NUM_PARENT_MESHES] = 1;

  if(n_blocks)mspa[BLOCK_PARENT_MESHES]   = new long long[n_blocks];
  if(num_ns)  mspa[NODESET_PARENT_MESHES] = new long long[num_ns];
  if(num_ss)  mspa[SIDESET_PARENT_MESHES] = new long long[num_ss];



  if (msia[NUM_NODES] && msia[DIM])
    coord = new double[msia[NUM_NODES] * msia[DIM]];

  if (msia[NUM_ELEMENTS])
    mspa[ELEM_ORDER_MAP] = new long long[msia[NUM_ELEMENTS]];

  if (msia[NUM_BLOCKS] > 0)
  {
    mspa[BLOCK_ID]                 = new long long[msia[NUM_BLOCKS]];
    mspa[ELEMENTS_IN_BLOCK]                 = new long long[msia[NUM_BLOCKS]];
    mspa[NODES_PER_ELEMENT]        = new long long[msia[NUM_BLOCKS]];
    mspa[ELEMENT_ATTRIBUTES]       = new long long[msia[NUM_BLOCKS]];
    block_element_type       = new Element_Type[msia[NUM_BLOCKS]];
    mspsa[ELEMENT_TYPES]            = new std::string[msia[NUM_BLOCKS]];

    msppa[ELMT_NODE_LINKAGE]        = new long long*[msia[NUM_BLOCKS]];
    msppda[ATTRIBUTES]               = new double*[msia[NUM_BLOCKS]];
  }

// set indexes double dimensioned arrays to NULL

  for (long long b = 0; b < msia[NUM_BLOCKS]; ++b)
  {
    mspa[BLOCK_ID][b]           = 0;
    mspa[BLOCK_PARENT_MESHES][b] = 0;
    mspa[ELEMENTS_IN_BLOCK][b]           = 0;
    mspa[NODES_PER_ELEMENT][b]  = 0;
    mspa[ELEMENT_ATTRIBUTES][b] = 0;
    block_element_type[b] = UNKNOWN_ELEMENT;

    msppa[ELMT_NODE_LINKAGE][b]  = NULL;
    msppda[ATTRIBUTES][b]         = NULL;
  }

  if (msia[NUM_NODE_SETS] > 0)
  {
    mspa[NODE_SET_ID]              = new long long[msia[NUM_NODE_SETS]];
    mspa[NUM_NODES_IN_NODE_SET]    = new long long[msia[NUM_NODE_SETS]];
    mspa[NUM_DF_IN_NODE_SET]       = new long long[msia[NUM_NODE_SETS]];

    msppa[NODE_SET_NODES]           = new long long*[msia[NUM_NODE_SETS]];
    msppda[NODE_SET_DF]              = new double*[msia[NUM_NODE_SETS]];

    for(i = 0; i < msia[NUM_NODE_SETS]; i++) {
      mspa[NODESET_PARENT_MESHES][i] = 0;
      mspa[NODE_SET_ID][i]           = 0;
      mspa[NUM_NODES_IN_NODE_SET][i] = 0;
      mspa[NUM_DF_IN_NODE_SET][i]    = 0;
      msppa[NODE_SET_NODES][i]        = NULL;
      msppda[NODE_SET_DF][i]           = NULL;
    }
  }

  if (msia[NUM_SIDE_SETS] > 0)
  {
    mspa[SIDE_SET_ID]              = new long long[msia[NUM_SIDE_SETS]];
    mspa[NUM_ELEMENTS_IN_SIDE_SET] = new long long[msia[NUM_SIDE_SETS]];
    mspa[NUM_NODES_IN_SIDE_SET]    = new long long[msia[NUM_SIDE_SETS]];
    mspa[NUM_DF_IN_SIDE_SET]       = new long long[msia[NUM_SIDE_SETS]];

    msppa[SIDE_SET_ELEMENTS]        = new long long*[msia[NUM_SIDE_SETS]];
    msppa[SIDE_SET_NODE_COUNTER]    = new long long*[msia[NUM_SIDE_SETS]];
    msppa[SIDE_SET_FACES]           = new long long*[msia[NUM_SIDE_SETS]];
    msppa[SIDE_SET_NODES]           = new long long*[msia[NUM_SIDE_SETS]];
    msppda[SIDE_SET_DF]              = new double*[msia[NUM_SIDE_SETS]];

    for (i = 0; i < msia[NUM_SIDE_SETS]; ++i)
    {
      mspa[SIDESET_PARENT_MESHES][i]    = 0;
      mspa[SIDE_SET_ID][i]              = 0;
      mspa[NUM_ELEMENTS_IN_SIDE_SET][i] = 0;
      mspa[NUM_NODES_IN_SIDE_SET][i]    = 0;
      mspa[NUM_DF_IN_SIDE_SET][i]       = 0;
      msppa[SIDE_SET_ELEMENTS][i]        = NULL;
      msppa[SIDE_SET_NODE_COUNTER][i]    = NULL;
      msppa[SIDE_SET_FACES][i]           = NULL;
      msppa[SIDE_SET_NODES][i]           = NULL;
      msppda[SIDE_SET_DF][i]              = NULL;
    }
  }

  if (msia[NUM_ELEMENTS] > 0) mspa[GLOBAL_ELEMENT_NUMBERS]  = new long long[msia[NUM_ELEMENTS]];
  if (msia[NUM_NODES] > 0)    mspa[GLOBAL_NODE_NUMBERS]     = new long long[msia[NUM_NODES]];

  mspsa[COORDINATE_NAMES] = new std::string[msia[DIM]];

  msia[NUM_QA_RECORDS] = num_qa;
  msia[NUM_INFO_RECORDS] = num_info;

  if (num_qa > 0)
    qa_strings = new std::string[num_qa][4];

  if (num_info > 0)
    mspsa[INFO_STRINGS] = new std::string[num_info];
}


/*****************************************************************************/
void Mesh_Specification::Specify_Block_Information(long long b,
                          long long id,
                          long long number_of_block_elements,
                          long long number_of_nodes_per_element,
                          long long number_of_element_attributes,
                          Element_Type type)
/*****************************************************************************/
{
  assert(mspa[BLOCK_ID] != 0);
  assert(mspa[ELEMENTS_IN_BLOCK] != 0);
  assert(mspa[NODES_PER_ELEMENT] != 0);
  assert(mspa[ELEMENT_ATTRIBUTES] != 0);
  assert(block_element_type != 0);
  assert(msppa[ELMT_NODE_LINKAGE] != 0);
  assert(msppda[ATTRIBUTES] != 0);

  mspa[BLOCK_ID][b]           = id;
  mspa[ELEMENTS_IN_BLOCK][b]           = number_of_block_elements;
  mspa[NODES_PER_ELEMENT][b]  = number_of_nodes_per_element;
  mspa[ELEMENT_ATTRIBUTES][b] = number_of_element_attributes;
  block_element_type[b] = type;

  msppa[ELMT_NODE_LINKAGE][b]  = new long long[mspa[ELEMENTS_IN_BLOCK][b]*mspa[NODES_PER_ELEMENT][b]];
  msppda[ATTRIBUTES][b]         = new double[mspa[ELEMENTS_IN_BLOCK][b]*mspa[ELEMENT_ATTRIBUTES][b]];
}


/*****************************************************************************/
void Mesh_Specification::
Specify_Node_Set_Information(long long i,
                             long long id,
                             long long number_of_nodes_in_node_set,
                             long long number_of_df_in_node_set)
/*****************************************************************************/
{
  assert(mspa[NODE_SET_ID]);
  assert(mspa[NUM_NODES_IN_NODE_SET]);
  assert(mspa[NUM_DF_IN_NODE_SET]);
  assert(msppa[NODE_SET_NODES]);
  assert(msppda[NODE_SET_DF]);

  mspa[NODE_SET_ID][i]           = id;
  mspa[NUM_NODES_IN_NODE_SET][i] = number_of_nodes_in_node_set;
  mspa[NUM_DF_IN_NODE_SET][i]    = number_of_df_in_node_set;
  msppa[NODE_SET_NODES][i]        = new long long[number_of_nodes_in_node_set];
  msppda[NODE_SET_DF][i]           = new double[number_of_df_in_node_set];
}

/*****************************************************************************/
void Mesh_Specification::Specify_Side_Set_Information(long long i,
                             long long id,
                             long long number_of_faces_in_side_set,
                             long long number_of_nodes_in_side_set,
                             long long number_of_df_in_side_set)
/*****************************************************************************/
{
  assert(mspa[SIDE_SET_ID]);
  assert(mspa[NUM_ELEMENTS_IN_SIDE_SET]);
  assert(mspa[NUM_NODES_IN_SIDE_SET]);
  assert(mspa[NUM_DF_IN_SIDE_SET]);
  assert(msppa[SIDE_SET_ELEMENTS]);
  assert(msppa[SIDE_SET_NODE_COUNTER]);
  assert(msppa[SIDE_SET_FACES]);
  assert(msppa[SIDE_SET_NODES]);
  assert(msppda[SIDE_SET_DF]);

  mspa[SIDE_SET_ID][i]              = id;
  mspa[NUM_ELEMENTS_IN_SIDE_SET][i] = number_of_faces_in_side_set;
  mspa[NUM_NODES_IN_SIDE_SET][i]    = number_of_nodes_in_side_set;
  mspa[NUM_DF_IN_SIDE_SET][i]       = number_of_df_in_side_set;
  msppa[SIDE_SET_ELEMENTS][i]        = new long long[number_of_faces_in_side_set];
  msppa[SIDE_SET_NODE_COUNTER][i]    = new long long[number_of_faces_in_side_set];
  msppa[SIDE_SET_FACES][i]           = new long long[number_of_faces_in_side_set];
  msppa[SIDE_SET_NODES][i]           = new long long[number_of_nodes_in_side_set];
  msppda[SIDE_SET_DF][i]              = new double[number_of_df_in_side_set];
}

/*****************************************************************************/
void Mesh_Specification::Free_NonTransient_Storage()
/*****************************************************************************/
//  Frees all storage except that which contains data required for an EXODUS
//  time step dump.  The data items retained are:
//      long long *mspa[BLOCK_ID];
//      long long *mspa[ELEMENTS_IN_BLOCK];
{
  long long i, b;

  delete[] coord;
  delete[] mspa[ELEM_ORDER_MAP];
  delete[] mspa[NODES_PER_ELEMENT];
  delete[] mspa[ELEMENT_ATTRIBUTES];
  delete[] block_element_type;

  delete[] mspa[BLOCK_PARENT_MESHES];
  delete[] mspa[NODESET_PARENT_MESHES];
  delete[] mspa[SIDESET_PARENT_MESHES];


  if (msppa[ELMT_NODE_LINKAGE]){
    for (b = 0; b < msia[NUM_BLOCKS]; b++) {
      if (msppa[ELMT_NODE_LINKAGE][b]) delete[] msppa[ELMT_NODE_LINKAGE][b];
    }
    delete[] msppa[ELMT_NODE_LINKAGE];
  }

  if (msppda[ATTRIBUTES]){
    for (b = 0; b < msia[NUM_BLOCKS]; b++) {
      if (msppda[ATTRIBUTES][b]) delete[] msppda[ATTRIBUTES][b];
    }
    delete[] msppda[ATTRIBUTES];
  }


  delete[] mspa[NODE_SET_ID];
  delete[] mspa[NUM_NODES_IN_NODE_SET];
  delete[] mspa[NUM_DF_IN_NODE_SET];
  if (msppa[NODE_SET_NODES]){
    for (i=0; i<msia[NUM_NODE_SETS]; i++){
      delete[] msppa[NODE_SET_NODES][i];
      delete[] msppda[NODE_SET_DF][i];
    }
    delete[] msppa[NODE_SET_NODES];
    delete[] msppda[NODE_SET_DF];
  }

  delete[] mspa[SIDE_SET_ID];
  delete[] mspa[NUM_ELEMENTS_IN_SIDE_SET];
  delete[] mspa[NUM_NODES_IN_SIDE_SET];
  delete[] mspa[NUM_DF_IN_SIDE_SET];

  if (msppa[SIDE_SET_ELEMENTS]){
    for (i=0; i<msia[NUM_SIDE_SETS]; i++){
      delete[] msppa[SIDE_SET_ELEMENTS][i];
      delete[] msppa[SIDE_SET_NODE_COUNTER][i];
      delete[] msppa[SIDE_SET_NODES][i];
      delete[] msppda[SIDE_SET_DF][i];
      delete[] msppa[SIDE_SET_FACES][i];
    }
  }

  delete[] msppa[SIDE_SET_ELEMENTS];
  delete[] msppa[SIDE_SET_NODE_COUNTER];
  delete[] msppa[SIDE_SET_FACES];
  delete[] msppa[SIDE_SET_NODES];
  delete[] msppda[SIDE_SET_DF];

  delete[] qa_strings;

  delete[] mspsa[INFO_STRINGS];

  delete[] mspsa[COORDINATE_NAMES];

  delete[] mspsa[ELEMENT_TYPES];

  delete[] mspa[GLOBAL_ELEMENT_NUMBERS];
  delete[] mspa[GLOBAL_NODE_NUMBERS];


  if (mspa[NBR_PROC_LIST])  delete [] mspa[NBR_PROC_LIST];

  //Parallel data
  Free_Parallel_Data();
  Free_Locational_Data();
  Free_Global_Data();


  msia[NUM_NODE_SETS]            = 0;

  msia[NUM_SIDE_SETS]            = 0;


  coord                    = NULL;
  mspa[ELEM_ORDER_MAP]           = NULL;

  mspa[NODES_PER_ELEMENT]        = NULL;
  mspa[ELEMENT_ATTRIBUTES]       = NULL;
  block_element_type       = NULL;

  msppa[ELMT_NODE_LINKAGE]        = NULL;
  msppda[ATTRIBUTES]               = NULL;

  mspa[NODE_SET_ID]              = NULL;
  mspa[NUM_NODES_IN_NODE_SET]    = NULL;
  mspa[NUM_DF_IN_NODE_SET]       = NULL;
  msppa[NODE_SET_NODES]           = NULL;
  msppda[NODE_SET_DF]              = NULL;

  mspa[SIDE_SET_ID]              = NULL;
  mspa[NUM_ELEMENTS_IN_SIDE_SET] = NULL;
  mspa[NUM_NODES_IN_SIDE_SET]    = NULL;
  mspa[NUM_DF_IN_SIDE_SET]       = NULL;
  msppa[SIDE_SET_ELEMENTS]        = NULL;
  msppa[SIDE_SET_NODE_COUNTER]    = NULL;
  msppa[SIDE_SET_FACES]           = NULL;
  msppa[SIDE_SET_NODES]           = NULL;
  msppda[SIDE_SET_DF]              = NULL;

  msia[NUM_QA_RECORDS]           = 0;
  msia[NUM_INFO_RECORDS]         = 0;

  qa_strings               = NULL;
  mspsa[INFO_STRINGS]             = NULL;
  mspsa[COORDINATE_NAMES]         = NULL;
  mspsa[ELEMENT_TYPES]            = NULL;

  mspa[GLOBAL_ELEMENT_NUMBERS]   = NULL;
  mspa[GLOBAL_NODE_NUMBERS]      = NULL;

  msia[NUM_EDGES]                = 0;
  msia[NUM_FACES]                = 0;

  msia[NUM_NBR_PROCS]            = 0;
  mspa[NBR_PROC_LIST]            = NULL;

  return;
}

/*****************************************************************************/
void Mesh_Specification::Zero_Set()
/*****************************************************************************/
// Set the mesh object to correspond to the empty mesh.
{

  for(long long i = 0; i < NUM_MSIA;i ++)msia[i] = 0;
  for(long long i = 0; i < NUM_MSPA;i ++)mspa[i] = NULL;
  for(long long i = 0; i < NUM_MSPPA;i ++)msppa[i] = NULL;
  for(long long i = 0; i < NUM_MSPPDA;i ++)msppda[i] = NULL;
  for(long long i = 0; i < NUM_MSPSA;i ++)mspsa[i] = NULL;

  next = NULL;

  coord                    = NULL;
  qa_strings               = NULL;
  block_element_type       = NULL;

}

/*****************************************************************************/
void Mesh_Specification::Free()
/*****************************************************************************/
{
  Free_NonTransient_Storage();
  if (mspa[BLOCK_ID] != 0) delete [] mspa[BLOCK_ID];
  if (mspa[ELEMENTS_IN_BLOCK] != 0) delete [] mspa[ELEMENTS_IN_BLOCK];

  Zero_Set();
  return;
}

/*****************************************************************************/
void Mesh_Specification::Parallel_Data_Size(long long i_num_internal_nodes,
                                  long long i_num_border_nodes,
                                  long long i_num_external_nodes,
                                  long long i_num_internal_elems,
                                  long long i_num_border_elems,
                                  long long i_num_node_comm_maps,
                                  long long i_num_elem_comm_maps)
/*****************************************************************************/
{
  msia[NUM_INTERNAL_NODES] = i_num_internal_nodes;
  msia[NUM_BORDER_NODES] = i_num_border_nodes;
  msia[NUM_EXTERNAL_NODES] = i_num_external_nodes;
  msia[NUM_INTERNAL_ELEMS] = i_num_internal_elems;
  msia[NUM_BORDER_ELEMS] = i_num_border_elems;
  msia[NUM_NODE_COMM_MAPS] = i_num_node_comm_maps;
  msia[NUM_ELEM_COMM_MAPS] = i_num_elem_comm_maps;
  Allocate_Locational_Data();
  Allocate_LoadBal_Data();
}


/*****************************************************************************/
void Mesh_Specification::Allocate_Locational_Data()
/*****************************************************************************/
{
  if(msia[NUM_INTERNAL_ELEMS]) mspa[INTERNAL_ELEMENTS] = new long long[msia[NUM_INTERNAL_ELEMS]];
  if(msia[NUM_BORDER_ELEMS])   mspa[BORDER_ELEMENTS]   = new long long[msia[NUM_BORDER_ELEMS]];

  if(msia[NUM_INTERNAL_NODES]) mspa[INTERNAL_NODES]    = new long long[msia[NUM_INTERNAL_NODES]];
  if(msia[NUM_BORDER_NODES])   mspa[BORDER_NODES]      = new long long[msia[NUM_BORDER_NODES]];
  if(msia[NUM_EXTERNAL_NODES]) mspa[EXTERNAL_NODES]    = new long long[msia[NUM_EXTERNAL_NODES]];
}

/*****************************************************************************/
void Mesh_Specification::Free_Locational_Data()
/*****************************************************************************/
{
  if(msia[NUM_INTERNAL_ELEMS]) delete [] mspa[INTERNAL_ELEMENTS];// = new int[msia[NUM_INTERNAL_ELEMS]];
  if(msia[NUM_BORDER_ELEMS])   delete [] mspa[BORDER_ELEMENTS];//   = new int[msia[NUM_BORDER_ELEMS]];

  if(msia[NUM_INTERNAL_NODES]) delete [] mspa[INTERNAL_NODES];//    = new int[msia[NUM_INTERNAL_NODES]];
  if(msia[NUM_BORDER_NODES])   delete [] mspa[BORDER_NODES];//      = new int[msia[NUM_BORDER_NODES]];
  if(msia[NUM_EXTERNAL_NODES]) delete [] mspa[EXTERNAL_NODES];//    = new int[msia[NUM_EXTERNAL_NODES]];
}

/*****************************************************************************/
void Mesh_Specification::Allocate_LoadBal_Data()
/*****************************************************************************/
{
  if(msia[NUM_NODE_COMM_MAPS]) {
    mspa[NODE_CMAP_NODE_CNTS]  = new long long[msia[NUM_NODE_COMM_MAPS]];
    mspa[NODE_CMAP_IDS]        = new long long[msia[NUM_NODE_COMM_MAPS]];
    msppa[COMM_NODE_IDS]        = new long long*[msia[NUM_NODE_COMM_MAPS]];
    msppa[COMM_NODE_PROC_IDS]   = new long long*[msia[NUM_NODE_COMM_MAPS]];
  }

  if(msia[NUM_ELEM_COMM_MAPS]) {
    mspa[ELEM_CMAP_ELEM_CNTS]  = new long long[msia[NUM_ELEM_COMM_MAPS]];
    mspa[ELEM_CMAP_IDS]        = new long long[msia[NUM_ELEM_COMM_MAPS]];
    msppa[COMM_ELEM_IDS]        = new long long*[msia[NUM_ELEM_COMM_MAPS]];
    msppa[COMM_SIDE_IDS]        = new long long*[msia[NUM_ELEM_COMM_MAPS]];
    msppa[COMM_ELEM_PROC_IDS]   = new long long*[msia[NUM_ELEM_COMM_MAPS]];
  }
}

/*****************************************************************************/
void Mesh_Specification::Allocate_Global_Data()
/*****************************************************************************/
{
  if( msia[NUM_ELM_BLKS_GLOBAL] ){
    mspa[ELEM_BLK_IDS_GLOBAL] = new long long[msia[NUM_ELM_BLKS_GLOBAL]];
    mspa[ELEM_BLK_CNTS_GLOBAL] = new long long[msia[NUM_ELM_BLKS_GLOBAL]];
  }
  if( msia[NUM_NODE_SETS_GLOBAL] ){
    mspa[NS_IDS_GLOBAL] = new long long[msia[NUM_NODE_SETS_GLOBAL]];
    mspa[NS_CNTS_GLOBAL] = new long long[msia[NUM_NODE_SETS_GLOBAL]];
    mspa[NS_DF_CNTS_GLOBAL] = new long long[msia[NUM_NODE_SETS_GLOBAL]];
    for(long long i = 0; i < msia[NUM_NODE_SETS_GLOBAL]; i++){
      mspa[NS_IDS_GLOBAL][i] = 0;
      mspa[NS_DF_CNTS_GLOBAL][i] = 0;
      mspa[NS_DF_CNTS_GLOBAL][i] = 0;
    }
  }
  if( msia[NUM_SIDE_SETS_GLOBAL] ){
    mspa[SS_IDS_GLOBAL] = new long long[msia[NUM_SIDE_SETS_GLOBAL]];
    mspa[SS_CNTS_GLOBAL] = new long long[msia[NUM_SIDE_SETS_GLOBAL]];
    mspa[SS_DF_CNTS_GLOBAL] = new long long[msia[NUM_SIDE_SETS_GLOBAL]];
    for(long long i = 0; i < msia[NUM_SIDE_SETS_GLOBAL];i ++){
      mspa[SS_IDS_GLOBAL][i] = 0;
      mspa[SS_CNTS_GLOBAL][i] = 0;
      mspa[SS_DF_CNTS_GLOBAL][i] = 0;
    }
  }
}

/*****************************************************************************/
void Mesh_Specification::Free_Global_Data()
/*****************************************************************************/
{
  if( msia[NUM_ELM_BLKS_GLOBAL] ){
    delete [] mspa[ELEM_BLK_IDS_GLOBAL];// = new int[msia[NUM_ELM_BLKS_GLOBAL]];
    delete [] mspa[ELEM_BLK_CNTS_GLOBAL];// = new int[msia[NUM_ELM_BLKS_GLOBAL]];
  }
  if( msia[NUM_NODE_SETS_GLOBAL] ){
    delete [] mspa[NS_IDS_GLOBAL];// = new int[msia[NUM_NODE_SETS_GLOBAL]];
    delete [] mspa[NS_CNTS_GLOBAL];// = new int[msia[NUM_NODE_SETS_GLOBAL]];
    delete [] mspa[NS_DF_CNTS_GLOBAL];// = new int[msia[NUM_NODE_SETS_GLOBAL]];
  }
  if( msia[NUM_SIDE_SETS_GLOBAL] ){
    delete [] mspa[SS_IDS_GLOBAL];// = new int[msia[NUM_SIDE_SETS_GLOBAL]];
    delete [] mspa[SS_CNTS_GLOBAL];// = new int[msia[NUM_SIDE_SETS_GLOBAL]];
    delete [] mspa[SS_DF_CNTS_GLOBAL];// = new int[msia[NUM_SIDE_SETS_GLOBAL]];
  }
}

//DMH the Rank here is redundant...
/*****************************************************************************/
void Mesh_Specification::Global_Data_Size( long long Num_Nodes_Global, long long Num_Elems_Global,
                                  long long Num_Elem_Blks, long long Num_Node_Sets,
                                  long long Num_Side_Sets, long long Num_Total_Procs,
                                  long long Rank)
/*****************************************************************************/
{
  msia[NUM_NODES_GLOBAL] = Num_Nodes_Global;
  msia[NUM_ELEMS_GLOBAL] = Num_Elems_Global;
  msia[NUM_ELM_BLKS_GLOBAL] = Num_Elem_Blks;
  msia[NUM_NODE_SETS_GLOBAL] = Num_Node_Sets;
  msia[NUM_SIDE_SETS_GLOBAL] = Num_Side_Sets;
  msia[NUM_TOTAL_PROC] = Num_Total_Procs;
  msia[PROC_ID] = Rank;
  msia[NUM_PROC_IN_FILE] = 1;
  std::string fileType("p");
  file_type = fileType;
  Allocate_Global_Data();
}


/*****************************************************************************/
void Mesh_Specification::Allocate_Parallel_Data()
/*****************************************************************************/
{
  for(long long i = 0; i < msia[NUM_NODE_COMM_MAPS]; i++) {
    if(mspa[NODE_CMAP_NODE_CNTS][i]) {
      msppa[COMM_NODE_IDS][i] = new long long[mspa[NODE_CMAP_NODE_CNTS][i]];
      msppa[COMM_NODE_PROC_IDS][i] = new long long[mspa[NODE_CMAP_NODE_CNTS][i]];
    }
  }

  for(long long i = 0; i < msia[NUM_ELEM_COMM_MAPS]; i++) {
    if(mspa[ELEM_CMAP_ELEM_CNTS][i]) {
      msppa[COMM_ELEM_IDS][i] = new long long[mspa[ELEM_CMAP_ELEM_CNTS][i]];
      msppa[COMM_SIDE_IDS][i] = new long long[mspa[ELEM_CMAP_ELEM_CNTS][i]];
      msppa[COMM_ELEM_PROC_IDS][i] = new long long[mspa[ELEM_CMAP_ELEM_CNTS][i]];
    }
  }

  return;
}

/*****************************************************************************/
void Mesh_Specification::Free_Parallel_Data()
/*****************************************************************************/
{
  for(long long i = 0; i < msia[NUM_NODE_COMM_MAPS]; i++) {
    if(mspa[NODE_CMAP_NODE_CNTS][i]) {
      delete [] msppa[COMM_NODE_IDS][i];// = new int[mspa[NODE_CMAP_NODE_CNTS][i]];
      delete [] msppa[COMM_NODE_PROC_IDS][i];// = new int[mspa[NODE_CMAP_NODE_CNTS][i]];
    }
  }
  delete [] msppa[COMM_NODE_IDS];
  delete [] msppa[COMM_NODE_PROC_IDS];

  for(long long i = 0; i < msia[NUM_ELEM_COMM_MAPS]; i++) {
    if(mspa[ELEM_CMAP_ELEM_CNTS][i]) {
      delete [] msppa[COMM_ELEM_IDS][i];// = new int[mspa[ELEM_CMAP_ELEM_CNTS][i]];
      delete [] msppa[COMM_SIDE_IDS][i];// = new int[mspa[ELEM_CMAP_ELEM_CNTS][i]];
      delete [] msppa[COMM_ELEM_PROC_IDS][i];// = new int[mspa[ELEM_CMAP_ELEM_CNTS][i]];
    }
  }
  delete [] msppa[COMM_ELEM_IDS];
  delete [] msppa[COMM_SIDE_IDS];
  delete [] msppa[COMM_ELEM_PROC_IDS];

  if(msia[NUM_NODE_COMM_MAPS]) {
    delete [] mspa[NODE_CMAP_NODE_CNTS];//  = new int[msia[NUM_NODE_COMM_MAPS]];
    delete [] mspa[NODE_CMAP_IDS];//        = new int[msia[NUM_NODE_COMM_MAPS]];
  }

  if(msia[NUM_ELEM_COMM_MAPS]) {
    delete [] mspa[ELEM_CMAP_ELEM_CNTS];//  = new int[msia[NUM_ELEM_COMM_MAPS]];
    delete [] mspa[ELEM_CMAP_IDS];//        = new int[msia[NUM_ELEM_COMM_MAPS]];
  }

  return;
}

/*****************************************************************************/
Mesh_Specification * Mesh_Specification::consolidateMS()
/*****************************************************************************/
{
  if(this->Next() == NULL)return this;
  ms_lt::Mesh_Specification * ndb = new ms_lt::Mesh_Specification();
  if(!ndb)return ndb;
  /*global info*/
  {
    long long nn = 0;
    long long ne = 0;
    long long nns = 0;
    long long nss = 0;
    long long nb = 0;
    long long submesh_count = 0;
    Mesh_Specification * ms = this;
    while(ms){
      nn += ms->getMSI(NUM_NODES);
      ne += ms->getMSI(NUM_ELEMENTS);
      nb += ms->getMSI(NUM_BLOCKS);
      nns += ms->getMSI(NUM_NODE_SETS);
      nss += ms->getMSI(NUM_SIDE_SETS);
      ms = ms->Next();
    }
    ndb->Specify_Global_Information(std::string ("PAMGEN Consolidated Mesh"),
                                    this->getMSI(DIM),
                                    nn,//num nodes
                                    ne,//num elements
                                    nb,//num blocks
                                    nns,//num node sets
                                    nss,//num side sets
                                    1,//num qa recs
                                    0//num info recs
                                    );

    long long *  bp  = ndb->getMSP(ms_lt::Mesh_Specification::BLOCK_PARENT_MESHES);
    long long *  nsp = ndb->getMSP(ms_lt::Mesh_Specification::NODESET_PARENT_MESHES);
    long long *  ssp = ndb->getMSP(ms_lt::Mesh_Specification::SIDESET_PARENT_MESHES);

    ms = this;
    nb = 0;
    nns = 0;
    nss = 0;
    while(ms){
      for(int nbct = nb;nbct < nb+ ms->getMSI(NUM_BLOCKS);nbct++)      bp[nbct] = submesh_count;
      for(int nsct = nns;nsct < nns+ ms->getMSI(NUM_NODE_SETS);nsct++)nsp[nsct] = submesh_count;
      for(int ssct = nss;ssct < nss+ ms->getMSI(NUM_SIDE_SETS);ssct++)ssp[ssct] = submesh_count;
      nb += ms->getMSI(NUM_BLOCKS);
      nns += ms->getMSI(NUM_NODE_SETS);
      nss += ms->getMSI(NUM_SIDE_SETS);
      submesh_count ++;
      ms = ms->Next();
    }
    msia[ms_lt::Mesh_Specification::NUM_PARENT_MESHES] = submesh_count;
  }


  time_t tim = time(NULL);
  char * s = ctime(&tim);
  s[strlen(s)-1]=0;

  typedef std::string QA_Record[4];
  QA_Record * qa_recs = ndb->QA_Records();
  qa_recs[0][0] = "PAMGEN";
  qa_recs[0][1] = "PArallel Mesh GENerator";
  qa_recs[0][2] = s;
  qa_recs[0][3] = s;


  std::string *coord_names = ndb->getMSPSA(ms_lt::Mesh_Specification::COORDINATE_NAMES);
  coord_names[0] = "X";
  coord_names[1] = "Y";
  if(this->getMSI(DIM) ==3){
    coord_names[2] = "Z";
  }

  /*copy coordinates*/
  {
    double * dest = ndb->Coord();
    int dim = this->getMSI(DIM);
    Mesh_Specification * ms = this;
    long long total_nodes = ndb->getMSI(NUM_NODES);
    long long tct = 0;
    while (ms){
      double * source = ms->Coord();
      long long sub_total_nodes =  ms->getMSI(NUM_NODES);
      for(long long ict = 0; ict < sub_total_nodes;ict ++,tct ++){
        for(int dct = 0; dct < dim; dct ++){
          dest[tct + dct * total_nodes] = source[ict + dct * sub_total_nodes];
        }
      }
      ms = ms->next;
    }
  }
  //Side set information
  {
    long long nsct = 0;
    long long the_num_side_set_nodes = 0;
    Mesh_Specification * ms = this;
    while(ms){
      long long * ne_in_ss = ms->getMSP(NUM_ELEMENTS_IN_SIDE_SET);
      long long * nn_in_ss = ms->getMSP(NUM_NODES_IN_SIDE_SET);
      long long * ss_id = ms->getMSP(SIDE_SET_ID);

      for(int ssct = 0; ssct <  ms->getMSI(NUM_SIDE_SETS);ssct ++,nsct ++){
        ndb->Specify_Side_Set_Information(nsct,//the index
                                          ss_id[ssct],//the id
                                          ne_in_ss[ssct],
                                          nn_in_ss[ssct],
                                          0);
        the_num_side_set_nodes +=  nn_in_ss[ssct];
      }
      ms = ms->Next();
    }
    ndb->setMSI(ms_lt::Mesh_Specification::NUM_SIDE_SET_NODES,the_num_side_set_nodes);
  }

  {
    long long * const * dside_set_elements     = ndb->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_ELEMENTS);
    long long * const * dside_set_faces        = ndb->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_FACES);
    long long * const * dside_set_nodes        = ndb->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_NODES);
    long long * const * dside_set_node_counter = ndb->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_NODE_COUNTER);
    Mesh_Specification * ms = this;
    int nsct = 0;
    int node_offset = 0;
    int element_offset = 0;
    while(ms){
      long long * ne_in_ss = ms->getMSP(NUM_ELEMENTS_IN_SIDE_SET);
      long long * nn_in_ss = ms->getMSP(NUM_NODES_IN_SIDE_SET);

      long long * const * sside_set_elements     = ms->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_ELEMENTS);
      long long * const * sside_set_faces        = ms->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_FACES);
      long long * const * sside_set_nodes        = ms->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_NODES);
      long long * const * sside_set_node_counter = ms->getMSPP(ms_lt::Mesh_Specification::SIDE_SET_NODE_COUNTER);
      for(int ssct = 0; ssct <  ms->getMSI(NUM_SIDE_SETS);ssct ++,nsct ++){
        long long * sthe_elements = sside_set_elements[ssct];
        long long * sthe_faces = sside_set_faces[ssct];
        long long * sthe_nodes = sside_set_nodes[ssct];
        long long * sthe_node_counter = sside_set_node_counter[ssct];

        long long * dthe_elements = dside_set_elements[nsct];
        long long * dthe_faces = dside_set_faces[nsct];
        long long * dthe_nodes = dside_set_nodes[nsct];
        long long * dthe_node_counter = dside_set_node_counter[nsct];

        for(int tct = 0; tct < ne_in_ss[ssct];tct ++){
          dthe_elements[tct] = sthe_elements[tct]+element_offset;
          dthe_faces[tct] = sthe_faces[tct];
          dthe_node_counter[tct] = sthe_node_counter[tct];
        }
        for(int nct = 0; nct < nn_in_ss[ssct];nct++){
          dthe_nodes[nct] = sthe_nodes[nct]+node_offset;
        }

      }
      node_offset += ms->getMSI(NUM_NODES);
      element_offset += ms->getMSI(NUM_ELEMENTS);
      ms = ms->Next();
    }
  }

 //Node set information
  {
    long long nsct = 0;
    Mesh_Specification * ms = this;
    while(ms){
      long long * nn_in_ns = ms->getMSP(NUM_NODES_IN_NODE_SET);
      long long * ns_id = ms->getMSP(NODE_SET_ID);

      for(int ct = 0; ct <  ms->getMSI(NUM_NODE_SETS);ct ++,nsct ++){
        ndb->Specify_Node_Set_Information(nsct,//the index
                                          ns_id[ct],//the id
                                          nn_in_ns[ct],
                                          0);
      }
      ms = ms->Next();
    }
  }
  // populate nodeset info
  {
    Mesh_Specification * ms = this;
    int nsct = 0;
    long long node_offset = 0;
    long long * const * dnode_set_nodes = ndb->getMSPP(ms_lt::Mesh_Specification::NODE_SET_NODES);
    while(ms){
      long long * const * snode_set_nodes = ms->getMSPP(ms_lt::Mesh_Specification::NODE_SET_NODES);
      long long * num_nodes_in_set =  ms->getMSP(ms_lt::Mesh_Specification::NUM_NODES_IN_NODE_SET);

      for(int ct = 0; ct <  ms->getMSI(NUM_NODE_SETS);ct ++,nsct ++){
        long long * dthe_nodes = dnode_set_nodes[nsct];
        long long * sthe_nodes = snode_set_nodes[ct];
        for(int nct = 0;nct < num_nodes_in_set[ct];nct ++){
          dthe_nodes[nct] = sthe_nodes[nct]+node_offset;
        }
      }
      node_offset += ms->getMSI(NUM_NODES);
      ms = ms->Next();
    }
  }
  //connectivity
  {
    int dim = this->getMSI(DIM);
    Element_Type the_element_type = HEX8;
    if(dim == 2){
      the_element_type = QUAD4;
    }

    int block_count = 0;
    int block_offset = 0;
    int tblock_offset = 0;

    Mesh_Specification * ms = this;
    while(ms){
      long long snb = ms->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);
      long long * sbids =  ms->getMSP(ms_lt::Mesh_Specification::BLOCK_ID);
      long long * seib =  ms->getMSP(ms_lt::Mesh_Specification::ELEMENTS_IN_BLOCK);
      long long * snpe =  ms->getMSP(ms_lt::Mesh_Specification::NODES_PER_ELEMENT);
      for(int bct = 0; bct < snb; bct ++,block_count ++){
        ndb->Specify_Block_Information(block_count,// 0 based index
                                       sbids[bct] + block_offset,// 1 base count
                                       seib[bct],// element in block
                                       snpe[bct],
                                       0,// attributes
                                       the_element_type);
	tblock_offset = sbids[bct]+block_offset;
      }
      block_offset = tblock_offset;
      ms = ms->Next();
    }
  }
  //populate connectivity
  {
    int dim = this->getMSI(DIM);
    int num_nodes_per_element = 8;
    if(dim == 2){
      num_nodes_per_element = 4;
    }

    int block_count = 0;
    int node_offset = 0;
    long long * const * dlinks = ndb->getMSPP(ms_lt::Mesh_Specification::ELMT_NODE_LINKAGE);
    Mesh_Specification * ms = this;
    while(ms){
      long long * const * slinks = ms->getMSPP(ms_lt::Mesh_Specification::ELMT_NODE_LINKAGE);
      long long snb = ms->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);
      for(int bct = 0; bct < snb; bct ++,block_count ++){
        long long * sconn = slinks[bct];
        long long * dconn = dlinks[block_count];
        long long * elem_in_block = ms->getMSP(ms_lt::Mesh_Specification::ELEMENTS_IN_BLOCK);
        for(int elct = 0; elct < elem_in_block[bct];elct ++){
          for(int nct = 0; nct < num_nodes_per_element; nct ++){
            dconn[elct*num_nodes_per_element+nct] = sconn[elct*num_nodes_per_element+nct]+node_offset;
          }
        }
      }
      node_offset += ms->getMSI(NUM_NODES);
      ms = ms->Next();
    }
  }

  //element types
  {
    int dim = this->getMSI(DIM);
    std::string *el_types = ndb->getMSPSA(ms_lt::Mesh_Specification::ELEMENT_TYPES);
    for(long long bct = 0; bct < ndb->getMSI(ms_lt::Mesh_Specification::NUM_BLOCKS);bct ++ ){
      el_types[bct] = "QUAD";
      if(dim==3){
        el_types[bct] = "HEX";
      }
    }
  }

  //Global element numbers OFFSET IS WRONG SHOULD BE ACROSS ALL PROCS
  {
    long long * the_map = ndb->getMSP(ms_lt::Mesh_Specification::ELEM_ORDER_MAP);
    long long * dglobal_element_numbers = ndb->getMSP(ms_lt::Mesh_Specification::GLOBAL_ELEMENT_NUMBERS);
    long long element_offset = 0;
    long long tict = 0;
    Mesh_Specification * ms = this;
    while(ms){
      long long * sthe_map = ms->getMSP(ms_lt::Mesh_Specification::ELEM_ORDER_MAP);
      long long * sglobal_element_numbers = ms->getMSP(ms_lt::Mesh_Specification::GLOBAL_ELEMENT_NUMBERS);
      long long nels =  ms->getMSI(NUM_ELEMENTS);
      for(int ict = 0; ict < nels; ict++,tict++){
        the_map[tict] = sthe_map[ict] + element_offset;
        dglobal_element_numbers[tict] = sglobal_element_numbers[ict] + element_offset;
      }

      element_offset += ms->getMSI(NUM_ELEMS_GLOBAL);
      ms = ms->Next();
    }
  }
  //global node numbers OFFSET is WRONG SHOULD BE ACROSS ALL PROCESSORS
  {
    long long * dglobal_node_numbers = ndb->getMSP(ms_lt::Mesh_Specification::GLOBAL_NODE_NUMBERS);
    long long node_offset = 0;
    int tict = 0;
    Mesh_Specification * ms = this;
    while(ms){
      long long * sglobal_node_numbers = ms->getMSP(ms_lt::Mesh_Specification::GLOBAL_NODE_NUMBERS);
      long long nnodes =  ms->getMSI(NUM_NODES);
      for(int ict = 0; ict < nnodes; ict++,tict++){
        dglobal_node_numbers[tict] = sglobal_node_numbers[ict] + node_offset;
      }

      node_offset += ms->getMSI(NUM_NODES_GLOBAL);
      ms = ms->Next();
    }
  }


  //global information
  {
    long long gnodes = 0;
    long long gels = 0;
    long long gblks = 0;
    long long gns = 0;
    long long gss = 0;
    Mesh_Specification * ms = this;
    while(ms){
      gnodes += ms->getMSI(NUM_NODES_GLOBAL);
      gels   += ms->getMSI(NUM_ELEMS_GLOBAL);
      gblks  += ms->getMSI(NUM_ELM_BLKS_GLOBAL);
      gns    += ms->getMSI(NUM_NODE_SETS_GLOBAL);
      gss    += ms->getMSI(NUM_SIDE_SETS_GLOBAL);

      ms = ms->Next();
    }
    ndb->Global_Data_Size(gnodes,
                          gels,
                          gblks,
                          gns,
                          gss,
                          this->getMSI(NUM_TOTAL_PROC),
                          this->getMSI(PROC_ID));

    long long* elem_blk_ids_global =    ndb->getMSP(ms_lt::Mesh_Specification::ELEM_BLK_IDS_GLOBAL);//Elem_Blk_Ids_Global();
    for(long long bct = 0; bct <  gblks; bct++)elem_blk_ids_global[bct] = bct+1;// add 1 for block index

  }

  // element block counts
  {
    long long* delem_blk_cnts_global =   ndb->getMSP(ms_lt::Mesh_Specification::ELEM_BLK_CNTS_GLOBAL);
    int blk_cnt = 0;
    Mesh_Specification * ms = this;
    while(ms){
      long long* sdelem_blk_cnts_global =   ms->getMSP(ms_lt::Mesh_Specification::ELEM_BLK_CNTS_GLOBAL);
      int nblk = ms->getMSI(NUM_ELM_BLKS_GLOBAL);
      for(int ict = 0; ict < nblk; ict ++,blk_cnt++){
        delem_blk_cnts_global[blk_cnt] = sdelem_blk_cnts_global[ict];
      }

      ms = ms->Next();
    }

  }

  // ns global info
  {
    long long* dns_ids_global =          ndb->getMSP(ms_lt::Mesh_Specification::NS_IDS_GLOBAL);
    long long* dns_cnts_global =         ndb->getMSP(ms_lt::Mesh_Specification::NS_CNTS_GLOBAL);
    long long* dns_df_cnts_global =      ndb->getMSP(ms_lt::Mesh_Specification::NS_DF_CNTS_GLOBAL);
    int nsct = 0;
    Mesh_Specification * ms = this;
    while(ms){
      long long* sns_ids_global =          ms->getMSP(ms_lt::Mesh_Specification::NS_IDS_GLOBAL);
      long long* sns_cnts_global =         ms->getMSP(ms_lt::Mesh_Specification::NS_CNTS_GLOBAL);
      long long* sns_df_cnts_global =      ms->getMSP(ms_lt::Mesh_Specification::NS_DF_CNTS_GLOBAL);
      int nns = ms->getMSI(ms_lt::Mesh_Specification::NUM_NODE_SETS_GLOBAL);
      for(int ict = 0; ict < nns; ict ++,nsct ++){
        dns_ids_global[nsct]     = sns_ids_global[ict];
        dns_cnts_global[nsct]    = sns_cnts_global[ict];
        dns_df_cnts_global[nsct] = sns_df_cnts_global[ict];
      }
      ms = ms->Next();
    }
  }

  //ss global info
  {
    long long* dss_ids_global =          ndb->getMSP(ms_lt::Mesh_Specification::SS_IDS_GLOBAL);
    long long* dss_cnts_global =         ndb->getMSP(ms_lt::Mesh_Specification::SS_CNTS_GLOBAL);
    long long* dss_df_cnts_global =      ndb->getMSP(ms_lt::Mesh_Specification::SS_DF_CNTS_GLOBAL);
    int ssct = 0;
    Mesh_Specification * ms = this;
    while(ms){
      long long* sss_ids_global =          ms->getMSP(ms_lt::Mesh_Specification::SS_IDS_GLOBAL);
      long long* sss_cnts_global =         ms->getMSP(ms_lt::Mesh_Specification::SS_CNTS_GLOBAL);
      long long* sss_df_cnts_global =      ms->getMSP(ms_lt::Mesh_Specification::SS_DF_CNTS_GLOBAL);
      int nss = ms->getMSI(ms_lt::Mesh_Specification::NUM_SIDE_SETS_GLOBAL);
      for(int ict = 0; ict < nss; ict ++,ssct ++){
        dss_ids_global[ssct]     = sss_ids_global[ict];
        dss_cnts_global[ssct]    = sss_cnts_global[ict];
        dss_df_cnts_global[ssct] = sss_df_cnts_global[ict];
      }
      ms = ms->Next();
    }
  }

  //parallel data
  {
    int internal_node_list_size = 0;
    int border_node_list_size = 0;
    int internal_element_list_size = 0;
    int border_element_list_size = 0;
    std::map<long,std::vector < long > > node_proc_ids;
    std::map<long,std::vector < std::pair <long/*element*/,long/*face*/ > > > element_proc_ids;
    std::vector<long>  elem_cmap_elem_cnts;

    Mesh_Specification * ms = this;
    int node_offset = 0;
    int element_offset = 0;
    while(ms){

      internal_node_list_size    += ms->getMSI(NUM_INTERNAL_NODES);
      border_node_list_size      += ms->getMSI(NUM_BORDER_NODES);
      internal_element_list_size += ms->getMSI(NUM_INTERNAL_ELEMS);
      border_element_list_size   += ms->getMSI(NUM_BORDER_ELEMS);

      long long * const * scnid    = ms->getMSPP(ms_lt::Mesh_Specification::COMM_NODE_IDS);
      long long * const * sceid    = ms->getMSPP(ms_lt::Mesh_Specification::COMM_ELEM_IDS);
      long long * const * scsid    = ms->getMSPP(ms_lt::Mesh_Specification::COMM_SIDE_IDS);

      {
        int nncm = ms->getMSI(NUM_NODE_COMM_MAPS);
        long long * ncmap_ids = ms->getMSP(NODE_CMAP_IDS);
        long long * ncmap_cnts = ms->getMSP(NODE_CMAP_NODE_CNTS);
        for(int ict = 0; ict < nncm; ict ++){
          long long proc_id = ncmap_ids[ict];
          long long nnodes =  ncmap_cnts[ict];

	  for (long long nctr = 0; nctr < nnodes; nctr ++){
	    (node_proc_ids[proc_id]).push_back(scnid[ict][nctr]+node_offset);
	    assert(proc_id == ms->getMSPP(ms_lt::Mesh_Specification::COMM_NODE_PROC_IDS)[ict][nctr]);
          }
        }
      }

      
      {
        int necm = ms->getMSI(NUM_ELEM_COMM_MAPS);
        long long * ecmap_ids = ms->getMSP(ELEM_CMAP_IDS);
        long long * ecmap_cnts = ms->getMSP(ELEM_CMAP_ELEM_CNTS);
        for(int ict = 0; ict < necm; ict ++){
          long long proc_id = ecmap_ids[ict];
          long long nels =  ecmap_cnts[ict];
	  for(long long ectr = 0; ectr < nels;ectr++){
	    (element_proc_ids[proc_id]).push_back(std::pair<long,long>(sceid[ict][ectr]+element_offset,scsid[ict][ectr]));
	    assert(proc_id == ms->getMSPP(ms_lt::Mesh_Specification::COMM_ELEM_PROC_IDS)[ict][ectr]);
          }
        }
      }

      node_offset += ms->getMSI(NUM_NODES);
      element_offset += ms->getMSI(NUM_ELEMENTS);
      ms = ms->Next();
    }  /* while(ms)*/

    ndb->Parallel_Data_Size(internal_node_list_size,
                            border_node_list_size,
                            0,//num_external_nodes,
                            internal_element_list_size,
                            border_element_list_size,
                            node_proc_ids.size(),
                            element_proc_ids.size()
                            );

    {
      long long * ncnc = ndb->getMSP(ms_lt::Mesh_Specification::NODE_CMAP_NODE_CNTS);
      long long * nci  = ndb->getMSP(ms_lt::Mesh_Specification::NODE_CMAP_IDS);
      long long * ecec = ndb->getMSP(ms_lt::Mesh_Specification::ELEM_CMAP_ELEM_CNTS);
      long long * eci  = ndb->getMSP(ms_lt::Mesh_Specification::ELEM_CMAP_IDS);
      std::map<long,std::vector<long>>::iterator mit;
      {
        int count = 0;
	for(mit = node_proc_ids.begin();mit != node_proc_ids.end();mit ++,count ++){
	  ncnc[count] = (*mit).second.size();//node_cmap_node_cnts[count];
	  nci[count]  = (*mit).first;
        }
      }
      std::map<long,std::vector < std::pair < long,long > > >::iterator it;

      {
        int count = 0;
        for(it = element_proc_ids.begin();it != element_proc_ids.end();it ++,count ++){
	  ecec[count] = (*it).second.size();//elem_cmap_elem_cnts[count];
          eci[count]  = (*it).first;
        }
      }
    }


    ndb->Allocate_Parallel_Data();
    {

      long long * const * dcnid    = ndb->getMSPP(ms_lt::Mesh_Specification::COMM_NODE_IDS);
      long long * const * dcnpid   = ndb->getMSPP(ms_lt::Mesh_Specification::COMM_NODE_PROC_IDS);
      long long * const * dceid    = ndb->getMSPP(ms_lt::Mesh_Specification::COMM_ELEM_IDS);
      long long * const * dcsid    = ndb->getMSPP(ms_lt::Mesh_Specification::COMM_SIDE_IDS);
      long long * const * dcepid   = ndb->getMSPP(ms_lt::Mesh_Specification::COMM_ELEM_PROC_IDS);

        {
	std::map<long,std::vector<long>>::iterator mit;
	int lcount = 0;
	for(mit = node_proc_ids.begin();mit != node_proc_ids.end();mit ++,lcount ++){
	  std::vector<long>::iterator vit;
	  int count = 0;
	  for(vit = (*mit).second.begin();vit != (*mit).second.end();vit ++,count ++){
	    dcnid[lcount][count] = *vit;
	    dcnpid[lcount][count] = (*mit).first;
            }
          }
        }

        {
	std::map< long, std::vector < std::pair < long, long > > > :: iterator mit;
	int lcount = 0;
	for(mit = element_proc_ids.begin();mit != element_proc_ids.end();mit ++,lcount ++){
	  std::vector< std::pair < long, long >>::iterator vit;
	  int count = 0;
	  for(vit = (*mit).second.begin();vit != (*mit).second.end();vit ++,count ++){
	    dceid[lcount][count] = (*vit).first;
	    dcsid[lcount][count] = (*vit).second;
	    dcepid[lcount][count] = (*mit).first;
          }
        }
      }
    }
  }

  //populate border nodes elements
  {
    long long * dint_elems    = ndb->getMSP(ms_lt::Mesh_Specification::INTERNAL_ELEMENTS);
    long long * dint_nodes    = ndb->getMSP(ms_lt::Mesh_Specification::INTERNAL_NODES);
    long long * dborder_elems = ndb->getMSP(ms_lt::Mesh_Specification::BORDER_ELEMENTS);
    long long * dborder_nodes = ndb->getMSP(ms_lt::Mesh_Specification::BORDER_NODES);

    int nin = 0;
    int nbn = 0;
    int nie = 0;
    int nbe = 0;

    Mesh_Specification * ms = this;
    int node_offset = 0;
    int element_offset = 0;
    while(ms){
      long long * sint_elems    = ms->getMSP(ms_lt::Mesh_Specification::INTERNAL_ELEMENTS);
      long long * sint_nodes    = ms->getMSP(ms_lt::Mesh_Specification::INTERNAL_NODES);
      long long * sborder_elems = ms->getMSP(ms_lt::Mesh_Specification::BORDER_ELEMENTS);
      long long * sborder_nodes = ms->getMSP(ms_lt::Mesh_Specification::BORDER_NODES);

      int internal_node_list_size    = ms->getMSI(NUM_INTERNAL_NODES);
      int border_node_list_size      = ms->getMSI(NUM_BORDER_NODES);
      int internal_element_list_size = ms->getMSI(NUM_INTERNAL_ELEMS);
      int border_element_list_size   = ms->getMSI(NUM_BORDER_ELEMS);

      for(int ict = 0; ict < internal_node_list_size;    ict ++, nin++){dint_nodes[nin]    = sint_nodes[ict]+node_offset;}
      for(int ict = 0; ict < border_node_list_size;      ict ++, nbn++){dborder_nodes[nbn] = sborder_nodes[ict]+node_offset;}
      for(int ict = 0; ict < internal_element_list_size; ict ++, nie++){dint_elems[nie]    = sint_elems[ict]+element_offset;}
      for(int ict = 0; ict < border_element_list_size;   ict ++, nbe++){dborder_elems[nbe] = sborder_elems[ict]+element_offset;}

      node_offset    += ms->getMSI(NUM_NODES);
      element_offset += ms->getMSI(NUM_ELEMENTS);

      ms = ms->Next();
    }
  }




  return ndb;
}


}
