// $Id$

#include <ctype.h>
#include <assert.h>

#include "../asrc/element_dictionary.h"
#include "pamgen_mesh_specification.h"
#include <string.h>

namespace ms_lt{

  Mesh_Specification * ms_lt::Mesh_Specification::static_storage = NULL;


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
void Mesh_Specification::
Specify_Side_Set_Information(long long i,
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
  register long long i;

  for(i = 0; i < msia[NUM_NODE_COMM_MAPS]; i++) {
    if(mspa[NODE_CMAP_NODE_CNTS][i]) {
      msppa[COMM_NODE_IDS][i] = new long long[mspa[NODE_CMAP_NODE_CNTS][i]];
      msppa[COMM_NODE_PROC_IDS][i] = new long long[mspa[NODE_CMAP_NODE_CNTS][i]];
    }
  }
   
  for(i = 0; i < msia[NUM_ELEM_COMM_MAPS]; i++) {
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
  register long long i;

  for(i = 0; i < msia[NUM_NODE_COMM_MAPS]; i++) {
    if(mspa[NODE_CMAP_NODE_CNTS][i]) {
      delete [] msppa[COMM_NODE_IDS][i];// = new int[mspa[NODE_CMAP_NODE_CNTS][i]];
      delete [] msppa[COMM_NODE_PROC_IDS][i];// = new int[mspa[NODE_CMAP_NODE_CNTS][i]];
    }
  }   
  delete [] msppa[COMM_NODE_IDS];
  delete [] msppa[COMM_NODE_PROC_IDS];

  for(i = 0; i < msia[NUM_ELEM_COMM_MAPS]; i++) {
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

}
