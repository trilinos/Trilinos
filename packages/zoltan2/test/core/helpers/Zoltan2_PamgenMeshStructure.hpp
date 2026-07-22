// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* \file Zoltan2_PamgenMeshStructure.hpp
* \brief A class for copying the state of a Pamgen mesh for easy access 
*/

#ifndef PAMGEN_MESH_STRUCTURE
#define PAMGEN_MESH_STRUCTURE

#ifdef HAVE_ZOLTAN2_PAMGEN

#include <Pamgen_config.h>
#include <create_inline_mesh.h>
#include <limits.h>
#include <pamgen_im_exodusII.h>
#include <pamgen_im_ne_nemesisI.h>


class PamgenMesh{
public:
  
  /* \brief Destructor */
  ~PamgenMesh(); // free memory
  
  /* \bried Method for calling the Pamgen Mesh constructor
   * \param[in] file_data is the input file describing the pamgen mesh read into a char array
   * \param[in] dimension is the dimension of the problem, i.e. 1,2, or 3
   * \param[in] comm is the process communicator
   */
  void createMesh(char * file_data, int dimension, const RCP<const Comm<int>> &comm);

  /* \bried Method for copying the state of a Pamgen mesh once it has been created
   */
  void storeMesh(); // read mesh to memory
  
  /* Metod for computing the coordinates of the element centers */
  void computeElementCoordinates();

  int num_dim;
  int num_nodes;
  int num_elem;
  int num_elem_blk;
  int num_node_sets;
  int num_side_sets;
  int num_node_set_nodes;
  int num_node_set_dfs;
  int num_side_set_elements;
  int num_side_set_nodes;
  int num_side_set_dfs;
  int num_block_properties;
  int num_node_set_properties;
  int num_side_set_properties;
  
  char title[MAX_STR_LENGTH];
  
  int version;
  double version_number;
  double * coord;
  double * element_coord;
  
  char buffer[3][MAX_STR_LENGTH + 1];
  char *bptr[3];
  
  int * element_order_map ;
  
  int * global_element_numbers ;
  int * global_node_numbers ;
  
  /*block info*/
  int * block_id ;
  char ** element_types ;
  int *   elements ;
  int *   nodes_per_element ;
  int *   element_attributes ;
  int ** elmt_node_linkage ;
  
  /*side sets*/
  int * side_set_id ;
  int * num_elements_in_side_set ;
  int * num_df_in_side_set ;
  int **side_set_elements ;
  int **side_set_faces ;
  
  /*node sets*/
  int * node_set_id ;
  int * num_nodes_in_node_set ;
  int * num_df_in_node_set ;
  int **node_set_nodes  ;
  
  /*qa*/
  int num_qa_records;
  int num_info_records;
  char* qaRecord[100][4];
  char** info_records ;
  
  /*nemesis data*/
  int num_nodes_global;
  int num_elems_global;
  int num_elm_blks_global;
  int num_node_sets_global;
  int num_side_sets_global;
  int num_total_proc;
  int num_proc_in_file;
  char type[2];
  
  /*nemesis data
   global info*/
  
  int * elem_blk_ids_global ;
  int * elem_blk_cnts_global  ;
  
  int * ns_ids_global ;
  int * ns_cnts_global ;
  int * ns_df_cnts_global ;
  int * ss_ids_global ;
  int * ss_cnts_global ;
  int * ss_df_cnts_global ;
  
  /*parallel info*/
  int num_internal_nodes;
  int num_border_nodes;
  int num_external_nodes;
  int num_internal_elems;
  int num_border_elems;
  int num_node_comm_maps;
  int num_elem_comm_maps;
  
  int * internal_elements ;
  int * border_elements ;
  int * internal_nodes ;
  int * border_nodes ;
  int * external_nodes ;
  
  int * node_cmap_node_cnts ;
  int * node_cmap_ids       ;
  int * elem_cmap_elem_cnts ;
  int * elem_cmap_ids       ;
  
  int ** comm_node_ids       ;
  int ** comm_node_proc_ids  ;
  int ** comm_elem_ids       ;
  int ** comm_side_ids       ;
  int ** comm_elem_proc_ids  ;
  
};


PamgenMesh::~PamgenMesh()
{
  // free mesh
  Delete_Pamgen_Mesh();
  
  // free storage
  int i;
  int j;
  int b;
  for( i = 0; i < 100; i++){
    for( j=0; j<4; j++){
      free(this->qaRecord[i][j]);
    }
  }
  
  free(this->coord); // free vertex coords
  free(this->element_coord);
  
  
  if (this->num_elem){
    free(this->element_order_map);
    
    if (this->num_elem){
      free(this->global_element_numbers);
    }
    
    if (this->num_nodes){
      free(this->global_node_numbers);
    }
    
    
    /*block info*/
    free(this->block_id);
    free(this->nodes_per_element);
    free(this->element_attributes);
    free(this->elements);
    
    
    for(i = 0; i < this->num_elem_blk; i ++){
      free(this->element_types[i]);
    }
    
    /*connectivity*/
    for(b = 0; b < this->num_elem_blk; b++){
      free(this->elmt_node_linkage[b]);
    }
    free(this->element_types);
    free(this->elmt_node_linkage);
    
    
    if(this->num_node_sets){
      
      for(i = 0; i < this->num_node_sets; i ++){
        if(this->num_nodes_in_node_set[i]) {
          free(this->node_set_nodes[i]);
        }
      }
      free(this->node_set_id);
      free(this->num_nodes_in_node_set);
      free(this->node_set_nodes);
      free(this->num_df_in_node_set);
      
    }
    
    /*side sets*/
    if(this->num_side_sets){
      
      for(i = 0; i < this->num_side_sets; i ++){
        
        free(this->side_set_elements[i]);
        free(this->side_set_faces[i]);
        
      }
      free(this->side_set_id);
      free(this->num_elements_in_side_set);
      free(this->num_df_in_side_set);
      free(this->side_set_elements);
      free(this->side_set_faces);
    }
    
    
    if(this->num_info_records) {
      for(i = 0; i < this->num_info_records; i ++){
        free(this->info_records[i]);
      }
      free(this->info_records);
    }
    
    
    /*nemesis data
     global info */
    free(this->elem_blk_ids_global);
    free(this->elem_blk_cnts_global);
    
    free(this->ns_ids_global);
    free(this->ns_cnts_global);
    free(this->ns_df_cnts_global);
    free(this->ss_ids_global);
    free(this->ss_cnts_global);
    free(this->ss_df_cnts_global);
    
    free(this->internal_elements);
    free(this->border_elements);
    free(this->internal_nodes);
    free(this->border_nodes);
    free(this->external_nodes);
    
    if(this->num_node_comm_maps > 0){
      
      for(j = 0; j < this->num_node_comm_maps; j++) {
        free(this->comm_node_ids[j]);
        free(this->comm_node_proc_ids[j]);
        
      }
      
      
      for(j = 0; j < this->num_elem_comm_maps; j++) {
        free(this->comm_elem_ids[j]);
        free(this->comm_side_ids[j]);
        free(this->comm_elem_proc_ids[j]);
        
      }/*loop over num_elem_co*/
    }
    
    if(this->num_node_comm_maps > 0)
    {
      free(this->node_cmap_node_cnts);
      free(this->node_cmap_ids);
      free(this->comm_node_ids);
      free(this->comm_node_proc_ids);
      
      free(this->elem_cmap_elem_cnts);
      free(this->elem_cmap_ids);
      free(this->comm_elem_ids);
      free(this->comm_side_ids);
      free(this->comm_elem_proc_ids);
    }
    
  }
  
}

void PamgenMesh::storeMesh()
{
  
  int idum = 0;
  float fdum;
  int i;
  int j;
  int b;
  char * cdum = NULL;
  int error = 0;
  int id = 0;
  
  this->bptr[0] = this->buffer[0];
  this->bptr[1] = this->buffer[1];
  this->bptr[2] = this->buffer[2];
  
  for(i = 0; i < 100; i++)
    for(j=0; j<4; j++) this->qaRecord[i][j] = (char*)malloc(MAX_STR_LENGTH+1) ;
  
  
  error += im_ex_get_init (  id,
                           this->title,
                           &this->num_dim,
                           &(this->num_nodes),
                           &this->num_elem,
                           &this->num_elem_blk,
                           &this->num_node_sets,
                           &this->num_side_sets);
  
  
  error += im_ex_inquire(id, IM_EX_INQ_NS_NODE_LEN, (int*)&this->num_node_set_nodes,
                         &fdum, cdum);
  error += im_ex_inquire
  (id, IM_EX_INQ_NS_DF_LEN,   (int*)&this->num_node_set_dfs,
   &fdum, cdum);
  error += im_ex_inquire(id, IM_EX_INQ_SS_ELEM_LEN, (int*)&this->num_side_set_elements,
                         &fdum, cdum);
  error += im_ex_inquire(id, IM_EX_INQ_SS_NODE_LEN, (int*)&this->num_side_set_nodes,
                         &fdum, cdum);
  error += im_ex_inquire(id, IM_EX_INQ_SS_DF_LEN,   (int*)&this->num_side_set_dfs,
                         &fdum, cdum);
  
  /* get version number */
  
  error += im_ex_inquire(id, IM_EX_INQ_API_VERS, &idum, &fdum, cdum);
  
  this->version_number = (double) fdum;
  
  this->version = (int) this->version_number;
  
  /* get genesis-II parameters */
  
  error += im_ex_inquire(id, IM_EX_INQ_EB_PROP, (int*)&this->num_block_properties, &fdum, cdum);
  
  error += im_ex_inquire(id, IM_EX_INQ_NS_PROP, (int*)&this->num_node_set_properties,
                         &fdum, cdum);
  
  error += im_ex_inquire(id, IM_EX_INQ_SS_PROP, (int*)&this->num_side_set_properties,
                         &fdum, cdum);
  
  this->coord = (double *)malloc(this->num_nodes*this->num_dim*sizeof(double));
  
  error += im_ex_get_coord(id,this->coord,this->coord+this->num_nodes,this->coord+2*this->num_nodes);
  
  
  error += im_ex_get_coord_names (id, this->bptr);
  
  if (this->num_elem){
    this->element_order_map = (int *)malloc(this->num_elem * sizeof(int));
    error += im_ex_get_map(id, this->element_order_map);
    
    if (this->num_elem){
      this->global_element_numbers = (int *)malloc(this->num_elem*sizeof(int));
      error += im_ex_get_elem_num_map(id, this->global_element_numbers);
    }
    
    if (this->num_nodes){
      this->global_node_numbers = (int *)malloc(this->num_nodes * sizeof(int));
      error += im_ex_get_node_num_map(id, this->global_node_numbers);
    }
    
    
    /*block info*/
    
    this->block_id           = (int *)malloc(this->num_elem_blk*sizeof(int));
    this->nodes_per_element  = (int *)malloc(this->num_elem_blk*sizeof(int));
    this->element_attributes = (int *)malloc(this->num_elem_blk*sizeof(int));
    this->elements           = (int *)malloc(this->num_elem_blk*sizeof(int));
    this->element_types      = (char **)malloc(this->num_elem_blk*sizeof(char *));
    this->elmt_node_linkage  = (int **)malloc(this->num_elem_blk*sizeof(int*));
    
    error += im_ex_get_elem_blk_ids(id, this->block_id);
    
    for(i = 0; i < this->num_elem_blk; i ++){
      this->element_types[i] = (char *)malloc((MAX_STR_LENGTH + 1)*sizeof(char));
      error += im_ex_get_elem_block(id,
                                    this->block_id[i],
                                    this->element_types[i],
                                    (int*)&(this->elements[i]),
                                    (int*)&(this->nodes_per_element[i]),
                                    (int*)&(this->element_attributes[i]));
    }
    
    /*connectivity*/
    for(b = 0; b < this->num_elem_blk; b++){
      this->elmt_node_linkage[b] = (int*)malloc(this->nodes_per_element[b]*this->elements[b]*sizeof(int));
      error += im_ex_get_elem_conn(id,this->block_id[b],this->elmt_node_linkage[b]);
    }
    
    
    if(this->num_node_sets){
      this->node_set_id           = (int *) malloc(this->num_node_sets*sizeof(int));
      this->num_nodes_in_node_set = (int *) malloc(this->num_node_sets*sizeof(int));
      this->node_set_nodes        = (int **)malloc(this->num_node_sets*sizeof(int*));
      this->num_df_in_node_set    = (int *) malloc(this->num_node_sets*sizeof(int*));
      
      error += im_ex_get_node_set_ids(id, this->node_set_id);
      
      
      for(i = 0; i < this->num_node_sets; i ++){
        error += im_ex_get_node_set_param(id, this->node_set_id[i],
                                          (int*)&this->num_nodes_in_node_set[i],
                                          (int*)&this->num_df_in_node_set[i]);
        
        this->node_set_nodes[i] = NULL;
        
        if(this->num_nodes_in_node_set[i]) {
          this->node_set_nodes[i] = (int *)malloc(this->num_nodes_in_node_set[i]*sizeof(int));
          error += im_ex_get_node_set(id, this->node_set_id[i], this->node_set_nodes[i]);
        }
      }
    }
    
    /*side sets*/
    if(this->num_side_sets){
      this->side_set_id = (int*)malloc(this->num_side_sets*sizeof(int));
      this->num_elements_in_side_set = (int*)malloc(this->num_side_sets*sizeof(int));
      this->num_df_in_side_set = (int*)malloc(this->num_side_sets*sizeof(int));
      this->side_set_elements = (int**)malloc(this->num_side_sets*sizeof(int *));
      this->side_set_faces = (int **)malloc(this->num_side_sets*sizeof(int*));
      
      error += im_ex_get_side_set_ids(id, this->side_set_id);
      for(i = 0; i < this->num_side_sets; i ++){
        int ne = 0;
        error += im_ex_get_side_set_param(id, this->side_set_id[i],
                                          (int*)&this->num_elements_in_side_set[i],
                                          (int*)&this->num_df_in_side_set[i]);
        
        ne = this->num_elements_in_side_set[i];
        this->side_set_elements[i] = (int*)malloc(ne*sizeof(int));
        this->side_set_faces[i] = (int*)malloc(ne*sizeof(int));
        if(ne){
          error += im_ex_get_side_set(id, this->side_set_id[i],
                                      this->side_set_elements[i],
                                      this->side_set_faces[i]);
          
        }
      }
    }
    
    error += im_ex_inquire(id, IM_EX_INQ_QA, (int*)&this->num_qa_records, &fdum, cdum);
    
    if(this->num_qa_records)error +=  im_ex_get_qa(id,this->qaRecord);
    
    
    error += im_ex_inquire(id, IM_EX_INQ_INFO, (int*)&this->num_info_records, &fdum, cdum);
    if(this->num_info_records) {
      this->info_records = (char **)malloc(this->num_info_records*sizeof(char *));/*new std::string[num_info_records];*/
      for(i = 0; i < this->num_info_records; i ++){
        this->info_records[i] = (char *)malloc(MAX_STR_LENGTH+1);
      }
      error += im_ex_get_info(id, this->info_records);
    }
    
    
    /*nemesis data
     global info*/
    if ( im_ne_get_init_global(id, &this->num_nodes_global, &this->num_elems_global,
                               &this->num_elm_blks_global, &this->num_node_sets_global,
                               &this->num_side_sets_global) < 0 )
      ++error;
    
    
    
    if ( im_ne_get_init_info(id, &this->num_total_proc, &this->num_proc_in_file, this->type) < 0 )
      ++error;
    
    this->elem_blk_ids_global = (int*)malloc(this->num_elm_blks_global*sizeof(int));
    this->elem_blk_cnts_global = (int*)malloc(this->num_elm_blks_global*sizeof(int));
    
    if ( im_ne_get_eb_info_global(id,this->elem_blk_ids_global,this->elem_blk_cnts_global) < 0 )
      ++error;
    
    this->ns_ids_global = (int *)malloc(this->num_node_sets_global*sizeof(int));
    this->ns_cnts_global = (int *)malloc(this->num_node_sets_global*sizeof(int));
    this->ns_df_cnts_global = (int *)malloc(this->num_node_sets_global*sizeof(int));
    this->ss_ids_global = (int *)malloc(this->num_side_sets_global*sizeof(int));
    this->ss_cnts_global = (int *)malloc(this->num_side_sets_global*sizeof(int));
    this->ss_df_cnts_global = (int *)malloc(this->num_side_sets_global*sizeof(int));
    
    
    if ( this->num_node_sets_global > 0 ) {
      if ( im_ne_get_ns_param_global(id,this->ns_ids_global,this->ns_cnts_global,
                                     this->ns_df_cnts_global) < 0 )++error;
    }
    
    if ( this->num_side_sets_global > 0 ) {
      if ( im_ne_get_ss_param_global(id,this->ss_ids_global,this->ss_cnts_global,
                                     this->ss_df_cnts_global) < 0 )  ++error;
    }
    
    /*parallel info*/
    if ( im_ne_get_loadbal_param( id,
                                 &this->num_internal_nodes,
                                 &this->num_border_nodes,
                                 &this->num_external_nodes,
                                 &this->num_internal_elems,
                                 &this->num_border_elems,
                                 &this->num_node_comm_maps,
                                 &this->num_elem_comm_maps,
                                 0/*unused*/ ) < 0 )++error;
    
    this->internal_elements = (int *)malloc(this->num_internal_elems*sizeof(int));
    this->border_elements   = (int *)malloc(this->num_border_elems*sizeof(int));
    this->internal_nodes    = (int *)malloc(this->num_internal_nodes*sizeof(int));
    this->border_nodes      = (int *)malloc(this->num_border_nodes*sizeof(int));
    this->external_nodes    = (int *)malloc(this->num_external_nodes*sizeof(int));
    
    if ( im_ne_get_elem_map( id,
                            this->internal_elements,
                            this->border_elements,
                            0/* not used proc_id*/ ) < 0 )++error;
    
    if ( im_ne_get_node_map( id,
                            this->internal_nodes,
                            this->border_nodes,
                            this->external_nodes,
                            0/* not used proc_id*/ ) < 0 )++error;
    
    
    if(this->num_node_comm_maps > 0){
      
      this->node_cmap_node_cnts = (int*) malloc(this->num_node_comm_maps*sizeof(int));
      this->node_cmap_ids       = (int*) malloc(this->num_node_comm_maps*sizeof(int));
      this->comm_node_ids       = (int**)malloc(this->num_node_comm_maps*sizeof(int*));
      this->comm_node_proc_ids  = (int**)malloc(this->num_node_comm_maps*sizeof(int*));
      
      this->elem_cmap_elem_cnts = (int*) malloc(this->num_elem_comm_maps*sizeof(int));
      this->elem_cmap_ids       = (int*) malloc(this->num_elem_comm_maps*sizeof(int));
      this->comm_elem_ids       = (int**)malloc(this->num_elem_comm_maps*sizeof(int*));
      this->comm_side_ids       = (int**)malloc(this->num_elem_comm_maps*sizeof(int*));
      this->comm_elem_proc_ids  = (int**)malloc(this->num_elem_comm_maps*sizeof(int*));
      
      if ( im_ne_get_cmap_params( id,
                                 this->node_cmap_ids,
                                 (int*)this->node_cmap_node_cnts,
                                 this->elem_cmap_ids,
                                 (int*)this->elem_cmap_elem_cnts,
                                 0/*not used proc_id*/ ) < 0 )++error;
      
      for(j = 0; j < this->num_node_comm_maps; j++) {
        this->comm_node_ids[j]       = (int *)malloc(this->node_cmap_node_cnts[j]*sizeof(int));
        this->comm_node_proc_ids[j]  = (int *)malloc(this->node_cmap_node_cnts[j]*sizeof(int));
        if ( im_ne_get_node_cmap( id,
                                 this->node_cmap_ids[j],
                                 this->comm_node_ids[j],
                                 this->comm_node_proc_ids[j],
                                 0/*not used proc_id*/ ) < 0 )++error;
        
      }
      
      
      
      for(j = 0; j < this->num_elem_comm_maps; j++) {
        this->comm_elem_ids[j]       = (int *)malloc(this->elem_cmap_elem_cnts[j]*sizeof(int));
        this->comm_side_ids[j]       = (int *)malloc(this->elem_cmap_elem_cnts[j]*sizeof(int));
        this->comm_elem_proc_ids[j]  = (int *)malloc(this->elem_cmap_elem_cnts[j]*sizeof(int));
        if ( im_ne_get_elem_cmap( id,
                                 this->elem_cmap_ids[j],
                                 this->comm_elem_ids[j],
                                 this->comm_side_ids[j],
                                 this->comm_elem_proc_ids[j],
                                 0 /*not used proc_id*/ ) < 0 )++error;
        
        
      }/*loop over num_elem_co*/
    }
  }
  
  // compute element center coordinates
  this->computeElementCoordinates();
}

void PamgenMesh::computeElementCoordinates()
{
  this->element_coord = (double * )malloc(this->num_dim * this->num_elem * sizeof(double));
  memset(this->element_coord, 0, this->num_dim * this->num_elem * sizeof(double));
  
  // loop over elements
  int n_id = 0;
  
  int el_count = 0;
  for(int i = 0; i <  this->num_elem_blk; i++)
  {
    int els = this->elements[i];
    int nperel = this->nodes_per_element[i];
    // nodes for el
    int * connect = this->elmt_node_linkage[i]; // get node con
    
    for(int j = 0; j < els; j++)
    {
      
      // sum
      for(int k = 0; k < nperel; k++)
      {
        n_id = connect[j*nperel + k]-1;
        for(int l = 0; l < this->num_dim; l++)
          element_coord[el_count + l * this->num_elem] += this->coord[n_id + l *this->num_nodes];
      }
      
      // complete average
      for(int k = 0; k < this->num_dim; k++)
        element_coord[el_count + k*this->num_elem] /= nperel;
      
      el_count ++;
    }
  }
  
}

void PamgenMesh::createMesh(char * file_data, int dimension, const RCP<const Comm<int>> &comm)
{
  int rank = comm->getRank();
  int nproc = comm->getSize();
  long long cr_result = Create_Pamgen_Mesh(file_data, dimension, rank, nproc, INT_MAX);
  
  if (cr_result == ERROR_PARSING_DEFINITION){
    long long essz = getPamgenEchoStreamSize();
    char * echo_char_array = (char *)malloc(essz+1);
    printf("PARSE ERROR\n");
    echo_char_array[essz] = '\0';
    echo_char_array = getPamgenEchoStream(echo_char_array);
    if(echo_char_array)printf("%s",echo_char_array);
    if(cr_result == ERROR_CREATING_IMD)printf("ERROR Failure to create Inline_Mesh_Desc creation\n");
    if(echo_char_array)free(echo_char_array);
  }
  
  if(cr_result == ERROR_CREATING_MS){
    long long essz = getPamgenErrorStreamSize();
    char * error_char_array = (char *)malloc(essz+1);
    error_char_array[essz] = '\0';
    error_char_array = getPamgenErrorStream(error_char_array);
    if(error_char_array)printf("%s",error_char_array);
    printf("\nERROR Failure to create Mesh_Specification\n");
    if(error_char_array)free(error_char_array);
  }
  
  
  long long wssz = getPamgenWarningStreamSize();
  if(wssz){
    char * warning_char_array = (char *)malloc(wssz+1);
    warning_char_array[wssz] = '\0';
    warning_char_array = getPamgenWarningStream(warning_char_array);
    printf("WARNING Records\n");
    printf("%s",warning_char_array);
    free(warning_char_array);
  }
  
  
//  int issz = getPamgenInfoStreamSize();
//  if(issz){
//    char * info_char_array = (char *)malloc(issz+1);
//    info_char_array[issz] = '\0';
//    info_char_array = getPamgenInfoStream(info_char_array);
//    printf("INFO Records\n");
//    printf("%s",info_char_array);
//    free(info_char_array);
//  }
  
}

#endif

#endif
