// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

#include <stdlib.h> 
#include <stdio.h>
#include "Teuchos_CommandLineProcessor.hpp"

#include "create_inline_mesh.h"
#include "pamgen_im_exodusII.h"
#include "pamgen_im_ne_nemesisI.h"
#include <string.h>
#include "limits.h"

void write_mesh_to_stdout();
void read_mesh_to_memory();
void free_memory();

void write_to_exodus(int proc_id,int num_procs,char * out_file_name);


struct mesh_storage_struct{
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

  /*nemesis data*/
  /* global info*/

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

}mss;




/*****************************************************************************/
int main(int argc, char** argv)
/*****************************************************************************/
{
  int issz;
  int start_rank = 0;
  int end_rank = 0;
  char * out_file_name = NULL;
  FILE * infile = NULL;
  long size;
  char *file_char_array = NULL;

  // Teuchos-based optins
  Teuchos::CommandLineProcessor clp(false);

  int rank =0;
  clp.setOption("rank",  &rank, "rank of processor");
  int num_procs=1;
  clp.setOption("num_procs",  &num_procs, "rank of processor");
  int all = 0;
  clp.setOption("all",  &all, "generate all meshes for number of processors");
  int dimension = 3;
  clp.setOption("dimension",  &dimension, "dimension");
  std::string inFileName;
  clp.setOption("file",  &inFileName, "mesh description file");  

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  const char * file_name = inFileName.c_str();
  
  if(rank < 0 || rank >= num_procs){
    printf("rtpamgen: rank must be positive and between 0 and the number of processors %i\n",rank);
    exit(1);
  }

  if(num_procs < 0){
    printf("rtpamgen: number of processors must be positive %i\n", num_procs);
    exit(1);
  }


  /*deal with all switch*/
  end_rank = num_procs;
  if(!all){
    start_rank = rank;
    end_rank = start_rank+1;
  }

  infile = fopen(file_name,"r");
  if(!infile){
    printf("unable to open input file ");
    return 1;
  }
  fseek(infile, 0, SEEK_END);
  size = ftell(infile);
  fseek(infile, 0, SEEK_SET);
  file_char_array = (char *)malloc(size + 1);
  file_char_array[size] = '\0';
  fread(file_char_array, sizeof(char), size, infile);
  fclose(infile);

  /*create the out_file_name*/
  out_file_name = (char*)malloc(MAX_STR_LENGTH+1);

  for( rank = start_rank; rank != end_rank; rank ++){
    int cr_result;
    sprintf(out_file_name,"%s.exo.%i.%i",file_name,num_procs,rank);    
    
    cr_result = Create_Pamgen_Mesh(file_char_array,
                                   dimension,
				   rank,
				   num_procs,
                                   INT_MAX);
    
    
    if (cr_result == ERROR_PARSING_DEFINITION){
      int essz = getPamgenEchoStreamSize();
      char * echo_char_array = (char *)malloc(essz+1);
      printf("PARSE ERROR\n");
      echo_char_array[essz] = '\0';
      echo_char_array = getPamgenEchoStream(echo_char_array);
      if(echo_char_array)printf(echo_char_array);
      if(cr_result == ERROR_CREATING_IMD)printf("ERROR Failure to create Inline_Mesh_Desc creation\n");
      if(echo_char_array)free(echo_char_array);
      return 1;
    }
    
    if(cr_result == ERROR_CREATING_MS){
      int essz = getPamgenErrorStreamSize();
      char * error_char_array = (char *)malloc(essz+1);
      error_char_array[essz] = '\0';
      error_char_array = getPamgenErrorStream(error_char_array);
      if(error_char_array)printf(error_char_array);
      printf("\nERROR Failure to create Mesh_Specification\n");
      if(error_char_array)free(error_char_array);
      return 1;
    }
    
    
    {
      int wssz = getPamgenWarningStreamSize();
      if(wssz){
	char * warning_char_array = (char *)malloc(wssz+1);
	warning_char_array[wssz] = '\0';
	warning_char_array = getPamgenWarningStream(warning_char_array);
	printf("WARNING Records\n");
	printf(warning_char_array);
	free(warning_char_array);
      }
    }
    
    
    {
      issz = getPamgenInfoStreamSize();
      if(issz){
	char * info_char_array = (char *)malloc(issz+1);
	info_char_array[issz] = '\0';
	info_char_array = getPamgenInfoStream(info_char_array);
	printf("INFO Records\n");
	printf(info_char_array);
	free(info_char_array);
      }
    }
    
    read_mesh_to_memory();

    write_mesh_to_stdout();

    write_to_exodus(rank,num_procs,out_file_name);
    
    Delete_Pamgen_Mesh();
    free_memory();
  }/* end loop over all output ranks*/

  if(file_char_array)free(file_char_array);
  if(out_file_name)free(out_file_name);
}

/*****************************************************************************/
void free_memory()
/*****************************************************************************/
{
  int i;
  int j; 
  int b;
  for( i = 0; i < 100; i++){
    for( j=0; j<4; j++){
      free(mss.qaRecord[i][j]);
    }
  }
  
  
  
  free(mss.coord);
  
  if (mss.num_elem){
    free(mss.element_order_map); 
    
    if (mss.num_elem){
      free(mss.global_element_numbers);
    }
    
    if (mss.num_nodes){
      free(mss.global_node_numbers);
    }
    
    
    /*block info*/
    
    free(mss.block_id);
    free(mss.nodes_per_element);
    free(mss.element_attributes);
    free(mss.elements);
    
    
    for(i = 0; i < mss.num_elem_blk; i ++){
      free(mss.element_types[i]);
    }
    
    /*connectivity*/
    for(b = 0; b < mss.num_elem_blk; b++){
      free(mss.elmt_node_linkage[b]);
    }
    free(mss.element_types);
    free(mss.elmt_node_linkage);
    
    
    if(mss.num_node_sets){
      
      for(i = 0; i < mss.num_node_sets; i ++){      
	if(mss.num_nodes_in_node_set[i]) {
	  free(mss.node_set_nodes[i]);
	}
      }
      free(mss.node_set_id);
      free(mss.num_nodes_in_node_set);
      free(mss.node_set_nodes);
      free(mss.num_df_in_node_set);
      
    }
    
    /*side sets*/
    if(mss.num_side_sets){
      
      for(i = 0; i < mss.num_side_sets; i ++){
	
	free(mss.side_set_elements[i]);
	free(mss.side_set_faces[i]);
	
      }
      free(mss.side_set_id);
      free(mss.num_elements_in_side_set);
      free(mss.num_df_in_side_set);
      free(mss.side_set_elements);
      free(mss.side_set_faces);
    }
    
    
    if(mss.num_info_records) { 
      for(i = 0; i < mss.num_info_records; i ++){
	free(mss.info_records[i]);
      }
      free(mss.info_records);
    }
    
    
    /*nemesis data*/
    /* global info*/
    
    free(mss.elem_blk_ids_global);
    free(mss.elem_blk_cnts_global);
    
    free(mss.ns_ids_global);
    free(mss.ns_cnts_global);
    free(mss.ns_df_cnts_global);
    free(mss.ss_ids_global);
    free(mss.ss_cnts_global);
    free(mss.ss_df_cnts_global);
    
    free(mss.internal_elements);
    free(mss.border_elements);
    free(mss.internal_nodes);
    free(mss.border_nodes);
    free(mss.external_nodes);
    
    
    if(mss.num_node_comm_maps > 0){
      
      for(j = 0; j < mss.num_node_comm_maps; j++) {
	free(mss.comm_node_ids[j]);
	free(mss.comm_node_proc_ids[j]);
	
      }
      
      
      for(j = 0; j < mss.num_elem_comm_maps; j++) {
	free(mss.comm_elem_ids[j]);
	free(mss.comm_side_ids[j]);
	free(mss.comm_elem_proc_ids[j]);
	
	
	
      }/*loop over num_elem_co*/
    }
    
    
    free(mss.node_cmap_node_cnts);
    free(mss.node_cmap_ids);
    free(mss.comm_node_ids);
    free(mss.comm_node_proc_ids);
    
    free(mss.elem_cmap_elem_cnts);
    free(mss.elem_cmap_ids);
    free(mss.comm_elem_ids);
    free(mss.comm_side_ids);
    free(mss.comm_elem_proc_ids);
    
  }
}


/*****************************************************************************/
void read_mesh_to_memory()
/*****************************************************************************/
{
  int idum = 0;
  float fdum;
  int i;
  int j;
  int b;
  char * cdum = NULL;
  int error = 0;
  int id = 0;

  mss.bptr[0] = mss.buffer[0];      
  mss.bptr[1] = mss.buffer[1];      
  mss.bptr[2] = mss.buffer[2];      
  
  for(i = 0; i < 100; i++)
    for(j=0; j<4; j++) mss.qaRecord[i][j] = (char*)malloc(MAX_STR_LENGTH+1) ;
  
  
  error += im_ex_get_init (  id,
			     mss.title,
			     &mss.num_dim,
			     &(mss.num_nodes),
			     &mss.num_elem, 
			     &mss.num_elem_blk,
			     &mss.num_node_sets,
			     &mss.num_side_sets);
  
  
  error += im_ex_inquire(id, IM_EX_INQ_NS_NODE_LEN, (int*)&mss.num_node_set_nodes, 
			 &fdum, cdum);
  error += im_ex_inquire(id, IM_EX_INQ_NS_DF_LEN,   (int*)&mss.num_node_set_dfs, 
			 &fdum, cdum);
  error += im_ex_inquire(id, IM_EX_INQ_SS_ELEM_LEN, (int*)&mss.num_side_set_elements,
			 &fdum, cdum);
  error += im_ex_inquire(id, IM_EX_INQ_SS_NODE_LEN, (int*)&mss.num_side_set_nodes,
			 &fdum, cdum);
  error += im_ex_inquire(id, IM_EX_INQ_SS_DF_LEN,   (int*)&mss.num_side_set_dfs, 
			 &fdum, cdum);
    
  /* get version number*/
    
  error += im_ex_inquire(id, IM_EX_INQ_API_VERS, &idum, &fdum, cdum);
    
  mss.version_number = (double) fdum;
    
  mss.version = (int) mss.version_number;
    
  /* get genesis-II parameters*/
    
  error += im_ex_inquire(id, IM_EX_INQ_EB_PROP, (int*)&mss.num_block_properties, &fdum, cdum);
    
  error += im_ex_inquire(id, IM_EX_INQ_NS_PROP, (int*)&mss.num_node_set_properties, 
			 &fdum, cdum);
    
  error += im_ex_inquire(id, IM_EX_INQ_SS_PROP, (int*)&mss.num_side_set_properties, 
			 &fdum, cdum);
    
  mss.coord = (double *)malloc(mss.num_nodes*mss.num_dim*sizeof(double));
    
  error += im_ex_get_coord(id,mss.coord,mss.coord+mss.num_nodes,mss.coord+2*mss.num_nodes);

    
  error += im_ex_get_coord_names (id, mss.bptr);
    
  if (mss.num_elem){
    mss.element_order_map = (int *)malloc(mss.num_elem * sizeof(int));
    error += im_ex_get_map(id, mss.element_order_map);
      
    if (mss.num_elem){
      mss.global_element_numbers = (int *)malloc(mss.num_elem*sizeof(int));
      error += im_ex_get_elem_num_map(id, mss.global_element_numbers);
    }
      
    if (mss.num_nodes){
      mss.global_node_numbers = (int *)malloc(mss.num_nodes * sizeof(int));
      error += im_ex_get_node_num_map(id, mss.global_node_numbers);
    }
      
     
    /*block info*/

    mss.block_id           = (int *)malloc(mss.num_elem_blk*sizeof(int));
    mss.nodes_per_element  = (int *)malloc(mss.num_elem_blk*sizeof(int));
    mss.element_attributes = (int *)malloc(mss.num_elem_blk*sizeof(int));
    mss.elements           = (int *)malloc(mss.num_elem_blk*sizeof(int));
    mss.element_types      = (char **)malloc(mss.num_elem_blk*sizeof(char *));
    mss.elmt_node_linkage  = (int **)malloc(mss.num_elem_blk*sizeof(int*));

    error += im_ex_get_elem_blk_ids(id, mss.block_id);

    for(i = 0; i < mss.num_elem_blk; i ++){
      mss.element_types[i] = (char *)malloc((MAX_STR_LENGTH + 1)*sizeof(char));
      error += im_ex_get_elem_block(id, 
				    mss.block_id[i], 
				    mss.element_types[i],
				    (int*)&(mss.elements[i]),
				    (int*)&(mss.nodes_per_element[i]), 
				    (int*)&(mss.element_attributes[i]));
    }
    
    /*connectivity*/
    for(b = 0; b < mss.num_elem_blk; b++){
      mss.elmt_node_linkage[b] = (int*)malloc(mss.nodes_per_element[b]*mss.elements[b]*sizeof(int));
      error += im_ex_get_elem_conn(id,mss.block_id[b],mss.elmt_node_linkage[b]);
    }

      
    if(mss.num_node_sets){
      mss.node_set_id           = (int *) malloc(mss.num_node_sets*sizeof(int));
      mss.num_nodes_in_node_set = (int *) malloc(mss.num_node_sets*sizeof(int));
      mss.node_set_nodes        = (int **)malloc(mss.num_node_sets*sizeof(int*));
      mss.num_df_in_node_set    = (int *) malloc(mss.num_node_sets*sizeof(int*));
      
      error += im_ex_get_node_set_ids(id, mss.node_set_id);


      for(i = 0; i < mss.num_node_sets; i ++){
	error += im_ex_get_node_set_param(id, mss.node_set_id[i],
					  (int*)&mss.num_nodes_in_node_set[i], 
					  (int*)&mss.num_df_in_node_set[i]);

      	mss.node_set_nodes[i] = NULL;

	if(mss.num_nodes_in_node_set[i]) {
	  mss.node_set_nodes[i] = (int *)malloc(mss.num_nodes_in_node_set[i]*sizeof(int));
	  error += im_ex_get_node_set(id, mss.node_set_id[i], mss.node_set_nodes[i]);
	}
      }
    }

    /*side sets*/
    if(mss.num_side_sets){
      mss.side_set_id = (int*)malloc(mss.num_side_sets*sizeof(int));
      mss.num_elements_in_side_set = (int*)malloc(mss.num_side_sets*sizeof(int));
      mss.num_df_in_side_set = (int*)malloc(mss.num_side_sets*sizeof(int));
      mss.side_set_elements = (int**)malloc(mss.num_side_sets*sizeof(int *));
      mss.side_set_faces = (int **)malloc(mss.num_side_sets*sizeof(int*));

      error += im_ex_get_side_set_ids(id, mss.side_set_id);
      for(i = 0; i < mss.num_side_sets; i ++){
	int ne = 0;  
	error += im_ex_get_side_set_param(id, mss.side_set_id[i], 
					  (int*)&mss.num_elements_in_side_set[i],
					  (int*)&mss.num_df_in_side_set[i]);

	ne = mss.num_elements_in_side_set[i];
	mss.side_set_elements[i] = (int*)malloc(ne*sizeof(int));
	mss.side_set_faces[i] = (int*)malloc(ne*sizeof(int));
	if(ne){
	  error += im_ex_get_side_set(id, mss.side_set_id[i], 
				      mss.side_set_elements[i], 
				      mss.side_set_faces[i]);

	}
      }
    }
      
    error += im_ex_inquire(id, IM_EX_INQ_QA, (int*)&mss.num_qa_records, &fdum, cdum);

    if(mss.num_qa_records)error +=  im_ex_get_qa(id,mss.qaRecord);


    error += im_ex_inquire(id, IM_EX_INQ_INFO, (int*)&mss.num_info_records, &fdum, cdum);
    if(mss.num_info_records) { 
      mss.info_records = (char **)malloc(mss.num_info_records*sizeof(char *));/*new std::string[num_info_records];*/
      for(i = 0; i < mss.num_info_records; i ++){
	mss.info_records[i] = (char *)malloc(MAX_STR_LENGTH+1);
      }
      error += im_ex_get_info(id, mss.info_records);
    }


    /*nemesis data*/
    /* global info*/
    if ( im_ne_get_init_global(id, &mss.num_nodes_global, &mss.num_elems_global,
			       &mss.num_elm_blks_global, &mss.num_node_sets_global,
			       &mss.num_side_sets_global) < 0 )
      ++error;

  

    if ( im_ne_get_init_info(id, &mss.num_total_proc, &mss.num_proc_in_file, mss.type) < 0 )
      ++error;

    mss.elem_blk_ids_global = (int*)malloc(mss.num_elm_blks_global*sizeof(int));
    mss.elem_blk_cnts_global = (int*)malloc(mss.num_elm_blks_global*sizeof(int));

    if ( im_ne_get_eb_info_global(id,mss.elem_blk_ids_global,mss.elem_blk_cnts_global) < 0 )
      ++error;

    mss.ns_ids_global = (int *)malloc(mss.num_node_sets_global*sizeof(int));
    mss.ns_cnts_global = (int *)malloc(mss.num_node_sets_global*sizeof(int));
    mss.ns_df_cnts_global = (int *)malloc(mss.num_node_sets_global*sizeof(int));
    mss.ss_ids_global = (int *)malloc(mss.num_side_sets_global*sizeof(int));
    mss.ss_cnts_global = (int *)malloc(mss.num_side_sets_global*sizeof(int));
    mss.ss_df_cnts_global = (int *)malloc(mss.num_side_sets_global*sizeof(int));

    
    if ( mss.num_node_sets_global > 0 ) {
      if ( im_ne_get_ns_param_global(id,mss.ns_ids_global,mss.ns_cnts_global,
				     mss.ns_df_cnts_global) < 0 )++error;
    }
    
    if ( mss.num_side_sets_global > 0 ) {
      if ( im_ne_get_ss_param_global(id,mss.ss_ids_global,mss.ss_cnts_global,
				     mss.ss_df_cnts_global) < 0 )  ++error;      
    }
    
    /*parallel info*/
    if ( im_ne_get_loadbal_param( id, 
				  &mss.num_internal_nodes,
				  &mss.num_border_nodes, 
				  &mss.num_external_nodes,
				  &mss.num_internal_elems, 
				  &mss.num_border_elems,
				  &mss.num_node_comm_maps,
				  &mss.num_elem_comm_maps,
				  0/*unused*/ ) < 0 )++error;
    
    mss.internal_elements = (int *)malloc(mss.num_internal_elems*sizeof(int));
    mss.border_elements   = (int *)malloc(mss.num_border_elems*sizeof(int));
    mss.internal_nodes    = (int *)malloc(mss.num_internal_nodes*sizeof(int));
    mss.border_nodes      = (int *)malloc(mss.num_border_nodes*sizeof(int));
    mss.external_nodes    = (int *)malloc(mss.num_external_nodes*sizeof(int));
   
    if ( im_ne_get_elem_map( id, 
			     mss.internal_elements, 
			     mss.border_elements, 
			     0/* not used proc_id*/ ) < 0 )++error;
    
    if ( im_ne_get_node_map( id, 
			     mss.internal_nodes, 
			     mss.border_nodes,
			     mss.external_nodes, 
			     0/* not used proc_id*/ ) < 0 )++error;


    if(mss.num_node_comm_maps > 0){

      mss.node_cmap_node_cnts = (int*) malloc(mss.num_node_comm_maps*sizeof(int));
      mss.node_cmap_ids       = (int*) malloc(mss.num_node_comm_maps*sizeof(int));
      mss.comm_node_ids       = (int**)malloc(mss.num_node_comm_maps*sizeof(int*));
      mss.comm_node_proc_ids  = (int**)malloc(mss.num_node_comm_maps*sizeof(int*));

      mss.elem_cmap_elem_cnts = (int*) malloc(mss.num_elem_comm_maps*sizeof(int));
      mss.elem_cmap_ids       = (int*) malloc(mss.num_elem_comm_maps*sizeof(int));
      mss.comm_elem_ids       = (int**)malloc(mss.num_elem_comm_maps*sizeof(int*));
      mss.comm_side_ids       = (int**)malloc(mss.num_elem_comm_maps*sizeof(int*));
      mss.comm_elem_proc_ids  = (int**)malloc(mss.num_elem_comm_maps*sizeof(int*));

      if ( im_ne_get_cmap_params( id, 
				  mss.node_cmap_ids,
				  (int*)mss.node_cmap_node_cnts, 
				  mss.elem_cmap_ids,
				  (int*)mss.elem_cmap_elem_cnts, 
				  0/*not used proc_id*/ ) < 0 )++error;
      
      for(j = 0; j < mss.num_node_comm_maps; j++) {
	mss.comm_node_ids[j]       = (int *)malloc(mss.node_cmap_node_cnts[j]*sizeof(int));
	mss.comm_node_proc_ids[j]  = (int *)malloc(mss.node_cmap_node_cnts[j]*sizeof(int));
	if ( im_ne_get_node_cmap( id, 
				  mss.node_cmap_ids[j], 
				  mss.comm_node_ids[j], 
				  mss.comm_node_proc_ids[j],
				  0/*not used proc_id*/ ) < 0 )++error;
	
      }
      
      
      
      for(j = 0; j < mss.num_elem_comm_maps; j++) {
	mss.comm_elem_ids[j]       = (int *)malloc(mss.elem_cmap_elem_cnts[j]*sizeof(int));
	mss.comm_side_ids[j]       = (int *)malloc(mss.elem_cmap_elem_cnts[j]*sizeof(int));
	mss.comm_elem_proc_ids[j]  = (int *)malloc(mss.elem_cmap_elem_cnts[j]*sizeof(int));
	if ( im_ne_get_elem_cmap( id, 
				  mss.elem_cmap_ids[j], 
				  mss.comm_elem_ids[j], 
				  mss.comm_side_ids[j],
				  mss.comm_elem_proc_ids[j],
				  0 /*not used proc_id*/ ) < 0 )++error;
	

      }/*loop over num_elem_co*/
    }
  }
}



/*****************************************************************************/
void write_mesh_to_stdout()
/*****************************************************************************/
{
int i;
int j;
int b;
int ict;
int nct;

    printf("\nExodus header info:\nTitle: %s\nDimension %i \nNumber of Nodes %i \nNumber of Elements %i \nNumber of Element Blocks %i \nNumber of Node Sets %i \nNumber of Side Sets %i \n\n",
	   mss.title,
	   mss.num_dim,
	   mss.num_nodes, 
	   mss.num_elem, 
	   mss.num_elem_blk, 
	   mss.num_node_sets, 
	   mss.num_side_sets);
  
  
  printf("num node set nodes %i\n",mss.num_node_set_nodes);
  printf("num node set dfs %i\n",mss.num_node_set_dfs);
 
  printf("num side set elements %i\n",mss.num_side_set_elements);
  printf("num side set nodes %i\n",mss.num_side_set_nodes);
  printf("num side set dfs %i\n",mss.num_side_set_dfs);
    
  /* get version number*/
  printf("num block properties %i\n",mss.num_block_properties);
  printf("num node set properties %i\n",mss.num_node_set_properties);
  printf("num side set properties %i\n",mss.num_side_set_properties);
    

  printf("A taste of coords\n");
  for(i = 0; i < mss.num_nodes && i < 10; i ++){
    if(mss.num_dim == 3) printf("X %f Y %f Z %f \n",mss.coord[i],mss.coord[i+mss.num_nodes],mss.coord[i+2*mss.num_nodes]);
    if(mss.num_dim == 2) printf("X %f Y %f \n",mss.coord[i],mss.coord[i+mss.num_nodes]);
  }
 
  for(i = 0; i < mss.num_dim; i ++)printf("coord name %i %s\n",i,mss.bptr[i]);
    
  if (mss.num_elem){
 
    printf ("A tast of map\n");
    for(i = 0; i < mss.num_elem && i <  10; i ++)printf("map i=%i, val=%i\n",i,mss.element_order_map[i]);
    if (mss.num_elem){
      printf ("A tast of global elem numbers\n");
      for(i = 0; i < mss.num_elem && i <  10; i ++)printf("global el num  i=%i, val=%i\n",i,mss.global_element_numbers[i]);
    }
      
    if (mss.num_nodes){
      printf ("A tast of global elem numbers\n");
      for(i = 0; i < mss.num_elem && i <  10; i ++)printf("global node num  i=%i, val=%i\n",i,mss.global_node_numbers[i]);
    }
      
     
    /*block info*/

    for(i = 0; i < mss.num_elem_blk; i ++){
      printf("block i = %i has id %i \n",i,mss.block_id[i]);
      printf("block i = %i\n",i);
      printf("block id %i\n",mss.block_id[i]);
      printf("element_type %s\n",mss.element_types[i]);
      printf("num elements %i\n",mss.elements[i]);
      printf("nodes per element %i\n",mss.nodes_per_element[i]);
      printf("element attributes %i\n",mss.element_attributes[i]);
    }
    
    /*connectivity*/
    for(b = 0; b < mss.num_elem_blk; b++){
      for(ict = 0; ict < mss.elements[b] && ict < 10;ict++){
	printf("block %i element %i connectivty ",mss.block_id[b],ict);
	for(nct = 0; nct < mss.nodes_per_element[b]; nct++){
	  printf("%i ",mss.elmt_node_linkage[b][nct+ict*mss.nodes_per_element[b]]);
	}
	printf("\n");
      }
    }

      
    if(mss.num_node_sets){
      for(i = 0; i < mss.num_node_sets; i ++){
	printf("Nodeset i = %i id = %i has %i nodes\n",i,mss.node_set_id[i],mss.num_nodes_in_node_set[i]);
	for(j = 0; j < mss.num_nodes_in_node_set[i] && j < 10; j ++){
	  printf("nodeset node i=%i = %i\n",j,mss.node_set_nodes[i][j]);
	}
      }
    }

    /*side sets*/
    if(mss.num_side_sets){
      for(i = 0; i < mss.num_side_sets; i ++){
        int ne = 0;
	printf("Side set index %i id %i has %i elements\n",i,mss.side_set_id[i],mss.num_elements_in_side_set[i]);
	ne = mss.num_elements_in_side_set[i];
	if(ne){
	  for(j = 0; j < ne && j < 10; j ++){
	    printf("element %i and face %i\n",mss.side_set_elements[i][j],mss.side_set_faces[i][j]);
	  }
	}
      }
    }
      
    printf("num qa records %i\n",mss.num_qa_records);

    for(i = 0; i < mss.num_qa_records; i ++){
      printf("\nQA Record %i\n %s\n%s\n",
	     i,
	     mss.qaRecord[i][0],
	     mss.qaRecord[i][1]);
    }

    printf("Num Info Records %i\n",mss.num_info_records);
    if(mss.num_info_records) { 
      printf("Info Records\n");
      for(i = 0; i < mss.num_info_records; i ++){
	printf(mss.info_records[i]);
      }
    }


    /*nemesis data*/
    /* global info*/

    printf("Nemesis data\n");
    printf("Num nodes global %i\n",mss.num_nodes_global);
    printf("Num elems global %i\n",mss.num_elems_global);
    printf("Num elm_blks global %i\n",mss.num_elm_blks_global);
    printf("Num node sets global %i\n",mss.num_node_sets_global);
    printf("Num side sets global %i\n",mss.num_side_sets_global);
    printf("Num total proc %i\n",mss.num_total_proc);
    printf("Num proc in file %i\n",mss.num_proc_in_file);

    for(i = 0; i < mss.num_elm_blks_global; i ++){
      printf("element block index %i has id %i and %i elements\n",i,mss.elem_blk_ids_global[i],mss.elem_blk_cnts_global[i]);
    }
    
    if ( mss.num_node_sets_global > 0 ) {
      for(i = 0; i < mss.num_node_sets_global;i ++){
	printf("global ns info for ns index %i id %i num_nodes = %i num_ns_df = %i\n",i,mss.ns_ids_global[i],mss.ns_cnts_global[i],mss.ns_df_cnts_global[i]);
      }
    }
    
    if ( mss.num_side_sets_global > 0 ) {
      for(i = 0; i < mss.num_side_sets_global;i ++){
	printf("global ss info for ss index %i id %i num_elements = %i num_ss_df = %i\n",i,mss.ss_ids_global[i],mss.ss_cnts_global[i],mss.ss_df_cnts_global[i]);
      }
      
    }
    
    /*parallel info*/

    printf("Loadbal params:\nnum_internal_nodes %i\nnum_border_nodes%i\nnum_external_nodes%i\nnum_internal_elems%i\nnum_border_elems%i\nnum_node_comm_maps%i\nnum_elem_comm_maps%i\n",
	   mss.num_internal_nodes,
	   mss.num_border_nodes, 
	   mss.num_external_nodes,
	   mss.num_internal_elems, 
	   mss.num_border_elems,
	   mss.num_node_comm_maps,
	   mss.num_elem_comm_maps);

    for(i = 0; i < mss.num_internal_nodes && i < 10; i ++)printf("internal node i=%i = %i\n",i,mss.internal_nodes[i]);
    for(i = 0; i < mss.num_border_nodes && i < 10;   i ++)printf("border node i=%i = %i\n",i,  mss.border_nodes[i]);
    for(i = 0; i < mss.num_external_nodes && i < 10; i ++)printf("external node i=%i = %i\n",i,mss.external_nodes[i]);
    for(i = 0; i < mss.num_internal_elems && i < 10; i ++)printf("internal elem i=%i = %i\n",i,mss.internal_elements[i]);
    for(i = 0; i < mss.num_border_elems && i < 10;   i ++)printf("border elem i=%i = %i\n",i,  mss.border_elements[i]);

    if(mss.num_node_comm_maps > 0){

      for(i =0; i < mss.num_node_comm_maps; i ++){
	printf("node_cmap_id i = %i node_cmap_id = %i node_cmap_node_cnts = %i\n",i,mss.node_cmap_ids[i],mss.node_cmap_node_cnts[i]);
      }
      for(i =0; i < mss.num_elem_comm_maps; i ++){
	printf("elem_cmap_id i = %i elem_cmap_id = %i elem_cmap_elem_cnts = %i\n",i,mss.elem_cmap_ids[i],mss.elem_cmap_elem_cnts[i]);
      }
      
      for(j = 0; j < mss.num_node_comm_maps; j++) {
	for(i = 0; i < mss.node_cmap_node_cnts[j] && i < 10;i ++){
	  printf("node_cmap_id i=%i = %i comm_node_ids = %i comm_node_proc_ids = %i\n",
		 i,
		 mss.node_cmap_ids[j],
		 mss.comm_node_ids[j][i],
		 mss.comm_node_proc_ids[j][i]);
	}
      }
      
      for(j = 0; j < mss.num_elem_comm_maps; j++) {
	for(i = 0; i < mss.elem_cmap_elem_cnts[j] && i < 10;i ++){
	  printf("elem_cmap_id i=%i = %i comm_elem_ids = %i comm_side_ids = %i comm_elem_proc_ids = %i\n",
		 i,
		 mss.elem_cmap_ids[j],
		 mss.comm_elem_ids[j][i],
		 mss.comm_side_ids[j][i],
		 mss.comm_elem_proc_ids[j][i]);
	}

      }/*loop over num_elem_co*/
    }
  }
}

// Because this clashes with pamgen's definition
#undef MAX_ERR_LENGTH

#include "exodusII.h"
#include "ne_nemesisI.h"

#define PERROR printf("error %i ec %i\n",error,ec);ec++;
#undef PERROR 
#define PERROR ;
/*****************************************************************************/
void write_to_exodus(int rank, int num_procs, char * out_file_name)
/*****************************************************************************/
{

  int exo_access = EX_CLOBBER;
  int cpu_word_size = sizeof(double);
  int io_word_size = sizeof(double);
  int out_id;
  int i; 
  int b;
  int error = 0;
  
  ex_init_params exinit;

  out_id = ex_create(out_file_name, exo_access, &cpu_word_size,
		 &io_word_size);

  if (out_id < 0){
    printf("error opening file");
  }

  strncpy( exinit.title, mss.title, MAX_LINE_LENGTH-1 );
  exinit.title[MAX_LINE_LENGTH-1] = 0;
  exinit.num_dim       = mss.num_dim;
  exinit.num_nodes     = mss.num_nodes;
  exinit.num_edge      = 0;
  exinit.num_edge_blk  = 0;
  exinit.num_face      = 0;
  exinit.num_face_blk  = 0;
  exinit.num_elem      = mss.num_elem;
  exinit.num_elem_blk  = mss.num_elem_blk;
  exinit.num_node_sets = mss.num_node_sets;
  exinit.num_edge_sets = 0;
  exinit.num_face_sets = 0;
  exinit.num_side_sets = mss.num_side_sets;
  exinit.num_elem_sets = 0;
  exinit.num_node_maps = 0;
  exinit.num_edge_maps = 0;
  exinit.num_face_maps = 0;
  exinit.num_elem_maps = 0;


  PERROR;


  if ( ex_put_init_ext(out_id, &exinit) < 0 )
    ++error;
  PERROR;
 
/* now write parallel global information*/

  if ( ne_put_init_global( out_id, 
			   mss.num_nodes_global, 
			   mss.num_elems_global,
                           mss.num_elm_blks_global, 
			   mss.num_node_sets_global,
                           mss.num_side_sets_global ) < 0 )
    ++error;
  PERROR;
  
  if ( ne_put_init_info( out_id, mss.num_total_proc, mss.num_proc_in_file, 
                         mss.type ) < 0 )
    ++error;
  PERROR;
  
  if ( ne_put_eb_info_global(out_id,mss.elem_blk_ids_global,mss.elem_blk_cnts_global) < 0 )
    ++error;
  PERROR;
  
  if ( mss.num_node_sets_global > 0 ) {
    if ( ne_put_ns_param_global( out_id,
				 mss.ns_ids_global,
				 mss.ns_cnts_global,
                                 mss.ns_df_cnts_global ) < 0 )
      ++error;
  }
  PERROR;
  
  if ( mss.num_side_sets_global > 0 ) {
    if ( ne_put_ss_param_global( out_id,
				 mss.ss_ids_global,
				 mss.ss_cnts_global,
                                 mss.ss_df_cnts_global ) < 0 )
      ++error;
  }
  PERROR;

  /*writingparallel info*/
  if ( ne_put_loadbal_param( out_id, 
                             mss.num_internal_nodes,
                             mss.num_border_nodes, 
                             mss.num_external_nodes,
                             mss.num_internal_elems, 
                             mss.num_border_elems,
                             mss.num_node_comm_maps,
                             mss.num_elem_comm_maps,
                             rank ) < 0 )
    ++error;
  PERROR;

  if ( ne_put_cmap_params( out_id, 
                           mss.node_cmap_ids,
                           (int*)mss.node_cmap_node_cnts, 
                           mss.elem_cmap_ids,
                           (int*)mss.elem_cmap_elem_cnts, 
                           rank ) < 0 )
    ++error;
  PERROR;

  if ( ne_put_elem_map( out_id, 
                        mss.internal_elements, 
                        mss.border_elements, 
                        rank ) < 0 )
    ++error;
  PERROR;


  if ( ne_put_node_map( out_id, 
                        mss.internal_nodes, 
                        mss.border_nodes,
                        mss.external_nodes, 
                        rank ) < 0 )
    ++error;
  PERROR;

  for (i = 0; i < mss.num_node_comm_maps; i++) {
    if ( ne_put_node_cmap( out_id, 
                           mss.node_cmap_ids[i], 
                           mss.comm_node_ids[i], 
                           mss.comm_node_proc_ids[i],
                           rank ) < 0 )
      ++error;
  }
  PERROR;


  for (i = 0; i < mss.num_elem_comm_maps; i++) {
    if ( ne_put_elem_cmap( out_id, 
                           mss.elem_cmap_ids[i], 
                           mss.comm_elem_ids[i], 
                           mss.comm_side_ids[i],
                           mss.comm_elem_proc_ids[i],
                           rank ) < 0 )
      ++error;
  }
  
  PERROR;

  /*coords*/
  error += ex_put_coord(out_id, mss.coord, (mss.coord)+mss.num_nodes, (mss.coord)+2*mss.num_nodes);
  PERROR;
  error += ex_put_coord_names(out_id, mss.bptr);
  PERROR;
  /*map*/
  error += ex_put_map(out_id, mss.element_order_map);
  PERROR;
  error += ex_put_elem_num_map(out_id, mss.global_element_numbers);
  PERROR;
  error += ex_put_node_num_map(out_id, mss.global_node_numbers);
  PERROR;



  /*block info*/
  for(b = 0; b < mss.num_elem_blk; b++)
  {
    int gpe = 0;
    int fpe = 0;
    error += ex_put_block( out_id,
                           EX_ELEM_BLOCK,
                           mss.block_id[b],
                           mss.element_types[b],
                           mss.elements[b],
                           mss.nodes_per_element[b],
                           gpe, fpe,
                           mss.element_attributes[b] );  /* num attr*/
    PERROR;
  }

/* write element connectivity information*/
  
  for (b = 0; b < mss.num_elem_blk; b++) {
    if ( mss.elements[b] > 0 ){
      error += ex_put_elem_conn(out_id,mss.block_id[b],mss.elmt_node_linkage[b]);
      PERROR;
    }
  }


/* write in nodal boundary sets for the body.*/

  for(i = 0; i < mss.num_node_sets; i++) {  
    error += ex_put_node_set_param(out_id, mss.node_set_id[i],
                                  mss.num_nodes_in_node_set[i], 
                                  mss.num_df_in_node_set[i]);
    PERROR;
    if(mss.num_nodes_in_node_set[i])
      error += ex_put_node_set(out_id, mss.node_set_id[i], mss.node_set_nodes[i]);
    PERROR;

  }

  for(i = 0; i < mss.num_side_sets; i++) {
    error += ex_put_side_set_param(out_id, mss.side_set_id[i], 
                                  mss.num_elements_in_side_set[i],
                                  mss.num_df_in_side_set[i]);

    PERROR;
    if(mss.num_elements_in_side_set[i]) 
      error += ex_put_side_set(out_id, mss.side_set_id[i], 
                              mss.side_set_elements[i], 
                              mss.side_set_faces[i]);
    PERROR;
  }

    error += ex_put_qa(out_id, mss.num_qa_records, mss.qaRecord);
    PERROR;
 
  ex_close(out_id);


}
