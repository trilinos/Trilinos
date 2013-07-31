#include "getopts.h"
#include <stdlib.h>
#include <stdio.h>
#include "create_inline_mesh.h"
#include "../mesh_spec_lt/im_exodusII_l.h"
#include "../mesh_spec_lt/im_ne_nemesisI_l.h"
#include <string.h>
#include <limits.h>

void write_mesh_to_stdout();
void read_mesh_to_memory();
void free_memory();

/*#define HAVE_EXODUS*/
#ifdef HAVE_EXODUS
void write_to_exodus(int proc_id,int num_procs,char * out_file_name);
#endif /*HAVE_EXODUS*/


struct mesh_storage_struct{
  long long num_dim;
  long long num_nodes;
  long long num_elem;
  long long num_elem_blk;
  long long num_node_sets;
  long long num_side_sets;
  long long num_node_set_nodes;
  long long num_node_set_dfs;
  long long num_side_set_elements;
  long long num_side_set_nodes;
  long long num_side_set_dfs;
  long long num_block_properties;
  long long num_node_set_properties;
  long long num_side_set_properties;

  char title[MAX_STR_LENGTH];

  long long version;
  double version_number;
  double * coord;

  char buffer[3][MAX_STR_LENGTH + 1];
  char *bptr[3];

  long long * element_order_map ;

  long long * global_element_numbers ;
  long long * global_node_numbers ;

  /*block info*/
  long long * block_id ;
  char ** element_types ;
  long long *   elements ;
  long long *   nodes_per_element ;
  long long *   element_attributes ;
  long long ** elmt_node_linkage ;

  /*side sets*/
  long long * side_set_id ;
  long long * num_elements_in_side_set ;
  long long * num_df_in_side_set ;
  long long **side_set_elements ;
  long long **side_set_faces ;

  /*node sets*/
  long long * node_set_id ;
  long long * num_nodes_in_node_set ;
  long long * num_df_in_node_set ;
  long long **node_set_nodes  ;

  /*qa*/
  long long num_qa_records;
  long long num_info_records;
  char* qaRecord[100][4];
  char** info_records ;

  /*nemesis data*/
  long long num_nodes_global;
  long long num_elems_global;
  long long num_elm_blks_global;
  long long num_node_sets_global;
  long long num_side_sets_global;
  long long num_total_proc;
  long long num_proc_in_file;
  char type[2];

  /*nemesis data
    global info*/

  long long * elem_blk_ids_global ;
  long long * elem_blk_cnts_global  ;

  long long * ns_ids_global ;
  long long * ns_cnts_global ;
  long long * ns_df_cnts_global ;
  long long * ss_ids_global ;
  long long * ss_cnts_global ;
  long long * ss_df_cnts_global ;

  /*parallel info*/
  long long num_internal_nodes;
  long long num_border_nodes;
  long long num_external_nodes;
  long long num_internal_elems;
  long long num_border_elems;
  long long num_node_comm_maps;
  long long num_elem_comm_maps;

  long long * internal_elements ;
  long long * border_elements ;
  long long * internal_nodes ;
  long long * border_nodes ;
  long long * external_nodes ;

  long long * node_cmap_node_cnts ;
  long long * node_cmap_ids       ;
  long long * elem_cmap_elem_cnts ;
  long long * elem_cmap_ids       ;

  long long ** comm_node_ids       ;
  long long ** comm_node_proc_ids  ;
  long long ** comm_elem_ids       ;
  long long ** comm_side_ids       ;
  long long ** comm_elem_proc_ids  ;

}mss;




/*****************************************************************************/
int main(int argc, char** argv)
/*****************************************************************************/
{

  struct options opts[] =
  {
    { 1, "rank",   "rank of processor (0) ", "r", 1 },
    { 2, "num_procs", "number of processors (1)", "n", 1 },
    { 3, "all",      "generate all meshes for number of processors (false)",      "a", 0 },
    { 4, "dimension",    "dimension (3D)",       "d", 1 },
    { 5, "file",    "mesh description file (no default)",       "f", 1 },
    { 0, NULL,      NULL,                      NULL, 0 }
  };

  char *args;
  long long c;
  long long all = FALSE;
  long long rank = 0;
  long long issz;
  long long num_procs = 1;
  long long start_rank = 0;
  long long end_rank = 0;
  long long dimension = 3;
  char * file_name = NULL;
  char * out_file_name = NULL;
  FILE * infile = NULL;
  long size;
  char *file_char_array = NULL;

  if(argc == 1){
    getopts_usage(argv[0],opts);
    return 1;
  }

  while ((c = getopts(argc, argv, opts, &args)) != 0) { switch(c)
        {
	  /* Special Case: Recognize options that we didn't set above. */
	case -2:
	  printf("Unknown Option: \n");/* <<  args << std::endl;*/
	  getopts_usage(argv[0],opts);
	  return 1;
	  break;
	  /* Special Case: getopts() can't allocate memory.
	     Should probably exit() here */
	case -1:
	  printf("Unabled to allocate memory from getopts().\n");
	  break;
	case 1:
	  rank = atoi(args);
	  break;
        case 2:
	  num_procs = atoi(args);
	  break;
        case 3:
	  all = TRUE;
	  break;
        case 4:
	  dimension = atoi(args);
	  break;
        case 5:
	  {
	    int sl = strlen(args);
	    file_name = (char *)malloc(sl+1);
	    file_name[sl] = '\0';
	    strcpy(file_name,args);
	  }
	  break;
	default:
	  break;
        }
/*   This free() is required since getopts() automagically allocates space */
/*     for "args" everytime it's called.  */
      free(args);
    }

  if(rank < 0 || rank >= num_procs){
    getopts_usage(argv[0],opts);
    printf("rank must be positive and between 0 and  the number of processors %lli\n",rank);
  }

  if(num_procs < 0){
    getopts_usage(argv[0],opts);
    printf("number of processors must be positive %lli\n", num_procs);
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
    sprintf(out_file_name,"%s.exo.%lli.%lli",file_name,num_procs,rank);

    cr_result = Create_Pamgen_Mesh(file_char_array,
				   dimension,
				   rank,
				   num_procs,
				   9223372036854775807LL);


    if (cr_result == ERROR_PARSING_DEFINITION){
      int essz = getPamgenEchoStreamSize();
      char * echo_char_array = (char *)malloc(essz+1);
      printf("PARSE ERROR\n");
      echo_char_array[essz] = '\0';
      echo_char_array = getPamgenEchoStream(echo_char_array);
      if(echo_char_array)printf("%s",echo_char_array);
      if(cr_result == ERROR_CREATING_IMD)printf("ERROR Failure to create Inline_Mesh_Desc creation\n");
      if(echo_char_array)free(echo_char_array);
      return 1;
    }

    if(cr_result == ERROR_CREATING_MS){
      int essz = getPamgenErrorStreamSize();
      char * error_char_array = (char *)malloc(essz+1);
      error_char_array[essz] = '\0';
      error_char_array = getPamgenErrorStream(error_char_array);
      if(error_char_array)printf("%s",error_char_array);
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
	printf("%s",warning_char_array);
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
	printf("%s",info_char_array);
	free(info_char_array);
      }
    }

    read_mesh_to_memory();

    write_mesh_to_stdout();

#ifdef HAVE_EXODUS
    write_to_exodus(rank,num_procs,out_file_name);
#endif/* HAVE_EXODUS*/

    Delete_Pamgen_Mesh();
    free_memory();
  }/* end loop over all output ranks*/

  if(file_char_array)free(file_char_array);
  if(out_file_name)free(out_file_name);
  if(file_name)free(file_name);

  return 0;
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


    /*nemesis data
      global info */

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
  long long idum = 0;
  float fdum;
  long long i;
  long long j;
  long long b;
  char * cdum = NULL;
  int error = 0;
  int id = 0;

  mss.bptr[0] = mss.buffer[0];
  mss.bptr[1] = mss.buffer[1];
  mss.bptr[2] = mss.buffer[2];

  for(i = 0; i < 100; i++)
    for(j=0; j<4; j++) mss.qaRecord[i][j] = (char*)malloc(MAX_STR_LENGTH+1) ;


  error += im_ex_get_init_l (  id,
			     mss.title,
			     &mss.num_dim,
			     &(mss.num_nodes),
			     &mss.num_elem,
			     &mss.num_elem_blk,
			     &mss.num_node_sets,
			     &mss.num_side_sets);


  error += im_ex_inquire_l(id, IM_EX_INQ_NS_NODE_LEN, (long long*)&mss.num_node_set_nodes,
			 &fdum, cdum);
  error += im_ex_inquire_l(id, IM_EX_INQ_NS_DF_LEN,   (long long*)&mss.num_node_set_dfs,
			 &fdum, cdum);
  error += im_ex_inquire_l(id, IM_EX_INQ_SS_ELEM_LEN, (long long*)&mss.num_side_set_elements,
			 &fdum, cdum);
  error += im_ex_inquire_l(id, IM_EX_INQ_SS_NODE_LEN, (long long*)&mss.num_side_set_nodes,
			 &fdum, cdum);
  error += im_ex_inquire_l(id, IM_EX_INQ_SS_DF_LEN,   (long long*)&mss.num_side_set_dfs,
			 &fdum, cdum);

  /* get version number */

  error += im_ex_inquire_l(id, IM_EX_INQ_API_VERS, &idum, &fdum, cdum);

  mss.version_number = (double) fdum;

  mss.version = (long long) mss.version_number;

  /* get genesis-II parameters */

  error += im_ex_inquire_l(id, IM_EX_INQ_EB_PROP, (long long*)&mss.num_block_properties, &fdum, cdum);

  error += im_ex_inquire_l(id, IM_EX_INQ_NS_PROP, (long long*)&mss.num_node_set_properties,
			 &fdum, cdum);

  error += im_ex_inquire_l(id, IM_EX_INQ_SS_PROP, (long long*)&mss.num_side_set_properties,
			 &fdum, cdum);

  mss.coord = (double *)malloc(mss.num_nodes*mss.num_dim*sizeof(double));

  error += im_ex_get_coord_l(id,mss.coord,mss.coord+mss.num_nodes,mss.coord+2*mss.num_nodes);


  error += im_ex_get_coord_names_l (id, mss.bptr);

  if (mss.num_elem){
    mss.element_order_map = (long long *)malloc(mss.num_elem * sizeof(long long));
    error += im_ex_get_map_l(id, mss.element_order_map);

    if (mss.num_elem){
      mss.global_element_numbers = (long long *)malloc(mss.num_elem*sizeof(long long));
      error += im_ex_get_elem_num_map_l(id, mss.global_element_numbers);
    }

    if (mss.num_nodes){
      mss.global_node_numbers = (long long *)malloc(mss.num_nodes * sizeof(long long));
      error += im_ex_get_node_num_map_l(id, mss.global_node_numbers);
    }


    /*block info*/

    mss.block_id           = (long long *)malloc(mss.num_elem_blk*sizeof(long long));
    mss.nodes_per_element  = (long long *)malloc(mss.num_elem_blk*sizeof(long long));
    mss.element_attributes = (long long *)malloc(mss.num_elem_blk*sizeof(long long));
    mss.elements           = (long long *)malloc(mss.num_elem_blk*sizeof(long long));
    mss.element_types      = (char **)malloc(mss.num_elem_blk*sizeof(char *));
    mss.elmt_node_linkage  = (long long **)malloc(mss.num_elem_blk*sizeof(long long*));

    error += im_ex_get_elem_blk_ids_l(id, mss.block_id);

    for(i = 0; i < mss.num_elem_blk; i ++){
      mss.element_types[i] = (char *)malloc((MAX_STR_LENGTH + 1)*sizeof(char));
      error += im_ex_get_elem_block_l(id,
				    mss.block_id[i],
				    mss.element_types[i],
				    (long long*)&(mss.elements[i]),
				    (long long*)&(mss.nodes_per_element[i]),
				    (long long*)&(mss.element_attributes[i]));
    }

    /*connectivity*/
    for(b = 0; b < mss.num_elem_blk; b++){
      mss.elmt_node_linkage[b] = (long long*)malloc(mss.nodes_per_element[b]*mss.elements[b]*sizeof(long long));
      error += im_ex_get_elem_conn_l(id,mss.block_id[b],mss.elmt_node_linkage[b]);
    }


    if(mss.num_node_sets){
      mss.node_set_id           = (long long *) malloc(mss.num_node_sets*sizeof(long long));
      mss.num_nodes_in_node_set = (long long *) malloc(mss.num_node_sets*sizeof(long long));
      mss.node_set_nodes        = (long long **)malloc(mss.num_node_sets*sizeof(long long*));
      mss.num_df_in_node_set    = (long long *) malloc(mss.num_node_sets*sizeof(long long*));

      error += im_ex_get_node_set_ids_l(id, mss.node_set_id);


      for(i = 0; i < mss.num_node_sets; i ++){
	error += im_ex_get_node_set_param_l(id, mss.node_set_id[i],
					  (long long*)&mss.num_nodes_in_node_set[i],
					  (long long*)&mss.num_df_in_node_set[i]);

      	mss.node_set_nodes[i] = NULL;

	if(mss.num_nodes_in_node_set[i]) {
	  mss.node_set_nodes[i] = (long long *)malloc(mss.num_nodes_in_node_set[i]*sizeof(long long));
	  error += im_ex_get_node_set_l(id, mss.node_set_id[i], mss.node_set_nodes[i]);
	}
      }
    }

    /*side sets*/
    if(mss.num_side_sets){
      mss.side_set_id = (long long*)malloc(mss.num_side_sets*sizeof(long long));
      mss.num_elements_in_side_set = (long long*)malloc(mss.num_side_sets*sizeof(long long));
      mss.num_df_in_side_set = (long long*)malloc(mss.num_side_sets*sizeof(long long));
      mss.side_set_elements = (long long**)malloc(mss.num_side_sets*sizeof(long long *));
      mss.side_set_faces = (long long **)malloc(mss.num_side_sets*sizeof(long long*));

      error += im_ex_get_side_set_ids_l(id, mss.side_set_id);
      for(i = 0; i < mss.num_side_sets; i ++){
	long long ne = 0;
	error += im_ex_get_side_set_param_l(id, mss.side_set_id[i],
					  (long long*)&mss.num_elements_in_side_set[i],
					  (long long*)&mss.num_df_in_side_set[i]);

	ne = mss.num_elements_in_side_set[i];
	mss.side_set_elements[i] = (long long*)malloc(ne*sizeof(long long));
	mss.side_set_faces[i] = (long long*)malloc(ne*sizeof(long long));
	if(ne){
	  error += im_ex_get_side_set_l(id, mss.side_set_id[i],
				      mss.side_set_elements[i],
				      mss.side_set_faces[i]);

	}
      }
    }

    error += im_ex_inquire_l(id, IM_EX_INQ_QA, (long long*)&mss.num_qa_records, &fdum, cdum);

    if(mss.num_qa_records)error +=  im_ex_get_qa_l(id,mss.qaRecord);


    error += im_ex_inquire_l(id, IM_EX_INQ_INFO, (long long*)&mss.num_info_records, &fdum, cdum);
    if(mss.num_info_records) {
      mss.info_records = (char **)malloc(mss.num_info_records*sizeof(char *));/*new std::string[num_info_records];*/
      for(i = 0; i < mss.num_info_records; i ++){
	mss.info_records[i] = (char *)malloc(MAX_STR_LENGTH+1);
      }
      error += im_ex_get_info_l(id, mss.info_records);
    }


    /*nemesis data
      global info*/
    if ( im_ne_get_init_global_l(id, &mss.num_nodes_global, &mss.num_elems_global,
			       &mss.num_elm_blks_global, &mss.num_node_sets_global,
			       &mss.num_side_sets_global) < 0 )
      ++error;



    if ( im_ne_get_init_info_l(id, &mss.num_total_proc, &mss.num_proc_in_file, mss.type) < 0 )
      ++error;

    mss.elem_blk_ids_global = (long long*)malloc(mss.num_elm_blks_global*sizeof(long long));
    mss.elem_blk_cnts_global = (long long*)malloc(mss.num_elm_blks_global*sizeof(long long));

    if ( im_ne_get_eb_info_global_l(id,mss.elem_blk_ids_global,mss.elem_blk_cnts_global) < 0 )
      ++error;

    mss.ns_ids_global = (long long *)malloc(mss.num_node_sets_global*sizeof(long long));
    mss.ns_cnts_global = (long long *)malloc(mss.num_node_sets_global*sizeof(long long));
    mss.ns_df_cnts_global = (long long *)malloc(mss.num_node_sets_global*sizeof(long long));
    mss.ss_ids_global = (long long *)malloc(mss.num_side_sets_global*sizeof(long long));
    mss.ss_cnts_global = (long long *)malloc(mss.num_side_sets_global*sizeof(long long));
    mss.ss_df_cnts_global = (long long *)malloc(mss.num_side_sets_global*sizeof(long long));


    if ( mss.num_node_sets_global > 0 ) {
      if ( im_ne_get_ns_param_global_l(id,mss.ns_ids_global,mss.ns_cnts_global,
				     mss.ns_df_cnts_global) < 0 )++error;
    }

    if ( mss.num_side_sets_global > 0 ) {
      if ( im_ne_get_ss_param_global_l(id,mss.ss_ids_global,mss.ss_cnts_global,
				     mss.ss_df_cnts_global) < 0 )  ++error;
    }

    /*parallel info*/
    if ( im_ne_get_loadbal_param_l( id,
				  &mss.num_internal_nodes,
				  &mss.num_border_nodes,
				  &mss.num_external_nodes,
				  &mss.num_internal_elems,
				  &mss.num_border_elems,
				  &mss.num_node_comm_maps,
				  &mss.num_elem_comm_maps,
				  0/*unused*/ ) < 0 )++error;

    mss.internal_elements = (long long *)malloc(mss.num_internal_elems*sizeof(long long));
    mss.border_elements   = (long long *)malloc(mss.num_border_elems*sizeof(long long));
    mss.internal_nodes    = (long long *)malloc(mss.num_internal_nodes*sizeof(long long));
    mss.border_nodes      = (long long *)malloc(mss.num_border_nodes*sizeof(long long));
    mss.external_nodes    = (long long *)malloc(mss.num_external_nodes*sizeof(long long));

    if ( im_ne_get_elem_map_l( id,
			     mss.internal_elements,
			     mss.border_elements,
			     0/* not used proc_id*/ ) < 0 )++error;

    if ( im_ne_get_node_map_l( id,
			     mss.internal_nodes,
			     mss.border_nodes,
			     mss.external_nodes,
			     0/* not used proc_id*/ ) < 0 )++error;


    if(mss.num_node_comm_maps > 0){

      mss.node_cmap_node_cnts = (long long*) malloc(mss.num_node_comm_maps*sizeof(long long));
      mss.node_cmap_ids       = (long long*) malloc(mss.num_node_comm_maps*sizeof(long long));
      mss.comm_node_ids       = (long long**)malloc(mss.num_node_comm_maps*sizeof(long long*));
      mss.comm_node_proc_ids  = (long long**)malloc(mss.num_node_comm_maps*sizeof(long long*));

      mss.elem_cmap_elem_cnts = (long long*) malloc(mss.num_elem_comm_maps*sizeof(long long));
      mss.elem_cmap_ids       = (long long*) malloc(mss.num_elem_comm_maps*sizeof(long long));
      mss.comm_elem_ids       = (long long**)malloc(mss.num_elem_comm_maps*sizeof(long long*));
      mss.comm_side_ids       = (long long**)malloc(mss.num_elem_comm_maps*sizeof(long long*));
      mss.comm_elem_proc_ids  = (long long**)malloc(mss.num_elem_comm_maps*sizeof(long long*));

      if ( im_ne_get_cmap_params_l( id,
				  mss.node_cmap_ids,
				  (long long*)mss.node_cmap_node_cnts,
				  mss.elem_cmap_ids,
				  (long long*)mss.elem_cmap_elem_cnts,
				  0/*not used proc_id*/ ) < 0 )++error;

      for(j = 0; j < mss.num_node_comm_maps; j++) {
	mss.comm_node_ids[j]       = (long long *)malloc(mss.node_cmap_node_cnts[j]*sizeof(long long));
	mss.comm_node_proc_ids[j]  = (long long *)malloc(mss.node_cmap_node_cnts[j]*sizeof(long long));
	if ( im_ne_get_node_cmap_l( id,
				  mss.node_cmap_ids[j],
				  mss.comm_node_ids[j],
				  mss.comm_node_proc_ids[j],
				  0/*not used proc_id*/ ) < 0 )++error;

      }



      for(j = 0; j < mss.num_elem_comm_maps; j++) {
	mss.comm_elem_ids[j]       = (long long *)malloc(mss.elem_cmap_elem_cnts[j]*sizeof(long long));
	mss.comm_side_ids[j]       = (long long *)malloc(mss.elem_cmap_elem_cnts[j]*sizeof(long long));
	mss.comm_elem_proc_ids[j]  = (long long *)malloc(mss.elem_cmap_elem_cnts[j]*sizeof(long long));
	if ( im_ne_get_elem_cmap_l( id,
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
long long i;
long long j;
long long b;
long long ict;
long long nct;

    printf("\nExodus header info:\nTitle: %s\nDimension %lli \nNumber of Nodes %lli \nNumber of Elements %lli \nNumber of Element Blocks %lli \nNumber of Node Sets %lli \nNumber of Side Sets %lli \n\n",
	   mss.title,
	   mss.num_dim,
	   mss.num_nodes,
	   mss.num_elem,
	   mss.num_elem_blk,
	   mss.num_node_sets,
	   mss.num_side_sets);


  printf("num node set nodes %lli\n",mss.num_node_set_nodes);
  printf("num node set dfs %lli\n",mss.num_node_set_dfs);

  printf("num side set elements %lli\n",mss.num_side_set_elements);
  printf("num side set nodes %lli\n",mss.num_side_set_nodes);
  printf("num side set dfs %lli\n",mss.num_side_set_dfs);

  /*get version number*/
  printf("num block properties %lli\n",mss.num_block_properties);
  printf("num node set properties %lli\n",mss.num_node_set_properties);
  printf("num side set properties %lli\n",mss.num_side_set_properties);


  printf("A taste of coords\n");
  for(i = 0; i < mss.num_nodes && i < 10; i ++){
    if(mss.num_dim == 3) printf("X %f Y %f Z %f \n",mss.coord[i],mss.coord[i+mss.num_nodes],mss.coord[i+2*mss.num_nodes]);
    if(mss.num_dim == 2) printf("X %f Y %f \n",mss.coord[i],mss.coord[i+mss.num_nodes]);
  }

  for(i = 0; i < mss.num_dim; i ++)printf("coord name %lli %s\n",i,mss.bptr[i]);

  if (mss.num_elem){

    printf ("A tast of map\n");
    for(i = 0; i < mss.num_elem && i <  10; i ++)printf("map i=%lli, val=%lli\n",i,mss.element_order_map[i]);
    if (mss.num_elem){
      printf ("A tast of global elem numbers\n");
      for(i = 0; i < mss.num_elem && i <  10; i ++)printf("global el num  i=%lli, val=%lli\n",i,mss.global_element_numbers[i]);
    }

    if (mss.num_nodes){
      printf ("A tast of global elem numbers\n");
      for(i = 0; i < mss.num_elem && i <  10; i ++)printf("global node num  i=%lli, val=%lli\n",i,mss.global_node_numbers[i]);
    }


    /*block info*/

    for(i = 0; i < mss.num_elem_blk; i ++){
      printf("block i = %lli has id %lli \n",i,mss.block_id[i]);
      printf("block i = %lli\n",i);
      printf("block id %lli\n",mss.block_id[i]);
      printf("element_type %s\n",mss.element_types[i]);
      printf("num elements %lli\n",mss.elements[i]);
      printf("nodes per element %lli\n",mss.nodes_per_element[i]);
      printf("element attributes %lli\n",mss.element_attributes[i]);
    }

    /*connectivity*/
    for(b = 0; b < mss.num_elem_blk; b++){
      for(ict = 0; ict < mss.elements[b] && ict < 10;ict++){
	printf("block %lli element %lli connectivty ",mss.block_id[b],ict);
	for(nct = 0; nct < mss.nodes_per_element[b]; nct++){
	  printf("%lli ",mss.elmt_node_linkage[b][nct+ict*mss.nodes_per_element[b]]);
	}
	printf("\n");
      }
    }


    if(mss.num_node_sets){
      for(i = 0; i < mss.num_node_sets; i ++){
	printf("Nodeset i = %lli id = %lli has %lli nodes\n",i,mss.node_set_id[i],mss.num_nodes_in_node_set[i]);
	for(j = 0; j < mss.num_nodes_in_node_set[i] && j < 10; j ++){
	  printf("nodeset node i=%lli = %lli\n",j,mss.node_set_nodes[i][j]);
	}
      }
    }

    /*side sets*/
    if(mss.num_side_sets){
      for(i = 0; i < mss.num_side_sets; i ++){
        int ne = 0;
	printf("Side set index %lli id %lli has %lli elements\n",i,mss.side_set_id[i],mss.num_elements_in_side_set[i]);
	ne = mss.num_elements_in_side_set[i];
	if(ne){
	  for(j = 0; j < ne && j < 10; j ++){
	    printf("element %lli and face %lli\n",mss.side_set_elements[i][j],mss.side_set_faces[i][j]);
	  }
	}
      }
    }

    printf("num qa records %lli\n",mss.num_qa_records);

    for(i = 0; i < mss.num_qa_records; i ++){
      printf("\nQA Record %lli\n %s\n%s\n",
	     i,
	     mss.qaRecord[i][0],
	     mss.qaRecord[i][1]);
    }

    printf("Num Info Records %lli\n",mss.num_info_records);
    if(mss.num_info_records) {
      printf("Info Records\n");
      for(i = 0; i < mss.num_info_records; i ++){
	printf("%s",mss.info_records[i]);
      }
    }


    /*nemesis data
      global info*/

    printf("Nemesis data\n");
    printf("Num nodes global %lli\n",mss.num_nodes_global);
    printf("Num elems global %lli\n",mss.num_elems_global);
    printf("Num elm_blks global %lli\n",mss.num_elm_blks_global);
    printf("Num node sets global %lli\n",mss.num_node_sets_global);
    printf("Num side sets global %lli\n",mss.num_side_sets_global);
    printf("Num total proc %lli\n",mss.num_total_proc);
    printf("Num proc in file %lli\n",mss.num_proc_in_file);

    for(i = 0; i < mss.num_elm_blks_global; i ++){
      printf("element block index %lli has id %lli and %lli elements\n",i,mss.elem_blk_ids_global[i],mss.elem_blk_cnts_global[i]);
    }

    if ( mss.num_node_sets_global > 0 ) {
      for(i = 0; i < mss.num_node_sets_global;i ++){
	printf("global ns info for ns index %lli id %lli num_nodes = %lli num_ns_df = %lli\n",i,mss.ns_ids_global[i],mss.ns_cnts_global[i],mss.ns_df_cnts_global[i]);
      }
    }

    if ( mss.num_side_sets_global > 0 ) {
      for(i = 0; i < mss.num_side_sets_global;i ++){
	printf("global ss info for ss index %lli id %lli num_elements = %lli num_ss_df = %lli\n",i,mss.ss_ids_global[i],mss.ss_cnts_global[i],mss.ss_df_cnts_global[i]);
      }

    }

    /*parallel info*/

    printf("Loadbal params:\nnum_internal_nodes %lli\nnum_border_nodes%lli\nnum_external_nodes%lli\nnum_internal_elems%lli\nnum_border_elems%lli\nnum_node_comm_maps%lli\nnum_elem_comm_maps%lli\n",
	   mss.num_internal_nodes,
	   mss.num_border_nodes,
	   mss.num_external_nodes,
	   mss.num_internal_elems,
	   mss.num_border_elems,
	   mss.num_node_comm_maps,
	   mss.num_elem_comm_maps);

    for(i = 0; i < mss.num_internal_nodes && i < 10; i ++)printf("internal node i=%lli = %lli\n",i,mss.internal_nodes[i]);
    for(i = 0; i < mss.num_border_nodes && i < 10;   i ++)printf("border node i=%lli = %lli\n",i,  mss.border_nodes[i]);
    for(i = 0; i < mss.num_external_nodes && i < 10; i ++)printf("external node i=%lli = %lli\n",i,mss.external_nodes[i]);
    for(i = 0; i < mss.num_internal_elems && i < 10; i ++)printf("internal elem i=%lli = %lli\n",i,mss.internal_elements[i]);
    for(i = 0; i < mss.num_border_elems && i < 10;   i ++)printf("border elem i=%lli = %lli\n",i,  mss.border_elements[i]);

    if(mss.num_node_comm_maps > 0){

      for(i =0; i < mss.num_node_comm_maps; i ++){
	printf("node_cmap_id i = %lli node_cmap_id = %lli node_cmap_node_cnts = %lli\n",i,mss.node_cmap_ids[i],mss.node_cmap_node_cnts[i]);
      }
      for(i =0; i < mss.num_elem_comm_maps; i ++){
	printf("elem_cmap_id i = %lli elem_cmap_id = %lli elem_cmap_elem_cnts = %lli\n",i,mss.elem_cmap_ids[i],mss.elem_cmap_elem_cnts[i]);
      }

      for(j = 0; j < mss.num_node_comm_maps; j++) {
	for(i = 0; i < mss.node_cmap_node_cnts[j] && i < 10;i ++){
	  printf("node_cmap_id i=%lli = %lli comm_node_ids = %lli comm_node_proc_ids = %lli\n",
		 i,
		 mss.node_cmap_ids[j],
		 mss.comm_node_ids[j][i],
		 mss.comm_node_proc_ids[j][i]);
	}
      }

      for(j = 0; j < mss.num_elem_comm_maps; j++) {
	for(i = 0; i < mss.elem_cmap_elem_cnts[j] && i < 10;i ++){
	  printf("elem_cmap_id i=%lli = %lli comm_elem_ids = %lli comm_side_ids = %lli comm_elem_proc_ids = %lli\n",
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

#ifdef HAVE_EXODUS
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
  int io_word_size = sizeof(float);
  int out_id;
  int i;
  int b;


  out_id = ex_create(out_file_name, exo_access, &cpu_word_size,
		 &io_word_size);

  if (out_id < 0){
    printf("error opening file");
  }

  ex_init_params exinit;
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

  int error = 0;

  PERROR;


  if ( ex_put_init_ext(out_id, &exinit) < 0 )
    ++error;
  PERROR;

/*now write parallel global information*/

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
  if ( ne_put_loadbal_param_l( out_id,
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
                           mss.element_attributes[b] );  /* num attr */
    PERROR;
  }

/* write element connectivity information */

  for (b = 0; b < mss.num_elem_blk; b++) {
    if ( mss.elements[b] > 0 ){
      error += ex_put_elem_conn(out_id,mss.block_id[b],mss.elmt_node_linkage[b]);
      PERROR;
    }
  }


/* write in nodal boundary sets for the body. */

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
#endif
