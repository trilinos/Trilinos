/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.  
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */


#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGE_FILES 1

#if defined(__LIBCATAMOUNT__)
#include <catamount/dclock.h>
#endif

#ifdef HAVE_PARALLEL
#include <mpi.h>
#else
#include <string.h>
#include <time.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>

#include "exodusII.h"

#define DEFAULT_NUM_FIELDS     0
#define DEFAULT_FILE_NAME      "mesh"
#define DEFAULT_NUM_ITERATIONS 1
#define EXODUS_FILE_TYPE       "e"
#define MBYTES                 (1024*1024)
#define MAX_STRING_LEN         128
#define NUM_NODES_PER_ELEM     8
#define QUIT                   FALSE
#define WRITE_FILE_TYPE        "new"

/*
 *	Prototypes
 */

typedef double realtyp;

void get_file_name(const char* base, const char* ext, int rank, int nprocs,
		   const char* other,char *output);

int parse_input(
	int	argc,
        char	*argv[],
	int	*exodus,
	int	*close_files,
	char	*file_name,
        int	*num_nodal_fields,
	int     *num_global_fields,
	int     *num_element_fields,
	int     *files_per_domain,
        int	*num_iterations,
	int     *sleep_time
);

int read_exo_mesh ( char *file_name, int rank, int *ndim, int num_domains, int *num_nodal_fields,
		    int *num_global_fields, int *num_element_fields,
		    int *num_timesteps, 
		    int sleep_time, int num_iterations, 
                    int *num_nodes, int **node_map, 
                    int *num_elems, int **elem_map,
                    realtyp **x_coords, realtyp **y_coords, realtyp **z_coords, int **loc_connect );

int write_exo_mesh ( char *file_name, int rank, int ndim, int num_domains, 
		     int num_nodal_fields, int num_global_fields, int num_element_fields,
		     int num_timesteps, int files_per_domain,
		     int sleep_time, int num_iterations,
                     int loc_num_nodes, int *node_map, 
                     int loc_num_elems, int *elem_map,
                     realtyp *x_coords, realtyp *y_coords, realtyp *z_coords, int *loc_connect, int close_files);

double my_timer()
{
  double t1 = 0.0;

#if !defined(__LIBCATAMOUNT__)
#ifdef HAVE_PARALLEL
  t1 = MPI_Wtime();
#else
  clock_t ctime = clock();
  t1 = ctime / (double)CLOCKS_PER_SEC;
#endif
#else
  return dclock();
#endif
  return t1;
}

/***********************************************************************
 *
 *  Main function
 * 
 ***********************************************************************/

int main( int argc, char **argv )
{
  int rank, num_domains;
  int quit=FALSE;
  int loc_num_nodes, loc_num_elems;
  int *loc_connect = NULL;

#ifdef HAVE_PARALLEL
  MPI_Info     mpi_info_object = MPI_INFO_NULL;				/* Copy of MPI Info object.		*/
#endif
  int         *elem_map                   = NULL;
  int          exodus =                     TRUE;			/* TRUE, perform EXODUS benchmark; FALSE don't */
  int          close_files = FALSE;
  char         file_name[MAX_STRING_LEN] =  DEFAULT_FILE_NAME;		/* Input file name.				*/
  /* object, FALSE otherwise. Should always be	*/
  /* TRUE in the current implementation.		*/
  int          num_nodal_fields =           DEFAULT_NUM_FIELDS;
  int          num_global_fields =          DEFAULT_NUM_FIELDS;
  int          num_element_fields =         DEFAULT_NUM_FIELDS;
  int          num_timesteps = 0;
  int          sleep_time = 0;           
  int          files_per_domain = 1;
  int          num_iterations =             DEFAULT_NUM_ITERATIONS;
#ifdef HAVE_PARALLEL
  static const char *hints[] = {					/* List of MPI Info hints that if defined in	*/
    "cb_buffer_size",							/* the environment process 0, will be used to	*/
    "cb_nodes",								/* set key/value pairs in the MPI	*/
    "ind_rd_buffer_size",						/* Info object.					*/
    "ind_wr_buffer_size",
    "cb_config_list",
    "romio_cb_read",
    "romio_cb_write",
    "romio_ds_read",
    "romio_ds_write",
    "romio_no_indep_rw"
  };
  char         key_name[MAX_STRING_LEN];				/* MPI Info object key name.			*/
  int          key;							/* MPI Info object key index.			*/
  int          key_exists;						/* TRUE, if the key exists in the MPI Info	*/
  const int    nhints =                     10;                     	/* Number of items in hints list.		*/
  int          nkeys;							/* Number of keys in a MPI Info object.		*/
  char         value[MAX_STRING_LEN];					/* Value of a key/value pair in a MPI Info	*/
#endif
  /* object.					*/
  realtyp       *x_coords = NULL;
  realtyp       *y_coords = NULL;
  realtyp       *z_coords = NULL;
  int         ndim;
#ifdef HAVE_PARALLEL
  MPI_Info    new_mpi_info_object;	
#endif
  
  /*
   *	Initialize Stuff
   */

   ex_opts (EX_VERBOSE | EX_ABORT);

#if defined(__LIBCATAMOUNT__)
  setlinebuf(stderr);
#endif

#ifdef HAVE_PARALLEL
  MPI_Init     ( &argc, &argv                 );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank        );
  MPI_Comm_size( MPI_COMM_WORLD, &num_domains );
#else
  rank = 0;
  num_domains = 1;
#endif	
  /*
   *	Processor 0: parse the command line arguments.
   */

  if ( rank == 0 ) {
    quit = (1 == parse_input(argc, argv,
			    &exodus,	       &close_files,	       file_name,
			    &num_nodal_fields, &num_global_fields, &num_element_fields,
			     &files_per_domain, &num_iterations,    &sleep_time));
  }


  /*
   *	Broadcast Input
   */
         
#ifdef HAVE_PARALLEL
  MPI_Bcast ( &quit,         	1,              MPI_INT,  0, MPI_COMM_WORLD );
#endif
  
  if (quit) {
#ifdef HAVE_PARALLEL
    MPI_Finalize();
#endif
    exit(0);
  }

#ifdef HAVE_PARALLEL
  MPI_Bcast( &exodus,		 1,              MPI_INT,  0, MPI_COMM_WORLD );
  MPI_Bcast( &close_files,	 1,              MPI_INT,  0, MPI_COMM_WORLD );
  MPI_Bcast( file_name,		 MAX_STRING_LEN, MPI_CHAR, 0, MPI_COMM_WORLD );
  MPI_Bcast( &num_nodal_fields,	 1,              MPI_INT,  0, MPI_COMM_WORLD );
  MPI_Bcast( &num_global_fields, 1,              MPI_INT,  0, MPI_COMM_WORLD );
  MPI_Bcast( &num_element_fields,1,              MPI_INT,  0, MPI_COMM_WORLD );
  MPI_Bcast( &num_iterations,	 1,              MPI_INT,  0, MPI_COMM_WORLD );
  MPI_Bcast( &files_per_domain,	 1,              MPI_INT,  0, MPI_COMM_WORLD );

  /* env_mpi_hints( nhints, hints ); */
  {
    char            *env;                           /* Contents of environmental variable.  */
    int              hint;                          /* ROMIO hint index.                    */
    char             hint_value[MAX_STRING_LEN];    /* ROMIO hint value.                    */
    int              rank;                          /* MPI process rank.                    */

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    /*      The "value" of the hint is obtained from the environment of
     *      processor 0 only. The value is broadcast to the other
     *      processors.
     */
    for ( hint = 0; hint < nhints; hint++ ) {
      if ( rank == 0 ) {
	env = getenv( hints[hint] );
	if ( env != NULL )
	  strcpy( hint_value, env );
	else
	  hint_value[0] = 0;
      }
      MPI_Bcast( hint_value, MAX_STRING_LEN, MPI_CHAR, 0, MPI_COMM_WORLD );
      if ( hint_value[0] ) {
	if ( mpi_info_object == MPI_INFO_NULL )
	  MPI_Info_create( &mpi_info_object );
	MPI_Info_set( mpi_info_object, hints[hint], hint_value );
      }
    }
  }
#endif

  if (rank == 0) {
    fprintf(stderr, "\nEXODUSII 3D Benchmark\n\n"              );
    fprintf(stderr, "   Number of Domains\t\t%8d\n",      num_domains );
    fprintf(stderr, "   Number of Files/Domain\t%8d\n",   files_per_domain );
    fprintf(stderr, "   Number of Iterations\t\t%8d\n",   num_iterations );
                 
#ifdef HAVE_PARALLEL
    if ( mpi_info_object != MPI_INFO_NULL ) {
      fprintf(stderr, "   MPI Hint Status\n" );

      MPI_Info_get_nkeys( mpi_info_object, &nkeys );

      for ( key = 0; key < nkeys; key++ ) {
	MPI_Info_get_nthkey( mpi_info_object, key, key_name                                );
	MPI_Info_get       ( mpi_info_object, key_name, MAX_STRING_LEN, value, &key_exists );

	fprintf(stderr, "      %s\t\t\t%s\n", key_name, value );
      }

      MPI_Info_free( &mpi_info_object );
    }
    else
      fprintf(stderr, "   MPI Hint Status\tMPI_INFO_NULL\n" );     		          
#endif
  }

  if ( exodus ) {
    int *node_map = NULL;
    if (0 == read_exo_mesh(file_name, rank, &ndim, num_domains, 
			   &num_nodal_fields, &num_global_fields, &num_element_fields,
			   &num_timesteps, sleep_time, num_iterations,
			   &loc_num_nodes, &node_map, &loc_num_elems, &elem_map,
			   &x_coords, &y_coords, &z_coords, &loc_connect)
	) {
      write_exo_mesh(file_name, rank, ndim, num_domains, 
		     num_nodal_fields, num_global_fields, num_element_fields,
		     num_timesteps, files_per_domain, sleep_time, num_iterations,
		     loc_num_nodes, node_map,  loc_num_elems, elem_map, 
		     x_coords, y_coords, z_coords, loc_connect, close_files);
			
      free( elem_map    );
      free( loc_connect );
      free( node_map    );
      free( x_coords    );
      free( y_coords    );
      free( z_coords    );
    }
}
#ifdef HAVE_PARALLEL
  MPI_Finalize();
#endif
  return(0);
}

/***********************************************************************
 *
 *  Parse Input
 * 
 ***********************************************************************/

int parse_input (int	argc,
		 char	*argv[],
		 int	*exodus,
		 int    *close_files,
		 char	*file_name,
		 int	*num_nodal_fields,
		 int	*num_global_fields,
		 int	*num_element_fields,
		 int    *files_per_domain,
		 int	*num_iterations,
		 int    *sleep_time) {
  int arg = 0;	/* Argument index.	*/
        
  while ( ++arg < argc ) {
    if ( strcmp( "-c", argv[arg] ) == 0 ) {
      if ( ++arg < argc ) *num_nodal_fields = atoi( argv[arg] );
    }
    else if ( strcmp( "-nv", argv[arg] ) == 0 ) {
      if ( ++arg < argc ) *num_nodal_fields = atoi( argv[arg] );
    }
    else if ( strcmp( "-gv", argv[arg] ) == 0 ) {
      if ( ++arg < argc ) *num_global_fields = atoi( argv[arg] );
    }
    else if ( strcmp( "-ev", argv[arg] ) == 0 ) {
      if ( ++arg < argc ) *num_element_fields = atoi( argv[arg] );
    }
    else if ( strcmp( "-f", argv[arg] ) == 0 ) { 
      if ( ++arg < argc ) strcpy( file_name, argv[arg] );
    }
    else if ( strcmp( "-M", argv[arg] ) == 0 ) { 
      if ( ++arg < argc ) *files_per_domain = atoi( argv[arg] );
    }
    else if ( strcmp( "-i", argv[arg] ) == 0 ) { 
      if ( ++arg < argc ) *num_iterations = atoi( argv[arg] );
    }
    else if ( strcmp( "-w", argv[arg] ) == 0 ) { 
      if ( ++arg < argc ) *sleep_time = atoi( argv[arg] );
    }
    else if ( strcmp( "-x", argv[arg] ) == 0 ) {
      *exodus = TRUE;
    }
    else if ( strcmp( "-C", argv[arg] ) == 0 ) {
      *close_files = TRUE;
    }
    else if ( (strcmp( "-h", argv[arg] ) == 0) || (strcmp( "-u", argv[arg] ) == 0) ) {
      fprintf(stderr, "                                                                \n" );
      fprintf(stderr, "NAME                                                            \n" );
      fprintf(stderr, "                                                                \n" );
      fprintf(stderr, "rd_wt_mesh - reads and writes a mesh in parallel for performance\n" );
      fprintf(stderr, "             benchmarking.                                      \n" );
      fprintf(stderr, "                                                                \n" );
      fprintf(stderr, "SYNOPSIS                                                        \n" );
      fprintf(stderr, "                                                                \n" );
      fprintf(stderr, "rd_wt_mesh [-S] [-c fields] [-f file_name] [-h] [-i iterations] \n" );
      fprintf(stderr, "           [-s] [-u] [-x] [-w]                                  \n" );
      fprintf(stderr, "                                                                \n" );
      fprintf(stderr, "DESCRIPTION                                                     \n" );
      fprintf(stderr, "                                                                \n" );
      fprintf(stderr, "This program reads and writes a mesh in parallel for performance\n" );
      fprintf(stderr, "benchmarking. The first Exodus database file read by \n" );
      fprintf(stderr, "this program is created by create_mesh. Performance summaries   \n" );
      fprintf(stderr, "are written to stdout.                                          \n" );
      fprintf(stderr, "                                                                \n" );
      fprintf(stderr, "OPTIONS                                                         \n" );
      fprintf(stderr, "                                                                \n" );
      fprintf(stderr, "-c fields      number of fields. Default: %d                    \n", DEFAULT_NUM_FIELDS     );
      fprintf(stderr, "-f file_name   file name prefix for all read files:             \n" );
      fprintf(stderr, "                                                                \n" );
      fprintf(stderr, "               'file_name'%s.nproc.rank [EXODUS II file]        \n", EXODUS_FILE_TYPE       );
      fprintf(stderr, "                                                                \n" );
      fprintf(stderr, "               Default: %s                                      \n", DEFAULT_FILE_NAME      );
      fprintf(stderr, "-M files       number of files/domain. Default: 1               \n" );
      fprintf(stderr, "-h             display help/usage information                   \n" );
      fprintf(stderr, "-i iterations  number of iterations. Default: %d                \n", DEFAULT_NUM_ITERATIONS );
      fprintf(stderr, "-C             minimize open files.                             \n" );
      fprintf(stderr, "-u             display help/usage information                   \n" );
      fprintf(stderr, "-w time        wait (sleep) specified time between timesteps.\n" );
      return(1);
    }
    else {
      fprintf(stderr, "Unknown option: %s\n", argv[arg]                         );
      fprintf(stderr, "Enter rd_wt_mesh -h for description of valid options.\n" );
                        
      return(1);
    }
  }
  return(0);
}

/***********************************************************************
 ***********************************************************************/
int read_exo_mesh (char *file_name, int rank, int *num_dim, int num_domains,
		   int *num_nodal_fields, int *num_global_fields, int *num_element_fields,
		   int *num_timesteps, 
		   int sleep_time, int num_iterations,
		   int *num_nodes, int **node_map,
		   int *num_elems, int **elem_map,
		   realtyp **x_coords, realtyp **y_coords, realtyp **z_coords, int **connect )

{
  int CPU_word_size=sizeof(realtyp);
  int IO_word_size=0;
  int exoid, err, num_elem_blk, num_node_sets, num_side_sets;
  int num_nodes_per_elem, num_attrs, num_vars, i, iter;
  size_t len_connect;
  size_t file_size;
  double fdum;
  char *cdum = 0;
  
  struct stat  file_status;  
  size_t glob_file_size;
  size_t glob_raw_data_vol;
  size_t raw_data_vol = 0;

  float version;
  realtyp *globals = NULL;
  
  double tstart, tend, t_tmp1, t_tmp2;
  double raw_read_time, max_raw_read_time=0.0, min_raw_read_time=DBL_MAX;
  double cum_raw_read_time=0.0;
  double total_time, max_total_time=0.0, min_total_time=DBL_MAX;
  double raw_sleep_time = 0.0;
  double cum_total_time=0.0;

  char tmp_name[MAX_STRING_LEN], title[MAX_STRING_LEN+1]; 
  char type[MAX_STRING_LEN+1];

  for (iter=0; iter<num_iterations; iter++) {
    /* open the EXODUS file */
    get_file_name(file_name, EXODUS_FILE_TYPE, rank, num_domains, NULL, tmp_name);
    exoid = ex_open (tmp_name, EX_READ, &CPU_word_size, &IO_word_size, &version);
      
    if (exoid < 0) {
      printf ("after ex_open\n");
      return(1);
    }

    raw_read_time = 0.0;
    raw_data_vol = 0;

    tstart = my_timer();
    

    err = ex_get_init (exoid, title, num_dim, 
		       num_nodes, num_elems, &num_elem_blk, 
		       &num_node_sets, &num_side_sets);
  
    if (err) {
      printf ("after ex_get_init, error = %d\n", err);
      ex_close (exoid);
      return(1);
    }
  
    len_connect = (size_t)NUM_NODES_PER_ELEM * (*num_elems);
  
    /* malloc things we need */
    if( iter == 0) {
      *elem_map = malloc(sizeof(int) * (*num_elems));
      assert(elem_map);
      *connect =  malloc(sizeof(int) * len_connect );
      assert(connect);
      *node_map = malloc(sizeof(int) * (*num_nodes));
      assert(node_map);
      
      *x_coords = malloc(sizeof(realtyp) * (*num_nodes));
      assert(x_coords);
      *y_coords = malloc(sizeof(realtyp) * (*num_nodes));
      assert(y_coords);
      *z_coords = malloc(sizeof(realtyp) * (*num_nodes));
      assert(z_coords);
    }
  
    t_tmp1 = my_timer();
  
    err = ex_get_coord (exoid, *x_coords, *y_coords, *z_coords);
  
    t_tmp2 = my_timer();
  
    raw_read_time += t_tmp2-t_tmp1;
    raw_data_vol += sizeof(realtyp) * 3 * (*num_nodes);
  
    if (err) {
      printf ("after ex_get_coord, error = %d\n", err);
      ex_close (exoid);
      return (1 );
    }
  
    err = ex_get_block 
      (exoid, EX_ELEM_BLOCK, 10, type, num_elems, &num_nodes_per_elem, 0, 0, &num_attrs);
  
    if (err) {
      printf ("after ex_get_elem_block, error = %d\n", err);
      ex_close (exoid);
      return( 1 );
    }
  
    t_tmp1 = my_timer();
  
    err = ex_get_conn (exoid, EX_ELEM_BLOCK, 10, *connect, 0, 0);
  
    t_tmp2 = my_timer();
  
    raw_read_time += t_tmp2-t_tmp1;
    raw_data_vol += sizeof(int) * len_connect;
  
    if (err) {
      printf ("after ex_get_elem_conn, error = %d\n", err);
      ex_close (exoid);
      return( 1 );
    }
  
    /* read element and node maps */
    t_tmp1 = my_timer();
  
    err = ex_get_id_map (exoid, EX_NODE_MAP, *node_map);
    t_tmp2 = my_timer();
    raw_data_vol += sizeof(int) * (*num_nodes);
    raw_read_time += t_tmp2-t_tmp1;
  
    t_tmp1 = my_timer();
    err = ex_get_id_map (exoid, EX_ELEM_MAP, *elem_map);
    t_tmp2 = my_timer();
    raw_data_vol += sizeof(int)*(*num_elems);
    raw_read_time += t_tmp2-t_tmp1;
  
  
    /* read results variables */
    err = ex_get_variable_param (exoid, EX_NODAL, &num_vars);
    
    if (err) {
      printf ("after ex_get_var_param, error = %d\n", err);
      ex_close (exoid);
      return( 1 );
    }
    *num_nodal_fields = num_vars;
    
    err = ex_get_variable_param (exoid, EX_GLOBAL, &num_vars);
    
    if (err) {
      printf ("after ex_get_var_param, error = %d\n", err);
      ex_close (exoid);
      return( 1 );
    }
    *num_global_fields = num_vars;
    if (*num_global_fields > 0) {
      globals = malloc(*num_global_fields * sizeof(realtyp));
      assert(globals);
    }
    
    err = ex_get_variable_param (exoid, EX_ELEM_BLOCK, &num_vars);
    
    if (err) {
      printf ("after ex_get_var_param, error = %d\n", err);
      ex_close (exoid);
      return( 1 );
    }
    *num_element_fields = num_vars;
    
    /* read number of timesteps */
    ex_inquire (exoid, EX_INQ_TIME, num_timesteps, &fdum, cdum);
    if (rank == 0) {
      fprintf(stderr, "   Number of Elements\t\t%8d (per domain)\n",   *num_elems );
      fprintf(stderr, "   Number of Nodes\t\t%8d (per domain)\n",      *num_nodes );
      fprintf(stderr, "   Number of Global  Fields\t%8d\n", *num_global_fields);
      fprintf(stderr, "   Number of Nodal   Fields\t%8d\n", *num_nodal_fields);
      fprintf(stderr, "   Number of Element Fields\t%8d\n", *num_element_fields);
      fprintf(stderr, "   Number of Timesteps\t\t%8d\n", *num_timesteps);
      if (sleep_time > 0) {
	fprintf(stderr, "   Timestep Sleep time\t%8d\n", sleep_time);
      }
    }		
    
    if (*num_nodal_fields+ *num_global_fields+ *num_element_fields > 0) {
      int t;
  
      if (rank == 0) {
	fprintf (stderr, "\nReading Timestep: ");
      }
      for (t = 0; t < *num_timesteps; t++) {
	if (sleep_time > 0) {
	  t_tmp1 = my_timer();
	  sleep(sleep_time);
	  t_tmp2 = my_timer();
	  raw_sleep_time += t_tmp2-t_tmp1;
	}
	if (rank == 0) {
	  fprintf(stderr," %d,", t+1);
	}

	for (i=1; i <= *num_nodal_fields; i++) {
	  t_tmp1 = my_timer();
	  err = ex_get_var (exoid, t+1, EX_NODAL, i, 0, *num_nodes, *x_coords);
	  t_tmp2 = my_timer();
  
	  raw_read_time += t_tmp2-t_tmp1;
	  raw_data_vol += sizeof(realtyp)*(*num_nodes);
  
	  if (err) {
	    printf ("after ex_get_nodal_var, error = %d\n", err);
	    ex_close (exoid);
	    return( 1 );
	  }
	}

	t_tmp1 = my_timer();
	err = ex_get_var (exoid, t+1, EX_GLOBAL, 0, 0, *num_global_fields, globals);
	t_tmp2 = my_timer();
  
	raw_read_time += t_tmp2-t_tmp1;
	raw_data_vol += sizeof(realtyp) * *num_global_fields;
  
	if (err) {
	  printf ("after ex_get_glob_vars, error = %d\n", err);
	  ex_close (exoid);
	  return( 1 );
	}

	for (i=1; i <= *num_element_fields; i++) {
	  t_tmp1 = my_timer();
	  err = ex_get_var (exoid, t+1, EX_ELEM_BLOCK, i, 10, *num_elems, *x_coords);
	  t_tmp2 = my_timer();
  
	  raw_read_time += t_tmp2-t_tmp1;
	  raw_data_vol += sizeof(realtyp) * (*num_elems);
  
	  if (err) {
	    printf ("after ex_get_elem_var, error = %d\n", err);
	    ex_close (exoid);
	    return( 1 );
	  }
	}
      }
	
      if (rank == 0) {
	fprintf(stderr,"\n");
      }
    }
  
    err = ex_close (exoid);
    if (err) {
      printf ("after ex_close, error = %d\n", err);
      return ( 1 );
    }
    tend = my_timer();
  
    total_time = tend - tstart - raw_sleep_time;
    if (total_time > max_total_time) max_total_time = total_time;
    if (total_time < min_total_time) min_total_time = total_time;
    cum_total_time += total_time;

    if (raw_read_time > max_raw_read_time) max_raw_read_time = raw_read_time;
    if (raw_read_time < min_raw_read_time) min_raw_read_time = raw_read_time;
    cum_raw_read_time += raw_read_time;
  }  /* end of for (iter...) */

#ifdef HAVE_PARALLEL
  MPI_Allreduce (&raw_data_vol, &glob_raw_data_vol, 1, MPI_LONG, MPI_SUM,
		 MPI_COMM_WORLD);
#else
  glob_raw_data_vol = raw_data_vol;
#endif
  /*
   *	Get File Sizes
   *
   *	Note: On ASCI Red, a specialized "stat", named "estat", was added to
   *	accommodate file sizes upto 16GB.                          3/27/2002
   */

  if (stat( tmp_name, &file_status)) {
    if ( rank == 0 ) 
      fprintf ( stderr, "Exodus Read: cannot get %s file size.\n", tmp_name );
    return( 1 );
  }
  else
    file_size = file_status.st_size;

#ifdef HAVE_PARALLEL
  MPI_Allreduce(&file_size,
		&glob_file_size,
		1,
		MPI_LONG,
		MPI_SUM,
                MPI_COMM_WORLD);
#else
  glob_file_size = file_size;
#endif
  
  if ( rank == 0 ) {
    fprintf(stderr, "                                                                \n" );
    fprintf(stderr, "   Exodus Read Results                                          \n" );
    fprintf(stderr, "                                                                \n" );
    fprintf(stderr, "      Sizes (bytes)                                             \n" );
    fprintf(stderr, "         File         %14ld                                   \n", (long)glob_file_size);
    fprintf(stderr, "         Raw Data     %14ld                                   \n", (long)glob_raw_data_vol);
    fprintf(stderr, "         Difference   %14ld (%5.2f%%)                         \n",
	    (long)(glob_file_size - glob_raw_data_vol),
	    (1.0*glob_file_size-1.0*glob_raw_data_vol)/(0.01*glob_file_size) );
    fprintf(stderr, "                                                                \n" );
    fprintf(stderr, "      Times (sec)                 \t Minimum\t Maximum\t Average\n" );
    fprintf(stderr, "         Raw Data Read (sec)     \t%8.4g\t%8.4g\t%8.4g         \n",
	    min_raw_read_time, max_raw_read_time, cum_raw_read_time / num_iterations );
    fprintf(stderr, "         All Other Read (sec)    \t        \t        \t%8.4g   \n",
	    (cum_total_time - cum_raw_read_time) / num_iterations );
    fprintf(stderr, "         Total Read (sec)        \t        \t        \t%8.4g   \n",
	    cum_total_time / num_iterations );
    fprintf(stderr, "                                                                \n" );
    fprintf(stderr, "      Input Bandwidths (MiB/sec)   \t Minimum\t Maximum\t Average\n" );
    fprintf(stderr, "         Raw Data Read (MiB/sec)   \t%8.4g\t%8.4g\t%8.4g         \n",
	    (double)glob_raw_data_vol / max_raw_read_time / MBYTES,
	    (double)glob_raw_data_vol / min_raw_read_time / MBYTES,
	    (double)glob_raw_data_vol / cum_raw_read_time / MBYTES * num_iterations);
    fprintf(stderr, "         Raw + Meta Data Read (MiB/sec)\t        \t        \t%8.4g   \n",
	    (double)glob_file_size / cum_total_time / MBYTES * num_iterations );
  }
  if (*num_global_fields > 0) 
    free(globals);
  return(0);
}

/***********************************************************************
 ***********************************************************************/
int write_exo_mesh (char *file_name, int rank, int num_dim, int num_domains, int num_nodal_fields,
		    int num_global_fields, int num_element_fields, int num_timesteps,
		    int files_per_domain, int sleep_time, int num_iterations,
		    int num_nodes, int *node_map,
		    int num_elems, int *elem_map,
		    realtyp *x_coords, realtyp *y_coords, realtyp *z_coords, int *connect, int close_files )

{
  int CPU_word_size = sizeof(realtyp);
  int IO_word_size  = sizeof(realtyp);
  int j, t, npd, err, num_elem_blk, num_node_sets, num_side_sets;
  int iter;
  int *elem_var_tab = NULL;
  size_t file_size;

  struct stat  file_status; 
  size_t glob_file_size;
  size_t glob_raw_data_vol;
  size_t raw_data_vol = 0;

  realtyp *globals = NULL;

  double raw_open_close_time = 0.0;
  double cum_open_close_time = 0.0;
  double min_raw_open_close_time = DBL_MAX;
  double max_raw_open_close_time = 0.0;
  double tstart, tend, t_tmp1, t_tmp2;
  double raw_write_time, max_raw_write_time=0.0, min_raw_write_time=DBL_MAX;
  double cum_raw_write_time=0.0;
  double total_time, max_total_time=0.0, min_total_time=DBL_MAX;
  double cum_total_time=0.0;
  double raw_sleep_time=0.0;
  double put_time_time=0.0;
  
  char tmp_name[MAX_STRING_LEN];
  char base_name[MAX_STRING_LEN];
  char **gvar_name = NULL;
  char **nvar_name = NULL;
  char **evar_name = NULL;

  int *exoid = NULL;
  exoid = malloc(files_per_domain * sizeof(int));

  raw_open_close_time = 0.0;
  for (iter=0; iter<num_iterations; iter++) {
    if (!close_files) {
      t_tmp1 = my_timer();
      for (npd = 0; npd < files_per_domain; npd++) {
	/* create the EXODUS file */
	
	sprintf(base_name, "%s_%d", file_name, npd);
	get_file_name(base_name, EXODUS_FILE_TYPE, rank, num_domains, WRITE_FILE_TYPE, tmp_name);
	exoid[npd] = ex_create (tmp_name, EX_CLOBBER, &CPU_word_size, &IO_word_size);
	
	if (exoid[npd] < 0) {
	  printf ("after ex_create\n");
	  
	  return( 1 );
	}
      }
      t_tmp2 = my_timer();
      raw_open_close_time += (t_tmp2 - t_tmp1);
    }
    raw_write_time = 0.0;
    raw_data_vol = 0;

    tstart = my_timer();
  
    num_elem_blk = 1;
    num_node_sets = 0;
    num_side_sets = 0;
  
    for (npd = 0; npd < files_per_domain; npd++) {

      if (close_files) {
	/* create the EXODUS file */
	sprintf(base_name, "%s_%d", file_name, npd);
	get_file_name(base_name, EXODUS_FILE_TYPE, rank, num_domains, WRITE_FILE_TYPE, tmp_name);
	exoid[npd] = ex_create (tmp_name, EX_CLOBBER, &CPU_word_size, &IO_word_size);
	
	if (exoid[npd] < 0) {
	  printf ("after ex_create\n");
	  return( 1 );
	}
      }
      err = ex_put_init (exoid[npd], "This is an EXODUSII performance test.", num_dim, 
			 num_nodes, num_elems, num_elem_blk, 
			 num_node_sets, num_side_sets);
  
      if (err) {
	printf ("after ex_put_init, error = %d\n", err);
	ex_close (exoid[npd]);

	return( 1);
      }
  
    
#if 0
      {
	int ids[1] = {10};
	int num_elem_per_block[1];
	char *names[1] = {"hex"};
	int num_node_per_elem[1];
	int num_attr_per_block[1];
	int write_map = num_domains > 1 ? TRUE : FALSE;
	num_elem_per_block[0] = num_elems;
	num_node_per_elem[0]  = NUM_NODES_PER_ELEM;
	num_attr_per_block[0] = 0;
	err = ex_put_concat_elem_block (exoid[npd], ids, names, num_elem_per_block,
					num_node_per_elem, num_attr_per_block, write_map);
      }
#else
      err = ex_put_block 
	(exoid[npd], EX_ELEM_BLOCK, 10, "hex", num_elems, NUM_NODES_PER_ELEM, 0, 0, 0);
#endif
      if (err) {
	printf ("after ex_put_elem_block, error = %d\n", err);
	ex_close (exoid[npd]);

	return( 1 );
      }
  
      t_tmp1 = my_timer();
      err = ex_put_coord (exoid[npd], x_coords, y_coords, z_coords);
      t_tmp2 = my_timer();
  
      raw_write_time += t_tmp2-t_tmp1;
      raw_data_vol += sizeof(realtyp)*num_dim * num_nodes;
  
      if (err) {
	printf ("after ex_put_coord, error = %d\n", err);
	ex_close (exoid[npd]);

	return( 1 );
      }
  
      t_tmp1 = my_timer();
      err = ex_put_elem_conn (exoid[npd], 10, connect);
      t_tmp2 = my_timer();
  
      raw_write_time += t_tmp2-t_tmp1;
      raw_data_vol += sizeof(int) * num_elems * NUM_NODES_PER_ELEM;
  
      if (err) {
	printf ("after ex_put_elem_conn, error = %d\n", err);
	ex_close (exoid[npd]);

	return( 1);
      }
  
      /* write out element and node maps */
      t_tmp1 = my_timer();
      err = ex_put_id_map (exoid[npd], EX_NODE_MAP, node_map);
      t_tmp2 = my_timer();
  
      raw_write_time += t_tmp2-t_tmp1;
      raw_data_vol += sizeof(int)*num_nodes;
  
      if (err) {
	printf ("after ex_put_id_map, error = %d\n", err);
	ex_close (exoid[npd]);

	return( 1 );
      }
  
      t_tmp1 = my_timer();
  
      err = ex_put_id_map (exoid[npd], EX_ELEM_MAP, elem_map);
  
      t_tmp2 = my_timer();
  
      raw_write_time += t_tmp2-t_tmp1;
      raw_data_vol += sizeof(int)*num_elems;
  
      if (err) {
	printf ("after ex_put_id_map, error = %d\n", err);
	ex_close (exoid[npd]);

	return( 1 );
      }

      /* write out simulated results fields;
	 we'll just write out the x coordinate field 'num_element_fields' times */
      if (num_element_fields > 0) {
	if (npd == 0) {
	  elem_var_tab = malloc(num_element_fields * sizeof(int));
	  assert(elem_var_tab);
	  for (j = 0; j < num_element_fields; j++)
	    elem_var_tab[j] = 1;
	}
      } else {
	elem_var_tab = 0;
      }
      err = ex_put_all_var_param(exoid[npd],
				 num_global_fields,
				 num_nodal_fields,
				 num_element_fields,
				 elem_var_tab,
				 0, 0, 0, 0);
      if (err) {
	fprintf(stderr, "after ex_put_all_var_param, error = %d\n", err);
	ex_close (exoid[npd]);
	exit(1);
      }

      if (num_nodal_fields > 0) {

	if (npd == 0) {
	  nvar_name = malloc (num_nodal_fields * sizeof(char *));
	  assert(nvar_name);
	  for (j=0; j<num_nodal_fields; j++) {
	    nvar_name[j] = malloc ((MAX_STRING_LEN+1) * sizeof (char));
	    sprintf (nvar_name[j], "node_field_%d", j+1);
	  }
	}
	err = ex_put_variable_names (exoid[npd], EX_NODAL, num_nodal_fields, nvar_name);
	if (npd == files_per_domain-1) {
	  for (j=0; j<num_nodal_fields; j++) {
	    free(nvar_name[j]);
	  }
	  free(nvar_name);
	}
      }

      if (num_global_fields > 0) {
	if (npd == 0) {
	  globals = malloc(num_global_fields * sizeof(realtyp));
	  gvar_name = malloc (num_global_fields * sizeof(char *));
	  for (j=0; j<num_global_fields; j++) {
	    gvar_name[j] = malloc ((MAX_STRING_LEN+1) * sizeof (char));
	    sprintf (gvar_name[j], "global_field_%d", j+1);
	    globals[j] = j;
	  }
	}

	err = ex_put_variable_names (exoid[npd], EX_GLOBAL, num_global_fields, gvar_name);
	
	if (npd == files_per_domain-1) {
	  for (j=0; j<num_global_fields; j++) {
	    free(gvar_name[j]);
	  }
	  free(gvar_name);
	}
      }

      if (num_element_fields > 0) {
	if (npd == 0) {
	  evar_name = malloc (num_element_fields * sizeof(char *));
	  for (j=0; j<num_element_fields; j++) {
	    evar_name[j] = malloc ((MAX_STRING_LEN+1) * sizeof (char));
	    sprintf (evar_name[j], "element_field_%d", j+1);
	  }
	}

	err = ex_put_variable_names (exoid[npd], EX_ELEM_BLOCK, num_element_fields, evar_name);

	if (npd == files_per_domain-1) {
	  free(elem_var_tab);
	  for (j=0; j<num_element_fields; j++) {
	    free(evar_name[j]);
	  }
	  free(evar_name);
	}
      }
      if (close_files) {
	t_tmp1 = my_timer();
	ex_close(exoid[npd]);
	t_tmp2 = my_timer();
	raw_open_close_time += (t_tmp2 - t_tmp1);
      }
    }
    
    if (num_nodal_fields+num_global_fields+num_element_fields > 0) {
      if ( rank == 0 ) {
	fprintf (stderr, "\nWriting Timestep: ");
      }
      for (t=0; t<num_timesteps; t++) {
	for (npd = 0; npd < files_per_domain; npd++) {
	  realtyp time = t;

	  if (close_files) {
	    float version;
	    sprintf(base_name, "%s_%d", file_name, npd);
	    t_tmp1 = my_timer();
	    get_file_name(base_name, EXODUS_FILE_TYPE, rank, num_domains, WRITE_FILE_TYPE, tmp_name);
	    exoid[npd] = ex_open(tmp_name, EX_WRITE, &CPU_word_size, &IO_word_size, &version);
	    t_tmp2 = my_timer();
	    raw_open_close_time += (t_tmp2 - t_tmp1);
	  }

	  if (sleep_time > 0) {
	    t_tmp1 = my_timer();
	    sleep(sleep_time);
	    t_tmp2 = my_timer();
	    raw_sleep_time += t_tmp2-t_tmp1;
	  }
	  t_tmp1 = my_timer();
	  ex_put_time(exoid[npd], t+1, &time);
	  t_tmp2 = my_timer();
	  put_time_time  += t_tmp2-t_tmp1;
	  if ( rank == 0) {
	    if (npd == 0) {
	      fprintf(stderr, " %d", t+1);
	    } else {
	      fprintf(stderr, ".");
	    }
	  }
	  if (num_global_fields > 0) {
	    t_tmp1 = my_timer();
	    err = ex_put_var (exoid[npd], t+1, EX_GLOBAL, 1, 0, num_global_fields, globals);
	    t_tmp2 = my_timer();
	    raw_write_time += t_tmp2-t_tmp1;
	    if (err) {
	      fprintf(stderr, "after ex_put_global_var, error = %d\n", err);
	      ex_close (exoid[npd]);
	      exit(1);
	    }
	  }
	  for (j=0; j<num_nodal_fields; j++) {
	    t_tmp1 = my_timer();
	    err = ex_put_var (exoid[npd], t+1, EX_NODAL, j+1, 0, num_nodes, x_coords);
	    t_tmp2 = my_timer();
	    raw_write_time += t_tmp2-t_tmp1;
	    if (err) {
	      fprintf(stderr, "after ex_put_nodal_var, error = %d\n", err);
	      ex_close (exoid[npd]);
	      exit(1);
	    }
	  }
	  for (j=0; j<num_element_fields; j++) {
	    t_tmp1 = my_timer();
	    err = ex_put_var (exoid[npd], t+1, EX_ELEM_BLOCK, j+1, 10, num_elems, x_coords);
	    t_tmp2 = my_timer();
	    raw_write_time += t_tmp2-t_tmp1;
	    if (err) {
	      fprintf(stderr, "after ex_put_element_var, error = %d\n", err);
	      ex_close (exoid[npd]);
	      exit(1);
	    }
	  }
	  if (close_files) {
	    t_tmp1 = my_timer();
	    ex_close(exoid[npd]);
	    t_tmp2 = my_timer();
	    raw_open_close_time += (t_tmp2 - t_tmp1);
	  }
	}
      }
      
      if ( rank == 0 ) {
	fprintf(stderr, "\n");
      }
    }
    t_tmp1 = my_timer();
    if (!close_files) {
      if ( rank == 0 ) {fprintf(stderr, "Closing database...");}
      for (npd = 0; npd < files_per_domain; npd++) {
	err = ex_close (exoid[npd]);
	if (err) {
	  printf ("after ex_close, error = %d\n", err);
	  return( 1);
	}
      }
      if ( rank == 0 ) {fprintf(stderr, "\n");}
    }
    t_tmp2 = my_timer();
    raw_open_close_time += (t_tmp2 - t_tmp1);
  
    if ( rank == 0 ) {
      fprintf(stderr, "ex_put_time = %5.2f seconds\n",put_time_time);
    }
    tend = my_timer();
  
    total_time = tend - tstart - raw_sleep_time;
		
    if (total_time > max_total_time) max_total_time = total_time;
    if (total_time < min_total_time) min_total_time = total_time;

    cum_total_time += total_time;

    if (raw_write_time > max_raw_write_time) max_raw_write_time = raw_write_time;
    if (raw_write_time < min_raw_write_time) min_raw_write_time = raw_write_time;

    cum_raw_write_time += raw_write_time;

    if (raw_open_close_time > max_raw_open_close_time) max_raw_open_close_time = raw_open_close_time;
    if (raw_open_close_time < min_raw_open_close_time) min_raw_open_close_time = raw_open_close_time;

    cum_open_close_time += raw_open_close_time;

    
    raw_data_vol += sizeof(realtyp)*(((size_t)num_nodes*num_nodal_fields) +
		     (num_elems * num_element_fields) +
		     (num_global_fields)) * num_timesteps*files_per_domain;
  
  }  /* end of for (iter...) */

#ifdef HAVE_PARALLEL
  MPI_Allreduce(
		&raw_data_vol,
		&glob_raw_data_vol,
		1,
		MPI_LONG,
		MPI_SUM,
               	MPI_COMM_WORLD
		);
#else
  glob_raw_data_vol = raw_data_vol;
#endif
  /*
   *	Get File Sizes
   *
   *	Note: On ASCI Red, a specialized "stat", named "estat", was added to
   *	accommodate file sizes upto 16GB.                          3/27/2002
   */

  if (  stat( tmp_name, &file_status ) ) {
    if ( rank == 0 ) 
      fprintf ( stderr, "Exodus Write: cannot get %s file size.\n", tmp_name );
  
    return( 1 );
  }
  else
    file_size = file_status.st_size * files_per_domain;

#ifdef HAVE_PARALLEL
  MPI_Allreduce (&file_size, &glob_file_size, 1, MPI_LONG, MPI_SUM,
                 MPI_COMM_WORLD);
#else
  glob_file_size = file_size;
#endif
  
  if ( rank == 0 ) {
    fprintf(stderr, "                                                                \n" );
    fprintf(stderr, "   Exodus Write Results                                         \n" );
    fprintf(stderr, "                                                                \n" );
    fprintf(stderr, "      Sizes (bytes)                                             \n" );
    fprintf(stderr, "         File         %14ld                                   \n", (long)glob_file_size);
    fprintf(stderr, "         Raw Data     %14ld                                   \n", (long)glob_raw_data_vol);
    fprintf(stderr, "         Difference   %14ld (%5.2f%%)                         \n",
	    (long)(glob_file_size - glob_raw_data_vol),
	    (1.0*glob_file_size-1.0*glob_raw_data_vol)/(0.01*glob_file_size) );
    fprintf(stderr, "                                                                \n" );
    fprintf(stderr, "      Times                       \t Minimum\t Maximum\t Average\n" );
    fprintf(stderr, "         Raw Data Write (sec)      \t%8.4g\t%8.4g\t%8.4g         \n",
	    min_raw_write_time, max_raw_write_time, cum_raw_write_time / num_iterations );
    fprintf(stderr, "         All Other Write (sec)     \t        \t        \t%8.4g   \n",
	    (cum_total_time - cum_raw_write_time) / num_iterations );
    fprintf(stderr, "         Open/Close Time (sec)     \t        \t        \t%8.4g   \n",
	    cum_open_close_time / num_iterations );
    fprintf(stderr, "         Total Write (sec)         \t        \t        \t%8.4g   \n",
	    cum_total_time / num_iterations );
    fprintf(stderr, "                                                                \n" );
    fprintf(stderr, "      Output Bandwidths            \t Minimum\t Maximum\t Average\n" );
    fprintf(stderr, "         Raw Data Write (MiB/sec)  \t%8.4g\t%8.4g\t%8.4g         \n",
	    (double)glob_raw_data_vol / max_raw_write_time / MBYTES,
	    (double)glob_raw_data_vol / min_raw_write_time / MBYTES,
	    (double)glob_raw_data_vol / cum_raw_write_time / MBYTES * num_iterations);
    fprintf(stderr, "         Raw + Meta Data Write (MiB/sec)        \t        \t%8.4g\n",
	    (double)glob_file_size / cum_total_time / MBYTES * num_iterations );
  }
  free(exoid);
  if (num_global_fields > 0) 
    free(globals);
  return(0);
}

/*****************************************************************************/
void get_file_name(const char* base, const char* ext, int rank, int nprocs,
		   const char* other,char *output)
{
  int i1, iTemp1;
  int iMaxDigit=0, iMyDigit=0;
  char cTemp[128];

  output[0] = '\0';
  strcpy(output, base);
  strcat(output, ".");
  strcat(output, ext);
  if (other != NULL) {
    strcat(output, ".");
    strcat(output, other);
  }

  if (nprocs > 1) {
    /*
     * Find out the number of digits needed to specify the processor ID.
     * This allows numbers like 01-99, i.e., prepending zeros to the
     * name to preserve proper alphabetic sorting of the files.
     */

    iTemp1 = nprocs;
    do {
      iTemp1 /= 10;
      iMaxDigit++;
    }
    while(iTemp1 >= 1);
  
    iTemp1 = rank;
    do {
      iTemp1 /= 10;
      iMyDigit++;
    }
    while(iTemp1 >= 1);

    strcat(output, ".");
    sprintf(cTemp, "%d", nprocs);
    strcat(output, cTemp);
    strcat(output, ".");

    /*
     * Append the proper number of zeros to the filename.
     */
    for(i1=0; i1 < iMaxDigit-iMyDigit; i1++)
      strcat(output, "0");

    sprintf(cTemp, "%d", rank);
    strcat(output, cTemp);
  }
}
