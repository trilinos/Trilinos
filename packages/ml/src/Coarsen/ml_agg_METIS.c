/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create tentative prolongators using METIS to create the      */
/* aggregate.                                                                */
/* ************************************************************************* */
/* Author        : Marzio Sala (SNL)                                         */
/* Date          : October 2003                                              */
/* ************************************************************************* */
/* Local Function :                                                          */
/*    ML_Aggregate_CoarsenMETIS                                              */
/*    ML_DecomposeGraph_with_METIS                                           */
/* ************************************************************************* */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_utils.h"
#include "ml_agg_METIS.h"
#include "ml_viz_stats.h"

/* ********************************************************************** */
/* metis.h is required to properly define idxtype, and to declare the     */
/* function `METIS_PartGraphRecursive' or `METIS_PartGraphKway'.          */
/* By default, idxtype is defined as int, so in principle you can         */
/* compile also without including metis.h. However, to be sure            */
/* that use has specified include and library, I require ML_METIS         */
/* to be defined at compilation time.                                     */
/* NOTE: without metis this function compiles, but does not do any  job.  */
/* ********************************************************************** */

#ifdef HAVE_ML_METIS
#include "metis.h"
#else
#define idxtype int
#endif

static int ML_DecomposeGraph_with_METIS( ML_Operator *Amatrix,
					 int N_nonzeros,
					 int N_parts,
					 int graph_decomposition[],
					 char bdry_nodes[],
					 int local_or_global,
					 int offsets[],
					 int reorder_flag,
					 int current_level );
static int find_max(int length, int vector[] );
static int find_index( int key, int list[], int N );
static int ML_LocalReorder_with_METIS( int Nrows, int xadj[], int adjncy[] ,
					 int Nparts, idxtype part[], int level,
					 ML_Comm *comm );
/* ======================================================================== */
/*!
 \brief Set the number of nodes for each aggregate (for graph-based
 decompositions).

 This function set the number of nodes to be included in each
 aggregate. The use can specify the desired value for a given level, or
 set all levels to the same value (passing level = -1).

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_Set_NodesPerAggr( ML *ml, ML_Aggregate *ag, 
				  int level, int Nnodes_per_aggre  )
{

  int i;
  ML_Aggregate_Options *pointer = NULL;
  int diff_level;
  int Nlevels = ml->ML_num_levels;
  
  /* ********************************************************************** */
  /* control on the input parameters                                        */
  /* ********************************************************************** */

  if ( ag->ML_id != ML_ID_AGGRE ) {
      printf("ML_Aggregate_SetNumberLocal : wrong object. \n");
      exit(-1);
  }

  if( Nnodes_per_aggre <= 0 ) {
    fprintf( stderr,
	     "*ML*ERR* Nlocal has an invalid value (%d)\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     Nnodes_per_aggre,
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  /* ********************************************************************** */
  /* take the pointer from the ag object. If it is NULL, this is the first  */
  /* time that this function is called, so allocate memory for all levels   */
  /* ********************************************************************** */

  pointer = (ML_Aggregate_Options *)ag->aggr_options;
  
  if( pointer == NULL ) {
    ML_memory_alloc((void*)&pointer, sizeof(ML_Aggregate_Options)*Nlevels,
		    "Naggregates");
    if( pointer == NULL ) {
      fprintf( stderr,
	       "*ML*ERR* not enough space to allocate %d bytes\n"
	       "*ML*ERR* (file %s, line %d)\n",
	       sizeof(int)*Nlevels,
	       __FILE__,
	       __LINE__ );
      exit( EXIT_FAILURE );
    }

    /* ******************************************************************** */
    /* set to the default values                                            */
    /* ******************************************************************** */

    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].id = ML_AGGREGATE_OPTIONS_ID;
      pointer[i].Naggregates = -1;
      pointer[i].Nnodes_per_aggregate = -1;
      pointer[i].local_or_global = ML_LOCAL_INDICES;
      pointer[i].reordering_flag = ML_NO;
    }
    
    ag->aggr_options = (void *)pointer;
  }

  if( level >= 0 ) {
    pointer[level].Nnodes_per_aggregate = Nnodes_per_aggre;
    pointer[level].Naggregates = -1;
    pointer[level].local_or_global = ML_LOCAL_INDICES;
  } else {
    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].Nnodes_per_aggregate = Nnodes_per_aggre;
      pointer[i].Naggregates = -1;
      pointer[i].local_or_global = ML_LOCAL_INDICES;
    }
  }
  
  return 0;
  
} /* ML_Aggregate_Set_NodesPerAggr */
  
/* ======================================================================== */
/*!
 \brief stored the required number of local aggregates for the specified level

 This function is used to specify the local number of aggregates which
 will be created by \c ML_Aggregate_CoarsenMETIS. The input variable
 \c level should reflect the structure of the multilevel solver. 
 
 Parameter list:
 - ml : ML object
 - ag : ML_Aggregate object, where the number of local aggregates
 - level :
 - Nlocal : this is the number of partition into which METIS will decompose
            the local graph of the matrix (that is, locally on each process,
	    ignoring the intra-processes pattern)
*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_Set_LocalNumber( ML *ml, ML_Aggregate *ag, 
				  int level, int Nlocal  )
{

  int i;
  ML_Aggregate_Options *pointer = NULL;
  int diff_level;
  int Nlevels = ml->ML_num_levels;

  /* ********************************************************************** */
  /* control on the input parameters                                        */
  /* ********************************************************************** */

  if ( ag->ML_id != ML_ID_AGGRE ) {
      printf("ML_Aggregate_SetNumberLocal : wrong object. \n");
      exit(-1);
  }

  if( Nlocal <= 0 ) {
    fprintf( stderr,
	     "*ML*ERR* Nlocal has an invalid value (%d)\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     Nlocal,
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  /* ********************************************************************** */
  /* take the pointer from the ag object. If it is NULL, this is the first  */
  /* time that this function is called, so allocate memory for all levels   */
  /* ********************************************************************** */

  pointer = (ML_Aggregate_Options *)ag->aggr_options;
  
  if( pointer == NULL ) {
    ML_memory_alloc((void*)&pointer, sizeof(ML_Aggregate_Options)*Nlevels,
		    "Naggregates");
    if( pointer == NULL ) {
      fprintf( stderr,
	       "*ML*ERR* not enough space to allocate %d bytes\n"
	       "*ML*ERR* (file %s, line %d)\n",
	       sizeof(int)*Nlevels,
	       __FILE__,
	       __LINE__ );
      exit( EXIT_FAILURE );
    }

    /* ******************************************************************** */
    /* set to the default values                                            */
    /* ******************************************************************** */

    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].id = ML_AGGREGATE_OPTIONS_ID;
      pointer[i].Naggregates = -1;
      pointer[i].Nnodes_per_aggregate = -1;
      pointer[i].local_or_global = ML_LOCAL_INDICES;
      pointer[i].reordering_flag = ML_NO;
    }
    
    ag->aggr_options = (void *)pointer;
  }

  if( level >= 0 ) {
    pointer[level].Naggregates = Nlocal;
    pointer[i].Nnodes_per_aggregate = -1;
    pointer[level].local_or_global = ML_LOCAL_INDICES;
  } else {
    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].Naggregates = Nlocal;
      pointer[i].Nnodes_per_aggregate = -1;
      pointer[i].local_or_global = ML_LOCAL_INDICES;
    }
  }
  
  return 0;
  
} /* ML_Aggregate_Set_LocalNumber */

int ML_Aggregate_Set_GlobalNumber( ML *ml, ML_Aggregate *ag, 
				   int level, int Nglobal  )
{

  int i;
  ML_Aggregate_Options *pointer = NULL;
  int diff_level;
  int Nlevels = ml->ML_num_levels;

  /* ********************************************************************** */
  /* control on the input parameters                                        */
  /* ********************************************************************** */

  if ( ag->ML_id != ML_ID_AGGRE ) {
      printf("ML_Aggregate_SetGlobalNumber : wrong object. \n");
      exit(-1);
  }

  if( Nglobal <= 0 ) {
    fprintf( stderr,
	     "*ML*ERR* Nlocal has an invalid value (%d)\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     Nglobal,
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  /* ********************************************************************** */
  /* take the pointer from the ag object. If it is NULL, this is the first  */
  /* time that this function is called, so allocate memory for all levels   */
  /* ********************************************************************** */

  pointer = (ML_Aggregate_Options *)ag->aggr_options;
  
  if( pointer == NULL ) {
    ML_memory_alloc((void*)&pointer, sizeof(ML_Aggregate_Options)*Nlevels,
		    "aggr_options");
    if( pointer == NULL ) {
      fprintf( stderr,
	       "*ML*ERR* not enough space to allocate %d bytes\n"
	       "*ML*ERR* (file %s, line %d)\n",
	       sizeof(int)*Nlevels,
	       __FILE__,
	       __LINE__ );
      exit( EXIT_FAILURE );
    }

    /* ******************************************************************** */
    /* set to the default values                                            */
    /* ******************************************************************** */
    
    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].id = ML_AGGREGATE_OPTIONS_ID;
      pointer[i].Naggregates = -1;
      pointer[i].Nnodes_per_aggregate = -1;
      pointer[i].local_or_global = ML_LOCAL_INDICES;
      pointer[i].reordering_flag = ML_NO;
    }
    
    ag->aggr_options = (void *)pointer;
  }
  
  if( level >= 0 ) {
    pointer[level].Naggregates = Nglobal;
    pointer[level].local_or_global = ML_GLOBAL_INDICES;
  } else {
    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].Naggregates = Nglobal;
      pointer[i].local_or_global = ML_GLOBAL_INDICES;
    }
  }
   
  return 0;
  
} /* ML_Aggregate_SetGlobalNumber */

/* ======================================================================== */
/*!
 \brief Set the reordering flag for METIS decompositions.

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_Set_ReorderingFlag( ML *ml, ML_Aggregate *ag, 
				     int level, int reordering_flag  )
{

  int i;
  ML_Aggregate_Options *pointer = NULL;
  int diff_level;
  int Nlevels = ml->ML_num_levels;
  
  /* ********************************************************************** */
  /* control on the input parameters                                        */
  /* ********************************************************************** */

  if ( ag->ML_id != ML_ID_AGGRE ) {
      printf("ML_Aggregate_SetNumberLocal : wrong object. \n");
      exit(-1);
  }

  if( (reordering_flag != ML_YES) && (reordering_flag != ML_NO) ) {
    fprintf( stderr,
	     "*ML*ERR* reordering_flag has a wrong value (%d)\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     reordering_flag,
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  /* ********************************************************************** */
  /* take the pointer from the ag object. If it is NULL, this is the first  */
  /* time that this function is called, so allocate memory for all levels   */
  /* ********************************************************************** */

  pointer = (ML_Aggregate_Options *)ag->aggr_options;
  
  if( pointer == NULL ) {
    ML_memory_alloc((void*)&pointer, sizeof(ML_Aggregate_Options)*Nlevels,
		    "Naggregates");
    if( pointer == NULL ) {
      fprintf( stderr,
	       "*ML*ERR* not enough space to allocate %d bytes\n"
	       "*ML*ERR* (file %s, line %d)\n",
	       sizeof(int)*Nlevels,
	       __FILE__,
	       __LINE__ );
      exit( EXIT_FAILURE );
    }

    /* ******************************************************************** */
    /* set to the default values                                            */
    /* ******************************************************************** */

    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].id = ML_AGGREGATE_OPTIONS_ID;
      pointer[i].Naggregates = -1;
      pointer[i].Nnodes_per_aggregate = -1;
      pointer[i].local_or_global = ML_LOCAL_INDICES;
      pointer[i].reordering_flag = ML_NO;
    }
    
    ag->aggr_options = (void *)pointer;
  }

  if( level >= 0 ) {
    pointer[level].reordering_flag = reordering_flag;
  } else {
    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].reordering_flag = reordering_flag;
    }
  }

  return 0;
  
} /* ML_Aggregate_Set_ReorderingFlag */

/* ======================================================================== */
/*!
 \brief Reorder the local graph for the coarser level matrix using a METIS
 function

 This function builds the graph of the coarser level matrix (without
 filling it with numerical values), then call METIS_NodeND to compute a
 reordering which reduces the fill-in during LU factorizations. It is
 just a LOCAL reordering.
 
*/
/* ------------------------------------------------------------------------ */

static int ML_LocalReorder_with_METIS( int Nrows, int xadj[], int adjncy[] ,
				       int Nparts, idxtype part[], int level,
				       ML_Comm *comm )
{

  int i, j, k, count;
  int row_coarsest, col_coarsest, col;
  int MaxNnzRow;
  idxtype * xadj2 = NULL, * adjncy2 = NULL;
  idxtype * perm = NULL, * iperm = NULL;
  idxtype options[8];
  size_t used_mem;
  int bandwidth_orig, bandwidth_perm;
  int mypid = comm->ML_mypid;
  double t0;
  int nbytes, nbytes_max;
  ML_Aggregate_Options aggr_options;
#ifdef METIS_DEBUG
  FILE *fp;
  char filename[80];
#endif
  
  /* ------------------- execution begins --------------------------------- */

  t0 = GetClock();

  if ( mypid == 0 &&  ML_Get_PrintLevel() > 8 ) {
    
    printf("Entering METIS reordering (level %d)\n",
	   level );
    
  }
  
  if( Nparts < 3 ) return; /* do nothing if 2 nodes only... */
  
  /* ********************************************************************** */
  /* In order to compute the fill-in reducin ordering, I need to form the   */
  /* graph of the matrix on the next level. This matrix is still not formed */
  /* so I have to construct it. It is a LOCAL matrix, and numerical entries */
  /* are not relevant. I need some memory for this. As I don't know exactly */
  /* the number of nonzeros in the coarsert matrix, I use the one of the    */
  /* finest level (cannot be more)                                          */
  /* ********************************************************************** */

  MaxNnzRow = 0;
  for( i=0 ; i<Nrows ; i++ ) {
    if( (xadj[i+1] - xadj[i])> MaxNnzRow )
      MaxNnzRow = xadj[i+1] - xadj[i];
  }

  MaxNnzRow *= MaxNnzRow;
  
  xadj2 = (idxtype *)malloc( sizeof(idxtype) * (Nparts+1) );
  adjncy2 = (idxtype *)malloc( sizeof(idxtype) * MaxNnzRow * Nparts );
  
  if( xadj2 == NULL || adjncy2 == NULL ) {
    fprintf( stderr,
	     "*ML*ERR* not nough memory to allocated %d and %d bytex\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     sizeof(idxtype) * Nparts,
	     sizeof(idxtype) * MaxNnzRow,
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }

  /* ********************************************************************** */
  /* I cycle over the rows of the finest level. For each nonzero (an edge   */
  /* for the METIS format), I look for the corresponding aggregate for form */
  /* the local row, added to the coarset level matrix. This matrix will be  */
  /* passed to the METIS function to compute orderings, and the new ordering*/
  /* will be used to enumerate the local aggregates.                        */
  /* ********************************************************************** */

  for( i=0 ; i<Nparts+1 ; i++ ) xadj2[i] = 0;
  for( i=0 ; i<Nparts*MaxNnzRow ; i++ ) adjncy2[i] = -1;
  
  /* cycle over the rows of finest lever */

  for( i=0 ; i<Nrows ; i++ ) {

    row_coarsest = part[i];
    
    /* cycle over the columns of the finest level */
    for( j=xadj[i] ; j<xadj[i+1] ; j++ ) {
      
      col = adjncy[j];

      col_coarsest = part[col];

      /* discard element diagonal */
      if( col_coarsest != row_coarsest ) {
      
	for( k= 0 ; k<xadj2[row_coarsest+1] ; k++ ) {

	  if( adjncy2[row_coarsest*MaxNnzRow+k] == col_coarsest ) {
	    k = -1;
	    break;
	  }
	  if( adjncy2[row_coarsest*MaxNnzRow+k] == -1 ) break;

	}
	if( k == -1 ) continue;
	else if( k == MaxNnzRow ) {
	  fprintf( stderr,
		   "*ML*ERR* something went wrong: k=%d, MaxNnzRow=%d\n"
		   "*ML*ERR* (file %s, line %d)\n",
		   k, MaxNnzRow,
		   __FILE__,
		   __LINE__ ) ;
	  exit( EXIT_FAILURE );
	}
	adjncy2[row_coarsest*MaxNnzRow+k] = col_coarsest;
	xadj2[row_coarsest+1]++;
      }
      
    }
    
  }

  /* compress the crs format */

  for( i=1 ; i<Nparts+1 ; i++ ) xadj2[i] += xadj2[i-1];

  /* count is the position of the first available zero position */
  
  count = 0;
  
  while( adjncy2[count] != -1 && count < Nparts*MaxNnzRow ) {
    count++;
  }

  for( i=count+1 ; i<Nparts*MaxNnzRow ; i++ ) {
    /* find the next nonzero position */
    if( adjncy2[i] == -1 ) continue;
    adjncy2[count] = adjncy2[i];
    adjncy2[i] = -1;
    
    while( adjncy2[count] != -1 ) {
      count++;
    }
  }

#ifdef METIS_DEBUG
  sprintf( filename,
	   "before_reordering_level%d_proc%d.m",
	   level,
	   mypid );

  fp = fopen( filename, "w" );

  fprintf( fp,
	   "A = zeros(%d,%d);\n",
	   Nparts,Nparts );
  for( i=0 ; i<Nparts ; i++ ) {
    fprintf( fp,
	     "A(%d,%d) = 1;\n",
	     i+1, i+1 );
    for( j=xadj2[i] ; j<xadj2[i+1] ; j++ ) {
      fprintf( fp,
	       "A(%d,%d) = 1;\n",
	       i+1, adjncy2[j]+1 );
    }
  }

  fclose( fp );
#endif
  
  used_mem = sizeof(idxtype) * xadj2[Nparts];
  
  if( count != xadj2[Nparts] ) {
    fprintf( stderr,
	     "*ML*ERR* something went wrong: count=%d, xadj2[Nparts]=%d\n"
	     "*ML*ERR* (file %s line %d)\n",
	     count,
	     xadj2[Nparts],
	     __FILE__,
	     __LINE__ ) ;
    exit( EXIT_FAILURE );
  }
  
  adjncy2 = (idxtype *) realloc( adjncy2, used_mem );

  perm = (idxtype *) malloc( sizeof(idxtype) * Nparts );
  iperm = (idxtype *) malloc( sizeof(idxtype) * Nparts );
  
  if( perm == NULL || iperm == NULL ) {
    fprintf( stderr,
	     "*ML*ERR* not nough memory to allocated %d bytex\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     sizeof(idxtype) * Nparts,
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }

  options[0] = 0; /* default values */
  /* those are instead possible non-default values. No idea about how
     to pick up them...
  options[0] = 1;
  options[1] = 2; 
  options[2] = 2;
  options[3] = 1;
  options[4] = 0;
  options[5] = 1;
  options[6] = 0;
  options[7] = 1; 
  */

  /* ********************************************************************** */
  /* estimate memory required by reordering                                 */
  /* ********************************************************************** */

  j = 0;
  i = 4;
#ifdef HAVE_ML_METIS
  METIS_EstimateMemory( &Nrows, xadj, adjncy, &j,
			&i, &nbytes );
#endif
  nbytes /= 1024;
  
  nbytes_max = ML_gmax_int( nbytes, comm );
  nbytes = ML_gsum_int( nbytes, comm);
  
  if( mypid == 0 &&  ML_Get_PrintLevel() > 8 ) {
   
    printf("METIS reordering (level %d) estimated required mem = %d Kb\n"
	   "METIS reordering (level %d) max estimated mem = %d Kb\n",
	   level,
	   nbytes,
	   level,
	   nbytes_max );
  }
  
  /* ********************************************************************** */
  /* the METIS manual conseils to use NodeNS and not EdgeND. I have tested  */
  /* EdgeNS (in one scalar case), and it does not seem too different from   */
  /* NodeND. No idea about CPU times.                                       */
  /* iperm is a vector defined as: row i of unreordered matrix is row       */
  /* iperm[i] of reordered matrix                                           */
  /* ********************************************************************** */

  i = 0; /* offset, C or FORTRAN */
#ifdef HAVE_ML_METIS
  METIS_NodeND( &Nparts, xadj2, adjncy2, &i, options, perm, iperm );
#else
  fprintf( stderr,
	   "*ERR*ML* This function has been compiled without -DHAVE_ML_METIS\n"
	   "*ERR*ML* To use METIS, please add --with-ml_metis to your\n"
	   "*ERR*ML* configure script (recall to specify include dir and\n"
	   "*ERR*ML* location of the METIS lib; see configure --help\n"
	   "*ERR*ML* for more defailts).\n"
	   "*ERR*ML* (file %s, line %d)\n",
	   __FILE__,
	   __LINE__);
  exit( EXIT_FAILURE );  
#endif

  /* ********************************************************************** */
  /* replace previous ordering in part with this fill-in reducing one       */
  /* Also, compute the bandwidth before and after reordering                */
  /* It seesm to me that the bandwidth is INCREASED after reordering.       */
  /* However, the number of zero elements (as given by OCTAVE's LU) is      */
  /* increased --- In the small test I run, from 88% to 93%.                */
  /* ********************************************************************** */

  bandwidth_orig = 0;
  bandwidth_perm = 0;
  for( i=0 ; i<Nparts ; i++ ) {
    for( j=xadj2[i] ; j<xadj2[i+1] ; j++ ) {
      k = abs(adjncy2[j] - i);
      if( k>bandwidth_orig ) bandwidth_orig = k;
      k = abs(iperm[adjncy2[j]] - iperm[i]);
      if( k>bandwidth_perm ) bandwidth_perm = k;
    }
  }

  if ( mypid == 0 &&  ML_Get_PrintLevel() > 8 ) {
   
    printf("METIS reordering (level %d) Bandwidth before = %d, after = %d\n",
	   level,
	   bandwidth_orig, bandwidth_perm);
    
  }
  
  for( i=0 ; i<Nrows ; i++ ) {
    j = part[i];
    part[i] = iperm[j];
  }
  
#ifdef METIS_DEBUG
  sprintf( filename,
	   "after_reordering_level%d_proc%d.m",
	   level,
	   mypid );
  
  fp = fopen( filename, "w" );

  fprintf( fp,
	   "A = zeros(%d,%d);\n",
	   Nparts,Nparts );
  for( i=0 ; i<Nparts ; i++ ) {
    fprintf( fp,
	     "A(%d,%d) = 1;\n",
	     iperm[i]+1, iperm[i]+1 );
    for( j=xadj2[i] ; j<xadj2[i+1] ; j++ ) {
      fprintf( fp,
	       "A(%d,%d) = 1;\n",
	       iperm[i]+1, iperm[adjncy2[j]]+1 );
    }
  }

  fclose( fp );
#endif

  /* ------------------- that's all folks --------------------------------- */

  free( xadj2 ) ;
  free( adjncy2 );
  free( perm );
  free( iperm );

  t0 = GetClock() - t0;

  if ( mypid == 0 &&  ML_Get_PrintLevel() > 8 ) {
   
    printf("METIS reordering (level %d) Time required = %e\n",
	   level,
	   t0 );
    
  }
  
  return 0;
  
} /* ML_LocalReorder_with_METIS */

/* ************************************************************************* */
/* This function calls METIS to decompose the graph of the local matrix      */
/* (that is, ignoring any inter-domain connections)                          */
/* ************************************************************************* */

static int ML_DecomposeGraph_with_METIS( ML_Operator *Amatrix,
					 int N_nonzeros,
					 int N_parts,
					 int graph_decomposition[],
					 char bdry_nodes[],
					 int local_or_global,
					 int offsets[],
					 int reorder_flag,
					 int current_level )
{

  int i, j,jj,  count;
  int Nrows, Nghosts;
  int nnz, *wgtflag=NULL, numflag, *options=NULL, edgecut;
  idxtype *xadj=NULL, *adjncy=NULL, *vwgt=NULL, *adjwgt=NULL;
  idxtype *part=NULL;
  ML_Comm * comm;
  int allocated = 0;
  int * rowi_col = NULL;
  int rowi_N;
  double * rowi_val = NULL;
  int nbytes = 0, nbytes_max = 0;
  
  /* ------------------- execution begins --------------------------------- */

  comm = Amatrix->comm;
  
  /* dimension of the problem (NOTE: only local matrices) */

  Nrows = Amatrix->getrow->Nrows;
  Nghosts = Amatrix->getrow->pre_comm->total_rcv_length;
  
  /* construct the CSR graph information of the LOCAL matrix
     using the get_row function */

  wgtflag = (idxtype *) malloc (4*sizeof(idxtype));
  options = (int *)     malloc (4*sizeof(int));
  
  /* set parameters */
   
  wgtflag[0] = 0;    /* no weights */
  numflag    = 0;    /* C style */
  options[0] = 0;    /* default options */
   
  xadj    = (idxtype *) malloc ((Nrows+1)*sizeof(idxtype));
  adjncy  = (idxtype *) malloc (N_nonzeros  *sizeof(idxtype));
   
  if(  xadj==NULL || adjncy==NULL ) {
    fprintf( stderr,
	     "on proc %d, not enought space for %d bytes.\n"
	     "file %s, line %d\n",
	     comm->ML_mypid, N_nonzeros,
	     __FILE__,
	     __LINE__);
  }
   
  xadj[0] = 0;
  count = 0;
   
  for (i = 0; i < Nrows; i++) {

    if( bdry_nodes[i] == 'T' ) {
      xadj[i+1] = xadj[i];
      continue;
    }
    
    /* retrive one row */
     
    ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
		      &rowi_N, 0);
/*
      printf("row %d: ", i);
      for( j=0 ; j<rowi_N ; j++ ) {
      printf("%d ", rowi_col[j]);
      }
      puts("");
*/

    xadj[i+1] = xadj[i];
    for( j=0 ; j<rowi_N ; j++ ) {
      jj = rowi_col[j];
      if( jj<Nrows ) {
	adjncy[count++] = rowi_col[j];
	xadj[i+1]++;
      }
       
    }      
  }
  
  if( count > N_nonzeros ) {
    fprintf( stderr,
	     "Warning: on proc %d, count > N_nonzeros (%d>%d)\n"
	     "a buffer overflow has probably occurred...\n",
	     comm->ML_mypid, count, N_nonzeros );
  }

#ifdef METIS_DEBUG
  for( i=0 ; i<Nrows ; i++ )
    printf("%d %d\n", xadj[i], xadj[i+1] );
#endif
  
  /* idxtype is by default int, but on some architectures can be
     slightly different (for instance, a short int). */
   
  part = (idxtype *) malloc( sizeof(idxtype) * Nrows);

  /* ********************************************************************** */
  /* Before calling METIS, I verify that the two extreme situations are     */
  /* handled separately.                                                    */
  /* ********************************************************************** */
   
  if( N_parts == 1 ) {

    for( i=0 ; i<Nrows ; i++ )
      part[i] = 0;
    edgecut = 0;
    
  } else if( N_parts == Nrows ) {

    fprintf( stderr,
	     "*ML*WRN*: on proc %d, N_part == N_nonzeros (%d==%d)\n",
	     comm->ML_mypid, N_parts, Nrows );
 
    for( i=0 ; i<Nrows ; i++ )
      part[i] = i;
    edgecut = 0;
  
  } else {

    /* ******************************************************************** */
    /* Put -1 in part, so that I can verify that METIS has filled each pos  */
    /* ******************************************************************** */

    for( i=0 ; i<Nrows ; i++ ) part[i] = -1;
    
    /* ******************************************************************** */
    /* Estimate memory required by METIS. This memory will be dynamically   */
    /* allocated inside; however this is a good estimation of how METIS     */
    /* will cost in terms of memory.                                        */
    /* Then, call METIS.                                                    */
    /* ******************************************************************** */

    if( N_parts < 8 ) {
     
      i = 1; /* optype in the METIS manual */
      numflag = 0;
#ifdef HAVE_ML_METIS
      METIS_EstimateMemory( &Nrows, xadj, adjncy, &numflag,
			    &i, &nbytes );
			   
      METIS_PartGraphRecursive (&Nrows, xadj, adjncy, vwgt, adjwgt,
				wgtflag, &numflag, &N_parts, options,
				&edgecut, part);
#else
      fprintf( stderr,
	       "*ML*ERR* Compile with metis...\n");
      exit( EXIT_FAILURE );
#endif
    } else {
     
      i = 2;
      numflag = 0;
#ifdef HAVE_ML_METIS
      METIS_EstimateMemory( &Nrows, xadj, adjncy, &numflag,
			    &i, &nbytes );

      METIS_PartGraphKway (&Nrows, xadj, adjncy, vwgt, adjwgt,
			   wgtflag, &numflag, &N_parts, options,
			   &edgecut, part);
#else
    if( Amatrix->comm->ML_mypid == 0 ) {
      fprintf( stderr,
	       "*ERR*ML* This function has been compiled without -DHAVE_ML_METIS\n"
	       "*ERR*ML* To use METIS, please add --with-ml_metis to your\n"
	       "*ERR*ML* configure script (recall to specify include dir and\n"
	       "*ERR*ML* location of the METIS lib; see configure --help\n"
	       "*ERR*ML* for more defailts).\n"
	       "*ERR*ML* (file %s, line %d)\n",
	       __FILE__,
	       __LINE__);
    }
    exit( EXIT_FAILURE );
#endif
    }
    /* ******************************************************************** */
    /* I verify that METIS did a good job. Problems may arise when the # of */
    /* aggregates is too high. In this case, we may end up with some        */
    /* aggregates which contain no nodes. In this case, I introduce at      */
    /* least one node, in a random fashion.                                 */
    /* ******************************************************************** */

    ML_Aggregates_CheckAggregates( N_parts, Nrows, part,
				   Amatrix->comm->ML_mypid );
  }

  /* ********************************************************************** */
  /* Some fancy output for memory usage.                                    */
  /* ********************************************************************** */

  nbytes /= 1024;
  
  nbytes_max = ML_gmax_int( nbytes, comm );
  nbytes = ML_gsum_int( nbytes, comm);

  if( Amatrix->comm->ML_mypid == 0 &&  ML_Get_PrintLevel() > 8 ) {
   
    printf("METIS aggregation (level %d) estimated required mem = %d Kb\n"
	   "METIS aggregation (level %d) max estimated mem = %d Kb\n",
	   current_level,
	   nbytes,
	   current_level,
	   nbytes_max );
  }
  
  /* ********************************************************************** */
  /* reordering using METIS to minimize the fill-in during factorization    */
  /* ********************************************************************** */

  if( reorder_flag == ML_YES ) {
    
    ML_LocalReorder_with_METIS( Nrows, xadj, adjncy ,
				N_parts,  part, current_level, comm );
    
  }
  
  /* copy back part into aggr_index, and set to -1
     the aggr_index corresponding to ghost nodes */

  for( i=0 ; i<Nrows ; i++ )  graph_decomposition[i] = (int)part[i];
  for( i = Nrows; i < Nrows+Nghosts ; i++ )  graph_decomposition[i] = -1;

  /* if global indices are required, modify the entries
     of graph_decomposition (only the LOCAL entries) so that
     they correspond to global indices. Also, set the array
     offsets, defined so that the global indices assigned
     to processor i are
     offsets[i] <= indices_of_i < offsets[i+1]
     Note that I do not suppose that N_parts is the same
     value among all the processors */
     
  if( local_or_global == ML_GLOBAL_INDICES ) {
    ML_DecomposeGraph_BuildOffsets( N_parts, offsets, comm->ML_nprocs);
  }
  
  /* ------------------- that's all folks --------------------------------- */

  ML_free(rowi_col); ML_free(rowi_val);
  rowi_col = NULL; rowi_val = NULL;
  allocated = 0; 

  if( options != NULL ) free( (void *)options );
  if( wgtflag != NULL ) free( (void *)wgtflag );
  if( adjncy != NULL  ) free( (void *)adjncy  );
  if( xadj != NULL    ) free( (void *)xadj    );
  if( part != NULL    ) free( (void *)part    );
   
  /* some checks on aggr_index: verify that we have the
     required number of aggregates */

  return 0;
  
} /* ML_DecomposeGraph_with_METIS */

/* ======================================================================== */
/*!
 \brief create non-smoothed aggregates using METIS. In order to use this
 function, the user has to define the number of aggregate (or the # of
 nodes in each aggregate) using the functions \c ML_Aggregate_Set_LocalNumber
 or \c ML_Aggregate_Set

 \note this function is derived from ML_Aggregate_CoarsenMIS

*/
/* ------------------------------------------------------------------------ */


int ML_Aggregate_CoarsenMETIS( ML_Aggregate *ml_ag, ML_Operator *Amatrix, 
			       ML_Operator **Pmatrix, ML_Comm *comm)
{
  unsigned int nbytes, length;
   int     i, j, jj, k, m, Nrows, exp_Nrows;
   int     N_neighbors, *neighbors = NULL, diff_level;
   double  printflag;
   int     *Asqrd_rcvleng= NULL, *Asqrd_sndleng= NULL, *send_list = NULL;
   int     total_recv_leng = 0, total_send_leng = 0, offset, msgtype;
   int     aggr_count, index, mypid, num_PDE_eqns;
   int     *aggr_index = NULL, *itmp_array = NULL, nullspace_dim;
   int     *sendlist_proc = NULL, Ncoarse, count, *int_buf = NULL;
   int     *int_buf2, *aggr_stat = NULL, procnum;
   int     index4, *new_send_leng = NULL, new_N_send;
   int     *new_send_neighbors = NULL, *new_send_list = NULL;
   int     max_count, *new_ia = NULL, *new_ja = NULL, new_Nrows;
   int     *new_recv_leng=NULL, exp_Ncoarse, new_N_recv;
   int     *new_recv_neighbors=NULL, *aggr_cnt_array = NULL;
   int     level, index3, count3, *recv_list = NULL, max_agg_size;
   int     **rows_in_aggs = NULL, lwork, info;
   double  dcompare1, *new_val = NULL, epsilon;
   double  *dble_buf = NULL, *nullspace_vect = NULL, *qr_tmp = NULL;
   double  *tmp_vect = NULL, *work = NULL, *new_null = NULL, *comm_val = NULL;
   double  *dble_buf2;
   ML_SuperNode          *aggr_head = NULL, *aggr_curr, *supernode;
   struct ML_CSR_MSRdata *csr_data;
   ML_Aggregate_Comm     *aggr_comm;
   ML_GetrowFunc         *getrow_obj;
   USR_REQ               *request=NULL;
   ML_Operator           *Asqrd = NULL, *tmatrix;
   struct ML_CSR_MSRdata *temp;
   char                  *vtype, *state, *bdry, *unamalg_bdry;
   int                   nvertices, *vlist, Asqrd_ntotal, Asqrd_Nneigh;
   int                   *Asqrd_neigh = NULL, max_element, Nghost;
   int                   **Asqrd_sndbuf = NULL, **Asqrd_rcvbuf = NULL;
   int                   **Asqrd_sendlist = NULL, **Asqrd_recvlist = NULL;
   ML_CommInfoOP         *mat_comm;
   int                   allocated = 0, *rowi_col = NULL, rowi_N, current;
   double                *rowi_val = NULL, *dtemp;
   int                   *templist, **proclist, *temp_index;
   int                   *temp_leng, *tem2_index, *tempaggr_index = NULL;
   int                   *send_leng = NULL, *recv_leng = NULL;
   int                   total_nz = 0;
   int                   count2;
   ML_agg_indx_comm      agg_indx_comm;
   int Nnonzeros = 0;

   int reorder_flag;
   /*   int kk, old_upper, nnzs, count2, newptr; */
#ifdef CLEAN_DEBUG
   char tlabel[80];
#endif

#ifdef DDEBUG
   int curagg,myagg,*good,*bad, kk;
#endif

#if defined(OUTPUT_AGGREGATES) || defined(DDEBUG) || defined(INPUT_AGGREGATES) || (ML_AGGR_INAGGR) || (ML_AGGR_OUTAGGR) || (ML_AGGR_MARKINAGGR)
FILE *fp;
char fname[80];
static int level_count = 0;
double *d2temp;
int agg_offset, vertex_offset;
#endif
 int * graph_decomposition = NULL;
 ML_Aggregate_Viz_Stats * aggr_viz_and_stats;
 ML_Aggregate_Options * aggr_options;
 int mod, Nprocs;
 
   /* ============================================================= */
   /* get the machine information and matrix references             */
   /* ============================================================= */

   mypid                   = comm->ML_mypid;
   Nprocs                  = comm->ML_nprocs;
   epsilon                 = ml_ag->threshold;
   num_PDE_eqns            = ml_ag->num_PDE_eqns;
   nullspace_dim           = ml_ag->nullspace_dim;
   nullspace_vect          = ml_ag->nullspace_vect;
   Nrows                   = Amatrix->outvec_leng;
   printflag               = ml_ag->print_flag;

   /* ============================================================= */
   /* check the system size versus null dimension size              */
   /* ============================================================= */

   if ( Nrows % num_PDE_eqns != 0 )
   {
      printf("ML_Aggregate_CoarsenMIS ERROR : Nrows must be multiples");
      printf(" of num_PDE_eqns.\n");
      exit(EXIT_FAILURE);
   }
   diff_level = ml_ag->max_levels - ml_ag->cur_level - 1;
   if ( diff_level > 0 ) num_PDE_eqns = nullspace_dim; /* ## 12/20/99 */

   /* ============================================================= */
   /* Figure out where the Dirichlet points are on the fine grid of */ 
   /* the unamalgmated system.                                      */
   /* ============================================================= */

   /* GOAL: construct the char vector `unamalg_bdry', whose components
      are 'T' is the node is a Dirichlet node, or 'F' if not.
      The size of this vector is Nrows+Nghost, because I need the
      ghost node also to construct the tentative prolongator */
   /* also, I use this part to retrive an estimation of the number
      of nonzeros (I need it to allocate memory for METIS later on */

   Nnonzeros = 0;
   
   Nghost = Amatrix->getrow->pre_comm->total_rcv_length;
   unamalg_bdry = (char *) ML_allocate(sizeof(char)*(Nrows+Nghost+1));
   dtemp = (double *) ML_allocate(sizeof(double)*(Nrows+Nghost+1));
   if (dtemp == NULL) pr_error("ml_agg_METIS: out of space.\n");

   for (i = 0; i < Nrows+Nghost; i++) dtemp[i] = 0.;

   for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      count2 = 0;
      for (j = 0; j < rowi_N; j++) if (rowi_val[j] != 0.) count2++;
      if (count2 <= 1) dtemp[i] = 1.;
      Nnonzeros += rowi_N; /* used by METIS */
   }
   ML_free(rowi_col); ML_free(rowi_val);
   rowi_col = NULL; rowi_val = NULL;
   allocated = 0; 

   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm,Amatrix->outvec_leng,
                    comm, ML_OVERWRITE,NULL);
   for (i = 0; i < Nrows+Nghost; i++) {
      if (dtemp[i] == 1.) unamalg_bdry[i] = 'T';
      else unamalg_bdry[i] = 'F';
   }
   ML_free(dtemp);

   /* ============================================================= */
   /* set up the threshold for weight-based coarsening              */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   ml_ag->curr_threshold *= 0.5;

   if ( mypid == 0 && printflag < ML_Get_PrintLevel())
   {
      printf("ML_Aggregate_CoarsenMETIS : current level = %d\n",
                                            ml_ag->cur_level);
      printf("ML_Aggregate_CoarsenMETIS : current eps = %e\n",epsilon);
   }
/*
   epsilon = epsilon * epsilon;
*/

#ifdef CLEAN_DEBUG
   sprintf(tlabel,"before amalg %d",comm->ML_mypid);
   ML_CommInfoOP_Print(Amatrix->getrow->pre_comm, tlabel);
#endif
   ML_Operator_AmalgamateAndDropWeak(Amatrix, num_PDE_eqns, epsilon);
   Nrows /= num_PDE_eqns;
#ifdef CLEAN_DEBUG
   sprintf(tlabel,"after amalg %d",comm->ML_mypid);
   ML_CommInfoOP_Print(Amatrix->getrow->pre_comm, tlabel);
#endif

   nvertices = Amatrix->outvec_leng;

   /* Need to set up communication for the matrix A */

   getrow_obj   = Amatrix->getrow;
   N_neighbors  = getrow_obj->pre_comm->N_neighbors;
   total_recv_leng = 0;
   for (i = 0; i < N_neighbors; i++) {
      total_recv_leng += getrow_obj->pre_comm->neighbors[i].N_rcv;
   }
   recv_list   = (int *) ML_allocate(sizeof(int *)*total_recv_leng);
   max_element = nvertices - 1;
   count = 0;
   for (i = 0; i < N_neighbors; i++) {
      for (j = 0; j < getrow_obj->pre_comm->neighbors[i].N_rcv; j++) {
         recv_list[count] = getrow_obj->pre_comm->neighbors[i].rcv_list[j];
         if (recv_list[count] > max_element ) max_element = recv_list[count];
         count++;
      }
   }
   Nghost = max_element - nvertices + 1;
   exp_Nrows = nvertices + Nghost;

   /* record the Dirichlet boundary */

   bdry = (char *) ML_allocate(sizeof(char)*(exp_Nrows + 1));
   for (i = Nrows ; i < exp_Nrows; i++) bdry[i] = 'F';
   for (i = 0; i < Nrows; i++) {
      bdry[i] = 'T';
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      if (rowi_N > 1) bdry[i] = 'F';
   }

   /* communicate the boundary information */

   dtemp = (double *) ML_allocate(sizeof(double)*(exp_Nrows+1));
   for (i = nvertices; i < exp_Nrows; i++) dtemp[i] = 0;
   for (i = 0; i < nvertices; i++) {
      if (bdry[i] == 'T') dtemp[i] = 1.;
      else  dtemp[i] = 0.;
   }
   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm,nvertices,comm,
                    ML_OVERWRITE,NULL);
   for (i = nvertices; i < exp_Nrows; i++) {
      if (dtemp[i] == 1.) bdry[i] = 'T';
      else bdry[i] = 'F';
   }
   ML_free(dtemp);

   /* ********************************************************************** */
   /* retrive the pointer to the ML_Aggregate_Options, which contains the    */
   /* number of aggregates (or their size), as well as few options for the   */
   /* constructions.                                                         */
   /* ********************************************************************** */

   aggr_options = (ML_Aggregate_Options *)ml_ag->aggr_options;
   if( aggr_options[ml_ag->cur_level].id != ML_AGGREGATE_OPTIONS_ID ) {
     fprintf( stderr,
	      "*ML*ERR* `ML_Aggregate_CoarsenMETIS' : wrong object\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }

   /* ********************************************************************** */
   /* allocate memory for aggr_index and call METIS to decompose the local   */
   /* graph into the number of parts specified by the user with a call       */
   /* ML_Aggregate_Set_LocalNumber( ml, ag, level, Nparts)                   */
   /* ********************************************************************** */

   nbytes = (Nrows+Nghost) * sizeof(int);

   if ( nbytes > 0 ) ML_memory_alloc((void**) &aggr_index, nbytes, "ACJ");
   else              aggr_index = NULL;

   if( aggr_options == NULL ) {
     fprintf( stderr,
	      "*ERR*ML* on proc %d, aggr_options is NULL\n"
	      "*ERR*ML* Please set the aggregate properties using\n"
	      "*ERR*ML* ML_SetNumberLocalAggregates or ML_SetAggregatesDimension\n",
	      "*ERR*ML* (file %s, line %d)\n",
	      mypid,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }

   aggr_count = aggr_options[ml_ag->cur_level].Naggregates;

   /* ********************************************************************** */
   /* if global number of aggregates has been specified, compute the local   */
   /* one (evenly dividing this global number)                               */
   /* ********************************************************************** */

   if( aggr_options[ml_ag->cur_level].local_or_global == ML_GLOBAL_INDICES ) {

     mod = aggr_count % Nprocs;
     
     aggr_count /= Nprocs;
     if( mypid == 0 ) {
       aggr_count += mod;
     }

     if( aggr_count < 1 ) aggr_count = 1;
     
   }
   
   /* ********************************************************************** */
   /* if aggr_count is minus 1,  this means that the user has set the # of   */
   /* nodes in each aggregate and not the local number. So, I compute this   */
   /* latter value which is the one needed by METIS. Note that value of      */
   /* aggr_count too close to the number of vertices may result in bizarre   */
   /* behavior of METIS                                                      */
   /* ********************************************************************** */

#ifdef METIS_DEBUG
   printf("aggr_count = %d\n",aggr_count);
#endif

   if( aggr_count == -1 ) {

     i = aggr_options[ml_ag->cur_level].Nnodes_per_aggregate;

     if( aggr_options[ml_ag->cur_level].Nnodes_per_aggregate >= Nrows) {
       aggr_count = 1;
     } else {
       aggr_count = (Nrows/i);
     }

     if ( mypid == 0 && 8 < ML_Get_PrintLevel() )  {
       printf("ML_Aggregate_CoarsenMETIS (level %d, proc 0) : "
	      "%d nodes/aggr (%d aggregates)\n",
	      ml_ag->cur_level,
	      i,
	      aggr_count );
     }

     aggr_options[ml_ag->cur_level].Naggregates = aggr_count;
     
   } else {

     if ( mypid == 0 && 8 < ML_Get_PrintLevel() )  {
       printf("ML_Aggregate_CoarsenMETIS (level %d, proc 0) : "
	      "%d aggregates\n",
	      ml_ag->cur_level,
	      aggr_count );
     }
     
   }
   
   if( aggr_count<=0 ) {
     fprintf( stderr,
	      "*ML*ERR* on proc %d, value of aggr_count not correct (%d)\n"
	      "*ML*ERR* Set this value using ML_Aggregate_Set_LocalNumber\n"
	      "*ML*ERR* or ML_Aggregate_Set_NodesPerAggr\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      mypid,
	      aggr_count,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   

#ifdef METIS_DEBUG
   printf("aggr_count (after analysis) = %d\n",aggr_count);
#endif
   
   /* ********************************************************************** */
   /* to call METIS, we have to create a graph using CSR data format         */
   /* Essentially, this requires a matrix-vector product (on the Amalgamated */
   /* matrix, so with dropped elements)                                      */
   /* ********************************************************************** */

   reorder_flag = aggr_options[ml_ag->cur_level].reordering_flag;
   
   ML_DecomposeGraph_with_METIS( Amatrix,Nnonzeros,aggr_count,
				 aggr_index, bdry,
				 ML_LOCAL_INDICES, NULL,
				 reorder_flag, ml_ag->cur_level );

   /* ********************************************************************** */
   /* I allocate room to copy aggr_index and pass this value to the user,    */
   /* who will be able to analyze and visualize this after the construction  */
   /* of the levels. This way, the only price we have to pay for stats and   */
   /* viz is essentially a little bit of memory.                             */
   /* this memory will be cleaned with the object ML_Aggregate ml_ag.        */
   /* I set the pointers using the ML_Aggregate_Info structure. This is      */
   /* allocated using ML_Aggregate_Info_Setup(ml,MaxNumLevels)               */
   /* ********************************************************************** */

#ifdef METIS_DEBUG
   printf("ml_ag->aggr_viz_and_stats = %x\n",ml_ag->aggr_viz_and_stats);
#endif

   if( ml_ag->aggr_viz_and_stats != NULL ) {

     ML_memory_alloc((void*)&graph_decomposition,
		     sizeof(int)*Nrows,"graph decomposition");
     if( graph_decomposition == NULL ) {
       fprintf( stderr,
		"*ML*ERR* Not enough memory for %d bytes\n"
		"*ML*ERR* (file %s, line %d)\n",
		sizeof(int)*Nrows,
		__FILE__,
	      __LINE__ );
       exit( EXIT_FAILURE );
     }

     for( i=0 ; i<Nrows ; i++ ) graph_decomposition[i] = aggr_index[i];

     aggr_viz_and_stats = (ML_Aggregate_Viz_Stats *) (ml_ag->aggr_viz_and_stats);
     aggr_viz_and_stats[ml_ag->cur_level].graph_decomposition = graph_decomposition;
     aggr_viz_and_stats[ml_ag->cur_level].Nlocal = Nrows;
     aggr_viz_and_stats[ml_ag->cur_level].Naggregates = aggr_count;
     aggr_viz_and_stats[ml_ag->cur_level].local_or_global = ML_LOCAL_INDICES;
     aggr_viz_and_stats[ml_ag->cur_level].is_filled = ML_YES;
     
   }
      
#ifdef METIS_DEBUG
   printf("proc %d: aggr_count = %d\n", mypid, aggr_count);
   printf("proc %d: Nrows = %d, nvertices = %d, Nghost = %d\n",
	  mypid, Nrows, nvertices, Nghost);
#endif
#ifdef NOT_METIS
   tempaggr_index = (int *) ML_allocate(sizeof(int)* (exp_Nrows+1)*num_PDE_eqns);
   for (i = 0;         i < nvertices ; i++) tempaggr_index[i] = aggr_index[i];
   for (i = nvertices; i < exp_Nrows ; i++) tempaggr_index[i] = -1;
   ML_free(aggr_index);
   aggr_index = tempaggr_index;
#endif
   
   N_neighbors = getrow_obj->pre_comm->N_neighbors;
   nbytes = N_neighbors * sizeof( int );
   if ( nbytes > 0 ) {
         ML_memory_alloc((void**) &neighbors,  nbytes, "AGL");
         ML_memory_alloc((void**) &recv_leng,  nbytes, "AGM");
         ML_memory_alloc((void**) &send_leng,  nbytes, "AGN");
   } 
   else {
         neighbors = recv_leng = send_leng = NULL;
   }
   for ( i = 0; i < N_neighbors; i++ ) {
      neighbors[i] = getrow_obj->pre_comm->neighbors[i].ML_id;
      recv_leng[i] = getrow_obj->pre_comm->neighbors[i].N_rcv;
      send_leng[i] = getrow_obj->pre_comm->neighbors[i].N_send;
   }
   total_recv_leng = total_send_leng = 0;
   for ( i = 0; i < N_neighbors; i++ ) {
         total_recv_leng += recv_leng[i];
         total_send_leng += send_leng[i];
   }
   nbytes = total_send_leng * sizeof( int );
   if ( nbytes > 0 ) ML_memory_alloc((void**) &send_list,nbytes,"AGO");
   else              send_list = NULL;
   if ( total_recv_leng+Nrows != exp_Nrows ) {
         printf("%d : ML_Aggregate_CoarsenMETIS - internal error.\n",mypid);
         printf("     lengths = %d %d \n",total_recv_leng+Nrows,exp_Nrows);
         exit(-1);
   }
   count = 0;
   for ( i = 0; i < N_neighbors; i++ ) {
         for (j = 0; j < send_leng[i]; j++)
            send_list[count++] = 
               getrow_obj->pre_comm->neighbors[i].send_list[j];
   }
   if ( count > total_send_leng ) {
         printf("%d : CoarsenMEIIS ERROR : count < total_send_leng\n",mypid);
         exit(1);
   }

   /* ********************************************************************** */
   /* take the decomposition as created by METIS and form the aggregates     */
   /* ********************************************************************** */

   total_nz = 0;
   count    = 0;
   for (i = 0; i < nvertices; i++) {
      if (aggr_index[i] >= 0) {
         current = -aggr_index[i]-2;
         aggr_index[i] = current;
         ML_get_matrix_row(Amatrix,1,&i,&allocated,&rowi_col,&rowi_val,&rowi_N,0);
         for (j = 0; j < rowi_N; j++) aggr_index[rowi_col[j]] = current;
      }
      else {
         /* still get the rows to count nonzeros */
         ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                           &rowi_N, 0);
      }
      total_nz += rowi_N;
   }
   total_nz = ML_Comm_GsumInt( comm, total_nz);
   i = ML_Comm_GsumInt( comm, nvertices);

   if( rowi_col != NULL ) ML_free(rowi_col); rowi_col = NULL;
   if( rowi_val != NULL ) ML_free(rowi_val); rowi_val = NULL;
   allocated = 0;
   
   if ( mypid == 0 && 8 < ML_Get_PrintLevel())
     printf("Aggregation(METIS) : Total nonzeros = %d (Nrows=%d)\n",total_nz,i);
   
   if ( ml_ag->operator_complexity == 0.0 ) {
      ml_ag->fine_complexity = total_nz;
      ml_ag->operator_complexity = total_nz;
   }
   else ml_ag->operator_complexity += total_nz;

   for (i = 0; i < exp_Nrows; i++) 
      aggr_index[i] = -aggr_index[i]-2;

   /* make sure that boundary nodes are not in any aggregate */

   for (i = 0; i < exp_Nrows; i++)
      if (bdry[i] == 'T') aggr_index[i] = -1;

   /* communicate the aggregate information.  Use temp_index and tem2_index as */
   /* communication buffer arrays.                                             */
   temp_index = (int *) ML_allocate(sizeof(int)*(1+total_send_leng+total_recv_leng));
   tem2_index = (int *) ML_allocate(sizeof(int)*(1+total_send_leng+total_recv_leng));
   temp_leng  = (int *) ML_allocate(sizeof(int)*(N_neighbors+1));

   ML_aggr_index_communicate(N_neighbors, temp_leng, send_leng, recv_leng, send_list, 
			     recv_list, nvertices, comm, aggr_index, 1563, tem2_index, 
			     neighbors, temp_index,1);

   agg_indx_comm.N_neighbors = N_neighbors;
   agg_indx_comm.temp_leng = temp_leng;
   agg_indx_comm.send_leng = send_leng;
   agg_indx_comm.recv_leng = recv_leng;
   agg_indx_comm.send_list = send_list;
   agg_indx_comm.recv_list = recv_list;
   agg_indx_comm.tem2_index = tem2_index;
   agg_indx_comm.neighbors = neighbors;
   agg_indx_comm.temp_index = temp_index;

   /* I am not sure that this is the good position for this... */
   
   dtemp = (double *) malloc( sizeof(double) * (Nrows+Nghost));
   for( i=0 ; i<Nrows ; i++ ) dtemp[i] = (double)aggr_index[i];
   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm,Amatrix->outvec_leng,
                    comm, ML_OVERWRITE,NULL);

   for( i=Nrows ; i<Nrows+Nghost ; i++ ) aggr_index[i] = (int)dtemp[i];
   /*   
   for( i=0 ; i<Nrows+Nghost ; i++ )
     printf("aggr_index[%d] = %d\n", i, aggr_index[i] );
   */
   free( dtemp );
   
   
#ifdef ML_AGGR_MARKINAGGR

   for (i = 0; i < exp_Nrows; i++) aggr_index[i] = -1;
   sprintf(fname,"agg_%d",level_count); level_count++;
   fp = fopen(fname,"r");
   aggr_count = 0;
   for (i = 0; i <nvertices; i++) {
      fscanf(fp,"%d%d",&k,&j);
      aggr_index[j] = k;
      if (k >= aggr_count) aggr_count = k+1;
   }
   fclose(fp);
#endif /*ifdef ML_AGGR_MARKINAGGR*/

   /* make sure that boundary nodes are not in any aggregate */
   
   for (i = 0; i < exp_Nrows; i++) {
     if (bdry[i] == 'T') aggr_index[i] = -1;
     else if (aggr_index[i] == -1)
	 {
	   printf("for node %d, bndry[.] = %c\n",i,bdry[i]);
       printf("ML_agg_METIS: I'm not sure who takes care of this guy\n");
	   printf("probably need to use the other aggregate file...\n");
     }
   }
   ML_free(bdry);

   for (i = exp_Nrows - 1; i >= 0; i-- ) {
      for (j = num_PDE_eqns-1; j >= 0; j--) {
         aggr_index[i*num_PDE_eqns+j] = aggr_index[i];
      }
   }

   if ( mypid == 0 && 8 < ML_Get_PrintLevel())
   {
      printf("Calling ML_Operator_UnAmalgamateAndDropWeak\n");
      fflush(stdout);
   }

   ML_Operator_UnAmalgamateAndDropWeak(Amatrix, num_PDE_eqns, epsilon);

   Nrows      *= num_PDE_eqns;
   nvertices  *= num_PDE_eqns;
   exp_Nrows  *= num_PDE_eqns;
   getrow_obj  = Amatrix->getrow;
   N_neighbors = getrow_obj->pre_comm->N_neighbors;

#ifdef ML_AGGR_INAGGR

   for (i = 0; i < exp_Nrows; i++) aggr_index[i] = -1;
   sprintf(fname,"agg_%d",level_count); level_count++;
   fp = fopen(fname,"r");
   if (fp == NULL)
   {
      printf("Cannot open aggregate file %s for reading.\n",fname);
      printf("exp_Nrows = %d\n",exp_Nrows);
      exit(1);
   }
   else
   {
      printf("Reading aggregate for level %d from file %s\n",
             level_count-1,fname);
      fflush(stdout);
   }
   aggr_count = 0;
   for (i = 0; i <nvertices; i++) {
     if ( fscanf(fp,"%d%d",&j,&k) != 2) break;
      aggr_index[j] = k;
      if (k >= aggr_count) aggr_count = k+1;
   }
   fclose(fp);
   printf("Read in %d aggregates\n\n",aggr_count);
#endif /*ifdef ML_AGGR_INAGGR*/


   /* some memory free (I am not completely sure about this...) */

   if( recv_list != NULL ) ML_free(recv_list);
   if( temp_index != NULL ) ML_free(temp_index);
   if( tem2_index != NULL ) ML_free(tem2_index);
   if( temp_leng != NULL ) ML_free(temp_leng);

   /* I'm not sure if I need most of this 'if' code. I just took it from */
   /* Charles ... but I guess that the majority of it is not needed.     */

   if ( num_PDE_eqns != 1 )
   {
      ML_memory_free((void**) &neighbors);
      ML_memory_free((void**) &recv_leng);
      ML_memory_free((void**) &send_leng);
      ML_memory_free((void**) &send_list);
      /* ---------------------------------------------------------- */
      /* allocate storage for the communication information         */
      /* ---------------------------------------------------------- */

      N_neighbors = getrow_obj->pre_comm->N_neighbors;
      nbytes = N_neighbors * sizeof( int );
      if ( nbytes > 0 )
      {
         ML_memory_alloc((void**) &neighbors,  nbytes, "AGL");
         ML_memory_alloc((void**) &recv_leng,  nbytes, "AGM");
         ML_memory_alloc((void**) &send_leng,  nbytes, "AGN");
      }
      else
      {
	neighbors = recv_leng = send_leng = NULL;
      }
      for ( i = 0; i < N_neighbors; i++ )
      {
         neighbors[i] = getrow_obj->pre_comm->neighbors[i].ML_id;
         recv_leng[i] = getrow_obj->pre_comm->neighbors[i].N_rcv;
         send_leng[i] = getrow_obj->pre_comm->neighbors[i].N_send;
      }
      total_recv_leng = 0;
      for ( i = 0; i < N_neighbors; i++ ) total_recv_leng += recv_leng[i];
      total_send_leng = 0;
      for ( i = 0; i < N_neighbors; i++ ) total_send_leng += send_leng[i];
      nbytes = total_send_leng * num_PDE_eqns * sizeof( int );
      if ( nbytes > 0 ) ML_memory_alloc((void**) &send_list,nbytes,"AGO");
      else              send_list = NULL;
      nbytes = total_recv_leng * sizeof( int );
      if ( nbytes > 0 ) ML_memory_alloc((void**) &recv_list,nbytes,"AGP");
      else              recv_list = NULL;

      /* ---------------------------------------------------------- */
      /* set up true external indices to be shipped to receive      */
      /* processors (true in view of that num_PDE_eqns can be > 1)  */
      /* ---------------------------------------------------------- */

      nbytes = Nrows * sizeof( int );
      if ( nbytes > 0 ) ML_memory_alloc((void**) &itmp_array,nbytes,"AGQ");
      count = 0;
      for ( i = 0; i < N_neighbors; i++ )
      {
         for ( j = 0; j < Nrows; j++ ) itmp_array[j] = -1;
         count3 = 0;
         for (j = 0; j < send_leng[i]; j++)
         {
            index3 = getrow_obj->pre_comm->neighbors[i].send_list[j];
            index3 = index3 / num_PDE_eqns * num_PDE_eqns;
            for (k = 0; k < num_PDE_eqns; k++)
            {
               if ( itmp_array[index3+k] < 0 )
                  itmp_array[index3+k] = count3++;
            }
         }
         for (j = 0; j < send_leng[i]; j++)
         {
            send_list[count+j] =
               getrow_obj->pre_comm->neighbors[i].send_list[j];
         }
         for ( j = 0; j < send_leng[i]; j++ )
         {
            index = send_list[count+j];
            if (itmp_array[index] >= 0) send_list[count+j] = itmp_array[index];
         }
         count += send_leng[i];
      }
      ML_memory_free((void**) &itmp_array);

      /* ---------------------------------------------------------- */
      /* send the adjusted indices to the receive processors        */
      /* ---------------------------------------------------------- */

      if ( N_neighbors > 0 )
         request = (USR_REQ *) ML_allocate(N_neighbors*sizeof(USR_REQ));

      offset = 0;
      for (i = 0; i < N_neighbors; i++)
      {
         msgtype = 2000;
         length = recv_leng[i] * sizeof( int );
         procnum = neighbors[i];
         comm->USR_irecvbytes((void *) &(recv_list[offset]),length,&procnum,
                              &msgtype, comm->USR_comm, request+i);
         offset += recv_leng[i];
      }
      offset = 0;
      for (i = 0; i < N_neighbors; i++)
      {
         msgtype = 2000;
         length = send_leng[i] * sizeof( int );
         procnum = neighbors[i];
         comm->USR_sendbytes((void *) &(send_list[offset]),length,procnum,
                              msgtype, comm->USR_comm);
         offset += send_leng[i];
      }
      offset = 0;
      for (i = 0; i < N_neighbors; i++)
      {
         msgtype = 2000;
         length = recv_leng[i] * sizeof( int );
         procnum = neighbors[i];
         comm->USR_cheapwaitbytes((void *) &(recv_list[offset]),length,&procnum,
                              &msgtype, comm->USR_comm, request+i);
         for (j = 0; j < recv_leng[i]; j++) recv_list[offset+j] += offset;
         offset += recv_leng[i];
      }
      if ( N_neighbors > 0 ) ML_free( request );

      ML_memory_free((void**) &recv_list);

      /* ---------------------------------------------------------- */
      /* update the send_list and send_leng's in line with remap    */
      /* ---------------------------------------------------------- */

      nbytes = Nrows * sizeof( int );
      if (nbytes > 0) ML_memory_alloc((void**) &itmp_array,nbytes,"AGR");
      total_send_leng = 0;
      for ( i = 0; i < N_neighbors; i++ )
      {
         count = 0;
         for ( j = 0; j < Nrows; j++ ) itmp_array[j] = -1;
         for (j = 0; j < send_leng[i]; j++)
         {
            index3 = getrow_obj->pre_comm->neighbors[i].send_list[j];
            index3 = index3 / num_PDE_eqns * num_PDE_eqns;
            for (k = 0; k < num_PDE_eqns; k++)
               itmp_array[index3+k] = 0;
         }
         for ( j = 0; j < Nrows; j++ )
         {
            if ( itmp_array[j] == 0 ) send_list[total_send_leng+count++] = j;
         }
         send_leng[i] = count;
         total_send_leng += count;
      }
      total_send_leng = 0;
      for ( i = 0; i < N_neighbors; i++ ) total_send_leng += send_leng[i];

      ML_memory_free((void**) &itmp_array);

      /* ---------------------------------------------------------- */
      /* update other processors with the new communication pattern */
      /* ---------------------------------------------------------- */

      if ( N_neighbors > 0 )
         request = (USR_REQ *) ML_allocate(N_neighbors*sizeof(USR_REQ));

      for (i = 0; i < N_neighbors; i++)
      {
         msgtype = 2002;
         length = sizeof( int );
         procnum = neighbors[i];
         comm->USR_irecvbytes((void *) &(recv_leng[i]), length, &procnum,
                              &msgtype, comm->USR_comm, request+i);
      }
      for (i = 0; i < N_neighbors; i++)
      {
         msgtype = 2002;
         length = sizeof( int );
         procnum = neighbors[i];
         comm->USR_sendbytes((void *) &(send_leng[i]), length, procnum,
                              msgtype, comm->USR_comm);
      }
      for (i = 0; i < N_neighbors; i++)
      {
         msgtype = 2002;
         length = sizeof( int );
         procnum = neighbors[i];
         comm->USR_cheapwaitbytes((void *) &(recv_leng[i]), length, &procnum,
                              &msgtype, comm->USR_comm, request+i);
      }
      if ( N_neighbors > 0 ) ML_free( request );

      total_recv_leng = 0;
      for (i = 0; i < N_neighbors; i++) total_recv_leng += recv_leng[i];
      exp_Nrows = Nrows + total_recv_leng;
   }


   /* count the size of each aggregate */

   aggr_cnt_array = (int *) ML_allocate(sizeof(int)*aggr_count);
   for (i = 0; i < aggr_count ; i++) aggr_cnt_array[i] = 0;
   for (i = 0; i < exp_Nrows; i++) 
      if (aggr_index[i] >= 0) 
         aggr_cnt_array[aggr_index[i]]++;

   /*
     for (i = 0; i < exp_Nrows; i++) printf("%d: AGG_INDX %d %d\n",comm->ML_mypid, i,aggr_index[i]);
     for (i = 0; i < aggr_count ; i++) printf("counts %d %d\n",i,aggr_cnt_array[i]);*/

#ifdef ML_AGGR_OUTAGGR
   sprintf(fname,"agg%d_%d",comm->ML_mypid,level_count);
   fp = fopen(fname,"w");
   agg_offset = ML_gpartialsum_int(aggr_count, comm);
   vertex_offset = ML_gpartialsum_int(nvertices, comm);
   for (i = 0; i < nvertices ; i++)
   {
#ifndef MAXWELL
#ifdef ALEGRA
      if (level_count == 0) { j = i; k = i+vertex_offset;}
#else
      if (level_count == 0) { j = update_index[i]; k = update[i];}
#endif /*ifdef ALEGRA*/
#else
      if (level_count == 0) { j = reordered_glob_nodes[i]; k = global_node_inds[i];}
#endif /* ifndef MAXWELL */
      else                  { j = i              ; k = i+vertex_offset;}
      if (aggr_index[j] >= 0)
         fprintf(fp,"%5d %5d\n",k, aggr_index[j]+agg_offset);
   }

   dtemp = (double *) ML_allocate(sizeof(double)*(exp_Nrows+1));
   for (i = 0; i < nvertices; i++) dtemp[i] = (double) (i + vertex_offset);
   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm, nvertices, comm,
                    ML_OVERWRITE,NULL);
   for (i = 0; i < exp_Nrows-nvertices; i++)
   {
#ifndef MAXWELL
#ifdef ALEGRA
      if (level_count == 0) { j = i+nvertices; k =  (int) dtemp[i+nvertices];}
#else
      if (level_count == 0) { j = extern_index[i]; k = external[i];} 
#endif /*ifdef ALEGRA*/
#else
      if (level_count == 0) { j = reordered_node_externs[i]; k =
	  global_node_externs[i];}
#endif /* ifndef MAXWELL */

      else                 { j = i+nvertices    ; k = (int) dtemp[i+nvertices];}
      if (aggr_index[j] >= 0)
         fprintf(fp,"%5d %5d\n", k, aggr_index[j]+agg_offset);
   }
   fclose(fp);
   level_count++;
   ML_free(dtemp);
#else
#ifdef INPUT_AGGREGATES
   agg_offset = ML_gpartialsum_int(aggr_count, comm);
   vertex_offset = ML_gpartialsum_int(nvertices, comm);
#endif
#endif /*ifdef ML_AGGR_OUTAGGR*/

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */

   Ncoarse = aggr_count;

   /* ============================================================= */
   /* update aggr_index to find out which local fine node is mapped */
   /* to which coarse aggregate in remote processors                */
   /* ------------------------------------------------------------- */

   nbytes = total_send_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf, nbytes, "AGg");
   else              int_buf = NULL;
   nbytes = total_recv_leng * sizeof(int);
   if ( nbytes > 0 ) ML_memory_alloc((void**) &int_buf2, nbytes, "AGh");
   else              int_buf2 = NULL;

   /* ------------------------------------------------------------- */
   /* send the remote node index back to remote processors, with   */
   /* added information on which remote nodes have been aggregated */
   /* by the local aggregates (and also the aggregate numbers).    */
   /* ------------------------------------------------------------- */

   offset = 0;
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      for ( j = 0; j < recv_leng[i]; j++ ) 
      {
         if ( aggr_index[Nrows+offset+j] < 0 ) int_buf2[offset+j] = -1;
         else int_buf2[offset+j] = aggr_index[Nrows+offset+j];
      }
      offset += recv_leng[i];
   }
   msgtype = 15963;
   ML_Aggregate_ExchangeData((char*) int_buf, (char*) int_buf2,
      N_neighbors, neighbors, send_leng, recv_leng, msgtype, ML_INT, comm);

   if ( int_buf2 != NULL ) ML_memory_free((void**) &int_buf2);

   /* ------------------------------------------------------------- */
   /* if int_buf[i] > 0, this means that aggr_index[send_list[i]]   */ 
   /* has been aggregated by a remote processor                     */
   /* int_buf2 is used to tabulate how many distinct aggregates     */
   /* in remote processors are used.                                */
   /* ------------------------------------------------------------- */

   offset = 0;
   m      = 0; /* store the local index offset for remote processors */ 
   new_N_recv = 0;
   nbytes = N_neighbors * sizeof(int);
   if ( nbytes > 0 ) 
   {
      ML_memory_alloc((void**) &new_recv_leng, nbytes, "AGi");
      ML_memory_alloc((void**) &new_recv_neighbors, nbytes, "AGj");
   } 
   else 
   {
      new_recv_leng = new_recv_neighbors = NULL;
   }
   for ( i = 0; i < N_neighbors; i++ ) 
   {
      /* ---------------------------------------------------------- */
      /* find out how large an array to allocate for int_buf2       */
      /* ---------------------------------------------------------- */

      max_count = -1;
      for ( j = 0; j < send_leng[i]; j++ ) 
      {
         index = int_buf[offset+j];
         max_count = (index > max_count ) ? index : max_count;
      }
      nbytes = ( max_count + 2 ) * sizeof(int);
      if (nbytes > 0) ML_memory_alloc((void **) &int_buf2, nbytes, "AGk");

      /* ---------------------------------------------------------- */
      /* see how many distinct remote aggregates are referenced by  */
      /* local fine nodes in aggregation in proc i ==> count        */
      /* ---------------------------------------------------------- */

      for ( j = 0; j <= max_count; j++ ) int_buf2[j] = 0;
      for ( j = 0; j < send_leng[i]; j++ ) 
      {
         index = int_buf[offset+j];
         if ( index >= 0 ) int_buf2[index]++;
         if (index >= 0 && index > max_count) 
            {printf("int_buf2 error : maxcount\n");exit(1);}
      }
      count = 0;
      for ( j = 0; j <= max_count; j++ ) 
      {
         if (int_buf2[j] > 0) 
         {
            count++; int_buf2[j] = 1;
         }
      }
      for ( j = max_count; j > 0; j-- ) int_buf2[j] = int_buf2[j-1];
      int_buf2[0] = 0;
      for ( j = 0; j < max_count; j++ ) int_buf2[j+1] += int_buf2[j];

      if ( count > 0 ) 
      {
         new_recv_leng[new_N_recv] = count * nullspace_dim;
         new_recv_neighbors[new_N_recv] = neighbors[i];
         new_N_recv++;
      } 

      /* ---------------------------------------------------------- */
      /* now assign local aggregate indices to local nodes that are */
      /* aggregated by remote processors                            */
      /* ---------------------------------------------------------- */

      for ( j = 0; j < send_leng[i]; j++ ) 
      {
         index = send_list[offset+j];

         /* ------------------------------------------------------- */
         /* The first condition indicates that the local node has   */
         /* been registered to have been aggregated by remote       */
         /* aggregates.  The second condition is needed in case     */
         /* the local node is linked to more than 1 remote          */
         /* processor (but only to one aggregate though)            */
         /* int_buf2 contains local indices of remote aggregates    */
         /* ------------------------------------------------------- */

         if ( aggr_index[index] <= -100 && int_buf[offset+j] >= 0 ) 
         {
            k = int_buf[offset+j];
            aggr_index[index] = int_buf2[k] + Ncoarse + m;
         } 
      }
      if (nbytes > 0) ML_memory_free((void **) &int_buf2);
      m += count;
      offset += send_leng[i];
   }
   exp_Ncoarse = Ncoarse + m;
 
   if ( int_buf  != NULL ) ML_memory_free((void**) &int_buf);

   /* ============================================================= */
   /* check and copy aggr_index                                     */
   /* ------------------------------------------------------------- */

   level = ml_ag->cur_level;
   nbytes = Nrows * sizeof( int );
   ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), nbytes, "AGl");
   count = aggr_count;
   for ( i = 0; i < Nrows; i+=num_PDE_eqns ) 
   {
      if ( aggr_index[i] >= 0 )
      {
         for ( j = 0; j < num_PDE_eqns; j++ ) 
            ml_ag->aggr_info[level][i+j] = aggr_index[i];
         if (aggr_index[i] >= count) count = aggr_index[i] + 1;
      }
      /*else
       *{
       *   printf("%d : CoarsenMIS error : aggr_index[%d] < 0\n",
       *          mypid,i);
       *   exit(1);
       *}*/
   }
   ml_ag->aggr_count[level] = count; /* for relaxing boundary points */
  
   /* ============================================================= */
   /* find out how many local coarse aggregates are needed by       */
   /* remote processors for interpolation (to construct the         */
   /* communicator - send info - for P)                             */
   /* ------------------------------------------------------------- */
#ifdef METIS_DEBUG
   printf("N_neighbors = %d\n", N_neighbors);
#endif
   new_N_send = 0;
   if ( N_neighbors > 0 ) 
   {
      nbytes = N_neighbors * sizeof(int);
      ML_memory_alloc((void**) &int_buf, nbytes, "AGm");
      nbytes = Ncoarse * sizeof(int);
      ML_memory_alloc((void**) &int_buf2, nbytes, "AGn");
      for ( i = 0; i < N_neighbors; i++ ) int_buf[i] = 0;

      /* ---------------------------------------------------------- */
      /* count which remote fine nodes belong to local aggregates   */
      /* in order to generate the communication pattern for         */
      /* the interpolation operator.                                */
      /* ---------------------------------------------------------- */

      offset = Nrows; 
      for ( i = 0; i < N_neighbors; i++ ) 
      {
         for ( j = 0; j < Ncoarse; j++ ) int_buf2[j] = 0;
         for ( j = 0; j < recv_leng[i]; j++ ) 
         {
            index = aggr_index[offset++];
            if ( index >= 0 ) int_buf2[index]++;
         }
         count = 0;
         for ( j = 0; j < Ncoarse; j++ ) if ( int_buf2[j] > 0 ) count++;
         int_buf[i] = count * nullspace_dim;
         if ( int_buf[i] > 0 ) new_N_send++;
      }

      /* ---------------------------------------------------------- */
      /* now the number of neighbors for P has been found, the next */
      /* step is to find the send_list and send_leng for the matvec */
      /* function for interpolation                                 */
      /* ---------------------------------------------------------- */

      nbytes = new_N_send * sizeof(int);
      if ( nbytes > 0 ) 
      {
         ML_memory_alloc((void**) &new_send_leng, nbytes, "AGo");
         ML_memory_alloc((void**) &new_send_neighbors, nbytes, "AGp");
         new_N_send = 0;
         for ( i = 0; i < N_neighbors; i++ ) 
         {
            if ( int_buf[i] > 0 ) 
            {
               new_send_leng[new_N_send] = int_buf[i]; 
               new_send_neighbors[new_N_send] = neighbors[i];
               new_N_send++;
            }
         }
         count = 0;
         for ( i = 0; i < new_N_send; i++ ) count += new_send_leng[i];
         nbytes = count * sizeof(int);
         ML_memory_alloc((void**) &new_send_list, nbytes, "AGq");
         offset = Nrows;
         m = count;
         count = 0;
         for ( i = 0; i < N_neighbors; i++ ) 
         {
            for ( j = 0; j < Ncoarse; j++ ) int_buf2[j] = 0;
            for ( j = 0; j < recv_leng[i]; j++ ) 
            {
               index = aggr_index[offset++];
               if ( index >= 0 ) int_buf2[index]++;
            }
            for ( j = 0; j < Ncoarse; j++ ) 
            {
               if ( int_buf2[j] > 0 ) 
               {
                  for ( jj = 0; jj < nullspace_dim; jj++ ) 
                     new_send_list[count++] = j * nullspace_dim + jj;
               } 
            } 
         } 
         if ( m != count ) 
         {
            printf("ML_Aggregate_CoarsenMIS : internal error (1).\n");
            exit(-1);
         }
      } 
      else 
      {
         new_send_leng = NULL;
         new_send_neighbors = NULL;
         new_send_list = NULL;
      }  
      ML_memory_free((void**) &int_buf);
      ML_memory_free((void**) &int_buf2);
   } 
   else 
   {
      new_send_leng = NULL;
      new_send_neighbors = NULL;
      new_send_list = NULL;
   }

   /* ============================================================= */
   /* set up the new operator                                       */
   /* ------------------------------------------------------------- */

   new_Nrows = Nrows;
   for ( i = 0; i < new_Nrows; i++ ) 
   {
      if ( aggr_index[i] >= exp_Ncoarse ) 
      {
         printf("WARNING : index out of bound %d = %d(%d)\n", i, aggr_index[i], 
                exp_Ncoarse);
/*
         for ( j = 0; j < new_Nrows; j++ ) 
            if ( aggr_index[j] >= exp_Ncoarse )
               printf("%d : aggr_index[%5d] = %5d *\n", mypid, j, aggr_index[j]); 
            else
               printf("%d : aggr_index[%5d] = %5d \n", mypid, j, aggr_index[j]); 
         exit(1);
*/
      }
   }
   nbytes = ( new_Nrows + 1 ) * sizeof(int); 
   ML_memory_alloc((void**)&(new_ia), nbytes, "AIA");
   nbytes = new_Nrows * nullspace_dim * sizeof(int); 
   ML_memory_alloc((void**)&(new_ja), nbytes, "AJA");
   nbytes = new_Nrows * nullspace_dim * sizeof(double); 
   ML_memory_alloc((void**)&(new_val), nbytes, "AVA");
   for ( i = 0; i < new_Nrows*nullspace_dim; i++ ) new_val[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* set up the space for storing the new null space               */
   /* ------------------------------------------------------------- */

   nbytes = Ncoarse * nullspace_dim * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&(new_null),nbytes,"AGr");
   for (i = 0; i < Ncoarse*nullspace_dim*nullspace_dim; i++) 
      new_null[i] = 0.0;

   /* ------------------------------------------------------------- */
   /* initialize the row pointer for the CSR prolongation operator  */
   /* (each row will have at most nullspace_dim nonzero entries)    */
   /* ------------------------------------------------------------- */

   for (i = 0; i <= Nrows; i++) new_ia[i] = i * nullspace_dim;

   /* trying this when a Dirichlet row is taken out */
   j = 0;
   new_ia[0] = 0;
   for (i = 0; i < Nrows; i++) {
      if (aggr_index[i] != -1) j += nullspace_dim;
      new_ia[i+1] = j;
   }

   /* ------------------------------------------------------------- */
   /* generate an array to store which aggregate has which rows.Then*/
   /* loop through the rows of A checking which aggregate each row  */
   /* is in, and adding it to the appropriate spot in rows_in_aggs  */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**)&rows_in_aggs,aggr_count*sizeof(int*),"MLs");
   for (i = 0; i < aggr_count; i++) 
   {
      rows_in_aggs[i] = (int *) ML_allocate(aggr_cnt_array[i]*sizeof(int));
      aggr_cnt_array[i] = 0;
      if (rows_in_aggs[i] == NULL) 
      {
         printf("Error: couldn't allocate memory in CoarsenMETIS\n");
         exit(1);
      }
   }
   for (i = 0; i < exp_Nrows; i+=num_PDE_eqns) 
   {
      if ( aggr_index[i] >= 0 && aggr_index[i] < aggr_count)
      {
         for (j = 0; j < num_PDE_eqns; j++)
         {
            index = aggr_cnt_array[aggr_index[i]]++; 
            rows_in_aggs[aggr_index[i]][index] = i + j;
         }
      }
   }

   /* ------------------------------------------------------------- */
   /* allocate work arrays for QR factorization                     */
   /* work and lwork are needed for lapack's QR routine.  These     */
   /* settings seemed easiest since I don't quite understand        */
   /* what they do, but may want to do something better here later  */
   /* ------------------------------------------------------------- */

   nbytes = total_recv_leng * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&comm_val, nbytes, "AGt");
   for (i = 0; i < total_recv_leng*nullspace_dim; i++) comm_val[i] = 0.0; 
   max_agg_size = 0;
   for (i = 0; i < aggr_count; i++) 
   {
      if (aggr_cnt_array[i] > max_agg_size) max_agg_size = aggr_cnt_array[i];
   }
   nbytes = max_agg_size * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&qr_tmp, nbytes, "AGu");
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&tmp_vect, nbytes, "AGv");

   lwork  = nullspace_dim;
   nbytes = nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&work, nbytes, "AGw");

   /* ------------------------------------------------------------- */
   /* ship the null space information to other processors           */
   /* ------------------------------------------------------------- */
 
   if (nullspace_vect != NULL) 
   {
      nbytes = total_send_leng * nullspace_dim * sizeof(double);
      ML_memory_alloc((void**) &dble_buf, nbytes,"AG1");
      nbytes = total_recv_leng * nullspace_dim * sizeof(double);
      ML_memory_alloc((void**) &dble_buf2, nbytes,"AG2");
      length = total_send_leng * nullspace_dim;
      for ( i = 0; i < total_send_leng; i++ ) 
      {
         index = send_list[i];
         for ( j = 0; j < nullspace_dim; j++ ) 
            dble_buf[i*nullspace_dim+j] = nullspace_vect[j*Nrows+index];
      }
      msgtype = 12093;
      length = sizeof(double) * nullspace_dim;
      ML_Aggregate_ExchangeData((char*)dble_buf2,(char*) dble_buf,
            N_neighbors, neighbors, recv_leng, send_leng,msgtype,
				(int) length,comm);
      ML_memory_free((void**) &dble_buf);
   } 

   /* ------------------------------------------------------------- */
   /* perform block QR decomposition                                */
   /* ------------------------------------------------------------- */


   for (i = 0; i < aggr_count; i++) 
   {
      /* ---------------------------------------------------------- */
      /* set up the matrix we want to decompose into Q and R:       */
      /* ---------------------------------------------------------- */

      length = aggr_cnt_array[i];
      if (nullspace_vect == NULL) 
      {
         for (j = 0; j < length; j++)
         {
            index = rows_in_aggs[i][j];
            for (k = 0; k < nullspace_dim; k++)
            {
              if ( unamalg_bdry[index] == 'T') qr_tmp[k*length+j] = 0.;
               else
               {
                  if (index % num_PDE_eqns == k) qr_tmp[k*length+j] = 1.0;
                  else                           qr_tmp[k*length+j] = 0.0;
               }
            }
         }
      }
      else 
      {
         for (k = 0; k < nullspace_dim; k++)
         {
            for (j = 0; j < length; j++)
            {
               index = rows_in_aggs[i][j];
               if ( unamalg_bdry[index] == 'T') qr_tmp[k*length+j] = 0.;
               else {
                  if (index < Nrows) {
                     qr_tmp[k*length+j] = nullspace_vect[k*Nrows+index];
                  }
                  else {
                     qr_tmp[k*length+j] = 
                        dble_buf2[(index-Nrows)*nullspace_dim+k];
                  }
               }
            }
         }
      }

      /* ---------------------------------------------------------- */
      /* now calculate QR using an LAPACK routine                   */
      /* ---------------------------------------------------------- */
      if (aggr_cnt_array[i] >= nullspace_dim) {

	MLFORTRAN(dgeqrf)(&(aggr_cnt_array[i]), &nullspace_dim, qr_tmp, 
			  &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
	if (info != 0)
	  pr_error("ErrOr in CoarsenMIS : dgeqrf returned a non-zero %d %d\n",
		   aggr_cnt_array[i],i);

	if (work[0] > lwork) 
	  {
	    lwork=(int) work[0]; 
	    ML_memory_free((void**) &work);
	    ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGx");
	  }
	else lwork=work[0];
		 
	/* ---------------------------------------------------------- */
	/* the upper triangle of qr_tmp is now R, so copy that into   */
	/* the new nullspace                                          */
	/* ---------------------------------------------------------- */

	for (j = 0; j < nullspace_dim; j++)
	  for (k = j; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse*nullspace_dim] = 
	      qr_tmp[j+aggr_cnt_array[i]*k];
		 
	/* ---------------------------------------------------------- */
	/* to get this block of P, need to run qr_tmp through another */
	/* LAPACK function:                                           */
	/* ---------------------------------------------------------- */

	if ( aggr_cnt_array[i] < nullspace_dim ){
	  printf("Error in dorgqr on %d row (dims are %d, %d)\n",i,aggr_cnt_array[i],
                 nullspace_dim);
	  printf("ERROR : performing QR on a MxN matrix where M<N.\n");
	}
	MLFORTRAN(dorgqr)(&(aggr_cnt_array[i]), &nullspace_dim, &nullspace_dim, 
			  qr_tmp, &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
	if (info != 0) {
	  printf("Error in dorgqr on %d row (dims are %d, %d)\n",i,aggr_cnt_array[i],
                 nullspace_dim);
	  pr_error("Error in CoarsenMIS: dorgqr returned a non-zero\n");
	}

	if (work[0] > lwork) 
	  {
	    lwork=(int) work[0]; 
	    ML_memory_free((void**) &work);
	    ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGy");
	  }
	else lwork=work[0];

	/* ---------------------------------------------------------- */
	/* now copy Q over into the appropriate part of P:            */
	/* The rows of P get calculated out of order, so I assume the */
	/* Q is totally dense and use what I know of how big each Q   */
	/* will be to determine where in ia, ja, etc each nonzero in  */
	/* Q belongs.  If I did not assume this, I would have to keep */
	/* all of P in memory in order to determine where each entry  */
	/* should go                                                  */
	/* ---------------------------------------------------------- */

	for (j = 0; j < aggr_cnt_array[i]; j++)
	  {
	    index = rows_in_aggs[i][j];
	    if ( index < Nrows )
	      {
		index3 = new_ia[index];
		for (k = 0; k < nullspace_dim; k++) 
		  {
		    new_ja [index3+k] = i * nullspace_dim + k;
		    new_val[index3+k] = qr_tmp[ k*aggr_cnt_array[i]+j];
		  }
	      }
	    else 
	      {
		index3 = (index - Nrows) * nullspace_dim;
		for (k = 0; k < nullspace_dim; k++)
		  comm_val[index3+k] = qr_tmp[ k*aggr_cnt_array[i]+j];
	      }
	  }
      }
      else {
	/* We have a small aggregate such that the QR factorization can not */
	/* be performed. Instead let us copy the null space from the fine   */
        /* into the coarse grid nullspace and put the identity for the      */
	/* prolongator????                                                  */
	for (j = 0; j < nullspace_dim; j++)
	  for (k = 0; k < nullspace_dim; k++)
            new_null[i*nullspace_dim+j+k*Ncoarse*nullspace_dim] = 
	      qr_tmp[j+aggr_cnt_array[i]*k];
	for (j = 0; j < aggr_cnt_array[i]; j++) {
	  index = rows_in_aggs[i][j];
	  index3 = new_ia[index];
	  for (k = 0; k < nullspace_dim; k++) {
	    new_ja [index3+k] = i * nullspace_dim + k;
	    if (k == j) new_val[index3+k] = 1.;
	    else new_val[index3+k] = 0.;
	  }
	}
      }


   }
	 
   ML_Aggregate_Set_NullSpace(ml_ag, num_PDE_eqns, nullspace_dim, 
                              new_null, Ncoarse*nullspace_dim);
   ML_memory_free( (void **) &new_null);
   if (nullspace_vect != NULL) ML_memory_free( (void **) &dble_buf2);

   /* ------------------------------------------------------------- */
   /* send the P rows back to its parent processor                  */
   /* ------------------------------------------------------------- */
 
   nbytes = total_send_leng * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**) &dble_buf, nbytes,"AGz");
   msgtype = 24945;
   length = sizeof(double) * nullspace_dim;
   ML_Aggregate_ExchangeData((char*)dble_buf,(char*) comm_val,
         N_neighbors, neighbors, send_leng, recv_leng,msgtype,(int)length,comm);
   for ( i = 0; i < total_send_leng; i++ )
   {
      index = send_list[i];
      if ( aggr_index[index] >= aggr_count )
      {
         dcompare1 = 0.0;
         for ( j = 0; j < nullspace_dim; j++ )
         {
            index4 = i * nullspace_dim + j;
            dcompare1 += dble_buf[index4];
         }
         if ( dcompare1 != 0.0 )
         {
            index4 = i * nullspace_dim;
            k      = index * nullspace_dim;
            for ( j = 0; j < nullspace_dim; j++ )
            {
               /*new_val[k+j] = dble_buf[index4+j];*/
               new_val[ new_ia[index]+j] = dble_buf[index4+j];
               /* new_ja = column index */
               /*new_ja[k+j]  = aggr_index[index]*nullspace_dim+j;*/
               new_ja[ new_ia[index]+j ]  = aggr_index[index]*nullspace_dim+j;
            }
         }
      }
   }
   ML_memory_free( (void **) &comm_val);
   ML_memory_free( (void **) &dble_buf);
 
   /* ------------------------------------------------------------- */
   /* check P (row sum = 1)                                         */
   /* ------------------------------------------------------------- */

#ifdef DDEBUG
   /* count up the number of connections/edges that leave aggregate */
   /* versus those that are within the aggregate.                   */
   good = (int *) ML_allocate(Nrows*sizeof(int));
   bad  = (int *) ML_allocate(Nrows*sizeof(int));
   for (i = 0; i < Nrows; i++) { good[i] = 0; bad[i] = 0; }
   count = 0; 
   for (i = 0; i < Nrows; i++) {
      /* figure out my aggregate */
      myagg = -1;
      for (kk = new_ia[i]; kk < new_ia[i+1]; kk++) {
         if (myagg != new_ja[kk]/6) {
            if (myagg == -1) myagg = new_ja[kk]/6;
            else {
               printf("a:something is wrong %d %d in row %d\n",
                       myagg, new_ja[kk]/6, i);
               /*exit(1);*/
            } 
          }
      }

      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);
      count2 = 0;
      for (j = 0; j < rowi_N; j++) {
         if (rowi_col[j]/num_PDE_eqns != -500 - i/num_PDE_eqns) {
            /* for each column figure out what aggregate it corresponds */
	    /* to. If it is the same as myagg, add 1 to good, otherwise */
	    /* add 1 to bad */
            curagg = -1;
	    for (kk = new_ia[rowi_col[j]]; kk < new_ia[rowi_col[j]+1]; kk++) {
	      if (curagg != new_ja[kk]/6) {
		if (curagg == -1) curagg = new_ja[kk]/6;
		else {
		  printf("b:Something is wrong %d %d in row %d\n",
			 curagg, new_ja[kk]/6, rowi_col[j]);
		  /*exit(1);*/
		} 
	      }
	    }
            if ((curagg != -1) && (myagg != -1)) {
               if (curagg == myagg) good[myagg]++;
               else bad[myagg]++;
            }
            
         }
      }
   }
   myagg = 0;
   sprintf(fname,"goodbad%d",level_count);
   fp = fopen(fname,"w");
   for (i = 0; i < Nrows; i++) { 
      if ((good[i] != 0) || (bad[i] != 0)) {
         myagg += good[i]; 
         myagg += bad[i]; 
         fprintf(fp,"%d (%d,%d)\n",i,good[i],bad[i]);
      }
   }
   fclose(fp);
   printf("total number of connections counted is %d\n",myagg);
   ML_free(bad);
   ML_free(good);
   ML_free(rowi_col); rowi_col = NULL;
   ML_free(rowi_val); rowi_val = NULL;
   allocated = 0;
#endif
   /* ------------------------------------------------------------- */
   /* set up the csr_data data structure                            */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &csr_data, sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;

   ML_Operator_Set_ApplyFuncData( *Pmatrix, nullspace_dim*Ncoarse, Nrows, 
                                  ML_EMPTY, csr_data, Nrows, NULL, 0);
   (*Pmatrix)->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   ML_memory_alloc((void**) &aggr_comm, sizeof(ML_Aggregate_Comm),"ACO");
   aggr_comm->comm = comm;
   aggr_comm->N_send_neighbors = new_N_send;
   aggr_comm->N_recv_neighbors = new_N_recv;
   aggr_comm->send_neighbors = new_send_neighbors;
   aggr_comm->recv_neighbors = new_recv_neighbors;
   aggr_comm->send_leng = new_send_leng;
   aggr_comm->recv_leng = new_recv_leng;
   aggr_comm->send_list = new_send_list;
   aggr_comm->local_nrows = Ncoarse * nullspace_dim;
   
   m = exp_Ncoarse - Ncoarse;
   ML_CommInfoOP_Generate( &((*Pmatrix)->getrow->pre_comm), 
                           ML_Aggregate_ExchangeBdry, aggr_comm, 
                           comm, Ncoarse*nullspace_dim, m*nullspace_dim);
#ifdef LEASTSQ_SERIAL
   ML_Operator_Set_Getrow((*Pmatrix), ML_EXTERNAL, Nrows, CSR_get_ones_rows);
#else
   ML_Operator_Set_Getrow((*Pmatrix), ML_EXTERNAL, Nrows, CSR_getrows);
#endif
   ML_Operator_Set_ApplyFunc((*Pmatrix), ML_INTERNAL, CSR_matvec);
   (*Pmatrix)->max_nz_per_row = 1;

#ifdef METIS_DEBUG
   printf("nullspace_dim = %d\n", nullspace_dim);
   printf("Ncoarse = %d\n", Ncoarse);
   printf("exp_Ncoarse = %d\n", exp_Ncoarse);
   /*   ML_PrintOut_MLOperator( *Pmatrix, comm ); */
#endif
   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_free(unamalg_bdry);
   ML_memory_free((void**) &comm_val);
   ML_memory_free((void**) &neighbors);
   ML_memory_free((void**) &recv_leng);
   ML_memory_free((void**) &send_leng);
   ML_memory_free((void**) &send_list);
   ML_free(aggr_index);
   ML_memory_free((void**) &aggr_stat);
   ML_memory_free((void**) &sendlist_proc);
   ML_free(aggr_cnt_array);
   for (i = 0; i < aggr_count; i++) ML_free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   ML_memory_free((void**)&qr_tmp);
   ML_memory_free((void**)&tmp_vect);
   ML_memory_free((void**)&work);

   if ( new_N_send > 0 ) 
   {
      ML_memory_free((void**) &new_send_leng);
      ML_memory_free((void**) &new_send_list);
      ML_memory_free((void**) &new_send_neighbors);
   }
   if ( N_neighbors > 0 ) 
   {
      ML_memory_free((void**) &new_recv_leng);
      ML_memory_free((void**) &new_recv_neighbors);
   }
   ML_memory_free((void**) &aggr_comm);

   aggr_curr = aggr_head;
   while ( aggr_curr != NULL ) 
   {
      supernode = aggr_curr;
      aggr_curr = aggr_curr->next;
      if ( supernode->length > 0 ) ML_free( supernode->list );
      ML_free( supernode );
   }
#if defined(OUTPUT_AGGREGATES) || defined(INPUT_AGGREGATES)
   /* Print Pmatrix*v (where v is constructed using global indices) */

   dtemp = (double *) ML_allocate(sizeof(double)*(*Pmatrix)->invec_leng);
   d2temp = (double *) ML_allocate(sizeof(double)*(*Pmatrix)->outvec_leng);
   for (i = 0; i < (*Pmatrix)->outvec_leng; i++) d2temp[i] = 0.;
   for (i = 0; i < (*Pmatrix)->invec_leng; i++)
      dtemp[i] = (double) (i+agg_offset);

   sprintf(fname,"PP%d_%d",comm->ML_mypid,level_count);
   fp = fopen(fname,"w");
   ML_Operator_Apply(*Pmatrix, (*Pmatrix)->invec_leng, dtemp, 
                     (*Pmatrix)->outvec_leng, d2temp);
   for (i = 0; i < nvertices; i++) {
#ifndef MAXWELL
#ifdef ALEGRA
      if (level_count == 1) { j = i; k = i;} 
#else
      if (level_count == 1) { j = update_index[i]; k = update[i];} 
#endif
#else
      if (level_count == 1) { j = reordered_glob_nodes[i]; k = global_node_inds[i];}
#endif /* ifndef MAXWELL */
      else                  { j = i              ; k = i+vertex_offset;}
      fprintf(fp,"PP%d(%d) (loc=%d) = %e\n",level_count,k,j, d2temp[j]);
   }
   fclose(fp);
   ML_free(dtemp);
   ML_free(d2temp);

   /*
   csr_data = (struct ML_CSR_MSRdata *) (*Pmatrix)->data;
   if (comm->ML_mypid == 1)
   {
      printf("%d : row_ptr = %d\nrow_ptr = %d\ncol = %d\nval = %e\n\n\n",
             comm->ML_mypid,
             csr_data->rowptr[14], csr_data->rowptr[15],
             csr_data->columns[14], csr_data->values[14]);
   }
   sprintf(fname,"comm%d",comm->ML_mypid);
   ML_CommInfoOP_Print((*Pmatrix)->getrow->pre_comm, fname);
   fflush(stdout);
   */

#endif
   return Ncoarse*nullspace_dim;
   
} /* ML_Aggregate_CoarsenMETIS */

/* ======================================================================== */
/*!
 \brief compute the global ID of the first aggregate (assuming a linear
 decomposition among the processes)

*/
/* ------------------------------------------------------------------------ */

int ML_DecomposeGraph_BuildOffsets( int N_parts,
				    int offsets[],
				    int N_procs )
{

  int i;
  
  if( offsets == NULL ) {
    fprintf( stderr,
	     "*ML*ERR*: array offsets set to NULL,\n"
	     "file %s, line %d\n",
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }
#ifdef ML_MPI
  MPI_Allgather(&N_parts, 1, MPI_INT,
		offsets+1, 1, MPI_INT,
		MPI_COMM_WORLD);
  offsets[0] = 0;
  for( i=2 ; i<N_procs+1 ; i++ ) {
    offsets[i] += offsets[i-1];
  }
  
#else
  offsets[0] = 0;
  offsets[1] = N_parts;
#endif

  return 0;
  
} /* ML_DecomposeGraph_BuildOffsets */

/* ======================================================================== */
/*!
 \brief find position of the maximum element in an integer vector

*/
/* ------------------------------------------------------------------------ */

static int find_max(int length, int vector[] ) 
{
  int max = -1;
  int pos = -1, i;
  
  for( i=0 ; i<length ; i++ )
    if( vector[i]>max ) {
      max = vector[i];
      pos = i;
    }
  

  return pos;
  
} /* find_max */

/* ======================================================================== */
/*!
 \brief find \c key in an integer vector

 \c I cannot use AZ_find_index because my vectors are not necessarely
 ordered in ascending order.

*/
/* ------------------------------------------------------------------------ */

static int find_index( int key, int list[], int N )
{
  int i;
  for( i=0 ; i<N ; i++ )
    if( list[i] == key )
      return i;

  return -1;

} /* find_index */

/* ======================================================================== */
/*!
 \brief check that each aggregate as at least one node

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregates_CheckAggregates( int Naggregates, int N_rows,
				   int graph_decomposition[], int mypid)
{

  int i, j, divided_aggre;
  int check[Naggregates];
  
  /* ------------------- execution begins --------------------------------- */

  for( i=0 ; i<Naggregates ; i++ )  check[i] = 0;

  for( i=0 ; i<N_rows ; i++ ) {
    j = graph_decomposition[i];
    if( j >= 0 ) check[j]++;
  }

  for( i=0 ; i<Naggregates ; i++ ) {

    if( check[i] == 0 ) {

      /* find an aggre that can give us a node, the one with
         more nodes*/
      divided_aggre = find_max( Naggregates, check );
      
      if( divided_aggre >= 0 && check[divided_aggre]>0 ) {

	j = find_index( divided_aggre, graph_decomposition, N_rows );
	if( j == -1 ) {
	  fprintf( stderr,
		   "*ML*ERR* internal error while refining METIS partition\n");
	  exit( EXIT_FAILURE );
	}
	
	graph_decomposition[j] = i;
	check[divided_aggre]--;
	check[i]++;
	
	fprintf( stderr,
		 "*ML*WRN* On proc %d, aggregate %d has no nodes\n"
		 "*ML*WRN* so I take one node from  aggregate %d\n",
		 mypid, i, divided_aggre );
      } else {
	fprintf( stderr,
		 "*ML*WRN* On proc %d, aggregate %d has no nodes\n"
		 "*ML*WRN* and I cannot decompose any aggregate...\n",
		 mypid, i );
	
      }
    }
    
  }

#ifdef METIS_DEBUG
  for( i=0 ; i<Naggregates ; i++ )
    printf("agg[%d] = %d\n",
	   i, check[i] );
#endif
  
  return 0;
  
} /* ML_Aggregates_CheckAggregates */
