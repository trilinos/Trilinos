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
#include "ml_agg_info.h"

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
					 int N_parts,
					 int graph_decomposition[],
					 char bdry_nodes[],
					 int local_or_global,
					 int offsets[],
					 int reorder_flag,
					 int current_level, int *total_nz, int * );
/*
static int find_max(int length, int vector[] );
static int find_index( int key, int list[], int N );
static int ML_Aggregates_CheckAggregates( int Naggregates, int N_rows,
					  int graph_decomposition[],
					  int mypid);
*/
static int ML_LocalReorder_with_METIS( int Nrows, int xadj[], int adjncy[] ,
					 int Nparts, idxtype part[], int level,
					 ML_Comm *comm );

#define OPTIMAL_VALUE (27*27)

#ifdef EXTREME_DEBUGGING
static int MyPID_ = 0;
void set_print(int MyPID ) 
{
  MyPID_ = MyPID;
}
#include <stdarg.h>
void print(char * str, ...) 
{
  
  va_list ArgList;
  va_start(ArgList, str);
  printf("===%d=== ", MyPID_);
  vprintf(str,ArgList);
  return;
}
#endif

int ML_Aggregate_Options_Defaults( ML_Aggregate_Options * pointer,
				   int NumLevels )
{

  int i;
  
  for( i=0 ; i<NumLevels ; i++ ) {
    pointer[i].id = ML_AGGREGATE_OPTIONS_ID;
    pointer[i].Naggregates_local = -1;
    pointer[i].Naggregates_global = -1;
    pointer[i].Nnodes_per_aggregate = -1;
    pointer[i].choice = -1;
    pointer[i].reordering_flag = ML_NO;
    pointer[i].desired_aggre_per_proc = -1;
  }

  return 0;
  
}

/* ======================================================================== */
/*!
 \brief Used to set the flag (used in Decompose_with_METIS) to compute
 the radius (based on graph information only) of each aggregate.

*/
/* ------------------------------------------------------------------------ */

int COMPUTE_GRAPH_RADIUS = ML_NO;

int ML_Get_Compute_GraphRadiusFlag() 
{
  return COMPUTE_GRAPH_RADIUS ;
}

int ML_Set_Compute_GraphRadiusFlag(int i) 
{
  COMPUTE_GRAPH_RADIUS = i;
  return 0;
}

int USE_DROPPING = ML_YES;

int ML_Aggregate_Set_UseDropping(int i)
{
  USE_DROPPING = i;
  return 0;
}

int ML_Aggregate_Get_UseDropping() 
{
  return USE_DROPPING;
}


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
	     "*ML*WRN* Nlocal has an invalid value (%d). Set to default.\n ",
	     Nnodes_per_aggre);
    Nnodes_per_aggre = ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate();
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

    ML_Aggregate_Options_Defaults( pointer, Nlevels );
    
    ag->aggr_options = (void *)pointer;
  }

  if( level >= 0 ) {
    pointer[level].Nnodes_per_aggregate = Nnodes_per_aggre;
    pointer[level].choice = ML_NUM_NODES_PER_AGGREGATE;
  } else {
    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].Nnodes_per_aggregate = Nnodes_per_aggre;
      pointer[i].choice = ML_NUM_NODES_PER_AGGREGATE;
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
	     "*ML*WRN* Nlocal has an invalid value (%d). Set to 1.\n",
	     Nlocal );
    Nlocal = 1;
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

    ML_Aggregate_Options_Defaults( pointer, Nlevels );
    
    ag->aggr_options = (void *)pointer;
  }

  if( level >= 0 ) {
    pointer[level].Naggregates_local = Nlocal;
    pointer[level].choice = ML_NUM_LOCAL_AGGREGATES;
  } else {
    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].Naggregates_local = Nlocal;
      pointer[i].choice = ML_NUM_LOCAL_AGGREGATES;
    }
  }
  
  return 0;
  
} /* ML_Aggregate_Set_LocalNumber */

int ML_Aggregate_Set_GlobalNumber( ML *ml, ML_Aggregate *ag, 
				   int level, int Nglobal  )
{

  int i;
  ML_Aggregate_Options *pointer = NULL;
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

    ML_Aggregate_Options_Defaults( pointer, Nlevels );
    
    ag->aggr_options = (void *)pointer;
  }
  
  if( level >= 0 ) {
    pointer[level].Naggregates_global = Nglobal;
    pointer[level].choice = ML_NUM_GLOBAL_AGGREGATES;
  } else {
    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].Naggregates_global = Nglobal;
      pointer[i].choice = ML_NUM_GLOBAL_AGGREGATES;
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

    ML_Aggregate_Options_Defaults( pointer, Nlevels );
    
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
#ifdef METIS_DEBUG
  FILE *fp;
  char filename[80];
#endif
  
  /* ------------------- execution begins --------------------------------- */

  t0 = GetClock();

  if ( mypid == 0 &&  ML_Get_PrintLevel() > 7 ) {
    
    printf("Entering METIS reordering (level %d)\n",
	   level );
    
  }
  
  if( Nparts < 3 ) return Nparts; /* do nothing if 2 nodes only... */
  
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
     to pick them up...
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
  
  nbytes_max = ML_gmax_int( nbytes, comm );
  nbytes = ML_gsum_int( nbytes, comm);
  
  if( mypid == 0 &&  ML_Get_PrintLevel() > 7 ) {
   
    printf("METIS reordering (level %d) estimated required mem = %d Kb\n"
	   "METIS reordering (level %d) max estimated mem = %d Kb\n",
	   level,
	   nbytes,
	   level,
	   nbytes_max );
  }
#endif
  
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
  
  /* ------------------- that's all folks --------------------------------- */

  ML_free( xadj2 ) ;
  ML_free( adjncy2 );
  ML_free( perm );
  ML_free( iperm );

  t0 = GetClock() - t0;

  if ( mypid == 0 &&  ML_Get_PrintLevel() > 8 ) {
   
    printf("METIS (level %d) Time required for reordering = %e (s)\n",
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
					 int N_parts,
					 int graph_decomposition[],
					 char bdry_nodes[],
					 int local_or_global,
					 int offsets[],
					 int reorder_flag,
					 int current_level,
					 int *total_nz,
					 int * radius)
{

  int i, j,jj,  count, count2, col;
  int Nrows, Nrows_global,NrowsMETIS, N_nonzeros, N_bdry_nodes;
  int *wgtflag=NULL, numflag, *options=NULL, edgecut;
  idxtype *xadj=NULL, *adjncy=NULL, *vwgt=NULL, *adjwgt=NULL;
  idxtype *part=NULL;
  ML_Comm * comm;
  int allocated = 0;
  int * rowi_col = NULL;
  int rowi_N;
  double * rowi_val = NULL;
  int nbytes = 0, nbytes_max = 0;
  int ok = 0;
  int * nodes_per_aggre = NULL;
  double t0;
  int * dep = NULL;
  int NcenterNodes;
  int * perm = NULL;
  char str[80];
  struct amalg_drop * temp;
  double * scaled_diag;
  
  /* ------------------- execution begins --------------------------------- */
  
  t0 = GetClock();

  sprintf( str, "METIS (level %d) :", current_level );
  
  comm = Amatrix->comm;

  /* forget dropping for a moment */

  if( ML_Aggregate_Get_UseDropping() == ML_NO ) {
    
    if( comm->ML_mypid == 0 && 2 < ML_Get_PrintLevel() ) {
      printf( "%s Warning : Dropping is not used\n", str);
    }
    
    temp = (struct amalg_drop *) Amatrix->data;
    scaled_diag  = temp->scaled_diag;
    temp->scaled_diag = NULL;
    
  }
  
  /* dimension of the problem (NOTE: only local matrices) */
  
  Nrows = Amatrix->getrow->Nrows;
  perm = (int *) malloc( sizeof(int) * Nrows );

  /* for some Epetra_matrices, N_nonzeros is set to -1.
     In this case, get all rows to allocate memory for adjncy.
     Also, define the set of boundary nodes. NOTE: the computation of
     nonzero elements is not really needed (ML_Operator usually have
     this number already compuuted. However, I still need to
     define the boundary nodes, and to handle epetra matrices.) 
     Finally, I need to compute the number of rows to give in input to
     METIS. Those do not include Dirichlet rows. */
     
  N_nonzeros = 0;
  NrowsMETIS = 0;
  for (i = 0; i < Nrows; i++) {  
    ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
  	              &rowi_N, 0);
    
    if( rowi_N <= 1 ) {
      bdry_nodes[i] = 'T';
      perm[i] = -1;
    } else {
      perm[i] = NrowsMETIS++;
      bdry_nodes[i] = 'F';
      N_nonzeros += rowi_N;
    }
  }

  N_bdry_nodes = ML_Comm_GsumInt(comm, Nrows-NrowsMETIS);
  Nrows_global = ML_Comm_GsumInt(comm, Nrows);
  
  if( comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
    printf("%s # bdry (block) nodes = %d, # (block) nodes = %d\n",
	   str,
	   N_bdry_nodes, Nrows_global);
  }
  
  /* construct the CSR graph information of the LOCAL matrix
     using the get_row function */

  wgtflag = (idxtype *) malloc (4*sizeof(idxtype));
  options = (int *)     malloc (4*sizeof(int));
  
  /* set parameters */
   
  wgtflag[0] = 0;    /* no weights */
  numflag    = 0;    /* C style */
  options[0] = 0;    /* default options */
   
  xadj    = (idxtype *) malloc ((NrowsMETIS+1)*sizeof(idxtype));
  adjncy  = (idxtype *) malloc ((N_nonzeros)*sizeof(idxtype));
   
  if(  xadj==NULL || adjncy==NULL ) {
    fprintf( stderr,
	     "on proc %d, not enought space for %d bytes.\n"
	     "file %s, line %d\n",
	     comm->ML_mypid, N_nonzeros,
	     __FILE__,
	     __LINE__);
  }
   
  count = 0; count2 = 0; xadj[0] = 0;
  
  for (i = 0; i < Nrows; i++) {

    if( bdry_nodes[i] == 'F' ) {

      xadj[count2+1] = xadj[count2]; /* nonzeros in row i-1 */
    
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
		        &rowi_N, 0);

      /* need to avoid boundary nodes in METIS vectors. Skip them */
      /* (I am not pretty sure that rows with zero elements are   */
      /* well eated by METIS.) perm has been allocates of size    */
      /* Nrows, so columns corresponding to external nodes can not*/
      /* be given as input to perm                                */

      for( j=0 ; j<rowi_N ; j++ ) {
	jj = rowi_col[j];
	if( jj<Nrows ) {
	  if( jj != i && perm[jj] != -1 ) {
	    adjncy[count++] = perm[jj];
	    xadj[count2+1]++;
	  }
	}
      }
      count2++;
    }      
  }

  *total_nz = count;

  if( ML_Aggregate_Get_UseDropping() == ML_NO ) temp->scaled_diag = scaled_diag;
  
#ifdef DUMP_MATLAB_FILE
      sprintf( str, "METIS_proc%d.m", comm->ML_mypid);
      fp = fopen(str,"w");
      fprintf(fp,"NrowsMETIS = %d;\n", NrowsMETIS);
      fprintf(fp,"xadj = zeros(NrowsMETIS,1);\n");
      for( i=0 ; i<NrowsMETIS+1 ; i++ ) {
	fprintf(fp,"xadj(%d) = %d;\n", i+1, xadj[i]+1);
      }
      fprintf(fp,"Nonzeros = %d;\n", count);
      fprintf(fp,"adjncy = zeros(Nonzeros,1);\n");
      for( i=0 ; i<count ; i++ ) {
	fprintf(fp,"adjncy(%d) = %d;\n", i+1, adjncy[i]+1);
      }
      fprintf(fp,"A = zeros(%d,%d)\n", NrowsMETIS,NrowsMETIS);
      for( i=0 ; i<NrowsMETIS ; i++ ) {
	for( j=xadj[i] ; j<xadj[i+1] ; j++ ) {
	  fprintf(fp,"A(%d,%d) = 1;\n",
		  i+1,adjncy[j]+1);
	}
      }
      fclose(fp);
#endif

#ifdef DUMP_WEST
      sprintf( str, "METIS_proc%d.m", comm->ML_mypid);
      fp = fopen(str,"w");
      fprintf(fp,"Nrows = %d\n", NrowsMETIS);
      for( i=0 ; i<NrowsMETIS+1 ; i++ ) {
	fprintf(fp,"%d\n", xadj[i]);
      }
      fprintf(fp,"Nonzeros = %d\n", count);
      for( i=0 ; i<count ; i++ ) {
	fprintf(fp,"%d\n", adjncy[i]);
      }
      fclose(fp);
#endif
      
  if( count > N_nonzeros || count2 != NrowsMETIS ) {
    fprintf( stderr,
	     "*ML*WRN* On proc %d, count  > N_nonzeros (%d>%d)\n"
	     "*ML*WRN* and count2 != NrowsMETIS (%d>%d)\n"
	     "a buffer overflow has probably occurred...\n",
	     comm->ML_mypid, count, N_nonzeros, count2, NrowsMETIS );
  }

  /* idxtype is by default int, but on some architectures can be
     slightly different (for instance, a short int). */
   
  part = (idxtype *) malloc( sizeof(idxtype) * NrowsMETIS );
  nodes_per_aggre  = (int *) malloc( sizeof(int) * N_parts );

  /* ********************************************************************** */
  /* Before calling METIS, I verify that the two extreme situations are     */
  /* handled separately.                                                    */
  /* ********************************************************************** */
  
  if( N_parts == 1 ) {

    for( i=0 ; i<NrowsMETIS ; i++ ) part[i] = 0;
    edgecut = 0;
    
  } else if( N_parts == NrowsMETIS ) {

    fprintf( stderr,
	     "*ML*WRN*: on proc %d, N_part == N_rows_noDirichlet (%d==%d)\n",
	     comm->ML_mypid, N_parts, NrowsMETIS );
 
    for( i=0 ; i<NrowsMETIS ; i++ ) part[i] = i;
    edgecut = 0;
  
  } else {

    ok = 0;

    while( ok == 0 ) {
      
      /* ****************************************************************** */
      /* Put -1 in part, so I can verify that METIS has filled each pos    */
      /* ****************************************************************** */

      for( i=0 ; i<NrowsMETIS ; i++ ) part[i] = -1;
    
      /* ****************************************************************** */
      /* Estimate memory required by METIS. This memory will be dynamically */
      /* allocated inside; however this is a good estimation of how METIS   */
      /* will cost in terms of memory.                                      */
      /* Then, call METIS.                                                  */
      /* ****************************************************************** */

      if( N_parts < 8 ) {

	i = 1; /* optype in the METIS manual */
	numflag = 0;
#ifdef HAVE_ML_METIS
	METIS_EstimateMemory( &NrowsMETIS, xadj, adjncy, &numflag,
			      &i, &nbytes );
	
	METIS_PartGraphRecursive (&NrowsMETIS, xadj, adjncy, vwgt, adjwgt,
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
	METIS_EstimateMemory( &NrowsMETIS, xadj, adjncy, &numflag,
			      &i, &nbytes );

	METIS_PartGraphKway (&NrowsMETIS, xadj, adjncy, vwgt, adjwgt,
			     wgtflag, &numflag, &N_parts, options,
			     &edgecut, part);
#else
	if( Amatrix->comm->ML_mypid == 0 ) {
	  fprintf( stderr,
		   "*ML*WRN* This function has been compiled without the configure\n"
		   "*ML*WRN* option --with-ml_metis or --with-ml_metis\n"
		   "*ML*WRN* I will put all the nodes in the same aggregate, this time...\n"
		   "*ML*WRN* (file %s, line %d)\n",
		   __FILE__,
		   __LINE__);
	}
	for( i=0 ; i<NrowsMETIS ; i++ ) part[i] = 0;
	N_parts = 1;
#endif
      }
      
      /* **************************************************************** */
      /* perform some checks. If aggregates with zero assigned nodes      */
      /* exist, then recall METIS, asking for a smaller number of sub     */
      /* graphs. This is the role of the `ok' variable.                   */
      /* Also, if the part vector contains some junk, recall METIS        */
      /* **************************************************************** */

      ok = 1;
      
      for( i=0 ; i<N_parts ; i++ ) nodes_per_aggre[i] = 0;
      for( i=0 ; i<NrowsMETIS ; i++ ) {
	j = part[i];
	if( j<0 || j>= N_parts ) {
	  ok = 0;
	  break;
	} 
	else nodes_per_aggre[j]++;
      }
      
      for( i=0 ; i<N_parts ; i++ ) {
	if( nodes_per_aggre[i] == 0 ) {
	  ok = 0;
	  break;
	}
      }
      
      if( ok == 0 ) {
	if( comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	  printf( "*ML*WRN* input # of (block) aggregates (%d) does not assure "
		  "non-empty aggregates.\n"
		  "*ML*WRN* Now recalling METIS with # aggregates = %d\n",
		  N_parts, N_parts/2 );
	}
	N_parts = N_parts/2;
      }
      
      if( N_parts == 0 ) {
	if( comm->ML_mypid == 0 && 9 < ML_Get_PrintLevel()) {
	  fprintf( stderr,
		   "*ML*WRN* something went **VERY** wrong in calling METIS\n"
		   "*ML*WRN* try to ask for a smaller number of subdomains\n"
		   "*ML*WRN* I will put all the nodes into one aggregate...\n"
		   "*ML*WRN* (file %s, line %d)\n",
		   __FILE__,
		   __LINE__ );
	}
	N_parts = 1;
      }
      
      /* ************************************************************** */
      /* handle the case N_parts = 1 separately. Do not recall METIS    */
      /* in this case, simply put everything to zero and continue       */
      /* ************************************************************** */
      
      if( N_parts == 1 ) {
	for( i=0 ; i<NrowsMETIS ; i++ ) part[i] = 0;
	ok = 1;
      }
      
    } /* while( ok == 0 ) */
  
  } /* if( N_parts == 1 ) */

  /* ********************************************************************** */
  /* Some fancy output for memory usage.                                    */
  /* ********************************************************************** */

  nbytes /= 1024;
  
  nbytes_max = ML_gmax_int( nbytes, comm );
  nbytes = ML_gsum_int( nbytes, comm);

  if( Amatrix->comm->ML_mypid == 0 &&  ML_Get_PrintLevel() > 7 ) {
   
    printf("%s Estimated required mem for METIS = %d Kb\n"
	   "%s Max estimated mem for METIS = %d Kb\n",
	   str,
	   nbytes,
	   str,
	   nbytes_max );
  }
  
  /* ********************************************************************** */
  /* reordering using METIS to minimize the fill-in during factorization    */
  /* ********************************************************************** */

  if( reorder_flag == ML_YES ) {
    
    ML_LocalReorder_with_METIS( NrowsMETIS, xadj, adjncy ,
				N_parts,  part, current_level, comm );
    
  }

  /* copy back part into aggr_index, and set to -1
     the aggr_index corresponding to ghost nodes */

  for( i=0 ; i<Nrows ; i++ ) {
    j = perm[i];
    if( j != -1 ) 
      graph_decomposition[i] = (int)part[j];
    else
      graph_decomposition[i] = -1;
  }

  /* if global indices are required, modify the entries
     of graph_decomposition (only the LOCAL entries) so that
     they correspond to global indices. Also, set the array
     offsets, defined so that the global indices assigned
     to processor i are
     offsets[i] <= indices_of_i < offsets[i+1]
     Note that I do not suppose that N_parts is the same
     value among all the processors */
     
  if( local_or_global == ML_GLOBAL_INDICES ) {
    ML_DecomposeGraph_BuildOffsets( N_parts, offsets, comm->ML_nprocs,
				    Amatrix->comm->USR_comm );
  }

  /* ********************************************************************** */
  /* Compute the number of edges one has to walk though to move from the    */
  /* "center" of each aggregate up to the boundaries. This is done by the   */
  /* function `ML_Compute_AggregateGraphRadius'. This works on CSR format   */
  /* (as METIS). I look for nodes on the boundaries of all the aggregates,  */
  /* and mark them in the `dep' vector with 0. All the other nodes (in the  */
  /* interior) will be marked -7). I will give                              */
  /* the entire graph in input; in output, we will have the max radius.     */
  /* ********************************************************************** */

  if( ML_Get_Compute_GraphRadiusFlag() == ML_YES ) {
  
    dep = (int *) malloc(sizeof(int) * Nrows );
    for( i=0 ; i<NrowsMETIS ; i++ ) dep[i] = -7;
    
    for( i=0 ; i<NrowsMETIS ; i++ ) {
      jj = graph_decomposition[i];
      ok = 0;
      for( j=xadj[i] ; j<xadj[i+1] ; j++ ) {
	col = adjncy[j];
	if( graph_decomposition[col] != jj ) {
	  dep[i] = 0;
	  break;
	}
      }
    }
        
    ML_Compute_AggregateGraphRadius( NrowsMETIS, xadj, adjncy, dep,
				     radius, &NcenterNodes );
    
    ML_free( dep );
    
  }

  /* ------------------- that's all folks --------------------------------- */

  ML_free(rowi_col); ML_free(rowi_val);
  rowi_col = NULL; rowi_val = NULL;
  allocated = 0; 

  if( options != NULL ) ML_free( options );
  if( wgtflag != NULL ) ML_free( wgtflag );
  if( adjncy != NULL  ) ML_free( adjncy  );
  if( xadj != NULL    ) ML_free( xadj    );
  if( part != NULL    ) ML_free( part    );
  if( perm != NULL    ) ML_free( perm    );
  if( nodes_per_aggre != NULL ) ML_free( nodes_per_aggre );
  
  t0 = GetClock() - t0;

  if ( comm->ML_mypid == 0 &&  ML_Get_PrintLevel() > 7 ) {
   
    printf("%s Time to partition graph = %e (s)\n",
	   str,
	   t0 );
    
  }
  
  return N_parts;
  
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
   int     i, j,  k, Nrows, exp_Nrows;
   int     diff_level;
   double  printflag;
   int     aggr_count, index, mypid, num_PDE_eqns;
   int     *aggr_index = NULL, nullspace_dim;
   int     Ncoarse, count;
   int     *new_ia = NULL, *new_ja = NULL, new_Nrows;
   int     exp_Ncoarse;
   int     *aggr_cnt_array = NULL;
   int     level, index3, max_agg_size;
   int     **rows_in_aggs = NULL, lwork, info;
   double  *new_val = NULL, epsilon;
   double  *nullspace_vect = NULL, *qr_tmp = NULL;
   double  *tmp_vect = NULL, *work = NULL, *new_null = NULL;
   ML_SuperNode          *aggr_head = NULL, *aggr_curr, *supernode;
   struct ML_CSR_MSRdata *csr_data;
   int                   total_nz = 0;
   char str[80];
   
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
 int optimal_value;
 char * unamalg_bdry = NULL;
 int radius;
 
#ifdef EXTREME_DEBUGGING
 set_print(comm->ML_mypid);
#endif
 
 /* ------------------- execution begins --------------------------------- */

 sprintf( str, "METIS (level %d) :", ml_ag->cur_level );

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
 
 if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s num PDE eqns = %d\n",
	    str,
	    num_PDE_eqns);
 }
 
   /* ============================================================= */
   /* check the system size versus null dimension size              */
   /* ============================================================= */

#ifdef EXTREME_DEBUGGING
   print("# rows orig = %d, # PDE eqns = %d\n", Nrows, num_PDE_eqns);
#endif
 
   if ( Nrows % num_PDE_eqns != 0 )
   {
      printf("ML_Aggregate_CoarsenMIS ERROR : Nrows must be multiples");
      printf(" of num_PDE_eqns.\n");
      exit(EXIT_FAILURE);
   }
   diff_level = ml_ag->max_levels - ml_ag->cur_level - 1;
   if ( diff_level > 0 ) num_PDE_eqns = nullspace_dim; /* ## 12/20/99 */

   /* ============================================================= */
   /* set up the threshold for weight-based coarsening              */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   ml_ag->curr_threshold *= 0.5;

   if ( mypid == 0 && 7 < ML_Get_PrintLevel())
   {
      printf("%s current eps = %e\n",
	     str,
	     epsilon);
   }
/*
   epsilon = epsilon * epsilon;
*/

   ML_Operator_AmalgamateAndDropWeak(Amatrix, num_PDE_eqns, epsilon);
   Nrows /= num_PDE_eqns;
   
   exp_Nrows = Nrows;

   /* ********************************************************************** */
   /* allocate memory for aggr_index, which will contain the METIS decomp.   */
   /* ********************************************************************** */

   nbytes = (Nrows*num_PDE_eqns) * sizeof(int);

   if ( nbytes > 0 ) {
     ML_memory_alloc((void**) &aggr_index, nbytes, "ACJ");
     if( aggr_index == NULL ) {
       fprintf( stderr,
		"*ML*ERR* not enough memory for %d bytes\n"
		"*ML*ERR* (file %s, line %d)\n",
		nbytes,
		__FILE__,
		__LINE__ );
       exit( EXIT_FAILURE );
     }
   }
   else              aggr_index = NULL;

   for( i=0 ; i<Nrows*num_PDE_eqns ; i++ ) aggr_index[i] = -1;
   
#ifdef EXTREME_DEBUGGING
   print("aggr_index has size %d bytes\n", nbytes);
   print("aggr_options pointer is %x\n", aggr_options);
   print("ml_ag->aggr_viz_and_stats pointer is  %x\n",
	 ml_ag->aggr_viz_and_stats );
#endif

   /* ********************************************************************** */
   /* retrive the pointer to the ML_Aggregate_Options, which contains the    */
   /* number of aggregates (or their size), as well as few options for the   */
   /* constructions.                                                         */
   /* ********************************************************************** */

   aggr_options = (ML_Aggregate_Options *)ml_ag->aggr_options;

   if( aggr_options == NULL ) {

     if( mypid == 0 && 8 < ML_Get_PrintLevel() ) {
       printf("%s Using default values\n",
	      str);
     }

     /* this value is hardwired in ml_agg_METIS.c, and can be set by      */
     /* the user with `ML_Aggregate_Set_OptimalNumberOfNodesPerAggregate' */
     
     optimal_value = ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate();
     
     aggr_count = Nrows/optimal_value;
     if( aggr_count < 1 ) aggr_count = 1;

     reorder_flag = ML_NO;
     
   } else {

     if( aggr_options[ml_ag->cur_level].id != ML_AGGREGATE_OPTIONS_ID ) {
       fprintf( stderr,
		"*ML*ERR* `ML_Aggregate_CoarsenMETIS' : wrong object\n"
		"*ML*ERR* (file %s, line %d)\n",
		__FILE__,
		__LINE__ );
       exit( EXIT_FAILURE );
     }

     /* ******************************************************************** */
     /* Retrive the user's defined choice to define the number of aggregates */
     /* For local number, it is ok.                                          */
     /* If global number of aggregates or nodes per aggregate have been      */
     /* specified, compute the local one (evenly dividing this global number)*/
     /* For those two latter cases, I suppose that the input value is the    */
     /* same on all processes (but I don't check ... )                       */
     /* ******************************************************************** */

     switch( aggr_options[ml_ag->cur_level].choice ) {

     case ML_NUM_LOCAL_AGGREGATES:

       aggr_count = aggr_options[ml_ag->cur_level].Naggregates_local;
       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Objective : %d local (block) aggregates (on proc 0)\n",
		 str,
		 aggr_count );
       }
       break;

     case ML_NUM_GLOBAL_AGGREGATES:
       
       aggr_count = aggr_options[ml_ag->cur_level].Naggregates_global;
       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Objective : %d global aggregates\n",
		 str,
		 aggr_count );
       }
       
       if( aggr_count < Nprocs ) {
	 if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	   fprintf( stderr,
		    "*ML*WRN* In CoarsenMETIS, %d global (block) aggregates are required,\n"
		    "*ML*WRN* but you have only %d processes. METIS requires at\n"
		    "*ML*WRN* one aggregate per process. Otherwise, you can use ParMETIS\n"
		    "*ML*WRN* as coarsen scheme. Now proceeding with 1 local (block) aggregate\n"
		    "*ML*WRN* (file %s, line %d)\n",
		    aggr_count,
		    Nprocs,
		    __FILE__,
		    __LINE__ );
	 }
	 aggr_count = 1;

       } else { 
       
	 mod = aggr_count % Nprocs;
     
	 aggr_count /= Nprocs;
	 if( mypid == 0 ) {
	   aggr_count += mod;
	 }
	 
	 if( aggr_count < 1 ) {
	   fprintf( stderr,
		    "*ML*WRN* something weird happened... Check the code !!\n"
		    "*ML*WRN* (file %s, line %d)\n",
		    __FILE__,
		    __LINE__ );
	   aggr_count = 1;
	 }
       }
       
       break;
       
     case ML_NUM_NODES_PER_AGGREGATE:
       
       aggr_count = aggr_options[ml_ag->cur_level].Nnodes_per_aggregate;

       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Objective : %d nodes per aggregate\n",
		 str,
		 aggr_count );
       }
       
       if( aggr_count >= Nrows) {

	 i = aggr_count;
	 
	 aggr_count = Nrows/OPTIMAL_VALUE;
	 if( aggr_count == 0 ) aggr_count = 1;
	 
	 if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	   fprintf( stderr,
		    "*ML*WRN* # (block) nodes per (block) aggregate (%d) > # (block) nodes (%d)\n"
		    "*ML*WRN* (on proc 0). Now proceeding with aggr_count = %d\n",
		    i,
		    Nrows,
		    aggr_count);
	 }
	 
       } else {

	 aggr_count = (Nrows/aggr_count);

	 if( aggr_count == 0 ) aggr_count = 1;
	 
       }

#ifdef ML_MPI
       MPI_Reduce( &aggr_count, &i, 1, MPI_INT, MPI_SUM, 0,
		   comm->USR_comm);
#else
       i = aggr_count;
#endif
       
       if ( mypid == 0 && 7 < ML_Get_PrintLevel() )  {
	 printf("%s avg %f (block) aggr/process\n",
		str,
		1.0*i/Nprocs );
       }
       
       break;
       
     } /* switch */
       
     reorder_flag = aggr_options[ml_ag->cur_level].reordering_flag;
     
   } /* if( aggr_options == NULL )*/
   
   if( aggr_count<=0 ) {
     fprintf( stderr,
	      "*ML*ERR* on proc %d, value of aggr_count not correct (%d)\n"
	      "*ML*ERR* Set this value using ML_Aggregate_Set_LocalNumber\n"
	      "*ML*ERR* or ML_Aggregate_Set_NodesPerAggr or ML_Aggregate_Set_GlobalNumber\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      mypid,
	      aggr_count,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }

   /* ********************************************************************** */
   /* to call METIS, we have to create a graph using CSR data format         */
   /* Essentially, this requires a matrix-vector product (on the Amalgamated */
   /* matrix, so with dropped elements)                                      */
   /* ********************************************************************** */

   if( ml_ag->aggr_viz_and_stats != NULL ) {
     radius = -1;
     ML_Set_Compute_GraphRadiusFlag(ML_YES);     
   }
   
   unamalg_bdry = (char *) ML_allocate( sizeof(char) * (Nrows+1) );

   if( unamalg_bdry == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* on proc %d, not enough space for %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      mypid,
	      sizeof(char) * Nrows,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }

#ifdef EXTREME_DEBUGGING
   print("# requested local aggregates = %d, # rows_METIS = %d\n",
	 aggr_count,
	 Nrows );
#endif

   aggr_count = ML_DecomposeGraph_with_METIS( Amatrix,aggr_count,
					      aggr_index, unamalg_bdry,
					      ML_LOCAL_INDICES, NULL,
					      reorder_flag, ml_ag->cur_level,
					      &total_nz,&radius);


#ifdef ML_MPI
   MPI_Allreduce( &Nrows, &i, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
   MPI_Allreduce( &aggr_count, &j, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
#else
   i = Nrows;
   j = aggr_count;
#endif

   if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
     printf("%s Using %d (block) aggregates (globally)\n",
	    str,
	    j );
     printf("%s # (block) aggre/ # (block) rows = %8.5f %% ( = %d / %d)\n",
	    str,
	    100.0*j/i,
	    j, i);
   }
   /*
   if ( mypid == 0 && 8 < ML_Get_PrintLevel() )  {
     printf("ML_Aggregate_CoarsenMETIS (level %d) : "
	    "%d Aggregates (on proc 0)\n",
	    ml_ag->cur_level,
	    aggr_count );
   }
   */
   j = ML_gsum_int( aggr_count, comm );
   if ( mypid == 0 && 7 < ML_Get_PrintLevel() )  {
     printf("%s %d (block) aggregates (globally)\n",
	    str,
	    j );
   }   
   
   /* ********************************************************************** */
   /* I allocate room to copy aggr_index and pass this value to the user,    */
   /* who will be able to analyze and visualize this after the construction  */
   /* of the levels. This way, the only price we have to pay for stats and   */
   /* viz is essentially a little bit of memory.                             */
   /* this memory will be cleaned with the object ML_Aggregate ml_ag.        */
   /* I set the pointers using the ML_Aggregate_Info structure. This is      */
   /* allocated using ML_Aggregate_Info_Setup(ml,MaxNumLevels)               */
   /* ********************************************************************** */

   if( ml_ag->aggr_viz_and_stats != NULL ) {

     graph_decomposition = (int *)ML_allocate(sizeof(int)*(Nrows+1));
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
     aggr_viz_and_stats[ml_ag->cur_level].Amatrix = Amatrix;
     aggr_viz_and_stats[ml_ag->cur_level].graph_radius = radius;
   }

   /* ********************************************************************** */
   /* take the decomposition as created by METIS and form the aggregates     */
   /* ********************************************************************** */
   
   total_nz = ML_Comm_GsumInt( comm, total_nz);
   i = ML_Comm_GsumInt( comm, Nrows);

   if ( mypid == 0 && 7 < ML_Get_PrintLevel())
     printf("%s Total (block) nnz = %d ( = %5.2f/(block)row)\n",
	    str,
	    total_nz,1.0*total_nz/i);
   
   if ( ml_ag->operator_complexity == 0.0 ) {
      ml_ag->fine_complexity = total_nz;
      ml_ag->operator_complexity = total_nz;
   }
   else ml_ag->operator_complexity += total_nz;

   /* fix aggr_index for num_PDE_eqns > 1 */
   
   for (i = Nrows - 1; i >= 0; i-- ) {
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

#ifdef EXTREME_DEBUGGING
   print("After `ML_Operator_UnAmalgamateAndDropWeak'\n");
#endif   
   
   Nrows      *= num_PDE_eqns;
   exp_Nrows  *= num_PDE_eqns;

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
   for (i = 0; i <Nrows; i++) {
     if ( fscanf(fp,"%d%d",&j,&k) != 2) break;
      aggr_index[j] = k;
      if (k >= aggr_count) aggr_count = k+1;
   }
   fclose(fp);
   printf("Read in %d aggregates\n\n",aggr_count);
#endif /*ifdef ML_AGGR_INAGGR*/
   
   /* count the size of each aggregate */

   aggr_cnt_array = (int *) ML_allocate(sizeof(int)*(aggr_count+1));
   for (i = 0; i < aggr_count ; i++) aggr_cnt_array[i] = 0;
   for (i = 0; i < exp_Nrows; i++) {
     if (aggr_index[i] >= 0) {
       if( aggr_index[i] >= aggr_count ) {
	 fprintf( stderr,
		  "*ML*WRN* on process %d, something weird happened...\n"
		  "*ML*WRN* node %d belong to aggregate %d (#aggr = %d)\n"
		  "*ML*WRN* (file %s, line %d)\n",
		  comm->ML_mypid,
		  i,
		  aggr_index[i],
		  aggr_count,
		  __FILE__,
		  __LINE__ );
       } else {
	 aggr_cnt_array[aggr_index[i]]++;
       }
     }
   }

#ifdef ML_AGGR_OUTAGGR
   sprintf(fname,"agg%d_%d",comm->ML_mypid,level_count);
   fp = fopen(fname,"w");
   agg_offset = ML_gpartialsum_int(aggr_count, comm);
   vertex_offset = ML_gpartialsum_int(Nrows, comm);
   for (i = 0; i < Nrows ; i++)
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
   for (i = 0; i < Nrows; i++) dtemp[i] = (double) (i + vertex_offset);
   ML_exchange_bdry(dtemp,Amatrix->getrow->pre_comm, Nrows, comm,
                    ML_OVERWRITE,NULL);
   for (i = 0; i < exp_Nrows-Nrows; i++)
   {
#ifndef MAXWELL
#ifdef ALEGRA
      if (level_count == 0) { j = i+Nrows; k =  (int) dtemp[i+Nrows];}
#else
      if (level_count == 0) { j = extern_index[i]; k = external[i];} 
#endif /*ifdef ALEGRA*/
#else
      if (level_count == 0) { j = reordered_node_externs[i]; k =
	  global_node_externs[i];}
#endif /* ifndef MAXWELL */

      else                 { j = i+Nrows    ; k = (int) dtemp[i+Nrows];}
      if (aggr_index[j] >= 0)
         fprintf(fp,"%5d %5d\n", k, aggr_index[j]+agg_offset);
   }
   fclose(fp);
   level_count++;
   ML_free(dtemp);
#else
#ifdef INPUT_AGGREGATES
   agg_offset = ML_gpartialsum_int(aggr_count, comm);
   vertex_offset = ML_gpartialsum_int(Nrows, comm);
#endif
#endif /*ifdef ML_AGGR_OUTAGGR*/

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */

#ifdef EXTREME_DEBUGGING
   print("Form tentative prolongatoe\n");
#endif

   Ncoarse = aggr_count;
   
   /* ============================================================= */
   /* check and copy aggr_index                                     */
   /* ------------------------------------------------------------- */

   level = ml_ag->cur_level;
   nbytes = (Nrows+1) * sizeof( int );
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
   /* set up the new operator                                       */
   /* ------------------------------------------------------------- */

#ifdef EXTREME_DEBUGGING
   print("nullspace_dim = %d\n", nullspace_dim );
   print("Nrows = %d\n", Nrows );
#endif
   
   new_Nrows = Nrows;
   exp_Ncoarse = Nrows;
   
   for ( i = 0; i < new_Nrows; i++ ) 
   {
      if ( aggr_index[i] >= exp_Ncoarse ) 
      {
         printf("*ML*WRN* index out of bound %d = %d(%d)\n",
		i, aggr_index[i], 
                exp_Ncoarse);
      }
   }
   nbytes = ( new_Nrows+1 ) * sizeof(int); 
   ML_memory_alloc((void**)&(new_ia), nbytes, "AIA");
   nbytes = ( new_Nrows+1)  * nullspace_dim * sizeof(int); 
   ML_memory_alloc((void**)&(new_ja), nbytes, "AJA");
   nbytes = ( new_Nrows+1)  * nullspace_dim * sizeof(double); 
   ML_memory_alloc((void**)&(new_val), nbytes, "AVA");
   for ( i = 0; i < new_Nrows*nullspace_dim; i++ ) new_val[i] = 0.0;
   
   /* ------------------------------------------------------------- */
   /* set up the space for storing the new null space               */
   /* ------------------------------------------------------------- */

   nbytes = (Ncoarse+1) * nullspace_dim * nullspace_dim * sizeof(double);
   ML_memory_alloc((void**)&(new_null),nbytes,"AGr");
   if( new_null == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* on process %d, not enough memory for %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      mypid,
	      nbytes,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   
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
   for (i = 0; i < aggr_count; i++) {
     nbytes = aggr_cnt_array[i]+1;
     rows_in_aggs[i] = (int *) ML_allocate(nbytes*sizeof(int));
     aggr_cnt_array[i] = 0;
     if (rows_in_aggs[i] == NULL)  {
       printf("*ML*ERR* couldn't allocate memory in CoarsenMETIS\n");
       exit(1);
     }
   }
   for (i = 0; i < exp_Nrows; i+=num_PDE_eqns) {
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
              if ( unamalg_bdry[index/num_PDE_eqns] == 'T')
		qr_tmp[k*length+j] = 0.;
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
               if ( unamalg_bdry[index/num_PDE_eqns] == 'T')
		 qr_tmp[k*length+j] = 0.;
               else {
                  if (index < Nrows) {
                     qr_tmp[k*length+j] = nullspace_vect[k*Nrows+index];
                  }
                  else {
		    fprintf( stderr,
			     "*ML*ERR* in QR\n"
			     "*ML*ERR* (file %s, line %d)\n",
			     __FILE__,
			     __LINE__ );
		    exit( EXIT_FAILURE );
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
		fprintf( stderr,
			 "*ML*ERR* in QR: index out of bounds (%d - %d)\n",
			 index,
			 Nrows );
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

   /* ------------------------------------------------------------- */
   /* set up the csr_data data structure                            */
   /* ------------------------------------------------------------- */

#ifdef EXTREME_DEBUGGING
   print("set up the csr_data data structure\n");
#endif
   
   ML_memory_alloc((void**) &csr_data, sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;

   ML_Operator_Set_ApplyFuncData( *Pmatrix, nullspace_dim*Ncoarse, Nrows, 
                                  ML_EMPTY, csr_data, Nrows, NULL, 0);
   (*Pmatrix)->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;
   (*Pmatrix)->getrow->pre_comm = ML_CommInfoOP_Create();
   (*Pmatrix)->max_nz_per_row = 1;
   
#ifdef LEASTSQ_SERIAL
   ML_Operator_Set_Getrow((*Pmatrix), ML_EXTERNAL, Nrows, CSR_get_ones_rows);
#else
   ML_Operator_Set_Getrow((*Pmatrix), ML_EXTERNAL, Nrows, CSR_getrows);
#endif
   ML_Operator_Set_ApplyFunc((*Pmatrix), ML_INTERNAL, CSR_matvec);
   (*Pmatrix)->max_nz_per_row = 1;

   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */
#ifdef EXTREME_DEBUGGING
   print("clean up\n");
#endif
   
   ML_free(unamalg_bdry);
   ML_memory_free((void**)&aggr_index);
   ML_free(aggr_cnt_array);
   for (i = 0; i < aggr_count; i++) ML_free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   ML_memory_free((void**)&qr_tmp);
   ML_memory_free((void**)&tmp_vect);
   ML_memory_free((void**)&work);

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
   for (i = 0; i < Nrows; i++) {
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
   if (comm->ML_myp`id == 1)
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
				    int N_procs,
				    USR_COMM comm)
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
		comm);
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

#ifdef LATER

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

static int ML_Aggregates_CheckAggregates( int Naggregates, int N_rows,
					  int graph_decomposition[], int mypid)
{

  int i, j, divided_aggre;
  int* check = (int *)malloc( sizeof(int) * Naggregates );
  
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

  ML_free( check ); check=NULL;  
  return 0;
  
} /* ML_Aggregates_CheckAggregates */
#endif
