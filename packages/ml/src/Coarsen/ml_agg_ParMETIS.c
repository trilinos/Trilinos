/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to create tentative prolongators using ParMETIS                 */
/*                                                                           */
/* NOTE ABOUT METIS AND PARMETIS: you can compile using ParMETIS version 2.x */
/* or 3.x. In the former case, run configure using the option                */
/*   --with-ml_partmetis2x,                                                  */
/* otherwise using                                                           */
/*   --with-ml_parmetis3x.                                                   */
/* Note that the two versions of ParMETIS seems to require the METIS library */
/* included in the distribution of ParMETIS, so you may need to change the   */
/* library you are linking while switching from version 2.x to 3.x.          */
/* A typical configure script may read as follows:                           

./configure \
 --prefix=/home/msala/Trilinos/LINUX_MPI/ \
 --enable-mpi --with-mpi-compilers \
 --enable-ml_timing \
 --enable-ml_flops \
 --with-ml_superlu \
 --with-ml_metis \
 --with-ml_parmetis3x \
 --with-ldflags="-L/home/msala/lib -lsuperlu -lparmetis-3.1 -lmetis-4.0" \
 --with-incdirs="-I/home/msala/Trilinos3PL/DSuperLU/SRC -I/home/msala/include/"

                                                                             */

/* ************************************************************************* */
/* Author        : Marzio Sala (SNL)                                         */
/* Date          : November 2003                                             */
/* ************************************************************************* */
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
#include "ml_agg_ParMETIS.h"

static int ML_BuildReorderedDecomposition( int starting_decomposition[],
					   int reordered_decomposition[],
					   int Nrows, int Naggregates,
					   int nodes_per_aggre[],
					   USR_COMM comm );
static int ML_DecomposeGraph_with_ParMETIS( ML_Operator *Amatrix,
					    int N_parts,
					    int graph_decomposition[],
					    double bdry_nodes[],
					    int, int);
static int ML_CountNodesPerAggre(int Nrows, int GraphDecomposition[],
					int Naggre, int * NnodesPerAggre,
					USR_COMM Comm);
extern ML_Operator * ML_BuildQ( int StartingNumElements,
				int ReorderedNumElements,
				int num_PDE_eqns, int,
				int * reordered_decomposition,
				double * StartingNullSpace,
				double * ReorderedNullSpace,
				int ComputeNewNullSpace,
				double * StartingBdry, double * ReorderedBdry,
				USR_COMM mpi_communicator,
				ML_Comm *ml_communicator );
extern void ML_DestroyQ(void);

/* ********************************************************************** */
/* parmetis.h is required to properly define idxtype, and to declare the  */
/* used functions. By default, idxtype is defined as int, so in principle */
/* one can compile also without including metis.h.                        */
/* As ParMETIS requires MPI, I suppose that the user will never call this */
/* function without mpi enabled, Should he do this, I simply put all the  */
/* nodes in the same aggregate and print a warning message. The same is   */
/* done is parmetis has not been linked (with mpi enabled).               */
/* ********************************************************************** */

#if defined(HAVE_ML_PARMETIS_2x) || defined(HAVE_ML_PARMETIS_3x)
#include "parmetis.h"
#else
#define idxtype int
#endif

static int OPTIMAL_VALUE = 27*27; /* don't ask me why */

static int PARMETIS_DEBUG_LEVEL = 0;
static int OPTIMAL_LOCAL_COARSE_SIZE = 128;

/* ======================================================================== */
/*!

 \brief Set the debug level for ParMETIS aggregation.

 Set the debug level for ParMETTIS aggregation related functions.
 Possible values:
 - 0 : no debug checks
 - 1 : debug checks, output only for warnings;
 - 2 : debug checks, output for each debug action + warnings;
 - 3 : as level 2, plus output each time the code enters and exits
       an "important" function.

  For multi-process runs, it is suggested to set the output level to 0
  on all processes, except one.
  
*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_Set_ParMETISDebugLevel(int level)
{
  PARMETIS_DEBUG_LEVEL = level;
  return 0;
}


/* ======================================================================== */
/*!

 \brief Set the optimal number of nodes per aggregates (it is
 supposed to be constant for all levels).

 \note I perform no checks to vereify whether OPTIMAL_VALUE has the same
 value among all the processes or not. It is assumed to be.
*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_Set_OptimalNumberOfNodesPerAggregate( int optimal_value ) 
{

   if( PARMETIS_DEBUG_LEVEL == 3 ) {
     printf("*ML*DBG* Entering ML_Aggregate_Set_OptimalNumberOfNodesPerAggregate'\n");
     printf("*ML*DBG* with input value %d\n", optimal_value);
   }
  
  if( optimal_value <= 1 && 9 < ML_Get_PrintLevel() ) {
    fprintf( stderr,
	     "*ML*WRN* invalid parameter for OptimalValue (%d)\n"
	     "*ML*WRN* (must be greater than 1)\n",
	     optimal_value );
  }
  
  OPTIMAL_VALUE = optimal_value;

  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Exiting ML_Aggregate_Set_OptimalNumberOfNodesPerAggregate'\n");
  }
  
  return 0;
  
} /* ML_Aggregate_Set_OptimalNumberOfNodesPerAggregate */

int ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate( ) 
{
  return OPTIMAL_VALUE;
}

int ML_Aggregate_Set_OptimalReqLocalCoarseSize( int value ) 
{
  OPTIMAL_LOCAL_COARSE_SIZE = value;
  return 0;
}

int ML_Aggregate_Get_OptimalReqLocalCoarseSize( int value ) 
{
  return( OPTIMAL_LOCAL_COARSE_SIZE );
}

/* ======================================================================== */
/*!
 \brief Set the required number of aggregates to be assigned to each processor.
 
 This function is used to set the required number of aggregates to be
 assigned to each processor. The logic behind is to put \c
 desired_aggre_per_proc if possible, so that the workload on a subset of
 processors in "optimal."

 This function is used to redistributed the tentative prolongator with
 ParMETIS as coarsen scheme.

 \date Albuquerque, 12-Nov-03
 
*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_Set_ReqLocalCoarseSize( ML *ml, ML_Aggregate *ag, 
					 int level,
					 int desired_aggre_per_proc )
{

  int i;
  ML_Aggregate_Options *pointer = NULL;
  int Nlevels = ml->ML_num_levels;
  double debug_starting_time;
  
  /* ------------------- execution begins --------------------------------- */

  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Entering `ML_Aggregate_Set_ReqLocalCoarseSize'\n");
    printf("*ML*DBG* with input value level=%d, desired_aggre_per_proc=%d\n",
	   level, desired_aggre_per_proc);
    debug_starting_time = GetClock(); 
  }
  
  /* ********************************************************************** */
  /* control on the input parameters                                        */
  /* ********************************************************************** */

  if ( ag->ML_id != ML_ID_AGGRE ) {
      printf("ML_Aggregate_Set_ReqLocalCoarseSize : wrong object. \n");
      exit(EXIT_FAILURE);
  }

  if( desired_aggre_per_proc <= 0 ) {
    fprintf( stderr,
	     "*ML*ERR* Nlocal has an invalid value (%d)\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     desired_aggre_per_proc,
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
    ML_memory_alloc((void**)&pointer, sizeof(ML_Aggregate_Options)*Nlevels,
		    "aggr_options");
    if( pointer == NULL ) {
      fprintf( stderr,
	       "*ML*ERR* not enough space to allocate %d bytes\n"
	       "*ML*ERR* (file %s, line %d)\n",
	       (int)sizeof(int)*Nlevels,
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
    pointer[level].desired_aggre_per_proc = desired_aggre_per_proc;
  } else {
    for( i=0 ; i<Nlevels ; i++ ) {
      pointer[i].desired_aggre_per_proc = desired_aggre_per_proc;
    }
  }

  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Exiting `ML_Aggregate_Set_ReqLocalCoarseSize'\n");
    printf("*ML*DBG* Total time = %e\n",  GetClock() - debug_starting_time);
  }

  return 0;
  
} /* ML_Aggregate_Set_ReqLocalCoarseSize */

/* ********************************************************************** */
/* The goal of this function is to define a reordered offset, so that the */
/* nodes corresponding to the first `desired_aggre_per_proc' are assigned */
/* to process 0, the other `desired_aggre_per_proc' to process 1, and so  */
/* on. This means that reorderd_decomposition still refers to the finer   */
/* grid; however, after coarsening, we will have the desired distribution */
/* on the coarser level.                                                  */
/* ********************************************************************** */

static int ML_BuildReorderedOffset( int starting_offset[],
				    int desired_aggre_per_proc, int Nprocs,
				    int nodes_per_aggre[], int Naggregates,
				    int reordered_offset[], int mypid ) 
{

  int i, iaggre, aggre_owner;
  int mod;
  int local_aggre = 0;
  double debug_starting_time;
  
  /* ------------------- execution begins --------------------------------- */

  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Entering `ML_BuildReorderedOffset'\n");
    printf("*ML*DBG* with input desired_aggre_per_proc=%d\n",
	   desired_aggre_per_proc);
    debug_starting_time = GetClock(); 
  }

  /* ********************************************************************** */
  /* check that desired_aggre_per_proc times number of proc is no less than */
  /* the total number of aggregates. If so, rearrange desired_aggre_per_proc*/
  /* I allways suppose that desired_aggre_per_proc * Nprocs is at least     */
  /* equal to the total number of aggregates.                               */
  /* ********************************************************************** */

  if( mypid == 0 && 8 < ML_Get_PrintLevel() ) 
    printf( "ParMETIS : Next-level matrix will have %d rows per process\n",
	    Naggregates / Nprocs+1 );
  
  if( desired_aggre_per_proc * Nprocs < Naggregates ) {
    mod = Naggregates % Nprocs;
    if( mod != 0 ) mod = 1;
    desired_aggre_per_proc = Naggregates / Nprocs+mod;
  }
  
  for( i=0 ; i<Nprocs+1 ; i++ ) {
    reordered_offset[i] = 0;
  }

  for( iaggre=0 ; iaggre<Naggregates ; iaggre++ ) {
    aggre_owner = iaggre / desired_aggre_per_proc;
    if( aggre_owner > Nprocs ) {
      fprintf( stderr,
	       "*ML*ERR* not a processor owner for aggregate %d\n"
	       "*ML*ERR* owner is %d, while Nprocs =%d\n"
	       "*ML*ERR* (file %s, line %d)\n",
	       iaggre,
	       aggre_owner, Nprocs, 
	       __FILE__,
	       __LINE__ );
      exit( EXIT_FAILURE );
    }
    reordered_offset[aggre_owner+1] += nodes_per_aggre[iaggre];
    if( aggre_owner == mypid ) local_aggre++;
  }

  for( i=2 ; i<Nprocs+1 ; i++ ) {
    reordered_offset[i] += reordered_offset[i-1];
  }

  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Exiting `ML_BuildReorderedOffset'\n");
    printf("*ML*DBG* Returning value local_aggre=%d\n", local_aggre);
    printf("*ML*DBG* Total time = %e\n",  GetClock() - debug_starting_time);
  }

  return local_aggre;
  
} /* ML_BuildReorderedOffset */

static int ML_BuildReorderedDecomposition( int starting_decomposition[],
					    int reordered_decomposition[],
					    int Nrows, int Naggregates,
					    int nodes_per_aggre[],
					    USR_COMM comm )
{

  int i, iaggre;
  int Nprocs;
  int * count = NULL,  * offset = NULL, * offset2 = NULL;
  double debug_starting_time;
  
  /* ------------------- execution begins --------------------------------- */
  
  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Entering `ML_BuildReorderedDecomposition'\n");
    printf("*ML*DBG* with Nrows=%d, Naggregates=%d\n",
	   Nrows, Naggregates);
    debug_starting_time = GetClock(); 
  }

#ifdef ML_MPI
  MPI_Comm_rank( comm, &Nprocs );
#else
  Nprocs = 1;
#endif

  /* count how many nodes belong to each aggregate (locally)                */

  count = (int *) malloc( sizeof(int) * Naggregates );
  offset = (int *) malloc( sizeof(int) * Naggregates );
  offset2 = (int *) malloc( sizeof(int) * Naggregates );
  
  for( i=0 ; i<Naggregates ; i++ ) {
    count[i] = 0;
    offset[i] = 0;
    offset2[i] = 0;
  }
  
  
  for( i=0 ; i<Nrows ; i++ ) {
    iaggre = starting_decomposition[i];
    count[iaggre]++;
  }

  /* compute the ID of the first node belonging to an aggregate             */
#ifdef ML_MPI
  MPI_Scan( count, offset, Naggregates, MPI_INT, MPI_SUM, comm);
  for( i=0 ; i<Naggregates ; i++ ) offset[i] -= count[i];
#endif  

  for( i=0 ; i<Naggregates ; i++ ) count[i] = 0, offset2[i] = 0;

  for( i=1 ; i<Naggregates ; i++ )
    offset2[i] += nodes_per_aggre[i-1] + offset2[i-1];;

  for( i=0 ; i<Nrows ; i++ ) {
    iaggre = starting_decomposition[i];
    reordered_decomposition[i] = count[iaggre] +  offset[iaggre] +
      offset2[iaggre];
    count[iaggre]++;
  }

  /* ------------------- that's all folks --------------------------------- */

  (void)free( (void *)count );
  (void)free( (void *)offset );
  (void)free( (void *)offset2 );

  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Exiting `ML_BuildReorderedDecomposition'\n");
    printf("*ML*DBG* Total time = %e\n",  GetClock() - debug_starting_time);
  }

  return 0;

} /* ML_BuildReorderedDecomposition */

/* ======================================================================== */
/*!
 \brief Reorder the local graph for the coarser level matrix using a ParMETIS
 function

 This function builds the graph of the coarser level matrix (without
 filling it with numerical values), then call ParMETIS_NodeND to compute a
 reordering which reduces the fill-in during LU factorizations. It is
 just a LOCAL reordering.
 
*/
/* ------------------------------------------------------------------------ */

static int ML_DecomposeGraph_with_ParMETIS( ML_Operator *Amatrix,
					    int N_parts,
					    int graph_decomposition[],
					    double bdry_nodes[],
					    int N_nonzeros,
					    int current_level)
{

  int i, j,jj,  count;
  int Nrows;
  int *wgtflag=NULL, numflag, *options=NULL, edgecut;
  idxtype *xadj=NULL, *adjncy=NULL;
#if defined(ML_MPI)
#if defined(HAVE_ML_PARMETIS_2x) || defined(HAVE_ML_PARMETIS_3x)
  idxtype *vwgt=NULL, *adjwgt=NULL;
#endif
#endif
  idxtype *part=NULL;
  ML_Comm * comm = Amatrix->comm;
  int allocated = 0;
  int * rowi_col = NULL;
  int rowi_N;
  double * rowi_val = NULL;
  int *global_numbering = NULL;
  int N_procs = Amatrix->comm->ML_nprocs;
  int mypid = Amatrix->comm->ML_mypid;
  int * offsets = NULL;
  idxtype * vtxdist = NULL;
  int ncon = 1;
  float * tpwgts = NULL;
  float ubvec; /* size = ncon */
  int * proc_with_parmetis = NULL;
#ifdef ML_MPI
  MPI_Group orig_group, parmetis_group;
  MPI_Comm orig_comm;
  MPI_Comm ParMETISComm;
#else
  int orig_group, parmetis_group, orig_comm, ParMETISComm;
#endif
  int N_procs_with_parmetis;
  int ok;
  int * nodes_per_aggre = NULL, * nodes_per_aggre2 = NULL;
  int skip_check = 0;
  double t0;
  double debug_starting_time;
  
  /* ------------------- execution begins --------------------------------- */

  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Entering `ML_DecomposeGraph_with_ParMETIS'\n");
    printf("*ML*DBG* with N_parts=%d, N_nonzeros=%d\n",
	   N_parts, N_nonzeros);
    debug_starting_time = GetClock();
  }

  t0 = GetClock();
  
  /* ********************************************************************** */
  /* some general variables                                                 */
  /* ********************************************************************** */
  
  Nrows = Amatrix->getrow->Nrows;
  /* Nghosts = Amatrix->getrow->pre_comm->total_rcv_length; */
#ifdef ML_MPI
  orig_comm = Amatrix->comm->USR_comm;
#endif
  
#ifdef ML_MPI
  MPI_Allreduce( &Nrows, &j, 1, MPI_INT, MPI_SUM, orig_comm );
#else
  j = Nrows;
#endif
  
  if( N_parts <= 0 ) N_parts = 1;

#ifdef MAYBE
  if( (j/OPTIMAL_VALUE>1) && (N_parts > j/OPTIMAL_VALUE) && N_parts != 1 ) {
    i = N_parts;
    N_parts = j/OPTIMAL_VALUE;
    if( N_parts == 0 ) N_parts = 1;
    if( mypid == 0 && 9 < ML_Get_PrintLevel() ) {
      printf("*ML*WRN* Changing value of N_parts from %d to %d (Nglobal = %d)\n",
	     i,
	     N_parts,
	     j );
    }
  }
#endif
  
  /* ********************************************************************** */
  /* no need to call parmetis if only one aggregate is required.            */
  /* ********************************************************************** */

  if( N_parts == 1 ) {
    for( i=0 ; i<Nrows ; i++ ) {
      graph_decomposition[i] = 0;
    }
    return 1;
  }

  /* for ParMETIS I need the global column number. Here I make use
     of the fact that ML required the ML_Operator rows to be
     decomposed linearly */

  /* allocates memory for global_ordering */

  ML_build_global_numbering(Amatrix, comm, &global_numbering);

  offsets = (int     *) malloc( sizeof(int) * (N_procs+1) );
  vtxdist = (idxtype *) malloc( sizeof(idxtype) * (N_procs+1) );
  
  ML_DecomposeGraph_BuildOffsets( Nrows, offsets, N_procs,
				  Amatrix->comm->USR_comm );

  proc_with_parmetis = (int *) malloc( sizeof(int) * N_procs );
  
  N_procs_with_parmetis = 0;
  vtxdist[0] = 0;
  for( i=1 ; i<N_procs+1 ; i++ ) {
    j = offsets[i]-offsets[i-1];
    if( j>0 ) {
      proc_with_parmetis[N_procs_with_parmetis] = i-1;  
      vtxdist[1+N_procs_with_parmetis++] = (idxtype)offsets[i];
    }
  }

  if( PARMETIS_DEBUG_LEVEL > 2 ) {
    printf("*ML*DBG* Including %d processes of %d in ParMETISComm\n",
	   N_procs_with_parmetis, N_procs );
  }
    
  if( Nrows > 0 ) {
    
    /* construct the CSR graph information of the LOCAL matrix
       using the get_row function */
    
    xadj    = (idxtype *) malloc ((Nrows+1)*sizeof(idxtype));
    adjncy  = (idxtype *) malloc ((N_nonzeros+1)*sizeof(idxtype));
    
    if(  xadj==NULL || adjncy==NULL ) {
      fprintf( stderr,
	       "on proc %d, not enought space for %d bytes.\n"
	       "file %s, line %d\n",
	       comm->ML_mypid, N_nonzeros,
	       __FILE__,
	       __LINE__);
      exit( EXIT_FAILURE );
    
    }
    
    xadj[0] = 0;
    count = 0;
    
    for (i = 0; i < Nrows; i++ ) {

      if( bdry_nodes[i] != 1.0 ) {
	
	ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
			  &rowi_N, 0);
	
	xadj[i+1] = xadj[i];
	for( j=0 ; j<rowi_N ; j++ ) {
	  jj = rowi_col[j];
	  if( jj != i && rowi_val[j] != 0.0 && bdry_nodes[jj] != 1.0 ) {
	    adjncy[count++] = global_numbering[jj];
	    xadj[i+1]++;
	  }
	} 
      }   else {
	xadj[i+1] = xadj[i];
      }   
    }
    
    if( count > N_nonzeros ) {
      fprintf( stderr,
	       "*ML_WRN* on proc %d, count > N_nonzeros (%d>%d)\n"
	       "*ML_WRN* a buffer overflow has probably occurred...\n"
	       "*ML_WRN* (file %s, line %d)\n",
	       comm->ML_mypid, count, N_nonzeros,
	       __FILE__,
	       __LINE__ );
    }
  }

  /* ********************************************************************** */
  /* idxtype is by default int, but on some architectures can be slightly   */
  /* different (for instance, a short int). For the sake of generality, here*/
  /* I put idxtype. The plus one in the allocation is to avoid cases with   */
  /* Nrows == 0 (this can happen is nodes have already been redistributed   */
  /* to a subset of processors, but that was not the last level).           */
  /* ********************************************************************** */
  
  part = (idxtype *) malloc( sizeof(idxtype) * (Nrows+1));
  tpwgts = (float *) malloc( sizeof(float) * N_parts );
  
  if( N_parts == 1 ) {

    /* should never be used here, print out something */
    puts("ehhhhhhhhhhhhhhhhhhhhhhhhhhhhh? check me !!!!");
    for( i=0 ; i<Nrows ; i++ )
      part[i] = 0;
    edgecut = 0;
    
  } else {
      
    /* set parameters for ParMETIS */
      
    wgtflag = (idxtype *) malloc (4*sizeof(idxtype));
    options = (int *)     malloc (4*sizeof(int));
      
    wgtflag[0] = 0;    /* no weights */
    numflag    = 0;    /* C style */
    options[0] = 0;    /* default options */
    
    ncon = 1;          /* number of weights of each vertex */
      
    /* fraction of vertex weight that should be distributed to
       each subdomain for each balance constraint */
      
    for( i=0 ; i<N_parts ; i++ )
      tpwgts[i] = 1.0/N_parts;
      
    /* imbalance tolerance for each vertex weight, 1 is the perfect
       balance and N_parts being perfect imbalance */
      
    ubvec = 1.05;
      
  }

  /* ********************************************************************** */
  /* ora viene il difficile... ParMETIS non sembra digerire tanto bene casi */
  /* con Nrows == 0. Quindi, devo un po` prenderlo in giro, e create un     */
  /* comunicatore che viva solo sui processi con Nrows diverso da zero.     */
  /* La chiamata a ParMETIS sara` quindi effettuata solo da questo subset   */
  /* di processori. ParMETIS mi sembra pure piu` rognoso per quanto riguarda*/
  /* il numero di partizioni. E` abbastanza facile che non tutti gli aggre  */
  /* abbiano almeno un nodo. Ho adottato la tecnica seguente: se non tutti  */
  /* gli aggre hanno un nodo, divido N_parts per 2, e ricomincio. Costa un  */
  /* po` in termini di CPU, ma e` sicuro. Il caso N_parts == 1 e` trattato  */
  /* separatamente.                                                         */
  /* ********************************************************************** */

#ifdef ML_MPI
  MPI_Comm_group( orig_comm, &orig_group );
  
  MPI_Group_incl( orig_group, N_procs_with_parmetis,
		  proc_with_parmetis, &parmetis_group );
  
  MPI_Comm_create( orig_comm, parmetis_group, &ParMETISComm );
#endif
  
  /* Only processes beloning to ParMETISComm will enter this `if'           */

  if( Nrows > 0 ) {

    ok = 0;

    while( ok == 0 ) {

      for( i=0 ; i<Nrows ; i++ ) part[i] = -7;
      skip_check = 0; 

#if defined(HAVE_ML_PARMETIS_2x) && defined(ML_MPI)
      ParMETIS_PartKway(vtxdist, xadj, adjncy, vwgt, adjwgt,
			wgtflag, &numflag, &N_parts,
			options, &edgecut, part, &ParMETISComm);
#elif defined(HAVE_ML_PARMETIS_3x)  && defined(ML_MPI)
      /* never tested right now... */
      ParMETIS_V3_PartKway (vtxdist, xadj, adjncy, vwgt, adjwgt,
			    wgtflag, &numflag, &ncon, &N_parts, tpwgts,
			    &ubvec,options,
			    &edgecut, part, &ParMETISComm);
#else
      if( Amatrix->comm->ML_mypid == 0 ) {
	fprintf( stderr,
		 "*ML*WRN* This function has been compiled without the configure\n"
		 "*ML*WRN* option --with-ml_parmetis2x or --with-ml_parmetis3x\n"
		 "*ML*WRN* I will put all the nodes in the same aggregate, this time...\n"
		 "*ML*WRN* (file %s, line %d)\n",
		 __FILE__,
		 __LINE__);
      }
      for( i=0 ; i<Nrows ; i++ ) part[i] = 0;
      edgecut = 0;
      N_parts = 1;
      ok = 1;
      skip_check = 1;
#endif
      
      ok = 1;

#ifdef ML_MPI      
      if( skip_check == 0 ) {

	/* **************************************************************** */
	/* perform some checks. If aggregates with zero assigned nodes      */
	/* exist, then recall ParMETIS, asking for a smaller number of sub  */
	/* graphs. This is the role of the `ok' variable.                   */
	/* Also, if the part vector contains some junk, recall ParMETIS     */
	/* **************************************************************** */

	nodes_per_aggre  = (int *) malloc( sizeof(int) * N_parts );
	nodes_per_aggre2 = (int *) malloc( sizeof(int) * N_parts );
	for( i=0 ; i<N_parts ; i++ ) nodes_per_aggre[i] = 0;
	for( i=0 ; i<Nrows ; i++ ) {
	  j = part[i];
	  if( j<0 || j>= N_parts ) {
	    ok = 0;
	    break;
	  } 
	  else nodes_per_aggre[j]++;
	}
	
	MPI_Allreduce( nodes_per_aggre, nodes_per_aggre2,
		       N_parts, MPI_INT, MPI_SUM, ParMETISComm );
	for( i=0 ; i<N_parts ; i++ ) {
	  if( nodes_per_aggre2[i] == 0 ) {
	    ok = 0;
	    break;
	  }
	}
	
	if( ok == 0 ) {
	  if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
	    printf( "*ML*WRN* input value of N_parts (%d) does not assure "
		    "non-empty aggregates.\n"
		    "*ML*WRN* Now recalling ParMETIS with N_parts = %d\n",
		    N_parts, N_parts/2 );
	  }
	  N_parts = N_parts/2;
	}
	
	if( N_parts == 0 ) {
	  if( mypid == 0 && 9 < ML_Get_PrintLevel()) {
	    fprintf( stderr,
		     "*ML*WRN* something went **VERY** wrong in calling ParMETIS\n"
		     "*ML*WRN* try to ask for a smaller number of subdomains\n"
		     "*ML*WRN* I will put all the nodes into one aggregate...\n"
		     "*ML*WRN* (file %s, line %d)\n",
		     __FILE__,
		     __LINE__ );
	  }
	  N_parts = 1;
	}

	/* ************************************************************** */
	/* handle the case N_parts = 1 separately. Do not recall parMETIS */
	/* in this case, simply put everything to zero and continue       */
	/* ************************************************************** */
	
	if( N_parts == 1 ) {
	  for( i=0 ; i<Nrows ; i++ ) part[i] = 0;
	  ok = 1;
	}

	if( nodes_per_aggre != NULL ) {
	  (void)free( (void *)nodes_per_aggre  );
	  nodes_per_aggre = NULL;
	}
	if( nodes_per_aggre2 != NULL ) {
	  (void)free( (void *)nodes_per_aggre2 );
	  nodes_per_aggre2 = NULL;
	}
	
      } /* if( skip_check == 1 ) */
#endif
      
    } /* while( ok == 0 ) */
#ifdef ML_MPI
    MPI_Group_free( &parmetis_group );
    MPI_Comm_free( &ParMETISComm );      
#endif
  } /* if( Nrows>0 ) */

  for( i=0 ; i<Nrows ; i++ )  graph_decomposition[i] = (int)part[i];

  /* ------------------- that's all folks --------------------------------- */

  if( global_numbering != NULL ) ML_free( global_numbering ); 
  if( rowi_col != NULL ) ML_free(rowi_col); rowi_col = NULL; 
  if( rowi_val != NULL ) ML_free(rowi_val); rowi_val = NULL;
  allocated = 0; 
  if( part               != NULL ) (void)free( (void*)part );
  if( proc_with_parmetis != NULL ) (void)free( (void *)proc_with_parmetis );
  if( offsets            != NULL ) (void)free( (void *)offsets );
  if( vtxdist            != NULL ) (void)free( (void *)vtxdist );
  if( tpwgts             != NULL ) (void)free( (void *)tpwgts );
  if( wgtflag            != NULL ) (void)free( (void *)wgtflag );
  if( options            != NULL ) (void)free( (void *)options );
  if( xadj               != NULL ) (void)free( (void *)xadj );
  if( adjncy             != NULL ) (void)free( (void *)adjncy );
  
  t0 = GetClock() - t0;

  if ( mypid == 0 &&  ML_Get_PrintLevel() > 7 ) {
   
    printf("ParMETIS (level %d) : time required = %e\n",
	   current_level,
	   t0 );
    
  }
  
  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Exiting `ML_DecomposeGraph_with_ParMETIS'\n");
    printf("*ML*DBG* With N_parts = %d\n", N_parts);
    printf("*ML*DBG* Total time = %e\n",  GetClock() - debug_starting_time);
  }

  return N_parts;
  
} /* ML_DecomposeGraph_with_ParMETIS */

/* ======================================================================== */
/*!
 \brief create non-smoothed aggregates using METIS. In order to use this
 function, the user has to define the number of aggregate (or the # of
 nodes in each aggregate) using the functions \c ML_Aggregate_Set_LocalNumber
 or \c ML_Aggregate_Set

 \note this function is derived from ML_Aggregate_CoarsenMIS, but then
 consistently modified. In particular, all the communication in the
 setup of Ptentative have been erased. Note also that this function
 CANNOT work without Epetra.

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_CoarsenParMETIS( ML_Aggregate *ml_ag, ML_Operator *Amatrix, 
				  ML_Operator **Pmatrix, ML_Comm *comm)
{
   unsigned int nbytes, length;
   int     i, j, jj, k, Nrows, exp_Nrows,  N_bdry_nodes;
   int     diff_level, Nrows_global;
   int     aggr_count, index, mypid, num_PDE_eqns;
   int     *aggr_index = NULL, nullspace_dim;
   int     Ncoarse, count;
   int     *new_ia = NULL, *new_ja = NULL, new_Nrows;
   int     *aggr_cnt_array = NULL;
   int     level, index3, max_agg_size;
   int     **rows_in_aggs = NULL, lwork, info;
   double  *new_val = NULL, epsilon;
   double  *nullspace_vect = NULL, *qr_tmp = NULL;
   double  *tmp_vect = NULL, *work = NULL, *new_null = NULL;
   ML_SuperNode          *aggr_head = NULL, *aggr_curr, *supernode;
   struct ML_CSR_MSRdata *csr_data;
   double                *starting_amalg_bdry, *reordered_amalg_bdry;
   int                   Nghost;
   int                   allocated = 0, *rowi_col = NULL, rowi_N;
   double                *rowi_val = NULL;
   int Nnonzeros2 = 0;
   int optimal_value;
   ML_Operator * Pmatrix2 = NULL;
   
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
   ML_Aggregate_Viz_Stats * aggr_viz_and_stats;
   ML_Aggregate_Options * aggr_options;
   int Nprocs;
   int * starting_offset = NULL, * reordered_offset = NULL;
   int desired_aggre_per_proc;
   int * nodes_per_aggre = NULL;
   int * starting_decomposition = NULL;
   int * reordered_decomposition = NULL;
   ML_Operator * QQ = NULL;
   ML_Operator *Pstart = NULL;
   int starting_aggr_count;
   char str[80], * str2;
   double * new_nullspace_vect = NULL;
   int * graph_decomposition = NULL;
   double debug_starting_time;
   
   /* ------------------- execution begins --------------------------------- */

   if( PARMETIS_DEBUG_LEVEL == 3 ) {
     printf("*ML*DBG* Entering `ML_Aggregate_CoarsenParMETIS'\n");
     debug_starting_time = GetClock();
   }

   sprintf( str, "ParMETIS (level %d) :", ml_ag->cur_level );
   
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

   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s num_PDE_eqns = %d\n",
	    str,
	    num_PDE_eqns);
   }

#ifdef ML_MPI
   MPI_Allreduce( &Nrows, &Nrows_global, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
#else
   Nrows_global = Nrows;
#endif
   
   /* ============================================================= */
   /* check the system size versus null dimension size              */
   /* ============================================================= */
     
   if ( Nrows % num_PDE_eqns != 0 )
   {
      printf("*ML*ERR* : Nrows must be multiples of num_PDE_eqns.\n");
      exit(EXIT_FAILURE);
   }
   diff_level = ml_ag->max_levels - ml_ag->cur_level - 1;
   if ( diff_level > 0 ) num_PDE_eqns = nullspace_dim; /* ## 12/20/99 */

   Nghost = Amatrix->getrow->pre_comm->total_rcv_length;
   
   /* ============================================================= */
   /* set up the threshold for weight-based coarsening              */
   /* ============================================================= */

   diff_level = ml_ag->begin_level - ml_ag->cur_level;
   if ( diff_level == 0 ) ml_ag->curr_threshold = ml_ag->threshold;
   epsilon = ml_ag->curr_threshold;
   ml_ag->curr_threshold *= 0.5;

   if ( mypid == 0 && 8 < ML_Get_PrintLevel() )
   {
      printf("%s current eps = %e\n",
	     str,
	     epsilon);
      if( epsilon != 0.0 ) {
	fprintf( stderr,
		 "WARNING: ParMETIS may not work with dropping!\n"
		 "WARNING: Now proceeding -- with fingers crossed\n" );
      }
   }
   
   ML_Operator_AmalgamateAndDropWeak(Amatrix, num_PDE_eqns, epsilon);
   Nrows /= num_PDE_eqns;
   Nrows_global /= num_PDE_eqns;
   Nghost /= num_PDE_eqns;
   exp_Nrows = Nrows+Nghost;

   /* record the Dirichlet boundary. unamalg_bdry[i] == 1 ==> the nodes is */
   /* a boundary node, unamalg_bdry[i] == 0 ==> the node is not            */

   nbytes = sizeof(double)*(exp_Nrows + 1);
   starting_amalg_bdry = (double *) ML_allocate(nbytes);
   if( starting_amalg_bdry == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* Not enough space to allocate %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      nbytes,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   for (i = Nrows ; i < exp_Nrows; i++) starting_amalg_bdry[i] = 0.0;

   Nnonzeros2 = 0, N_bdry_nodes = 0;
   for (i = 0; i < Nrows; i++) {
      ML_get_matrix_row(Amatrix, 1, &i, &allocated, &rowi_col, &rowi_val,
                        &rowi_N, 0);

      if (rowi_N > 1) {
        starting_amalg_bdry[i] = 0.0;
        Nnonzeros2 += rowi_N;
      } else {
	starting_amalg_bdry[i] = 1.0;
	N_bdry_nodes++;
      }
   }

   if( rowi_col != NULL ) ML_free(rowi_col );
   if( rowi_val != NULL ) ML_free(rowi_val );
   
   i = ML_Comm_GsumInt(comm, N_bdry_nodes);
   
   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s # bdry (block) nodes = %d, # (block) nodes = %d\n",
	    str,
	    i, Nrows_global);
   }
   
   /* communicate the boundary information */
   
   ML_exchange_bdry(starting_amalg_bdry,Amatrix->getrow->pre_comm,
		    Nrows,comm, ML_OVERWRITE,NULL);

   /* ********************************************************************** */
   /* allocate memory for starting_decomposition and call ParMETIS to        */
   /* decompose the local                                                    */
   /* graph into the number of parts specified by the user with a call       */
   /* ML_Aggregate_Set_LocalNumber( ml, ag, level, Nparts)                   */
   /* ********************************************************************** */

   nbytes = (Nrows+Nghost) * sizeof(int);

   if ( nbytes > 0 ) starting_decomposition = (int *)malloc(nbytes);
   else              starting_decomposition = NULL;
   
   if( starting_decomposition == NULL && nbytes > 0 ) {
     
     fprintf( stderr,
	      "*ML*ERR* not enough memory to allocated %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      nbytes,
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }

   /* ********************************************************************** */
   /* Retrive the user's defined data. If not set, pick up default settings  */
   /* ********************************************************************** */
     
   aggr_options = (ML_Aggregate_Options *) ml_ag->aggr_options;
   
   if( aggr_options == NULL ) {

     if( mypid == 0 && 8 < ML_Get_PrintLevel() ) {
       printf("%s Using default values\n",
	      str );
     }
    
     optimal_value = ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate();
     starting_aggr_count = (int)1.0*Nrows_global/optimal_value;
     if( starting_aggr_count < 1 ) starting_aggr_count = 1;
     /* 'sto xxxcio e` un po` piu` difficile, non ne ho idea, piglio
	a caso... */
     desired_aggre_per_proc = ML_max( 128, Nrows );
       
   } else {

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
       i = aggr_options[ml_ag->cur_level].Naggregates_local;
       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Requested %d local aggregates (on proc 0)\n",
		 str,
		 i );
       }
#ifdef ML_MPI
       MPI_Allreduce(&i,&starting_aggr_count,1,MPI_INT,MPI_SUM,comm->USR_comm);
#else
       starting_aggr_count = i;
#endif
       break;
       
     case ML_NUM_GLOBAL_AGGREGATES:
       
       starting_aggr_count = aggr_options[ml_ag->cur_level].Naggregates_global;
       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Requested %d global aggregates\n",
		 str,
		 starting_aggr_count );
       }

       break;
       
     case ML_NUM_NODES_PER_AGGREGATE:

       starting_aggr_count = aggr_options[ml_ag->cur_level].Nnodes_per_aggregate;

       if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	 printf( "%s Requested %d nodes per aggregate\n",
		 str,
		 starting_aggr_count );
       }
       
       if( starting_aggr_count >= Nrows_global) {

	 i = starting_aggr_count;

	 starting_aggr_count = Nrows/OPTIMAL_VALUE;
	 if( starting_aggr_count == 0) starting_aggr_count = 1;
	 
	 if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
	   fprintf( stderr,
		    "%s WARNING : # nodes per aggregate (%d) > # nodes (%d)\n"
		    "%s WARNING : now proceeding with # aggregates = %d\n",
		    str,
		    i,
		    Nrows,
		    str,
		    starting_aggr_count);
	 }

       } else {

	 starting_aggr_count = Nrows_global/starting_aggr_count;
	 if( starting_aggr_count == 0 ) starting_aggr_count = 1;

       }
       
       break;
       
     } /* switch */

     /* reorder_flag = aggr_options[ml_ag->cur_level].reordering_flag; */

     desired_aggre_per_proc = aggr_options[ml_ag->cur_level].desired_aggre_per_proc;

     if( desired_aggre_per_proc <= 0 )
       desired_aggre_per_proc = OPTIMAL_LOCAL_COARSE_SIZE;
     
   } /* if( aggr_options == NULL )*/

   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf("%s Objective : %d aggregates over %d nodes\n",
	    str,
	    starting_aggr_count,
	    Nrows_global);
     printf("%s Objective : %d aggregates on each process\n",
	    str,
	    desired_aggre_per_proc );
   } 
   
   starting_aggr_count =
     ML_DecomposeGraph_with_ParMETIS( Amatrix, starting_aggr_count,
				      starting_decomposition,
				      starting_amalg_bdry,
				      Nnonzeros2, ml_ag->cur_level );
   
   if( starting_aggr_count <= 0 ) {
     fprintf( stderr,
	      "*ML*ERR* Something went *very* wrong in ParMETIS...\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   
   if( mypid == 0 && 7 < ML_Get_PrintLevel() ) 
     printf("%s Using %d aggregates (globally)\n",
	    str,
	    starting_aggr_count );
   
   if( mypid == 0 && 7 < ML_Get_PrintLevel() ) {
     printf("%s # aggre/ # (block) rows = %7.3f %%  (= %d/%d)\n",
	    str,
	    100.0*starting_aggr_count/Nrows_global,
	    starting_aggr_count,
	    Nrows_global);
   }
   
   /* ********************************************************************** */
   /* compute operator complexity                                            */
   /* ********************************************************************** */
   
   Nnonzeros2 = ML_Comm_GsumInt( comm, Nnonzeros2);

   if ( mypid == 0 && 7 < ML_Get_PrintLevel())
     printf("%s Total (block) nnz = %d ( = %5.2f/(block)row)\n",
	    str,
	    Nnonzeros2,1.0*Nnonzeros2/Nrows_global);
   
   if ( ml_ag->operator_complexity == 0.0 ) {
      ml_ag->fine_complexity = Nnonzeros2;
      ml_ag->operator_complexity = Nnonzeros2;
   }
   else ml_ag->operator_complexity += Nnonzeros2;

   /* FIXME: erase meeeeeeeee
      fix aggr_index for num_PDE_eqns > 1 
   
   for (i = Nrows - 1; i >= 0; i-- ) {
      for (j = num_PDE_eqns-1; j >= 0; j--) {
         aggr_index[i*num_PDE_eqns+j] = aggr_index[i];
      }
   }
   */
   
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

     graph_decomposition = (int *) ML_allocate(sizeof(int)*Nrows );

     if( graph_decomposition == NULL ) {
       fprintf( stderr,
		"*ML*ERR* Not enough memory for %d bytes\n"
		"*ML*ERR* (file %s, line %d)\n",
		(int)sizeof(int)*Nrows,
		__FILE__,
	      __LINE__ );
       exit( EXIT_FAILURE );
     }

     for( i=0 ; i<Nrows ; i++ )
       graph_decomposition[i] = starting_decomposition[i];

     aggr_viz_and_stats = (ML_Aggregate_Viz_Stats *) (ml_ag->aggr_viz_and_stats);
     aggr_viz_and_stats[ml_ag->cur_level].graph_decomposition = graph_decomposition;
     aggr_viz_and_stats[ml_ag->cur_level].Nlocal = Nrows;
     aggr_viz_and_stats[ml_ag->cur_level].Naggregates = starting_aggr_count;
     aggr_viz_and_stats[ml_ag->cur_level].local_or_global = ML_GLOBAL_INDICES;
     aggr_viz_and_stats[ml_ag->cur_level].is_filled = ML_YES;
     
   }

   /* ********************************************************************** */
   /* Compute the new distribution, so that `desired_aggre_per_proc' aggre   */
   /* are stored on each processor (up to the maximum number of aggregates). */
   /* - starting_offset : decomposition of the unknowns for the finer grid   */
   /*                     before redistribution                              */
   /* - reordered_offset : decomposition of the unknowns for the finer grid  */
   /*                      as will be after redistribution, as done by       */
   /*                      operator QQ                                       */
   /* - Nrows, new_Nrows : number of local rows for the finer grid before    */
   /*                      and after redistribution                          */
   /* ********************************************************************** */

   starting_offset  = (int *)malloc( sizeof(int) * (Nprocs+1));
   reordered_offset = (int *)malloc( sizeof(int) * (Nprocs+1));
   nodes_per_aggre = (int *) malloc( sizeof(int) * starting_aggr_count );

   if( starting_offset == NULL || reordered_offset == NULL
       || nodes_per_aggre == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* Not enough memory\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }
   
   ML_DecomposeGraph_BuildOffsets( Nrows, starting_offset, Nprocs,
				   Amatrix->comm->USR_comm);
   
   /* ********************************************************************** */
   /* Compute how many nodes are contained in each aggregate. This will be   */
   /* done for all the aggregates (so some communications will occur).       */
   /* ********************************************************************** */
   
   ML_CountNodesPerAggre( Nrows, starting_decomposition,
			  starting_aggr_count, nodes_per_aggre,
			  Amatrix->comm->USR_comm );

   /* ********************************************************************** */
   /* Compute how many aggregates will be stored on this process. This is    */
   /* based on the `desired_aggre_per_proc', so that the first processes will*/
   /* have about this number (and then maybe some processes will have none). */
   /* This is used to determine a reorderd offset, so that each processor    */
   /* will hold the rows of the matrix required to form the given aggregates */
   /* This new row decomposition is hold in `reordered_decomposition'        */
   /* ********************************************************************** */

   aggr_count = ML_BuildReorderedOffset( starting_offset,
					 desired_aggre_per_proc,
					 Nprocs, nodes_per_aggre,
					 starting_aggr_count,
					 reordered_offset, mypid );

   new_Nrows = reordered_offset[mypid+1] - reordered_offset[mypid];
   
   i = 0;
   if( new_Nrows > 0 ) i++;

#ifdef ML_MPI
   MPI_Reduce( &i, &j, 1, MPI_INT, MPI_SUM, 0, comm->USR_comm);
#else
   j = i;
#endif

   if( mypid == 0 && 5 < ML_Get_PrintLevel() ) {
     printf( "%s Processes with at least 1 row at next level = %d\n",
	     str,
	     j );
   } 
   
   reordered_decomposition = (int *) malloc( sizeof(int) * (Nrows+1) );
   if( reordered_decomposition == NULL ) {
     fprintf( stderr,
	      "*ML*ERR* Not enough memory to allocate %d bytes\n"
	      "*ML*ERR* (file %s, line %d)\n",
	      (int)sizeof(int) * (Nrows+1),
	      __FILE__,
	      __LINE__ );
     exit( EXIT_FAILURE );
   }

   ML_BuildReorderedDecomposition( starting_decomposition,
				   reordered_decomposition, Nrows,
				   starting_aggr_count, nodes_per_aggre,
				   Amatrix->comm->USR_comm );
   
   /* ********************************************************************** */
   /* Finally, built the operator QQ, moving from the reorderd decomposition */
   /* to the starting one. So, QQ will be applied to vectors (or matrices)   */
   /* whose row decomposition is the reordered one. I need QQ because we have*/
   /* P_tent = QQ * \hat{P}_tent, where P_tent is the tentative prolongator  */
   /* as used by ML (from the next level to this level, in the starting      */
   /* decomposition), and \hat{P}_tent the one based on the reordered dec.   */
   /* ********************************************************************** */

   if( nullspace_vect == NULL /*&& diff_level == 0*/ ) {
     new_nullspace_vect = nullspace_vect;
     i = 0;
   } else {
     nbytes = sizeof(double) * (new_Nrows * num_PDE_eqns * nullspace_dim );

     if( nbytes == 0 ) new_nullspace_vect = NULL;
     else {
       new_nullspace_vect = (double *) malloc( nbytes );
       if( new_nullspace_vect == NULL ) {
	 fprintf( stderr,
		  "*ML*ERR* Not enough memory to allocate %d bytes\n"
		  "*ML*ERR* (file %s, line %d)\n",
		  nbytes,
		  __FILE__,
		  __LINE__ );
	 exit( EXIT_FAILURE );
       }
     }
       
     i = 1;
   }

   reordered_amalg_bdry = (double *) ML_allocate(sizeof(double)*(new_Nrows+1));
   
#ifdef ML_WITH_EPETRA
   if( nullspace_dim != num_PDE_eqns ) {
     printf("---------> Never tested with nullspace_dim != num_PDE_eqns\n"
	    "---------> Memory allocation within  ML_BuildQ to be checked...\n" );
   }
   QQ = ML_BuildQ( Nrows, new_Nrows, num_PDE_eqns, nullspace_dim,
		   reordered_decomposition,
		   nullspace_vect, new_nullspace_vect, i,
		   starting_amalg_bdry, reordered_amalg_bdry,
		   Amatrix->comm->USR_comm,
		   comm );
#else
   if( mypid == 0 ) 
     fprintf( stderr,
	      "*ML*ERR* Sorry, you cannot redistribute matrices within the ParMETIS\n"
	      "*ML*ERR* aggregation without epetra. Please recompile using epetra...\n" );
   exit( EXIT_FAILURE );
#endif

   if( starting_decomposition != NULL ) {
     (void)free( (void *)starting_decomposition );
     starting_decomposition = NULL;
   }
   if( reordered_decomposition != NULL ) {
     (void)free( (void *)reordered_decomposition );
     reordered_decomposition = NULL;
   }
   if( starting_amalg_bdry != NULL ) {
     ML_free( starting_amalg_bdry );
     starting_amalg_bdry = NULL;
   }
   if( starting_offset != NULL ) {
     (void)free( (void *)starting_offset );
     starting_offset = NULL;
   }
   if( reordered_offset != NULL ) {
     (void)free( (void *)reordered_offset );
     reordered_offset = NULL;
   }
          
   /* ********************************************************************** */
   /* Now reallocating aggr_index so that we can build the prolongator       */
   /* as if all the aggregates were local. Need some reallocations here.     */
   /* ********************************************************************** */

   nbytes = sizeof(int) * new_Nrows * num_PDE_eqns;
   
   if ( nbytes > 0 ) ML_memory_alloc((void**) &aggr_index, nbytes, "ACJ");
   else              aggr_index = NULL;

   /* k is the ID of the first aggregate on this proc */
#ifdef ML_MPI
   MPI_Scan( &aggr_count, &k, 1, MPI_INT, MPI_SUM, Amatrix->comm->USR_comm );
   k -= aggr_count;
#else
   k = 0;
#endif

   /*
     if( mypid < 5 && 8 < ML_Get_PrintLevel() ) {
     printf("ML_Aggregate_CoarsenParMETIS (level %d) : Assigning %d aggregates to each process\n",
	    ml_ag->cur_level, 
	    

   }
   */
   
   if( new_Nrows != 0 ) {

     jj = 0;
     for( i=0 ; i<aggr_count ; i++ ) {
       for( j=0 ; j<nodes_per_aggre[i+k] ; j++ ) {
	 aggr_index[jj] = i;
	 jj++;
       }
     }

     if( new_Nrows != jj ) {
       fprintf( stderr,
		"*ML*ERR* something went wrong in coarsening with ParMETIS:\n"
		"*ML*ERR* new_Nrows = %d, jj = %d\n"
		"*ML*ERR* (file %s, line %d)\n",
		new_Nrows, jj,
		__FILE__,
		__LINE__ );
       exit( EXIT_FAILURE );
     }
   }
   
   if( nodes_per_aggre != NULL ) {
     (void)free( (void *)nodes_per_aggre );
     nodes_per_aggre = NULL;
   }

   /* *************************** */
   /* TO ADD: OPERATOR COMPLEXITY */
   /* *************************** */

   for (i = new_Nrows - 1; i >= 0; i-- ) {
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
   
   new_Nrows  *= num_PDE_eqns;
   Nrows_global  *= num_PDE_eqns;

   /* count the size of each aggregate. Now all aggregates are local */
   
   aggr_cnt_array = (int *) ML_allocate(sizeof(int)*aggr_count);
   for (i = 0; i < aggr_count ; i++) aggr_cnt_array[i] = 0;
   for (i = 0; i < new_Nrows; i++) 
      if (aggr_index[i] >= 0) 
         aggr_cnt_array[aggr_index[i]]++;

   /* ============================================================= */
   /* Form tentative prolongator                                    */
   /* ============================================================= */
   
   Ncoarse = aggr_count;
   
   /* ============================================================= */
   /* check and copy aggr_index                                     */
   /* ------------------------------------------------------------- */

   level = ml_ag->cur_level;
   nbytes = new_Nrows * sizeof( int );
   ML_memory_alloc((void**) &(ml_ag->aggr_info[level]), nbytes, "AGl");
   count = aggr_count;
   for ( i = 0; i < new_Nrows; i+=num_PDE_eqns ) 
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

   for ( i = 0; i < new_Nrows; i++ ) 
   {
      if ( aggr_index[i] >= Ncoarse ) 
      {
         printf("*ML*WRN* index out of bound %d = %d (%d)\n"
		"*ML*WRN* (file %s, line %d)\n",
		i, aggr_index[i], 
                Ncoarse,
		__FILE__,
		__LINE__ );
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
   if( nbytes != 0 ) {
     ML_memory_alloc((void**)&(new_null),nbytes,"AGr");
     for (i = 0; i < Ncoarse*nullspace_dim*nullspace_dim; i++) 
       new_null[i] = 0.0;
   } else {
     new_null = NULL;
   }
   
   /* ------------------------------------------------------------- */
   /* initialize the row pointer for the CSR prolongation operator  */
   /* (each row will have at most nullspace_dim nonzero entries)    */
   /* ------------------------------------------------------------- */

   for (i = 0; i <= new_Nrows; i++) new_ia[i] = i * nullspace_dim;

   /* trying this when a Dirichlet row is taken out */
   j = 0;
   new_ia[0] = 0;
   for (i = 0; i < new_Nrows; i++) {
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
         printf("*ML*ERR* couldn't allocate memory in CoarsenParMETIS\n");
         exit(1);
      }
   }
   for (i = 0; i < new_Nrows; i+=num_PDE_eqns) 
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

   max_agg_size = 0;
   for (i = 0; i < aggr_count; i++) 
   {
      if (aggr_cnt_array[i] > max_agg_size) max_agg_size = aggr_cnt_array[i];
   }
   nbytes = max_agg_size * nullspace_dim * sizeof(double);
   if( nbytes > 0 ) ML_memory_alloc((void**)&qr_tmp, nbytes, "AGu");
   else             qr_tmp = NULL;
   nbytes = nullspace_dim * sizeof(double);
   if( nbytes > 0 ) ML_memory_alloc((void**)&tmp_vect, nbytes, "AGv");
   else             tmp_vect = NULL;
   
   lwork  = nullspace_dim;
   nbytes = nullspace_dim * sizeof(double);
   if( nbytes > 0 ) ML_memory_alloc((void**)&work, nbytes, "AGw");
   else             work = NULL;
  
   /* ------------------------------------------------------------- */
   /* perform block QR decomposition                                */
   /* ------------------------------------------------------------- */
     
   for (i = 0; i < aggr_count; i++) 
   {
      /* ---------------------------------------------------------- */
      /* set up the matrix we want to decompose into Q and R:       */
      /* ---------------------------------------------------------- */

      length = aggr_cnt_array[i];

      if (new_nullspace_vect == NULL) 
      {
         for (j = 0; j < (int) length; j++)
         {
            index = rows_in_aggs[i][j];
            for (k = 0; k < nullspace_dim; k++)
            {
              if ( reordered_amalg_bdry[index/num_PDE_eqns] == 1.0) qr_tmp[k*length+j] = 0.;
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
	   
            for (j = 0; j < (int) length; j++)
            {
               index = rows_in_aggs[i][j];
	       
               if ( reordered_amalg_bdry[index/num_PDE_eqns] == 1.0) qr_tmp[k*length+j] = 0.;
               else {
                  if (index < new_Nrows) {
		    qr_tmp[k*length+j] = new_nullspace_vect[k*new_Nrows+index];
                  }
                  else {
		    fprintf( stderr,
			     "*ML*ERR* error in QR factorization within ParMETIS aggregation\n"
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

	DGEQRF_F77(&(aggr_cnt_array[i]), &nullspace_dim, qr_tmp, 
			  &(aggr_cnt_array[i]), tmp_vect, work, &lwork, &info);
	if (info != 0)
	  pr_error("ErrOr in CoarsenParMETIS : "
		   "dgeqrf returned a non-zero %d %d\n",
		   aggr_cnt_array[i],i);

	if (work[0] > lwork) 
	  {
	    lwork=(int) work[0]; 
	    ML_memory_free((void**) &work);
	    ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGx");
	  }
	else lwork=(int) work[0];
		 
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
	  printf("*ML*ERR* in dorgqr on %d row (dims are %d, %d)\n",
		 i,aggr_cnt_array[i],
                 nullspace_dim);
	  printf("*ML*ERR* performing QR on a MxN matrix where M<N.\n");
	}
	DORGQR_F77(&(aggr_cnt_array[i]), &nullspace_dim,
			  &nullspace_dim, qr_tmp, &(aggr_cnt_array[i]),
			  tmp_vect, work, &lwork, &info);
	if (info != 0) {
	  printf("*ML*ERR* in dorgqr on %d row (dims are %d, %d)\n",
		 i,aggr_cnt_array[i],
                 nullspace_dim);
	  pr_error("Error in CoarsenParMETIS: dorgqr returned a non-zero\n");
	}

	if (work[0] > lwork) 
	  {
	    lwork=(int) work[0]; 
	    ML_memory_free((void**) &work);
	    ML_memory_alloc((void**) &work, sizeof(double)*lwork, "AGy");
	  }
	else lwork=(int) work[0];

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
	    if ( index < new_Nrows )
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
			 "*ML*ERR* error in QR factorization within ParMETIS\n" );
		exit( EXIT_FAILURE );
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

   } /* for( i over aggregates ) */
   
   ML_Aggregate_Set_NullSpace(ml_ag, num_PDE_eqns, nullspace_dim, 
                              new_null, Ncoarse*nullspace_dim);
   if( new_null != NULL ) {
     ML_memory_free( (void **) &new_null);
     new_null = NULL;
   }
   
   /* ------------------------------------------------------------- */
   /* set up the csr_data data structure                            */
   /* ------------------------------------------------------------- */

   ML_memory_alloc((void**) &csr_data, sizeof(struct ML_CSR_MSRdata),"CSR");
   csr_data->rowptr  = new_ia;
   csr_data->columns = new_ja;
   csr_data->values  = new_val;

   Pstart = ML_Operator_Create( Amatrix->comm );
			       
   ML_Operator_Set_ApplyFuncData( Pstart, nullspace_dim*Ncoarse, new_Nrows, 
                                  csr_data, new_Nrows, NULL, 0);
   Pstart->data_destroy = ML_CSR_MSR_ML_memorydata_Destroy;

   Pstart->getrow->pre_comm = ML_CommInfoOP_Create();
   
   ML_Operator_Set_Getrow((Pstart), new_Nrows, CSR_getrow);
   ML_Operator_Set_ApplyFunc(Pstart, CSR_matvec);
   Pstart->max_nz_per_row = 1;

   Pmatrix2 = ML_Operator_Create( Amatrix->comm );
   
   ML_2matmult(QQ, Pstart, Pmatrix2, ML_CSR_MATRIX );
   
   ML_Operator_Set_1Levels(Pmatrix2, (*Pmatrix)->from, (*Pmatrix)->to);
   ML_Operator_Set_BdryPts(Pmatrix2, (*Pmatrix)->bc);
   str2 = (char *)ML_allocate(80*sizeof(char));
   sprintf(str2,"%s",(*Pmatrix)->label);
   ML_Operator_Set_Label( Pmatrix2,str2);
   
   ML_free(str2);

   ML_Operator_Clean( *Pmatrix );

   memcpy((void *) *Pmatrix, (void *)Pmatrix2, sizeof(ML_Operator));
   /* FIXME : am I ok  ????? */
   ML_free(Pmatrix2);
      
   /* ********************************************************************** */
   /* I have to destroy the tentative local matrix, and the redistribution   */
   /* matrix QQ. This is actually an ML_Operator on the top of an Epetra     */
   /* object. So, I call ML_DestroyQ, which is a CPP function, to delete the */
   /* memory interally used by Epetra.                                       */
   /* ********************************************************************** */

   ML_Operator_Destroy( &Pstart ); 
   ML_Operator_Destroy( &QQ );

#ifdef ML_WITH_EPETRA
   ML_DestroyQ( );
#endif
   
   /* ------------------------------------------------------------- */
   /* clean up                                                      */
   /* ------------------------------------------------------------- */

   ML_memory_free((void**) &aggr_index);
   if( aggr_cnt_array != NULL ) ML_free(aggr_cnt_array);
   for (i = 0; i < aggr_count; i++) ML_free(rows_in_aggs[i]);
   ML_memory_free((void**)&rows_in_aggs);
   if( qr_tmp != NULL ) {
     ML_memory_free((void**)&qr_tmp);
     qr_tmp = NULL;
   }
   if( tmp_vect != NULL ) {
     ML_memory_free((void**)&tmp_vect);
     tmp_vect = NULL;
   }
   if( work != NULL ) {
     ML_memory_free((void**)&work);
     work = NULL;
   }

   aggr_curr = aggr_head;
   while ( aggr_curr != NULL ) 
   {
      supernode = aggr_curr;
      aggr_curr = aggr_curr->next;
      if ( supernode->length > 0 ) ML_free( supernode->list );
      ML_free( supernode );
   }

   if( reordered_amalg_bdry != NULL ) {
     ML_free( reordered_amalg_bdry );
     reordered_amalg_bdry = NULL;
   }
   if( new_nullspace_vect != NULL ) {
     ML_free( new_nullspace_vect );
     new_nullspace_vect = NULL;
   }

   /* ------------------- that's all folks --------------------------------- */

   if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Exiting `ML_Aggregate_CoarsenParMETIS'\n");
    printf("*ML*DBG* Ncoarse = %d, nummspace_dim = %d\n",
	   Ncoarse, nullspace_dim);
    printf("*ML*DBG* Total time = %e\n",  GetClock() - debug_starting_time);
  }

   return Ncoarse*nullspace_dim;

} /* ML_Aggregate_CoarsenParMETIS */

/* ********************************************************************** */
/* Count the nodes contained in each global aggregate                     */
/* ********************************************************************** */

static int ML_CountNodesPerAggre(int Nrows, int GraphDecomposition[],
				 int Naggre, int * NnodesPerAggre,
				 USR_COMM Comm) 
{
  
  int i, iaggre;
  int * count = NULL;
  int mypid;
  double debug_starting_time;
  
  /* ------------------- execution begins --------------------------------- */

  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Entering `ML_CountNodesPerAggre'\n");
    debug_starting_time = GetClock();
  }

  count = (int *) malloc( sizeof(int) * (Naggre) );

  if( count == NULL ) {
    fprintf( stderr,
	     "*ML*ERR* Not enough memory to allocate %d bytes\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     (int)sizeof(int) * (Naggre),
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }

#ifdef ML_MPI
  MPI_Comm_rank( Comm, &mypid );
#else
  mypid = 0;
#endif
  
  for( i=0 ; i<Naggre; i++ ) count[i] = 0;
  
  for( i=0 ; i<Nrows ; i++ ) {
    iaggre = GraphDecomposition[i];
    if( iaggre > Naggre || iaggre < 0 ) {
      fprintf( stderr,
	       "*ML*ERR* something went wrong in counting the nodes per aggre\n"
	       "*ML*ERR* node %d is assigned to global aggregate %d, but you\n"
	       "*ML*ERR* have only %d aggregates. This is proc %d.\n",
	       i,
	       iaggre,
	       Naggre,
	       mypid );
    }
    if( iaggre >= 0 && iaggre < Naggre ) count[iaggre]++;
  }

#ifdef ML_MPI
  MPI_Allreduce( count, NnodesPerAggre,
		 Naggre, MPI_INT, MPI_SUM, Comm);
#else
  for( i=0 ; i<Naggre ; i++ )
    NnodesPerAggre[i] = count[i];
#endif

  /* ********************************************************************** */
  /* some checks on the counts. I print out a warning if an aggregate has   */
  /* zero nodes or just one node.                                           */
  /* ********************************************************************** */

  for( i=0 ; i<Naggre ; i++ ) {
    if( NnodesPerAggre[i] == 0 && 2 < ML_Get_PrintLevel() ) {
      fprintf( stderr,
	       "*ML*WRN* aggregate %d on proc %d has zero nodes\n",
	       i,
	       mypid );
    } else if( NnodesPerAggre[i] == 1 && 8 < ML_Get_PrintLevel() ) {
      fprintf( stderr,
	       "*ML*WRN* aggregate %d on proc %d has only one node\n",
	       i,
	       mypid );
    }
  }

  if( count != NULL ) (void)free( (void *)count ); count = NULL;

  /* ------------------- that's all folks --------------------------------- */

  if( PARMETIS_DEBUG_LEVEL == 3 ) {
    printf("*ML*DBG* Exiting `ML_CountNodesPerAggre'\n");
    printf("*ML*DBG* Total time = %e\n",  GetClock() - debug_starting_time);
  }

  return 0;
  
} /* ML_CountNodesPerAggre */
