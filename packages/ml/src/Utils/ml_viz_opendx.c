/* ************************************************************************* */
/* ************************************************************************* */
/* Functions to write OpenDX-readable files (about graph decomposition)      */
/* ************************************************************************* */
/* Author        : Marzio Sala (SNL)                                         */
/* Date          : October 2003                                              */
/* ************************************************************************* */
/* Local Function :                                                          */
/*    ML_Aggregate_VisualizeWithOpenDX                                       */
/* ************************************************************************* */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_utils.h"
#include "ml_viz_opendx.h"
#include "ml_agg_METIS.h"

int ML_DecomposeGraph_LocalToGlobal( ML_Comm *comm,
				     int N_rows, int N_parts,
				     int graph_decomposition[] )
{

  int i;
  int N_procs = comm->ML_nprocs;
  int *offsets = (int*)malloc(sizeof(int)*(N_procs+1));
  
  ML_DecomposeGraph_BuildOffsets( N_parts, offsets, comm->ML_nprocs, comm->USR_comm);

  for( i=0 ; i<N_rows ; i++ )
    graph_decomposition[i] += offsets[comm->ML_mypid];
  
  ML_free( offsets ); offsets=NULL;

  return 0;
  
} /* ML_DecomposeGraph_LocalToGlobal */

int ML_DecomposeGraph_GlobalToLocal( ML_Comm *comm,
				     int N_rows, int N_parts,
				     int graph_decomposition[] )
{

  int i;
  int N_procs = comm->ML_nprocs;
  int *offsets = (int*)malloc(sizeof(int)*(N_procs+1));
  
  ML_DecomposeGraph_BuildOffsets( N_parts, offsets, comm->ML_nprocs, comm->USR_comm);

  for( i=0 ; i<N_rows ; i++ )
    graph_decomposition[i] -= offsets[comm->ML_mypid];

  ML_free( offsets ); offsets=NULL;

  return 0;
  
} /* ML_DecomposeGraph_LocalToGlobal */

int ML_DecomposeGraph_ConvertToDouble(ML_Operator *Amatrix,
				      int local_or_global,
				      int N_local_parts,
				      int graph_decomposition[],
				      double *values)
{

  int i;
  int N_rows = Amatrix->getrow->Nrows;
  int convert_to_local = 0;

  /* ------------------- execution begins --------------------------------- */

  /* need aggregate numbering in local form */
  
  if( local_or_global == ML_LOCAL_INDICES ) {
    ML_DecomposeGraph_LocalToGlobal( Amatrix->comm,
				     N_rows, N_local_parts,
				     graph_decomposition );
     convert_to_local = 1;
  }

  /* put them in values as doubles */

  for( i=0 ; i<N_rows ; i++ )
    values[i] = (double)graph_decomposition[i];
  
  /* back to input data */
  
  if( convert_to_local == 1 ) {
    ML_DecomposeGraph_GlobalToLocal( Amatrix->comm,
				     N_rows, N_local_parts,
				     graph_decomposition );
  }

  return 0;
  
} /* ML_DecomposeGraph_ConvertToDouble */

/* ======================================================================== */
/*!
 \brief write graph decomposition of the current level in a graphical
 format readable by OpenDX

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_VisualizeWithOpenDX( ML_Aggregate_Viz_Stats info,
				      char base_filename[],
				      ML_Comm *comm )
{

  int i,j, irow;
  ML_Operator *Amatrix = (ML_Operator *)(info.Amatrix);
  double *x = info.x;
  double *y = info.y;
  double *z = info.z;
  int local_or_global = info.local_or_global;
  int *graph_decomposition = info.graph_decomposition;
  int Naggregates = info.Naggregates;
  int Nlocal = info.Nlocal;
  int mypid = comm->ML_mypid;
  int nprocs = comm->ML_nprocs;
  int Nrows = Amatrix->getrow->Nrows;
  int *Nnz_row = (int*)malloc(sizeof( int ) * Nrows);
  int N_edges = 0;
  char filename[FILENAME_MAX];
  FILE *fp;
  int allocated = 0;
  int * rowi_col = NULL;
  int rowi_N;
  double * rowi_val = NULL;
  int offset;
  int *values;
  int ok;
  int shuffle, * reorder = NULL;
  
  /* ------------------- execution begins --------------------------------- */

  if( Nlocal != Nrows ) {
    fprintf( stderr,
	     "*ML*ERR* number of rows and lenght of graph_decomposition\n"
	     "*ML*ERR* differs (%d - %d)\n"
	     "*ML*ERR* (file %s, liine %d)\n",
	     Nrows,
	     Nlocal,
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  /* define file name for this processor */
  
  sprintf( filename,
	   "%s%d",
	   base_filename,
	   Amatrix->comm->ML_mypid );

  if( (fp = fopen( filename, "w" )) == NULL ) {
    fprintf( stderr,
	     "*VIZ*ERR* cannot open file `%s'\n",
	     filename );
    exit( EXIT_FAILURE );
  }
 
  /* write on file the nodal coordinates */
  
  fprintf( fp,
	   "\nobject 1 class array type float rank 1"
	   " shape 3 items %d data follows\n" ,
	   Nrows );

  /* handle 1D, 2D and 3D cases */
  
  if( y == NULL ) 
    for( irow=0 ; irow<Nrows ; irow++ ) {
      fprintf( fp,
	       "%f 0 0\n",
	       x[irow] );
    }
  else if( z == NULL )
    for( irow=0 ; irow<Nrows ; irow++ ) {
      fprintf( fp,
	       "%f %f 0\n",
	       x[irow], y[irow] );
    }
  else
    for( irow=0 ; irow<Nrows ; irow++ ) {
      fprintf( fp,
	       "%f %f %f\n",
	       x[irow], y[irow], z[irow] );
    }

  /* need to know the number of edges first */

  for( irow=0 ; irow<Nrows ; irow++ ) {

    ML_get_matrix_row(Amatrix, 1, &irow, &allocated, &rowi_col, &rowi_val,
		      &rowi_N, 0);

    for( i=0 ; i<rowi_N ; i++ ) {
      if( rowi_col[i]<Nrows ) N_edges++;
      
    }
    
  }
  
  /* write on file the edges connectivity */

  fprintf( fp,
	   "\nobject 2 class array type int"
	   " rank 1 shape 2 items %d data follows\n",
	   N_edges );

  for( irow=0 ; irow<Nrows ; irow++ ) {

    ML_get_matrix_row(Amatrix, 1, &irow, &allocated, &rowi_col, &rowi_val,
		      &rowi_N, 0);

    for( i=0 ; i<rowi_N ; i++ ) {
      if( rowi_col[i]<Nrows )
	fprintf( fp,
		 "%d %d\n",
		 irow, rowi_col[i] );
      
    }
    
  }
  
  fprintf( fp,
	   "attribute \"element type\" string \"lines\"\n"
	   "attribute \"ref\" string \"positions\"\n" );
  
  /* how write the partition_index */
  
  fprintf( fp,
	   "\nobject 3 class array type float"
	   " rank 0 items %d data follows\n",
	   Nrows );

  /* ********************************************************************** */
  /* if graph_decomposition contains global numbers for the aggregates, as  */
  /* obtained for instance using ParMETIS, fine. Otherwise, I suppose that  */
  /* the first Naggregates are hold on process 0, the next Naggregates on   */
  /* proc 1, and so on. As Naggregate may vary from proc to proc, I use the */
  /* MPI function MPI_Scan to obtain the offset. offset is set to 0 for     */
  /* global indices.                                                        */
  /* ********************************************************************** */

  values = (int *) malloc( sizeof(int) * Nrows );
  if( values == NULL ) {
    fprintf( stderr,
	     "*ML*ERR* not enough memory for %d bytes\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     Nrows * (int)sizeof(int),
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  if( local_or_global == ML_LOCAL_INDICES ) {

#ifdef ML_MPI
    MPI_Scan ( &Naggregates, &offset, 1, MPI_INT, MPI_SUM,
	       MPI_COMM_WORLD );
    offset -= Naggregates;
#else
    offset = 0;
#endif
  }

  /* ******************************************************************** */
  /* natural ordering, with lots of aggregates results in colors too      */
  /*  close for aggregates of the same domain. So, here I rearrange the   */
  /* node numbering so that coloring is (hopefully) improved.             */
  /* Also, if `shuffle' == 1, then I create a temp vector, with random    */
  /* reodering, so that visualization will (probably) work also for the 1 */
  /* proc case.                                                           */
  /* ******************************************************************** */

  shuffle = 1;

  reorder = (int *) malloc( sizeof(int) * Naggregates );

  if( shuffle == 1 ) { /* now always reordering */
    
    if( reorder == NULL ) {
      fprintf( stderr,
	       "*ML*ERR* not enough memory for %d bytes\n"
	       "*ML*ERR* (file %s, line %d)\n",
	       (int)sizeof(int) * Naggregates,
	       __FILE__,
	       __LINE__ );
      exit( EXIT_FAILURE );
    }
    for( i=0 ; i<Naggregates ; ++i ) reorder[i] = -1;

    srand(0);
    
    for( i=0 ; i<Naggregates ; ++i ) {

      do {

	ok = 0;
	
	j = (int)(1.0*(Naggregates)*rand()/RAND_MAX);

	if( reorder[j] == -1 && j<Naggregates ) {
	  reorder[j] = i;
	  ok = 1;
	}
      } while( ok == 0 );

    }

  } else {

    for( i=0 ; i<Naggregates ; ++i ) reorder[i] = i;
    
  }

  if( local_or_global == ML_LOCAL_INDICES ) {
    
    /* max_Naggregates = ML_gmax_int( Naggregates, comm); */

    for( i=0 ; i<Nrows ; i++ ) {
      values[i] = mypid +  nprocs * reorder[graph_decomposition[i]];
    }
    
  } else {

    for( i=0 ; i<Nrows ; i++ ) {
      values[i] = reorder[graph_decomposition[i]];
    }
  }

  ML_free(reorder);
  
  for( irow=0 ; irow<Nrows ; irow++ ) 
    fprintf( fp,
	     "%f\n",
	     (float)values[irow] );

  
  /* still some stuff for OpenDX */
  
  fprintf( fp,
	   "attribute \"dep\" string \"positions\"\n"
	   "\nobject \"viz mamma\" class field\n"
	   "component \"positions\" value 1\n"
	   "component \"connections\" value 2\n"
	   "component \"data\" value 3\n"
	   "end\n" );

  fclose( fp );

  /* ------------------- that's all folks --------------------------------- */

  ML_free( values ); values=NULL;
  ML_free( Nnz_row ); Nnz_row=NULL;
  ML_free(rowi_col); ML_free(rowi_val);
  rowi_col = NULL; rowi_val = NULL;
  allocated = 0; 
  
  return 0;
  
} /* ML_VisualizeWithOpenDX */

/* ======================================================================== */
/*!
 \brief write graph decomposition of the current level in a graphical
 format readable by OpenDX

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_VisualizeXYZ( ML_Aggregate_Viz_Stats info,
			      char base_filename[],
			      ML_Comm *comm,
			      double * vector)
{

  int i,j, irow, ipid;
  ML_Operator *Amatrix = (ML_Operator *)(info.Amatrix);
  double *x = info.x;
  double *y = info.y;
  int local_or_global = info.local_or_global;
  int *graph_decomposition = info.graph_decomposition;
  int Naggregates = info.Naggregates;
  int mypid = comm->ML_mypid;
  int nprocs = comm->ML_nprocs;
  int Nrows = Amatrix->getrow->Nrows;
  FILE *fp;
  double val;
  int AggrToVisualize = -1;
  char * str;
  char filemode[2];
  int Nlocal = info.Nlocal;

  /* ------------------- execution begins --------------------------------- */

  if( Nlocal != Nrows ) {
    fprintf( stderr,
	     "*ML*ERR* number of rows and lenght of graph_decomposition\n"
	     "*ML*ERR* differs (%d - %d)\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     Nrows,
	     Nlocal,
	     __FILE__,
	     __LINE__ );
    exit( EXIT_FAILURE );
  }
  
  if( mypid == 0 ) filemode[0] = 'w';
  else             filemode[0] = 'a';
  filemode[1] = '\0';

  str = getenv("ML_VIZ_AGGR");
  if( str != NULL ) AggrToVisualize = atoi(str);

  for( ipid=0 ; ipid<nprocs ; ++ipid ) {
    if( ipid == mypid ) {
      if( (fp = fopen( base_filename, filemode )) == NULL ) {
	fprintf( stderr,
		"*VIZ*ERR* cannot open file `%s'\n",
		base_filename );
	exit( EXIT_FAILURE );
      }

      for( irow=0 ; irow<Nrows ; irow++ ) {
	if( vector != NULL ) val = vector[irow];
	else if( AggrToVisualize != -1 ) { 
	  if( graph_decomposition[irow] == AggrToVisualize ) val = 1.0;
	  else                                               val = 0.0;
	} else 
	  val = 1.0*graph_decomposition[irow];
	fprintf( fp,
		"%f %f %f\n",
		x[irow], y[irow], val );
      }
      fclose(fp);
    }
#ifdef ML_MPI
    /* FIXME: I can get the communicator from ML_Comm... */
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }
  
  return 0;
  
} /* ML_VisualizeWithXYZ */
