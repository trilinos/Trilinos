/* ======================================================================== */
/*!
 \file ml_agg_info.c

 \brief Various stats on aggregates

 \author Marzio Sala, SNL, 9214

 \date 23-October-2003
 
*/
/* ------------------------------------------------------------------------ */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ml_aggregate.h"
#include "ml_lapack.h"
#include "ml_utils.h"
#include "ml_viz_opendx.h"
#include "ml_viz_stats.h"
#include "ml_agg_info.h"

/* ======================================================================== */
/*!
 \brief allocate memory required to keep trace of the graph decomposition
 among the various levels

 MaxLevels should be the same of the one used to create the multigrid
 structure

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_VizAndStats_Setup( ML_Aggregate *ag, int MaxLevels )
{

  int i;
  ML_Aggregate_Viz_Stats * info;
   
  if ( ag->ML_id != ML_ID_AGGRE ) 
    {
      printf("ML_Set_NumberLocalAggregates_METIS : wrong object. \n");
      exit(-1);
    }

  info = (ML_Aggregate_Viz_Stats *) ML_allocate(sizeof(ML_Aggregate_Viz_Stats)*(MaxLevels+1) );
  
  if( info == NULL ) {
    fprintf( stderr,
	     "*ML*ERR* not enough memory for %d bytes\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     sizeof(ML_Aggregate_Viz_Stats)*MaxLevels,
	     __FILE__,
	     __LINE__ );
  }

  for( i=0 ; i<MaxLevels ; i++ ) {
    info[i].x = NULL;
    info[i].y = NULL;
    info[i].z = NULL;
    info[i].graph_decomposition = NULL;
    info[i].Nlocal = 0;
    info[i].Naggregates = 0;
    info[i].is_filled = ML_NO;
  }

  ag->aggr_viz_and_stats = (void *)info;

  return 0;
  
} /* ML_Aggregate_VizAndStats_Setup */

/* ======================================================================== */
/*!
 \brief free memory allocated by \c ML_Aggregate_VizAnsStats_Setup

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_VizAndStats_Clean( ML_Aggregate *ag, int MaxLevels )
{

  int i;
  ML_Aggregate_Viz_Stats *info;
  
  if ( ag->ML_id != ML_ID_AGGRE ) 
    {
      printf("ML_Aggregate_VizAndStats_Clean : wrong object. \n");
      exit(-1);
    }

  for( i=0 ; i<MaxLevels ; i++ ) {
    info = (ML_Aggregate_Viz_Stats *)(ag->aggr_viz_and_stats);
    if( info[i].x != NULL ) {
      ML_free( info[i].x );
      info[i].x = NULL;
    }
    if( info[i].y != NULL ) {
      ML_free( info[i].y );
      info[i].y = NULL;
    }
    if( info[i].z != NULL ) {
      ML_free( info[i].z );
      info[i].z = NULL;
    }
    if( info[i].graph_decomposition != NULL ) {
      ML_free( info[i].graph_decomposition );
      info[i].graph_decomposition = NULL;
    }
    info[i].Nlocal = -1;
    info[i].Naggregates = 0;
    info[i].is_filled = ML_NO;
  }

  ML_free( ag->aggr_viz_and_stats );  ag->aggr_viz_and_stats = NULL;

  return 0;
  
} /* ML_Aggregate_VizAndStats_Clean */

/* ======================================================================== */
/*!
 \brief compute the center of gravity of the aggregates.

*/
/* ------------------------------------------------------------------------ */
 
void ML_Aggregate_ComputeRadius( ML_Aggregate_Viz_Stats finer_level,
				 ML_Aggregate_Viz_Stats coarser_level,
				 double R[] )
{

  int i,iaggre;
  double Ri;
  int N_fine = finer_level.Nlocal;
  int N_aggregates = finer_level.Naggregates;
  int *graph_decomposition = finer_level.graph_decomposition;
  int local_or_global = finer_level.local_or_global;
  double *x = finer_level.x;
  double *y = finer_level.y;
  double *z = finer_level.z;
  double *xm = coarser_level.x;
  double *ym = coarser_level.y;
  double *zm = coarser_level.z;
  
  /* ------------------- execution begins --------------------------------- */

  /* compute the radius of a ball circumscribing  each aggregate */
  /* for each node, I compute the distance between this node and the
     center of gravity of its aggregate. R[iaggre] is the maximum of all
     those distances */

  for( i=0 ; i<N_aggregates ; i++ ) R[i] = 0.0;
  
  if( local_or_global == ML_LOCAL_INDICES ) {
    
    for( i=0 ; i<N_fine ; i++ ) {
      iaggre = graph_decomposition[i];
      Ri = pow(x[i] - xm[iaggre], 2);
      if( ym != NULL )  Ri += pow(y[i] - ym[iaggre], 2);
      if( zm != NULL )  Ri += pow(z[i] - zm[iaggre], 2);
      if( Ri > R[iaggre] ) R[iaggre]=sqrt(Ri);
    }
    
  } else if( local_or_global == ML_GLOBAL_INDICES ) {

    printf("To do...\n");

  } else {

    fprintf( stderr,
	     "*ML*ERR* input parameter 4 has an incorrect value (%d)\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     local_or_global,
	     __FILE__,
	     __LINE__ );
  }

  return;
  
} /* ML_aggregate_ComputeVolumne */

#include "float.h"
#include "limits.h"

/* ======================================================================== */
/*!
 \brief Compute the maximum dimension of a box circumscribing each
 aggregate.

*/
/* ------------------------------------------------------------------------ */

void ML_Aggregate_ComputeBox( ML_Aggregate_Viz_Stats finer_level,
			      int Naggregates_global, double R[],
			      int offset, ML_Comm * comm)
{

  int i,iaggre;
  int N_fine = finer_level.Nlocal;
  int *graph_decomposition = finer_level.graph_decomposition;
  double *x = finer_level.x;
  double *y = finer_level.y;
  double *z = finer_level.z;
  double* xmin;
  double* xmax;
  double* ymin;
  double* ymax;
  double* zmin;
  double* zmax;
  double * dtemp;
  
  /* ------------------- execution begins --------------------------------- */

  xmin = (double*)malloc(sizeof(double)*Naggregates_global);
  xmax = (double*)malloc(sizeof(double)*Naggregates_global);
  ymin = (double*)malloc(sizeof(double)*Naggregates_global);
  ymax = (double*)malloc(sizeof(double)*Naggregates_global);
  zmin = (double*)malloc(sizeof(double)*Naggregates_global);
  zmax = (double*)malloc(sizeof(double)*Naggregates_global);
  dtemp = (double*)malloc(sizeof(double)*Naggregates_global);
  
  for( i=0 ; i<Naggregates_global ; i++ ) R[i] = 0.0;
  
  for( i=0 ; i<Naggregates_global ; i++ ) {
    xmin[i] =  DBL_MAX; ymin[i] =  DBL_MAX, zmin[i] =  DBL_MAX;
    xmax[i] = -DBL_MAX; ymax[i] = -DBL_MAX, zmax[i] = -DBL_MAX;
  }

  for( i=0 ; i<N_fine ; i++ ) {
    iaggre = graph_decomposition[i]+offset;
    xmin[iaggre] = ML_min(xmin[iaggre], x[i]);
    xmax[iaggre] = ML_max(xmax[iaggre], x[i]);
    if( y != NULL ) {
      ymin[iaggre] = ML_min(ymin[iaggre], y[i]);
      ymax[iaggre] = ML_max(ymax[iaggre], y[i]);
    }
    if( z != NULL ) {
      zmin[iaggre] = ML_min(zmin[iaggre], z[i]);
      zmax[iaggre] = ML_max(zmax[iaggre], z[i]);
    }
  }

#ifdef ML_MPI
  MPI_Reduce(xmin,dtemp,Naggregates_global,MPI_DOUBLE,MPI_MIN,0,
	     comm->USR_comm);
  if( comm->ML_mypid == 0 )
    for( i=0 ; i<Naggregates_global ; ++i ) xmin[i] = dtemp[i];
  
  MPI_Reduce(xmax,dtemp,Naggregates_global,MPI_DOUBLE,MPI_MAX,0,
	     comm->USR_comm);
  if( comm->ML_mypid == 0 )
    for( i=0 ; i<Naggregates_global ; ++i ) xmax[i] = dtemp[i];

  if( y != NULL ) {
    MPI_Reduce(ymin,dtemp,Naggregates_global,MPI_DOUBLE,MPI_MIN,0,
	       comm->USR_comm);
    if( comm->ML_mypid == 0 )
      for( i=0 ; i<Naggregates_global ; ++i ) ymin[i] = dtemp[i];
    
    MPI_Reduce(ymax,dtemp,Naggregates_global,MPI_DOUBLE,MPI_MAX,0,
	       comm->USR_comm);
    if( comm->ML_mypid == 0 )
      for( i=0 ; i<Naggregates_global ; ++i ) ymax[i] = dtemp[i];
  }
  
  if( z != NULL ) {
    MPI_Reduce(zmin,dtemp,Naggregates_global,MPI_DOUBLE,MPI_MIN,0,
	       comm->USR_comm);
    if( comm->ML_mypid == 0 )
      for( i=0 ; i<Naggregates_global ; ++i ) zmin[i] = dtemp[i];
    
    MPI_Reduce(zmax,dtemp,Naggregates_global,MPI_DOUBLE,MPI_MAX,0,
	       comm->USR_comm);
    if( comm->ML_mypid == 0 )
      for( i=0 ; i<Naggregates_global ; ++i ) zmax[i] = dtemp[i];
  }
#endif

  if( comm->ML_mypid == 0 ){
    for( iaggre=0 ; iaggre<Naggregates_global ; iaggre++ ) {
      R[iaggre] = xmax[iaggre] - xmin[iaggre];
      if( y != NULL )
	R[iaggre] = ML_max(R[iaggre], ymax[iaggre] - ymin[iaggre]);
      if( z != NULL )
	R[iaggre] = ML_max(R[iaggre], zmax[iaggre] - zmin[iaggre]);
      
    }
  }
  
  ML_free ( xmin ); 
  ML_free ( xmax ); 
  ML_free ( ymin ); 
  ML_free ( ymax ); 
  ML_free ( zmin );
  ML_free ( zmax );
  ML_free( dtemp );
  
  return;
  
} /* ML_Aggregate_ComputeBox */

/* ======================================================================== */
/*!
 \brief Compute the center of gravity of each aggregate

 This function computes the center of gravity of the non-smoothed
 aggregates, encoded in the integer vector \c graph_decomposition'.

 For 2D problems, simply set \c z=NULL.
 For 1D problems, set \c y=NULL,z=NULL.
 
*/
/* ------------------------------------------------------------------------ */

void ML_Aggregate_ComputeCenterOfGravity( ML_Aggregate_Viz_Stats finer_level,
					  ML_Aggregate_Viz_Stats coarser_level ,
					  ML_Comm * comm)

{

  int i,j,iaggre;
  int N_coarser = finer_level.Naggregates, N_coarser_global;
  int N_finer = finer_level.Nlocal, N_finer_global;
  int *graph_decomposition = finer_level.graph_decomposition;
  double *x = finer_level.x;
  double *y = finer_level.y;
  double *z = finer_level.z;
  double *xm = coarser_level.x;
  double *ym = coarser_level.y;
  double *zm = coarser_level.z;
  int * count = NULL, * itemp = NULL;
  double * dtemp = NULL;
  int offset;
  double * xmtemp = NULL, * ymtemp = NULL, * zmtemp = NULL;
  
  /* ------------------- execution begins --------------------------------- */

  /* too expensive now for large matrices... To treat LOCAL_INDICES
     sepaately */
  
  switch( finer_level.local_or_global ) {
  case ML_LOCAL_INDICES:
    N_finer_global = ML_gsum_int(N_finer,comm); 
    N_coarser_global = ML_gsum_int(N_coarser,comm);
#ifdef ML_MPI
    MPI_Scan(&N_coarser, &offset, 1, MPI_INT, MPI_SUM, comm->USR_comm);
    offset -= N_coarser;
#else
    offset = 0;
#endif
    break;
  case ML_GLOBAL_INDICES:
    N_finer_global = N_finer;
    N_coarser_global = N_coarser;
    offset = 0;
    break;
  }

  count = (int *) malloc( sizeof(int)*N_coarser_global);
  xmtemp = (double *) malloc( sizeof(double)*N_coarser_global);
  if( ym != NULL ) ymtemp = (double *) malloc( sizeof(double)*N_coarser_global);
  if( zm != NULL ) zmtemp = (double *) malloc( sizeof(double)*N_coarser_global);
  
  for( i=0 ; i<N_coarser_global ; i++ ) {
    count[i] = 0;
    xmtemp[i] = 0.0;
    if( ymtemp != NULL ) ymtemp[i] = 0.0;
    if( zmtemp != NULL ) zmtemp[i] = 0.0;
  }

  for( i=0 ; i<N_finer ; i++ ) {
    iaggre = graph_decomposition[i]+offset;
    xmtemp[iaggre] += x[i];
    if( ymtemp != NULL )  ymtemp[iaggre] += y[i];
    if( zmtemp != NULL )  zmtemp[iaggre] += z[i];
    count[iaggre]++;
  }

#ifdef ML_MPI
  itemp = (int *) malloc( sizeof(int)*N_coarser_global);
  dtemp = (double *) malloc( sizeof(double)*N_coarser_global);
  MPI_Allreduce(count,itemp,N_coarser_global,MPI_INT,MPI_SUM,comm->USR_comm);
  for( i=0 ; i<N_coarser_global ; i++ ) count[i] = itemp[i];
  
  MPI_Allreduce(xmtemp,dtemp,N_coarser_global,MPI_DOUBLE,MPI_SUM,comm->USR_comm);
  for( i=0 ; i<N_coarser_global ; i++ ) xmtemp[i] = dtemp[i];

  if( ymtemp != NULL ) {
    MPI_Allreduce(ymtemp,dtemp,N_coarser_global,MPI_DOUBLE,MPI_SUM,comm->USR_comm);
    for( i=0 ; i<N_coarser_global ; i++ ) ymtemp[i] = dtemp[i];
  }
  
  if( zmtemp != NULL ) {
    MPI_Allreduce(zmtemp,dtemp,N_coarser_global,MPI_DOUBLE,MPI_SUM,comm->USR_comm);
    for( i=0 ; i<N_coarser_global ; i++ ) zmtemp[i] = dtemp[i];
  }
#endif
  
  for( i=0 ; i<N_coarser ; i++ ) {
    j = i+offset;
    if( count[j] != 0 ) {
      xm[i] = xmtemp[j]/count[j];
      if( ym != NULL )  ym[i] = ymtemp[j]/count[j];
      if( zm != NULL )  zm[i] = zmtemp[j]/count[j];
    }
  }
  
  ML_free( count );
#ifdef ML_MPI
  ML_free( itemp );
  ML_free( dtemp );
#endif
  ML_free( xmtemp );
  if( ym != NULL ) ML_free( ymtemp );
  if( zm != NULL ) ML_free( zmtemp );
  
  return;
  
} /* ML_Aggregate_ComputeCenterOfGravity */

/* ======================================================================== */
/*!
 \brief compute the volume of the aggregates.

*/
/* ------------------------------------------------------------------------ */

void ML_Aggregate_ComputeVolume( int N_fine,
				 int N_aggregates,
				 int graph_decomposition[],
				 int local_or_global,
				 double volume[],
				 double V[] )
{

  int i,iaggre;

  /* ------------------- execution begins --------------------------------- */

  if( local_or_global == ML_LOCAL_INDICES ) {

    for( i=0 ; i<N_aggregates ; i++ ) {
      V[i] = 0;
    }
    
    for( i=0 ; i<N_fine ; i++ ) {
      iaggre = graph_decomposition[i];
      V[iaggre] += volume[i];
    }

  } else if( local_or_global == ML_GLOBAL_INDICES ) {

    printf("To do...\n");

  } else {

    fprintf( stderr,
	     "*ML*ERR* input parameter 4 has an incorrect value (%d)\n"
	     "*ML*ERR* (file %s, line %d)\n",
	     local_or_global,
	     __FILE__,
	     __LINE__ );
    
  }
  
  return;
    
} /* ML_aggre_compute_volume */

/* ======================================================================== */
/*!
 \brief Some analysis on the decomposition into aggregates

*/
/* ------------------------------------------------------------------------ */

void ML_Aggregate_AnalyzeLocalGraphDec( int N_aggregates,
					int *nodes_per_aggregate,
					int *min_loc_aggre,
					int *max_loc_aggre,
					double *avg_loc_aggre,
					double *std_loc_aggre,
					ML_Comm *comm ) 
{
  int i, j, min, max, sum, gsum, N;
  double avg, std;

  /* ------------------- execution begins --------------------------------- */

  N = ML_gsum_int( N_aggregates, comm );
  
  min = INT_MAX;
  max = 0;
  sum = 0;

  /* max, min, and sum (to compute the average) */
  
  for( i=0 ; i<N_aggregates ; i++ ) {
    j = nodes_per_aggregate[i];
    sum += j;
    if( j>max ) max = j;
    if( j<min ) min = j;
  }
  
  min  = ML_gmin_int( min, comm  );
  sum  = ML_gsum_int( sum, comm );
  max  = ML_gmax_int( max, comm );
  gsum = ML_gsum_int( sum, comm );
  
  avg = 1.0*gsum/N;
  
  /* now standard deviation */

  std = 0.0;
  
  for( i=0 ; i<N_aggregates ; i++ ) {
    std += pow(1.0*nodes_per_aggregate[i] - avg,2);
  }

  /* compute the standard deviation as

                 \sum){i=1}^N ( n_p_a[i] - ave )^e
     std = sqrt( --------------------------------- )
                           N - 1
  */
  
  if( std > 0.00001 && N>1 ) {
  
    std = ML_gsum_int( std, comm );
    std = sqrt(std/(N-1));

  }
  
  *min_loc_aggre = min;
  *max_loc_aggre = max;
  *avg_loc_aggre = avg;
  *std_loc_aggre = std;

  return;
  
} 

/* ======================================================================== */
/*!
 \brief Some statistics on a double vector

 This function finds out the minimum, maximum, average values of a
 double vector. It also computes the standard deviation.

 Parameter list:
 - \c N_update : size of the local vector
 - \c vector : double vector of size \c N_update
 - \c min_dec,max_vec,_avg_vec,srd_vec : in output, results of the computations
 
*/
/* ------------------------------------------------------------------------ */

void ML_Aggregate_AnalyzeVector( int Nlocal,
				 double vector[],
				 double *min_vec,
				 double *max_vec,
				 double *avg_vec,
				 double *std_vec,
				 ML_Comm *comm ) 
{
  int i, N;
  double d, dmin, dmax, sum;
  double avg, std;

  /* ------------------- execution begins --------------------------------- */

  N = ML_gsum_int( Nlocal, comm );
  
  dmin = DBL_MAX;
  dmax = DBL_MIN;
  sum = 0;

  /* local max, min, and sum (to compute the average) */
  
  for( i=0 ; i<Nlocal ; i++ ) {
    d = vector[i];
    sum += d;
    if( d>dmax ) dmax = d;
    if( d<dmin ) dmin = d;
  }

  dmin = ML_gmin_double( dmin, comm );
  sum  = ML_gsum_double( sum, comm  );
  dmax = ML_gmax_double( dmax, comm );

  avg = sum/N;
  
  /* now standard deviation */

  std = 0.0;
  
  for( i=0 ; i<Nlocal ; i++ ) {
    std += pow(vector[i] - avg,2);
  }

  /* compute the standard deviation as

                 \sum){i=1}^N ( n_p_a[i] - ave )^e
     std = sqrt( --------------------------------- )
                           N - 1
  */
  
  std = ML_gsum_double( std, comm );
  
  if( std > 0.00001 && N>1 ) {
  
    std = sqrt(std/(N-1));

  }
  
  *min_vec = dmin;
  *max_vec = dmax;
  *avg_vec = avg;
  *std_vec = std;

  return;
  
} /*  ML_Aggregate_AnalyzeVector */
 
/* ======================================================== */
/* the folliwing is a very simple function to format output */
/* -------------------------------------------------------- */

void ML_print_dash_line( int mypid )
{
  if( mypid == 0 )
    printf( "*ML* ------------------------- : "
	    "------------------------------------------\n" );
}

void ML_Aggregate_CountLocal( int N_fine, int graph_decomposition[],
			      int N_aggregates, int nodes_per_aggregate[] )

{
  int i, j;  

  for( i=0 ; i<N_aggregates ; i++ ) nodes_per_aggregate[i] = 0;

  for( i=0 ; i<N_fine ; i++ ) {
    j = graph_decomposition[i];
    if( j < 0 ) {
      fprintf( stderr,
	       "*ML*ERR* Something went wrong in buildind up the \n"
	       "*ML*ERR* the aggregates! graph_decomposition[%d] has\n"
	       "*ML*ERR* been setted (value = %d)\n",
	       i, j );
      exit( EXIT_FAILURE );
    }
    nodes_per_aggregate[j]++;
  }
  
}

int ML_Aggregate_VizAndStats_SetUpLevel( ML_Aggregate_Viz_Stats finer_level,
					 ML_Aggregate_Viz_Stats *coarser_level,
					 int dim ) 
{

  int Nlocal;
  size_t size;
  
  Nlocal = finer_level.Naggregates;

  size = sizeof(double)*Nlocal;
  ML_memory_alloc((void*)&(coarser_level->x),size,"x for info");
  if( dim > 1 ) ML_memory_alloc((void*)&(coarser_level->y),size,"y for info");
  if( dim > 2 ) ML_memory_alloc((void*)&(coarser_level->z),size,"z for info");

  return 0;  
}

  

/* ======================================================================== */
/*!
 \brief visualize aggregates and compute some statistics on the aggregates,
 like their diameter and number of nodes (for all levels).

 This function performs a post-processing of data. To use this function, the
 user should do the following:
 - before calling \c ML_Gen_MGHierarchy_UsingAggregation, he/she has to
   initialize a ML_Aggregate_Viz_Stats list (internally stored in the
   ML_Aggregate object). This is done as
   "ML_Aggregate_Viz_Stats_Create( ml_ag, MaxLevels )".
 - then, function `ML_Aggregate_CoarsenMETIS' will recognize that an
   ML_Aggregate_Viz_Stats object has been created, and she will put a copy
   of the graph decomposition into the ML_Aggre_Info object of that level
   (note that levels moves from 0 <fine> to # levels )
 - after a call to this function, he/she need to free memory, using
   "ML_Aggregate_Viz_Stats_Destroy"

 \date Albuquerque, 23-October-03  

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_VizAndStats_Compute( ML *ml, ML_Aggregate *ag,
				      int MaxMgLevels,
				      double *x, double *y, double *z,int dummy,
				      char *base_filename )
{

  int i, ilevel, iaggre;
  int Nlocal, Naggregates;
  char graphfile[132];
  ML_Aggregate_Viz_Stats * info;
  int dim;
  double dmin, davg, dmax, dstd;
  int  imin, iavg, imax;
  ML_Comm *comm;
  int finest_level = ml->ML_finest_level;
  int coarsest_level = ml->ML_coarsest_level;
  int incr_or_decr;
  int num_PDE_eqns  = ag->num_PDE_eqns;
  double h, H;
  int begin, end, diff;
  int radius;
  int * itemp = NULL, * itemp2 = NULL;
  double * dtemp = NULL, dsum;
  int Nrows, Naggregates_global, Nrows_global, offset;
  
  /* ------------------- execution begins --------------------------------- */
  
  info = ag->aggr_viz_and_stats;
  comm = ml->comm;
  
  if( finest_level > coarsest_level ) incr_or_decr = ML_DECREASING;
  else                                incr_or_decr = ML_INCREASING;

  /* ********************************************************************** */
  /* check out the number of dimensions                                     */
  /* ********************************************************************** */

  if( x == NULL )      dim = 0;           /* this means only stats on graph */
  else if( y == NULL ) dim = 1;
  else if( z == NULL ) dim = 2;
  else dim = 3;

  /* ********************************************************************** */
  /* copy ML_Operator pointers so that they can be used to viz graphs       */
  /* ********************************************************************** */

  for( i=0 ; i<MaxMgLevels ; i++ ) {
    info[i].Amatrix = &(ml->Amat[i]);
    if( info[i].Amatrix != NULL )
      ML_Operator_AmalgamateAndDropWeak(info[i].Amatrix, num_PDE_eqns, 0.0);  
  }

  /* ********************************************************************** */
  /* find out how many levels have been used                                */
  /* is_filled == ML_YES means that we have the graph_decomposition for     */
  /* level. However, we still miss the nodal coordinates for level>0, which */
  /* must be defined.  As AMG is without grid, I consider as "node" of      */
  /* levels>0 the coordinates of the center of gravity.                     */
  /* NOTE: I suppose that the first level is 0 (fine grid), then increasing */
  /* NOTE2: level 0 is set separately from the others                       */
  /* ********************************************************************** */

  info[finest_level].x = x;
  info[finest_level].y = y;
  info[finest_level].z = z;

  if( incr_or_decr == ML_INCREASING ) {

    if( comm->ML_mypid == 0 )  {
      puts("------------------------------------------------------------------------");
      printf( "Stats : levels are increasing, finest level is %d\n",
	      finest_level);
    }
    
    begin = finest_level;
    end   = coarsest_level;
    diff  = +1;

  } else {

    if( comm->ML_mypid == 0 ) {
      puts("------------------------------------------------------------------------");
      printf( "Stats : levels are decreasing, finest level is %d\n",
	      finest_level);
    }
    
    begin = coarsest_level+1;
    end   = finest_level+1;
    diff  = -1;

  }

  if( dim > 0 ) {
    switch( incr_or_decr ) {
    case ML_INCREASING:
      for( ilevel=finest_level ; ilevel<coarsest_level ; ilevel++ ) {

	Naggregates = info[ilevel].Naggregates;
	switch( dim ) {
	case 3:
	  info[ilevel+1].z = (double *)malloc(sizeof(double)*Naggregates);
	case 2:
	  info[ilevel+1].y = (double *)malloc(sizeof(double)*Naggregates);
	case 1:
	  info[ilevel+1].x = (double *)malloc(sizeof(double)*Naggregates);
	}
	ML_Aggregate_ComputeCenterOfGravity( info[ilevel],
					     info[ilevel+1],
					     comm);
      }
      break;
    case ML_DECREASING:
      for( ilevel=finest_level ; ilevel>coarsest_level ; ilevel-- ) {

	Naggregates = info[ilevel].Naggregates;
	
	switch( dim ) {
	case 3:
	  info[ilevel-1].z = (double *)malloc(sizeof(double)*Naggregates);
	case 2:
	  info[ilevel-1].y = (double *)malloc(sizeof(double)*Naggregates);
	case 1:
	  info[ilevel-1].x = (double *)malloc(sizeof(double)*Naggregates);
	}
	ML_Aggregate_ComputeCenterOfGravity( info[ilevel],
					     info[ilevel-1],
					     comm);
      }
      
      break;
    }
    
  }

  /* ********************************************************************** */
  /* statistics about the decomposition into subdomains                     */
  /* ********************************************************************** */

  Nrows = ml->Amat[finest_level].outvec_leng;

  imin = ML_gmin_int(Nrows,comm);
  iavg = ML_gsum_int(Nrows,comm)/(comm->ML_nprocs);
  imax = ML_gmax_int(Nrows,comm);
  
  if(  comm->ML_mypid == 0 ) {
    puts("------------------------------------------------------------------------");
    printf( "Stats : rows per process (min) = %d\n", imin);
    printf( "Stats : rows per process (avg) = %d\n", iavg);
    printf( "Stats : rows per process (max) = %d\n", imax);
  }
  
  if( dim > 0 ) {
    
    ML_Info_DomainDecomp( info[finest_level], comm, &H, &h );
    
    ML_Aggregate_AnalyzeVector( 1, &H,
				&dmin, &dmax, &davg, &dstd, comm );

    if(  comm->ML_mypid == 0 ) {
      puts("------------------------------------------------------------------------");
      printf( "Stats : subdomain linear dimension (min) = %f\n",
	    dmin );
      printf( "Stats : subdomain linear dimension (avg) = %f\n",
	      davg );
      printf( "Stats : subdomain linear dimension (max) = %f\n",
	      dmax );
    }

    ML_Aggregate_AnalyzeVector( 1, &h,
				&dmin, &dmax, &davg, &dstd, comm );
      
    if(  comm->ML_mypid == 0 ) {
      puts("------------------------------------------------------------------------");
      printf( "Stats : element linear dimension (min) = %f\n",
	    dmin );
      printf( "Stats : element linear dimension (avg) = %f\n",
	      davg );
      printf( "Stats : element linear dimension (max) = %f\n",
	      dmax );
    }

  }

  /* ********************************************************************** */
  /* Statistics about the ratio aggregates/nodes for each level (globally)  */
  /* ********************************************************************** */

  for( ilevel=begin ; ilevel<end ; ilevel++ ) {

    if( info[ilevel].is_filled == ML_YES ) {

      	Nlocal = info[ilevel].Nlocal;
	Naggregates = info[ilevel].Naggregates;
	radius = info[ilevel].graph_radius;

	Nrows_global = ML_gsum_int(Nlocal,comm);
	radius = ML_gsum_int(radius,comm);

	switch( info[ilevel].local_or_global ) {
	case ML_LOCAL_INDICES:
	  Naggregates_global = ML_gsum_int(Naggregates,comm);
#ifdef ML_MPI
	  MPI_Scan(&offset, &Naggregates, 1, MPI_INT, MPI_SUM, comm->USR_comm);
	  offset -= Naggregates;
#else
	  offset = 0;
#endif
	  break;
	case ML_GLOBAL_INDICES:
	  Naggregates_global = Naggregates;
	  offset = 0;
	  break;
	}

	/* computes how many nodes are in each aggregate for this level,
	   globally */
	itemp = (int *) malloc( sizeof(int) * Naggregates_global );
	for( i=0 ; i<Naggregates_global ; i++ ) itemp[i] = 0;
	
	for( i=0 ; i<Nlocal ; i++ ) {

	  iaggre = info[ilevel].graph_decomposition[i]+offset;

	  if( iaggre >= Naggregates_global ) {
	    fprintf(stderr,
		    "*ML*ERR* error : %d >= %d\n"
		    "*ML*ERR* (file %s, line %d\n",
		    iaggre,
		    Naggregates_global,
		    __FILE__,
		    __LINE__ );
	    exit( EXIT_FAILURE );
	  }
	  itemp[iaggre]++;


	}
	
#ifdef ML_MPI
	itemp2 = (int *)malloc( sizeof(int)*Naggregates_global);
	MPI_Reduce(itemp,itemp2,Naggregates_global,MPI_INT,MPI_SUM,0,
		   comm->USR_comm);
#else
	itemp2 = itemp;
#endif

	if( comm->ML_mypid == 0 ) {

	  imin = INT_MAX;
	  imax = INT_MIN;
	  
	  for( i=0 ; i<Naggregates_global ; i++ ) {
	    if( itemp2[i] > imax ) imax =  itemp2[i];
	    if( itemp2[i] < imin ) imin =  itemp2[i];
	  }
		 
	  puts("------------------------------------------------------------------------");
    
	  printf( "Stats (level %d) : NumAggr = %5d, NumNodes = %d\n",
		  ilevel,
		  Naggregates_global,Nrows_global);
	  printf( "Stats (level %d) : NumAggr/NumNodes  (avg)   = %7.5f %%\n",
		  ilevel,
		  100.0*Naggregates_global/Nrows_global );
	  printf( "Stats (level %d) : NumNodes per aggr (min)   = %d\n",
		  ilevel,
		  imin);
	  printf( "Stats (level %d) : NumNodes per aggr (avg)   = %d\n",
		  ilevel,
		  Nrows_global/Naggregates_global);
	  printf( "Stats (level %d) : NumNodes per aggr (max)   = %d\n",
		  ilevel,
		  imax);
	  if( info[ilevel].local_or_global == ML_LOCAL_INDICES ) 
	    printf( "Stats (level %d) : graph radius      (avg)   = %f\n",
		    ilevel,
		    1.0*radius/(comm->ML_nprocs) );
	
	}
    }
	
    ML_free(itemp);
#ifdef ML_MPI
    ML_free(itemp2);
#endif
	     

  }

  /* ********************************************************************** */
  /* some statistics on the aggregates. Note that I need some nodal         */
  /* coordinates to perform this task. The nodal coordinates of the         */
  /* current level aggregates are stored in level+1. This means that the    */
  /* last level will hold only coordinates, and no graph_decomposition.     */
  /* I allocate/deallocate R and H so that their shape is right for the     */
  /* level we are considering (and not waste space)                         */
  /* RorH is a double vector which will contain the radius of each          */
  /* aggregate, and then its linear size.                                   */
  /* ********************************************************************** */

#ifdef COMPUTE_RADIUS
    
    /* this works only for ML_LOCAL INDICES */
    for( ilevel=begin ; ilevel<end ; ilevel++ ) {
      
      Nlocal = info[ilevel+diff].Nlocal;
      Naggregates = info[ilevel].Naggregates;

      RorH = (double *) malloc( sizeof(double) * Naggregates );
    
      if( RorH == NULL ) {
	fprintf( stderr,
		 "*ML*ERR* not enough memory for %d bytes\n"
		 "*ML*ERR* (file %s, line %d)\n",
		 sizeof(double) * Nlocal,
		 __FILE__,
		 __LINE__ );
	exit( EXIT_FAILURE );
      }
    
      ML_Aggregate_ComputeRadius( info[ilevel], info[ilevel+diff], RorH );
      
      ML_Aggregate_AnalyzeVector( Naggregates, RorH,
				  &dmin, &dmax, &davg, &dstd, comm );
      
      if( comm->ML_mypid == 0 ) {
	puts("------------------------------------------------------------------------");
	printf( "Stats : aggregate radius (min) = %f\n",
		dmin );
	printf( "Stats : aggregate radius (avg) = %f\n",
		davg );
	printf( "Stats : aggregate radius (max) = %f\n",
		dmax );
      }
      
    }    
#endif
  if( dim > 0 ) {

    for( ilevel=begin ; ilevel<end ; ilevel++ ) {

      if( info[ilevel].is_filled == ML_YES ) {

	Naggregates = info[ilevel].Naggregates;
	radius = info[ilevel].graph_radius;

	switch( info[ilevel].local_or_global ) {
	case ML_LOCAL_INDICES:
	  Naggregates_global = ML_gsum_int(Naggregates,comm);
#ifdef ML_MPI
	  MPI_Scan(&Naggregates,&offset, 1, MPI_INT, MPI_SUM, comm->USR_comm);
	  offset -= Naggregates;
#else
	  offset = 0;
#endif
	  break;
	case ML_GLOBAL_INDICES:
	  Naggregates_global = Naggregates;
	  offset = 0;
	  break;
	}

	dtemp = (double *) malloc( sizeof(double)*Naggregates_global);
	ML_Aggregate_ComputeBox( info[ilevel], Naggregates_global,
				 dtemp, offset,comm);
      
	if( comm->ML_mypid == 0 ) {
	  
	  dmin =  DBL_MAX;
	  dmax = -DBL_MAX;
	  dsum = 0;
	  
	  for( i=0 ; i<Naggregates_global ; i++ ) {
	    if( dtemp[i] > dmax ) dmax =  dtemp[i];
	    if( dtemp[i] < dmin ) dmin =  dtemp[i];
	    dsum += dtemp[i];
	  }
	  
	  puts("------------------------------------------------------------------------");
	  printf( "Stats (level %d) : aggregate linear dimension (min) = %f\n",
		  ilevel, dmin );
	  printf( "Stats (level %d) : aggregate linear dimension (avg) = %f\n",
		  ilevel,
		  dsum/Naggregates_global );
	  printf( "Stats (level %d) : aggregate linear dimension (max) = %f\n",
		  ilevel,
		  dmax );
	}
	ML_free(dtemp);
      }
    }
  }
  
  /* ********************************************************************** */
  /* output to OpenDX. This required the nodal coordinates on the starting  */
  /* grid, as provided using x, y and z (y and z can be NULL)               */
  /* ********************************************************************** */
  
  if( dim > 0 ) {
    
    for( ilevel=begin ; ilevel<end ; ilevel++ ) {
      
      if( info[ilevel].is_filled == ML_YES ) {
	if( base_filename == NULL ) {
	  sprintf( graphfile,
		   ".graph_level%d_proc",
		   ilevel );
	} else {
	  sprintf( graphfile,
		   "%s_level%d_proc",
		   base_filename,
		   ilevel );
	}
	
	if( comm->ML_mypid == 0 ) {
	  puts("------------------------------------------------------------------------");
	  printf("Viz (level %d) : Writing OpenDX file `%s'\n",
		 ilevel, graphfile );
	}
	ML_Aggregate_VisualizeWithOpenDX( info[ilevel], graphfile,
					  comm );
      }
      
    }
  }
  
  /* ********************************************************************** */
  /* Clear memory allocated while creating the center of gravity (but not   */
  /* for the finest level, as those arrays are provided by the user). So I  */
  /* put to NULL those pointers in the info array.                          */
  /* ********************************************************************** */

  info[finest_level].x = NULL;
  info[finest_level].y = NULL;
  info[finest_level].z = NULL;

  for( ilevel=0 ; ilevel<MaxMgLevels ; ilevel++ ) {

    if( info[ilevel].x != NULL ) ML_free( info[ilevel].x );
    if( info[ilevel].y != NULL ) ML_free( info[ilevel].y );
    if( info[ilevel].z != NULL ) ML_free( info[ilevel].z );

    if( info[ilevel].Amatrix != NULL )
      ML_Operator_UnAmalgamateAndDropWeak(info[ilevel].Amatrix, num_PDE_eqns,
					  0.0);
    
  }

  if( comm->ML_mypid == 0 )
    puts("------------------------------------------------------------------------");
    
  return 0;
  
} /* ML_Aggregate_Visualize */

/* ======================================================================== */
/*!
 \brief Information about the decomposition into subdomains

*/
/* ------------------------------------------------------------------------ */

int ML_Info_DomainDecomp( ML_Aggregate_Viz_Stats info,
			  ML_Comm *comm, double *H, double *h )
{

  int j, irow, col;  
  ML_Operator *Amatrix = (ML_Operator *)(info.Amatrix);
  int N_dimensions;
  double *x = info.x;
  double *y = info.y;
  double *z = info.z;
  int Nrows = Amatrix->getrow->Nrows;
  int allocated = 0;
  int * rowi_col = NULL;
  int rowi_N;
  double * rowi_val = NULL;
  double xmin, xmax, ymin, ymax, zmin, zmax, h_row;
  double x_row, x_col, y_row, y_col, z_row, z_col;
  
  /* ------------------- execution begins --------------------------------- */

  xmin = DBL_MAX;
  xmax = -DBL_MAX;
  ymin = DBL_MAX;
  ymax = - DBL_MAX;
  zmin = DBL_MAX;
  zmax = - DBL_MAX;
  *h = 0.0;
  
  if( y == NULL ) N_dimensions = 1;
  else if( z == NULL ) N_dimensions = 2;
  else N_dimensions = 3;

  z_col = 0.0;
  z_row = 0.0;
  y_col = 0.0;
  z_col = 0.0;
  
  for( irow=0 ; irow<Nrows ; irow++ ) {

    if( z != NULL ) z_row = z[irow];
    if( y != NULL ) y_row = y[irow];

    x_row = x[irow];
 
    ML_get_matrix_row(Amatrix, 1, &irow, &allocated, &rowi_col, &rowi_val,
		      &rowi_N, 0);

    for( j=0 ; j<rowi_N ; j++ ) {
      col = rowi_col[j];
      switch( N_dimensions ) {
      case 3:
	z_col = z[col];
	zmax = ML_max(zmax,z_col);
	zmin = ML_min(zmin,z_col);
      case 2:
	y_col = y[col];
	ymax = ML_max(ymax,y_col);
	ymin = ML_min(ymin,y_col);
      case 1:
	x_col = x[col];
	xmax = ML_max(xmin,x_col);
	xmin = ML_min(xmin,x_col);
      }

      h_row = sqrt( pow(x_row-x_col+0.00000001,2) +
		    pow(y_row-y_col+0.00000001,2) +
		    pow(z_row-z_col+0.00000001,2) );

      if( h_row > (*h) ) *h = h_row;
      
    }
    
  }

  /* ********************************************************************** */
  /* find out the linear dimension of the domain                            */
  /* ********************************************************************** */

  *H = 0.0;
  *H = ML_max( *H, xmax-xmin );
  if( N_dimensions >1 ) *H = ML_max( *H, ymax-ymin );
  if( N_dimensions >2 ) *H = ML_max( *H, zmax-zmin );
  
  /* ------------------- that's all folks --------------------------------- */

  ML_free(rowi_col); ML_free(rowi_val);
  rowi_col = NULL; rowi_val = NULL;
  allocated = 0; 
  
  return 0;
  
} /* ML_VisualizeWithOpenDX */

/* ********************************************************************** */
/* used in `ML_Compute_AggregateGraphRadius'                              */
/* ********************************************************************** */

static int get_max_dep_line( int Nrows, int ia[], int ja[],
			     int dep [] ) 
{
  int i, j, col, ok, contour;

  contour = 1;
  ok = 0;
  
  while( ok == 0 ) {
    
    ok = 1;
    for( i=0 ; i<Nrows ; i++ ) {
      if( dep[i] == contour-1 ) {
	for( j=ia[i] ; j<ia[i+1] ; j++ ) {
	  col = ja[j];
	  if( dep[col] == -7 ) {
	    dep[col] = contour;
	    ok = 0;
	  }
	}
      }
    }
    if( ok == 0 ) contour++;

  }

  return( --contour );

} 


/* ======================================================================== */
/*!
 \brief Compute the radius an aggregates, coded in a CSR matrix.

 This function computes the radius (in matrix-graph sense) of a set of
 nodes, coded in the CSR vectors ia and ja.
 \param dep in : this vector, of size Nrows, is defined as follows:
  - dep[i] == -7 for a boundary node
  - dep[i] ==  0 for non-boundary ndes
  
*/
/* ------------------------------------------------------------------------ */

int ML_Compute_AggregateGraphRadius( int Nrows, int ia[], int ja[],
				     int dep [],
				     int *pradius, int *pNcenter ) 
{
  int i, j, radius, Ncenter;
  int max_dep;
  int * center;
  int * orig_dep = (int*)malloc(sizeof(int)*Nrows);

  for( i=0 ; i<Nrows ; i++ )
    orig_dep[i] = dep[i];
  
  max_dep = get_max_dep_line( Nrows, ia, ja, dep );
  /* define the center nodes */

  Ncenter = 0;
  center = (int *) malloc( sizeof(int) * Nrows );
  
  for( i=0 ; i<Nrows ; i++ ) {
    if( dep[i] == max_dep ) {
      center[Ncenter] = i;
      Ncenter++;
    }
  }

  /* now for all center nodes compute the radius */

  radius = 0;
  for( i=0 ; i<Ncenter ; i++ ) {
    for( j=0 ; j<Nrows ; j++ ) {
      if( orig_dep[j] == 0 )  dep[j] = -1;
      else                    dep[j] = -7;
    }
    dep[center[i]] = 0;
    j = get_max_dep_line( Nrows, ia, ja, dep );
    
    if( j > radius ) radius = j;
  }

  if( radius < max_dep ) {
    fprintf( stderr,
	     "*ML*ERR* error in `ML_Compute_AggregateGraphRadius'\n"
	     "*ML*ERR* radius < max_dep ( %d - %d )\n",
	     radius,
	     max_dep );
  }
  
  *pradius = radius;
  *pNcenter = Ncenter;
  
  ML_free( center );
  ML_free( orig_dep ); orig_dep = NULL;

  return 0;
  
}
