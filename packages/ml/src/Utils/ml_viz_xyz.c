/* ======================================================================== */
/*!
 \file ml_viz_xyz.c

 \brief Prints out information in a simple XYZ format

 \author Marzio Sala, SNL, 9214

 \date 28-May-2003
 
*/
/* ------------------------------------------------------------------------ */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "ml_viz_stats.h"
#include "ml_viz_xyz.h"

/* ======================================================================== */
/*!
 \brief write graph decomposition of the current level in a graphical
 format readable by XD3D

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
  double *z = info.z;
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
  int * reorder = NULL, ok;

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

  /* may need to reshuffle the aggregate ordering */
  if( vector == NULL && str == NULL) {

    reorder = (int *) ML_allocate( sizeof(int) * Naggregates );

    if( reorder == NULL ) {
      fprintf( stderr,
	       "*ML*ERR* not enough memory for %d bytes\n"
	       "*ML*ERR* (file %s, line %d)\n",
	       (int)sizeof(int) * Naggregates,
	       __FILE__,
	       __LINE__ );
      exit( EXIT_FAILURE );
    }
    /* sometimes ML_SHUFFLE is a mess, may need to be
     * deactivated */
#define ML_SHUFFLE
#ifdef ML_SHUFFLE
    for( i=0 ; i<Naggregates ; ++i ) reorder[i] = -1;

    srand(0); /* Initialize random seed */
    
    for( i=0 ; i<Naggregates ; ++i ) {

      do {

	ok = 0;
	
	j = (int)(1.0*(Naggregates)*rand()/RAND_MAX);

	if( j >=0 && j<Naggregates ) {
	  if (reorder[j] == -1) {
	    reorder[j] = i;
	    ok = 1;
	  }
	}
      } while( ok == 0 );

    }

    for( i=0 ; i<Naggregates ; ++i ) {
      if( reorder[i] < 0 || reorder[i] >= Naggregates ) {
	fprintf( stderr,
		"*ML*ERR* reorder failed.\n"
		"*ML*ERR* (file %s, line %d)\n",
	       __FILE__,
	       __LINE__ );
	exit( EXIT_FAILURE );
      }
    }
#else
    for (i = 0 ; i < Naggregates ; ++i)
      reorder[i] = i;
#endif
  }

  /* cycle over all local rows, plot corresponding value on file */

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
	  val = ipid +  nprocs * 1.0*graph_decomposition[irow];

	/* XD3D does not work in 3D, but maybe other codes will */
	if( z == NULL ) 
	  fprintf( fp,
		  "%f %f %f\n",
		  x[irow], y[irow], val );
	else
	  fprintf( fp,
		  "%f %f %f %f\n",
		  x[irow], y[irow], z[irow], val );

      }
      fclose(fp);
    }
#ifdef ML_MPI
    MPI_Barrier(comm->USR_comm);
#endif
  }
  
  if( reorder != NULL ) free(reorder);

  return 0;
  
} /* ML_VisualizeWithXYZ */

/* ======================================================================== */
/* ------------------------------------------------------------------------ */

int ML_PlotXYZ(int Npoints, double* x, double* y, double* z,
	       char base_filename[],
	       USR_COMM comm, double * vector)
{

  int irow;
  int ipid;
  int mypid = 0;
  int nprocs = 1;
  FILE *fp;
  double val;
  char filemode[2];

  /* ------------------- execution begins --------------------------------- */

#ifdef ML_MPI
  MPI_Comm_rank(comm,&mypid);
  MPI_Comm_size(comm,&nprocs);
#endif

  if( mypid == 0 ) filemode[0] = 'w';
  else             filemode[0] = 'a';
  filemode[1] = '\0';

  /* cycle over all local rows, plot corresponding value on file */

  for (ipid = 0 ; ipid < nprocs ; ++ipid) {
    if (ipid == mypid) {
      if ((fp = fopen( base_filename, filemode )) == NULL) {
	fprintf( stderr,
		"*ML*ERR* cannot open file `%s'\n",
		base_filename );
	exit( EXIT_FAILURE );
      }

      for (irow = 0 ; irow < Npoints ; ++irow) {
	val = vector[irow];

	if( z == NULL ) 
	  fprintf( fp,
		  "%f %f %f\n",
		  x[irow], y[irow], val );
	else
	  fprintf( fp,
		  "%f %f %f %f\n",
		  x[irow], y[irow], z[irow], val );

      }
      fclose(fp);
    }
#ifdef ML_MPI
    MPI_Barrier(comm);
#endif
  }
  
  return 0;
  
} /* ML_PlotXYZ */


