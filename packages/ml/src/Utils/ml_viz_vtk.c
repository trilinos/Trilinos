/* ======================================================================== */
/*!
 \file ml_viz_vtk.c

 \brief Prints out information in VTK format, readable by Paraview.

 \author Jonathan Hu, SNL, 9214

 \date 8-March-2005
 
*/
/* ------------------------------------------------------------------------ */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "ml_viz_stats.h"
#include "ml_viz_vtk.h"

/* ======================================================================== */
/*!
 \brief write graph decomposition of the current level in an ASCII
 format readable by Paraview

*/
/* ------------------------------------------------------------------------ */

int ML_Aggregate_VisualizeVTK( ML_Aggregate_Viz_Stats info,
                  char base_filename[],
                  ML_Comm *comm,
                  double * vector)
{

/* FIXME JJH Paraview supports visualization from multiple files (i.e., files
   written out by a parallel application.  Right now, we're just writing to one
   file.
*/

  int i,j, k, irow, ipid;
  ML_Operator *Amatrix = (ML_Operator *)(info.Amatrix);
  double *x = info.x;
  double *y = info.y;
  double *z = info.z;
  int *graph_decomposition = info.graph_decomposition;
  int Naggregates = info.Naggregates;
  int Nglobaggregates;
  int mypid = comm->ML_mypid;
  int nprocs = comm->ML_nprocs;
  int Nrows = Amatrix->getrow->Nrows;
  int Nglobrows;
  FILE *fp;
  double val;
  int AggrToVisualize = -1;
  char * str;
  char filemode[2];
  char comment[257];
  int Nlocal = info.Nlocal;
  int * reorder = NULL, ok;
  int offset;
  int vertex_offset;
  int cell_type;
  int *aggregates;
  int maxAggSize = 600;
  int *index;
  int myLocalAgg;

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

  str = getenv("ML_VIZ_AGGR");
  if( str != NULL ) {
    AggrToVisualize = atoi(str);
    printf("\n\n\aWARNING:  you have asked to visualize only aggregate %d by setting ML_VIZ_AGGR\n\n",AggrToVisualize);
  }

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
#if 0
    j = 0;
    for( i=0 ; i<Naggregates ; ++i ) {
      if (reorder[i] == i)
        ++j;
    }
    printf("On processor %d, reorder effectiveness = %f\n",
           comm->ML_mypid, (double)j / Naggregates);
#endif
#else
    for (i = 0 ; i < Naggregates ; ++i)
      reorder[i] = i;
#endif
  }

  /* Calculate the global number of vertices */
  Nglobrows = Nrows;
  ML_gsum_scalar_int(&Nglobrows, &i, comm);

  /********************************/
  /* Write out header information */
  /********************************/
  sprintf(filemode,"w");
  if (mypid == 0) {
    if ((fp = fopen( base_filename, filemode )) == NULL ) {
      fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
      exit( EXIT_FAILURE );
    }
    fprintf(fp,"# vtk DataFile Version 2.0\n");
    /* comment line, 256 characters max */
    strncpy(comment,"unstructured grid with aggregates represented as poly-vertices (cell type 2)\n",(size_t) 256);
    fprintf(fp,"%s",comment);
    fprintf(fp,"ASCII\n");
    fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
    /* number of vertices*/
    fprintf(fp,"POINTS %d float\n",Nglobrows);
    fclose(fp);
  }

#ifdef ML_MPI
  MPI_Barrier(comm->USR_comm);
  MPI_Scan (&Naggregates, &offset, 1, MPI_INT, MPI_SUM, comm->USR_comm);
  offset -= Naggregates;
#else
  offset = 0;
#endif  

#ifdef ML_MPI
  MPI_Barrier(comm->USR_comm);
  MPI_Scan (&Nrows, &vertex_offset, 1, MPI_INT, MPI_SUM, comm->USR_comm);
  vertex_offset -= Nrows;
#else
  vertex_offset = 0;
#endif  

  aggregates = (int *) ML_allocate(maxAggSize * Naggregates * sizeof(int));
  index = (int *) ML_allocate(Naggregates * sizeof(int));

  for (i=0; i < maxAggSize*Naggregates; i++) aggregates[i] = -1;
  for (i=0; i < Naggregates; i++) index[i] = 0;

  sprintf(filemode,"a");
  /* cycle over all local rows, write coordinates to file */
  for (ipid=0 ; ipid < nprocs ; ++ipid) {
    if (ipid == mypid)
    {
      if ((fp = fopen( base_filename, filemode )) == NULL ) {
        fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
        exit( EXIT_FAILURE );
      }

      for (irow=0 ; irow<Nrows ; irow++ )
      {
        if( vector != NULL ) val = vector[irow];
        else if( AggrToVisualize != -1 ) {
          if( graph_decomposition[irow] == AggrToVisualize ) {
            val = 1.0;
            myLocalAgg = 1;
          }
          else {
            val = 0.0;
            myLocalAgg = 0;
          }  
        }
        else {
          val = (double)(reorder[graph_decomposition[irow]] + offset);
          myLocalAgg = reorder[graph_decomposition[irow]];
        }

        aggregates[maxAggSize*(myLocalAgg) + index[myLocalAgg]++] =
            irow +vertex_offset; 

        if( z == NULL ) 
          fprintf( fp, "%f %f 0.0\n", x[irow], y[irow]);
        else
          fprintf( fp, "%f %f %f\n", x[irow], y[irow], z[irow]);
      }
      fclose(fp);
    }
#ifdef ML_MPI
    MPI_Barrier(comm->USR_comm);
#endif
  }

  /* Calculate global number of aggregates. */
  Nglobaggregates = Naggregates;
  ML_gsum_scalar_int(&Nglobaggregates, &i, comm);

  /********************************************
     Write out aggregate information in form:

         N n1 n2 n3 ...

     where N is number of vertices in aggregate
     and n1,... are the vertex numbers.
  *********************************************/

  if (mypid == 0) {
    if ((fp = fopen( base_filename, filemode )) == NULL ) {
      fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
      exit( EXIT_FAILURE );
    }
    fprintf(fp,"CELLS %d %d\n",Nglobaggregates, Nglobaggregates+Nglobrows);
    fclose(fp);
  }
#ifdef ML_MPI
  MPI_Barrier(comm->USR_comm);
#endif

  for (ipid=0 ; ipid < nprocs ; ++ipid) {
    if (ipid == mypid)
    {
      if ((fp = fopen( base_filename, filemode )) == NULL ) {
        fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
        exit( EXIT_FAILURE );
      }

      for (i=0; i< Naggregates; i++) {
        fprintf(fp,"%d ",index[i]);
        for (j=0; j<index[i]; j++)
          fprintf(fp,"%d ",aggregates[i*maxAggSize + j]);
        fprintf(fp,"\n"); 
      }
      fclose(fp);
    }
#ifdef ML_MPI
    MPI_Barrier(comm->USR_comm);
#endif
  }

  ML_free(aggregates);
  ML_free(index);

  /**************************************************************/
  /* Write out cell types (for right now, just polyvertex type) */
  /**************************************************************/
  if (mypid == 0) {
    if ((fp = fopen( base_filename, filemode )) == NULL ) {
      fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
      exit( EXIT_FAILURE );
    }
    fprintf(fp,"CELL_TYPES %d\n",Nglobaggregates);
    /* TODO JJH we could put in support for more generic cell types */
    cell_type = 2; /*VTK code for a polyvertex*/
    for (i=0; i<Nglobaggregates; i++)
      fprintf(fp,"%d\n",cell_type);
    fclose(fp);
  }

  /*************************/
  /* Write out color table */
  /*************************/
  if (mypid == 0) {
    if ((fp = fopen( base_filename, filemode )) == NULL ) {
      fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
      exit( EXIT_FAILURE );
    }
    fprintf(fp,"CELL_DATA %d\n",Nglobaggregates);
    fprintf(fp,"SCALARS cell_scalars int 1\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (i=0; i<Nglobaggregates; i++)
      fprintf(fp,"%d\n",i);
    fclose(fp);
  }
  
  if( reorder != NULL ) ML_free(reorder);

  return 0;
  
} /* ML_VisualizeWithVTK */

/* ======================================================================== */
/* ------------------------------------------------------------------------ */

int ML_PlotVTK(int Npoints, double* x, double* y, double* z,
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
  
} /* ML_PlotVTK */


