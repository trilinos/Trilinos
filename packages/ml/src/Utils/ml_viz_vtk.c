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
   written out by a parallel application).  I believe this requires the use
   of xml.  Right now, I'm just writing to one plain, ASCII file (not xml).
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
  int maxAggSize = 1000;
  int *index;
  int myLocalAgg;
  int NPosConnections; /* number of edge connections (in matrix graph) from
                          given unknown to unknown with a larger index --
                          corresponds to nnz in upper triangular part of
                          local matrix */
  int NGlobPosConnections;
  int NCells, NumPtsInCells, NGlobCells;
  int *bindx, allocated, row_length;
  double *values;
  int *GlobRowOrdering;
  int NGlobPtsInCells;

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
    printf("\n********************************************************************************\n\aWARNING: you have asked to visualize only aggregate %d by setting ML_VIZ_AGGR\n********************************************************************************\n",AggrToVisualize);
  }

  /* may need to reshuffle the aggregate ordering */
/*
  if( vector == NULL && str == NULL) {
*/

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
    /* JJH I believe ML_SHUFFLE is intended to improve color contrast in
       visualization by separating close aggregates.  This improves the chance
       that they're colored differently. */
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
/*
  }
*/

  /* Calculate the global number of vertices */
  Nglobrows = Nrows;
  ML_gsum_scalar_int(&Nglobrows, &i, comm);

  /* Impose a global numbering on the matrix. */

  ML_build_global_numbering(Amatrix, comm, &GlobRowOrdering);

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
    strncpy(comment,"Unstructured grid.  Matrix connectivity is given by the vtk line cell type (cell type 3).  Aggregates are defined by the vtk poly-vertiex cell type (cell type 2).  Function values are given as point scalar values.\n",(size_t) 256);
    fprintf(fp,"%s",comment);
    fprintf(fp,"ASCII\n");
#ifdef UNSTRUCTURED
    fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");
#else
    fprintf(fp,"DATASET POLYDATA\n");
#endif
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
        if( vector != NULL ) {
          val = vector[irow];
          myLocalAgg = graph_decomposition[irow];
        }
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

  i = Amatrix->max_nz_per_row;
  allocated = (i > 0) ? i : 10;
  bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
  values = (double *)  ML_allocate( allocated*sizeof(double));

  NPosConnections = 0;
  for (i = 0 ; i < Nrows; i++) {
    ML_get_matrix_row(Amatrix, 1, &i, &allocated, &bindx, &values,
            &row_length, 0);
    for  (j = 0; j < row_length; j++) {
      if (GlobRowOrdering[i] < GlobRowOrdering[bindx[j]])
        NPosConnections++;
    }
  }

/*
  printf("\n\n******* NPosConnections = %d\n\n\n",NPosConnections);
*/
  NGlobPosConnections = NPosConnections;
  ML_gsum_scalar_int(&NGlobPosConnections , &NPosConnections, comm);
/*
  printf("\n\n******* NGlobPosConnections = %d\n\n\n",NGlobPosConnections);
  printf("\n\n******* NGlobAggregates = %d\n\n\n",Nglobaggregates);
*/
#ifdef UNSTRUCTURED
  NGlobCells = Nglobaggregates+NGlobPosConnections;
  NGlobPtsInCells = Nglobrows             /* rows are in disjoint aggregates */
               +   2*NGlobPosConnections; /* edges have two vertices each */
#else
  NGlobCells = NGlobPosConnections;
  NGlobPtsInCells = 2 * NGlobPosConnections;
#endif

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
#ifdef UNSTRUCTURED
    fprintf(fp,"CELLS %d %d\n",NGlobCells, NGlobCells + NGlobPtsInCells);
#else
    fprintf(fp,"LINES %d %d\n",NGlobCells, NGlobCells + NGlobPtsInCells);
#endif
    fclose(fp);
  }
#ifdef ML_MPI
  MPI_Barrier(comm->USR_comm);
#endif

  /********************************************
    Write out connectivity information. 
  *********************************************/

  for (ipid=0 ; ipid < nprocs ; ++ipid)
  {
    if (ipid == mypid)
    {
      if ((fp = fopen( base_filename, filemode )) == NULL ) {
        fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
        exit( EXIT_FAILURE );
      }
      for (i = 0 ; i < Nrows; i++) {
        ML_get_matrix_row(Amatrix, 1, &i, &allocated, &bindx, &values,
                &row_length, 0);
        for  (j = 0; j < row_length; j++) {
          if (GlobRowOrdering[i] < GlobRowOrdering[bindx[j]])
            fprintf(fp,"2  %d  %d\n",
                    GlobRowOrdering[i], GlobRowOrdering[bindx[j]]);
        }
      }
      fclose(fp);
    }
#ifdef ML_MPI
    MPI_Barrier(comm->USR_comm);
#endif
  }

  /****************************************************
    Write out polyvertex (aggregate) information. 
  ****************************************************/

#ifdef UNSTRUCTURED
  for (ipid=0 ; ipid < nprocs ; ++ipid)
  {
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
#endif /*ifdef UNSTRUCTURED*/

  ML_free(aggregates);
  ML_free(index);

  /**********************************************************************
    Write out cell types (for right now, just line and polyvertex types)
  **********************************************************************/
#ifdef UNSTRUCTURED
  if (mypid == 0) {
    if ((fp = fopen( base_filename, filemode )) == NULL ) {
      fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
      exit( EXIT_FAILURE );
    }
    fprintf(fp,"CELL_TYPES %d\n",NGlobCells);
    cell_type = 3; /*VTK code for a line*/
    for (i=0; i<NGlobPosConnections; i++)
      fprintf(fp,"%d\n",cell_type);
    cell_type = 2; /*VTK code for a polyvertex*/
    for (i=0; i<Nglobaggregates; i++)
      fprintf(fp,"%d\n",cell_type);
    fclose(fp);
  }
#endif

  /******************************************************/
  /* Write out vertex data (there can be multiple sets) */
  /******************************************************/

  /* which aggregate each node belongs to */
  if (mypid == 0) {
    if ((fp = fopen( base_filename, filemode )) == NULL ) {
      fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
      exit( EXIT_FAILURE );
    }
    fprintf(fp,"POINT_DATA %d\n",Nglobrows);
    fprintf(fp,"SCALARS aggregate_number int 1\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    fclose(fp);
  }

  for (ipid=0 ; ipid < nprocs ; ++ipid)
  {
    if (ipid == mypid)
    {
      if ((fp = fopen( base_filename, filemode )) == NULL ) {
        fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
        exit( EXIT_FAILURE );
      }
      for (irow=0 ; irow<Nrows ; irow++ )
      {
        /* the user can concentrate on just one aggregate */
        if( AggrToVisualize != -1 ) {
          if( graph_decomposition[irow] == AggrToVisualize )
            val = 1.0;
          else
            val = 0.0;
        }
        else 
          val = (double)(reorder[graph_decomposition[irow]] + offset);

        fprintf(fp,"%d\n",(int) val);   
      }
      fclose(fp);
    }
#ifdef ML_MPI
    MPI_Barrier(comm->USR_comm);
#endif
  }

  /* any scalar value on the nodes */
  if( vector != NULL )
  {
    if (mypid == 0) {
      if ((fp = fopen( base_filename, filemode )) == NULL ) {
        fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
        exit( EXIT_FAILURE );
      }
      fprintf(fp,"SCALARS residual float 1\n");
      fprintf(fp,"LOOKUP_TABLE default\n");
      fclose(fp);
    }
  
    for (ipid=0 ; ipid < nprocs ; ++ipid)
    {
      if (ipid == mypid)
      {
        if ((fp = fopen( base_filename, filemode )) == NULL ) {
          fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
          exit( EXIT_FAILURE );
        }
        for (irow=0 ; irow<Nrows ; irow++ )
          fprintf(fp,"%lf\n",vector[irow]);
        fclose(fp);
      }
#ifdef ML_MPI
      MPI_Barrier(comm->USR_comm);
#endif
    }
  }



#ifdef NOTYETIMPLEMENTED
  if (mypid == 0) {
    if ((fp = fopen( base_filename, filemode )) == NULL ) {
      fprintf( stderr, "*VIZ*ERR* cannot open file `%s'\n", base_filename );
      exit( EXIT_FAILURE );
    }
    fprintf(fp,"CELL_DATA %d\n",NglobPts);
    fprintf(fp,"SCALARS cell_scalars int 1\n");
    fprintf(fp,"LOOKUP_TABLE default\n");
    for (i=0; i<Nglobaggregates; i++)
      fprintf(fp,"%d\n",i);
    fclose(fp);
  }
#endif
  ML_free(GlobRowOrdering);
  ML_free(values);
  ML_free(bindx);
  
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


