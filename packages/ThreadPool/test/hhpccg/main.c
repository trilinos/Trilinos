/*------------------------------------------------------------------------*/
/*                    TPI: Thread Pool Interface                          */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/* Redistribution and use in source and binary forms, with or without     */
/* modification, are permitted provided that the following conditions are */
/* met:                                                                   */
/*                                                                        */
/* 1. Redistributions of source code must retain the above copyright      */
/* notice, this list of conditions and the following disclaimer.          */
/*                                                                        */
/* 2. Redistributions in binary form must reproduce the above copyright   */
/* notice, this list of conditions and the following disclaimer in the    */
/* documentation and/or other materials provided with the distribution.   */
/*                                                                        */
/* 3. Neither the name of the Corporation nor the names of the            */
/* contributors may be used to endorse or promote products derived from   */
/* this software without specific prior written permission.               */
/*                                                                        */
/* THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY        */
/* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE      */
/* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR     */
/* PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE    */
/* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,    */
/* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR     */
/* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF */
/* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING   */
/* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS     */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.           */
/*------------------------------------------------------------------------*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ThreadPool_config.h>
#include <TPI.h>
#include <BoxPartitionIB.h>
#include <dcrs_matrix.h>
#include <CGSolver.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*--------------------------------------------------------------------*/
static
void hpccg_alloc_and_fill( const int np ,
                           const int my_p ,
                           const int gbox[][2] ,
                           const int ghost ,
                           struct distributed_crs_matrix * const matrix );

/*--------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
  const int ghost = 1 ;
  const int max_cube = 20 ;
  int ncube[20] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
                    0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

  FILE * print_file = stdout ;
  int print_iter = 500 ;
  int max_iter = 50 ;
  int overlap_comm = 0 ;

  float tolerance = 0.0 ; /* Force max iterations */

  int gbox[3][2] = { { 0 , 16 } , { 0 , 16 } , { 0 , 16 } };
  int nt = 0 ;
  int trials = 6 ;
  int ntest ;
  int np = 1;
  int my_p = 0 ;

#ifdef HAVE_MPI
  MPI_Init( & argc , & argv );
  MPI_Comm_size( MPI_COMM_WORLD , & np );
  MPI_Comm_rank( MPI_COMM_WORLD , & my_p );
#endif

  if ( ! my_p ) {
    const char arg_threads[] = "threads=" ;
    const char arg_cube[] = "cube=" ;
    const char arg_box[] = "box=" ;
    const char arg_max[] = "max_iter=" ;
    const char arg_trials[] = "trials=" ;
    const char arg_print[] = "print_iter=" ;
    const char arg_file[] = "print_file=" ;
    const char arg_comm[] = "overlap_comm=" ;
    const char arg_tolerance[] = "tolerance=" ;
    int i ;
    for ( i = 1 ; i < argc ; ++i ) {
      if ( ! strncmp(argv[i],arg_threads,strlen(arg_threads)) ) {
        sscanf(argv[i]+strlen(arg_threads),"%d",&nt);
      }
      else if ( ! strncmp(argv[i],arg_box,strlen(arg_box)) ) {
        sscanf(argv[i]+strlen(arg_box),"%d%*[x]%d%*[x]%d",
               & gbox[0][1] , & gbox[1][1] , & gbox[2][1] );
      }
      else if ( ! strncmp(argv[i],arg_cube,strlen(arg_cube)) ) {
        sscanf(argv[i]+strlen(arg_cube),
               "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
               ncube+0, ncube+1, ncube+2, ncube+3, ncube+4,
               ncube+5, ncube+6, ncube+7, ncube+8, ncube+9,
               ncube+10, ncube+11, ncube+12, ncube+13, ncube+14,
               ncube+15, ncube+16, ncube+17, ncube+18, ncube+19);
      }
      else if ( ! strncmp(argv[i],arg_max,strlen(arg_max)) ) {
        sscanf(argv[i]+strlen(arg_max),"%d",&max_iter);
      }
      else if ( ! strncmp(argv[i],arg_trials,strlen(arg_trials)) ) {
        sscanf(argv[i]+strlen(arg_trials),"%d",&trials);
      }
      else if ( ! strncmp(argv[i],arg_print,strlen(arg_print)) ) {
        sscanf(argv[i]+strlen(arg_print),"%d",&print_iter);
      }
      else if ( ! strncmp(argv[i],arg_comm,strlen(arg_comm)) ) {
        sscanf(argv[i]+strlen(arg_print),"%d",&overlap_comm);
      }
      else if ( ! strncmp(argv[i],arg_tolerance,strlen(arg_tolerance)) ) {
        sscanf(argv[i]+strlen(arg_print),"%f",&tolerance);
      }
      else if ( ! strncmp(argv[i],arg_file,strlen(arg_file)) ) {
        char buffer[256] ;
        sscanf(argv[i]+strlen(arg_file),"%s",buffer);
        print_file = fopen(buffer,"a");
      }
    }
  }

#ifdef HAVE_MPI
  {
    MPI_Bcast( & nt , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
    MPI_Bcast( & gbox[0][0] , 6 , MPI_INT , 0 , MPI_COMM_WORLD );
    MPI_Bcast( ncube , max_cube , MPI_INT , 0 , MPI_COMM_WORLD );
    MPI_Bcast( & overlap_comm , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
    MPI_Bcast( & max_iter , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
    MPI_Bcast( & print_iter , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
    MPI_Bcast( & trials , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
    MPI_Bcast( & tolerance , 1 , MPI_FLOAT , 0 , MPI_COMM_WORLD );
  }
#endif

  if ( nt ) {
    TPI_Init( nt );
    TPI_Block();
    TPI_Unblock();
  }

  if ( ! my_p ) {
    fprintf(print_file,"\"PROC\" , \"THREAD\" , \"EQUATION\" , \"NON-ZERO\" , \"FUSED-AVG\", \"FUSED-MAX\", \"BLAS-AVG\", \"BLAS-MAX\", \"FUSED\", \"BLAS\"  , \"Iter\"\n");
    fprintf(print_file,"\"COUNT\", \"COUNT\"  , \"COUNT\"    , \"COUNT\"    , \"Mflops\"   , \"Mflops\"   , \"Mflops\"  , \"Mflops\"  , \"error\", \"error\" , \"COUNT\"\n");
  }

  for ( ntest = 0 ; ! ntest || ( ntest < max_cube && ncube[ntest] ) ; ++ntest ) {
    struct distributed_crs_matrix matrix ;

    if ( ncube[ntest] ) {
      gbox[0][1] = gbox[1][1] = gbox[2][1] = ncube[ntest] ;
    }

    hpccg_alloc_and_fill( np, my_p, (const int (*)[2]) gbox, ghost, &matrix);

    {
      const int nRow = matrix.n_local_row ;

      double solve_dt[2] = { 0 , 0 };
      double solve_blas_dt[2] = { 0 , 0 };
      VECTOR_SCALAR norm_resid = 0.0 ;
      VECTOR_SCALAR norm_resid_blas = 0.0 ;
      int iter_count = 0 ;
      int iter_count_blas = 0 ;
      int k ;

      VECTOR_SCALAR * const b      = (VECTOR_SCALAR *) malloc( sizeof(VECTOR_SCALAR) * nRow );
      VECTOR_SCALAR * const x      = (VECTOR_SCALAR *) malloc( sizeof(VECTOR_SCALAR) * nRow );
      VECTOR_SCALAR * const x_blas = (VECTOR_SCALAR *) malloc( sizeof(VECTOR_SCALAR) * nRow );
      VECTOR_SCALAR * const xexact = (VECTOR_SCALAR *) malloc( sizeof(VECTOR_SCALAR) * nRow );

      {
        const VECTOR_SCALAR value = 1.0 /* 1.0 / 3.0 */ ;
        int i ;
        for ( i = 0 ; i < nRow ; ++i ) xexact[i] = value ;
      }

      for ( k = 0 ; k < trials ; ++k ) {
        double dt = 0 ;
        int i ;

        for ( i = 0 ; i < nRow ; ++i ) { x_blas[i] = 0.0 ; }

        cgsolve_set_lhs( & matrix , xexact , b );

        cgsolve_blas( & matrix, b, x_blas,
                      tolerance , max_iter , print_iter ,
                      & iter_count_blas, & norm_resid_blas, & dt );

        solve_blas_dt[0] += dt ;
        if ( ! k || dt < solve_blas_dt[1] ) { solve_blas_dt[1] = dt ; }
      }

      for ( k = 0 ; k < trials ; ++k ) {
        double dt = 0 ;
        int i ;

        for ( i = 0 ; i < nRow ; ++i ) { x[i] = 0.0 ; }

        cgsolve_set_lhs( & matrix , xexact , b );

        cgsolve( & matrix, b, x, overlap_comm,
                 tolerance , max_iter , print_iter ,
                 & iter_count, & norm_resid, & dt );

        solve_dt[0] += dt ;
        if ( ! k || dt < solve_dt[1] ) { solve_dt[1] = dt ; }
      }

      {
        int nnzGlobal = matrix.A_pc[ nRow ];
        double error[3] = { 0 , 0 , 0 };

        for ( k = 0 ; k < nRow ; ++k ) {
          error[0] += xexact[k] * xexact[k] ;
          error[1] += ( x[k] - xexact[k] ) * ( x[k] - xexact[k] );
          error[2] += ( x_blas[k] - xexact[k] ) * ( x_blas[k] - xexact[k] );
        }

#ifdef HAVE_MPI
        {
          double error_global[3] = { 0.0 , 0.0 , 0.0 };
          int nnz = nnzGlobal ;

          MPI_Allreduce( & nnz , & nnzGlobal , 1 , MPI_INT , MPI_SUM ,
                         MPI_COMM_WORLD );

          MPI_Allreduce( error , error_global , 3 , MPI_DOUBLE , MPI_SUM ,
                         MPI_COMM_WORLD );

          error[0] = error_global[0];
          error[1] = error_global[1];
          error[2] = error_global[2];
        }
#endif

        error[0] = sqrt( error[0] );
        error[1] = sqrt( error[1] );
        error[2] = sqrt( error[2] );

        if ( ! my_p ) {
          const int nRowGlobal = ( gbox[0][1] - gbox[0][0] ) *
                                 ( gbox[1][1] - gbox[1][0] ) *
                                 ( gbox[2][1] - gbox[2][0] );

          const double dt_mean_fuse_step = 1.0e6 * solve_dt[0]      / (double) trials ;
          const double dt_mean_blas_step = 1.0e6 * solve_blas_dt[0] / (double) trials ;
          const double dt_min_fuse_step  = 1.0e6 * solve_dt[1] ;
          const double dt_min_blas_step  = 1.0e6 * solve_blas_dt[1] ;

          const double Mflop_step = 2 * nnzGlobal 
                                  + 3 * 2 * nRowGlobal 
                                  + 2 * 2 * nRowGlobal ;

          const double Mflop_mean_fuse = Mflop_step * iter_count / dt_mean_fuse_step ;
          const double Mflop_mean_blas = Mflop_step * iter_count_blas / dt_mean_blas_step ;

          const double Mflop_max_fuse = Mflop_step * iter_count / dt_min_fuse_step ;
          const double Mflop_max_blas = Mflop_step * iter_count_blas / dt_min_blas_step ;

          fprintf(print_file,"%8d , %8d , %8d , %8d , %10g , %10g , %10g , %10g , %10g , %10g , %d\n",
                  np , nt , nRowGlobal , nnzGlobal ,
                  Mflop_mean_fuse , Mflop_max_fuse ,
                  Mflop_mean_blas , Mflop_max_blas ,
                  error[1] / error[0] , error[2] / error[0] , iter_count );
          fflush(print_file);
        }
      }

      free( xexact );
      free( x_blas );
      free( x );
      free( b );
    }
    free( matrix.A_a );
    free( matrix.A_ia );
    free( matrix.A_pc );
    free( matrix.p_recv_pc );
    free( matrix.p_send_pc );
    free( matrix.p_send_id );
  }

  if ( nt ) { TPI_Finalize(); }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0 ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static
void hpccg_alloc_and_fill( const int np ,
                           const int my_p ,
                           const int gbox[][2] ,
                           const int ghost ,
                           struct distributed_crs_matrix * const matrix )
{
  int (* const pbox)[3][2] = (int (*)[3][2]) malloc( sizeof(int)*np*3*2 );

  const int (* const my_box)[2] = (const int (*)[2]) pbox[my_p] ;

  int my_uses_box[3][2] ;
  int * map_local_ord = NULL;

  matrix->n_local_row     = 0 ;
  matrix->n_internal_row  = 0 ;
  matrix->A_pc            = NULL ;
  matrix->A_ia            = NULL ;
  matrix->A_a             = NULL ;

  matrix->p_size    = np ;
  matrix->p_rank    = my_p ;
  matrix->p_recv_pc = NULL  ;
  matrix->p_send_pc = NULL ;
  matrix->p_send_id = NULL ;

  /* Partition the global box */
  box_partition_rcb( np , gbox , pbox );

  /* Upper bound */
  map_local_ord = (int *) malloc( sizeof(int) *
                                  ( 2 * ghost + my_box[0][1]- my_box[0][0] ) *
                                  ( 2 * ghost + my_box[1][1]- my_box[1][0] ) *
                                  ( 2 * ghost + my_box[2][1]- my_box[2][0] ) );

  /* Generate local layout with ghosting. */
  box_partition_map( np, my_p, gbox,
                     (const int (* const)[3][2]) pbox,
                     ghost,
                     my_uses_box , map_local_ord ,
                     & matrix->n_internal_row ,
                     & matrix->n_local_row ,
                     & matrix->n_local_column ,
                     & matrix->p_recv_pc ,
                     & matrix->p_send_pc ,
                     & matrix->p_send_id );

  {
    const int nrow = matrix->n_local_row ;
    int * const pc = (int *) malloc( sizeof(int) * ( nrow + 1 ) );
    int * ia = NULL ;
    MATRIX_SCALAR * a = NULL ;

    int ix , iy , iz ;
    int sx , sy , sz ;

    /* Number of non zeros in each matrix row,
     * then prefix the array for offsets.
     */
    pc[0] = 0 ;

    for ( iz = my_box[2][0] ; iz < my_box[2][1] ; ++iz ) {
    for ( iy = my_box[1][0] ; iy < my_box[1][1] ; ++iy ) {
    for ( ix = my_box[0][0] ; ix < my_box[0][1] ; ++ix ) {
      const int irow = box_map_local( (const int (*const)[2]) my_uses_box, map_local_ord, ix, iy, iz );
      int count = 1 ; /* Count the diagonal */

      /* Count the off-diagonal terms to follow */
      for ( sz = -1 ; sz <= 1 ; ++sz ) {
      for ( sy = -1 ; sy <= 1 ; ++sy ) {
      for ( sx = -1 ; sx <= 1 ; ++sx ) {
        const int g_ix = ix + sx ;
        const int g_iy = iy + sy ;
        const int g_iz = iz + sz ;

        if ( my_uses_box[0][0] <= g_ix && g_ix < my_uses_box[0][1] &&
             my_uses_box[1][0] <= g_iy && g_iy < my_uses_box[1][1] &&
             my_uses_box[2][0] <= g_iz && g_iz < my_uses_box[2][1] &&
             ! ( sz == 0 && sy == 0 && sx == 0 ) ) {
          /* This column is within global bounds and is not a diagonal */
          ++count ;
        }
      }
      }
      }
      pc[ irow + 1 ] = count ;
    }
    }
    }

    for ( ix = 0 ; ix < nrow ; ++ix ) { pc[ix+1] += pc[ix] ; }

    ia = (int *)           malloc( sizeof(int)           * pc[ nrow ]  );
    a  = (MATRIX_SCALAR *) malloc( sizeof(MATRIX_SCALAR) * pc[ nrow ]  );

    for ( iz = my_box[2][0] ; iz < my_box[2][1] ; ++iz ) {
    for ( iy = my_box[1][0] ; iy < my_box[1][1] ; ++iy ) {
    for ( ix = my_box[0][0] ; ix < my_box[0][1] ; ++ix ) {
      const int irow = box_map_local( (const int (*const)[2]) my_uses_box, map_local_ord, ix, iy, iz );
      int ipc = pc[ irow ];

      /* Diagonal term first */
      ia[ ipc ] = irow ;
      a[  ipc ] = 27.0f ;
      ++ipc ;

      /* Off-diagonal terms to follow */
      for ( sz = -1 ; sz <= 1 ; ++sz ) {
      for ( sy = -1 ; sy <= 1 ; ++sy ) {
      for ( sx = -1 ; sx <= 1 ; ++sx ) {
        const int g_ix = ix + sx ;
        const int g_iy = iy + sy ;
        const int g_iz = iz + sz ;

        if ( my_uses_box[0][0] <= g_ix && g_ix < my_uses_box[0][1] &&
             my_uses_box[1][0] <= g_iy && g_iy < my_uses_box[1][1] &&
             my_uses_box[2][0] <= g_iz && g_iz < my_uses_box[2][1] &&
             ! ( sz == 0 && sy == 0 && sx == 0 ) ) {
          /* Column is within global bounds and is not a diagonal */
          /* 'icol' is mapped for communication */

          const int icol =
            box_map_local( (const int (*const)[2]) my_uses_box, map_local_ord, g_ix, g_iy, g_iz );

          if ( icol < 0 ) { abort(); }

          ia[ ipc ] = icol ;
          a[  ipc ] = -1.0f ;
          ++ipc ;
        }
      }
      }
      }
      if ( ipc != pc[ irow + 1 ] ) { abort(); }
    }
    }
    }

    matrix->A_pc = pc ;
    matrix->A_ia = ia ;
    matrix->A_a  = a ;
  }

  free( map_local_ord );
  free( pbox );
}

