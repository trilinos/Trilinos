
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ThreadPool_config.h>
#include <TPI.h>
#include <BoxPartition.h>
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

  VECTOR_SCALAR tolerance = 0.0 ; /* Force max iterations */

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
    MPI_Bcast( & max_iter , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
    MPI_Bcast( & print_iter , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
    MPI_Bcast( & trials , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
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

        cgsolve( & matrix, b, x,
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
  int (*pbox)[3][2] = NULL ;
  int * map_local_ord = NULL;

  matrix->n_local_row     = 0 ;
  matrix->n_internal_row  = 0 ;
  matrix->A_pc            = NULL ;
  matrix->A_ia            = NULL ;
  matrix->A_a             = NULL ;
  matrix->A_row_partition = NULL ;

  matrix->p_size    = np ;
  matrix->p_rank    = my_p ;
  matrix->p_recv_pc = NULL  ;
  matrix->p_send_pc = NULL ;
  matrix->p_send_id = NULL ;

  box_partition_rcb( np, my_p,
                     (const int (*)[2]) gbox, ghost,
                     & pbox ,
                     & map_local_ord ,
                     & matrix->p_recv_pc ,
                     & matrix->p_send_pc ,
                     & matrix->p_send_id );

  {
    const int (* const my_box)[2] = (const int (*)[2]) pbox[my_p] ;
    const int bx = my_box[0][0] ;
    const int by = my_box[1][0] ;
    const int bz = my_box[2][0] ;
    const int nx = my_box[0][1] - bx ;
    const int ny = my_box[1][1] - by ;
    const int nz = my_box[2][1] - bz ;
    const int n = nx * ny * nz ;
    const int nnz = 27 * n ; /* Upper bound */
    int    * const pc = (int *)   malloc( sizeof(int) * ( n + 1 ) );
    int    * const ia = (int *)   malloc( sizeof(int) * nnz );
    MATRIX_SCALAR  * const a =
      (MATRIX_SCALAR *) malloc( sizeof(MATRIX_SCALAR) * nnz );

    int irow = 0 ;
    int ipc  = 0 ;
    int ix , iy , iz ;
    int sx , sy , sz ;

    for ( iz = 0 ; iz < nz ; ++iz ) {
    for ( iy = 0 ; iy < ny ; ++iy ) {
    for ( ix = 0 ; ix < nx ; ++ix , ++irow ) {

      if ( irow != box_map_local( my_box, ghost, map_local_ord,ix,iy,iz) ) {
        fprintf(stderr,"P%d:  irow[%d] != box_map_local(%d,%d,%d) = %d\n",
                my_p,irow,ix,iy,iz,
                box_map_local( my_box, ghost, map_local_ord, ix, iy, iz) );
      }

      pc[ irow ] = ipc ;   /* Beginning of row coefficients */
      /* Diagonal term first */
      ia[ ipc ] = irow ;
      a[  ipc ] = 27.0f ;
      ++ipc ;

      /* Off-diagonal terms to follow */
      for ( sz = -1 ; sz <= 1 ; ++sz ) {
      for ( sy = -1 ; sy <= 1 ; ++sy ) {
      for ( sx = -1 ; sx <= 1 ; ++sx ) {
        const int dx = ix + sx ;
        const int dy = iy + sy ;
        const int dz = iz + sz ;
        const int global_x = dx + bx ;
        const int global_y = dy + by ;
        const int global_z = dz + bz ;

        if ( gbox[0][0] <= global_x && global_x < gbox[0][1] &&
             gbox[1][0] <= global_y && global_y < gbox[1][1] &&
             gbox[2][0] <= global_z && global_z < gbox[2][1] &&
             ! ( sz == 0 && sy == 0 && sx == 0 ) ) {
          /* Column is within global bounds and is not a diagonal */
          /* 'icol' is mapped for communication */

          const int icol =
            box_map_local(my_box,ghost,map_local_ord,dx,dy,dz);

          if ( icol < 0 ) {
            fprintf(stderr,"P%d : bad column at local (%d,%d,%d) global(%d,%d,%d)\n",
                    my_p, dx,dy,dz,global_x,global_y,global_z);
            fflush(stderr);
            abort();
          }

          ia[ ipc ] = icol ;
          a[  ipc ] = -1.0f ;
          ++ipc ;
        }
      }
      }
      }
    }
    }
    }

    pc[irow] = ipc ;

    matrix->n_local_row = irow ;
    matrix->A_pc = pc ;
    matrix->A_ia = ia ;
    matrix->A_a  = a ;
  }

  free( map_local_ord );
  free( pbox );

  /* Now partition rows into purely-internal span and
   * rows with with-off-process columns span.
   */
  {
    const int nrow = matrix->n_local_row ;
    const int * const pc = matrix->A_pc ;
    const int * const ia = matrix->A_ia ;
    int * const row_part = (int *) malloc( sizeof(int) * nrow );
    int irow = 0 , ipart = 0 ;

    for ( irow = 0 ; irow < matrix->n_local_row ; ++irow ) {
            int icol     = pc[ irow ];
      const int icol_end = pc[ irow + 1 ];
      for ( ; icol < icol_end && ia[ icol ] < nrow ; ++icol );
      if ( icol == icol_end ) {
        /* Entirely local */
        row_part[ ipart ] = irow ;
        ++ipart ;
      }
    }
    matrix->n_internal_row = ipart ;

    for ( irow = 0 ; irow < matrix->n_local_row ; ++irow ) {
            int icol     = pc[ irow ];
      const int icol_end = pc[ irow + 1 ];
      for ( ; icol < icol_end && ia[ icol ] < nrow ; ++icol );
      if ( icol < icol_end ) {
        /* Has off-processor */
        row_part[ ipart ] = irow ;
        ++ipart ;
      }
    }
    matrix->A_row_partition = row_part ;
  }
}

