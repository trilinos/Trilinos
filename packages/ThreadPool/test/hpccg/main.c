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
#include <BoxPartition.h>
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
                           struct cgsolve_data * const data )
{
  int (*pbox)[3][2] = NULL ;
  int * map_local_ord = NULL;

  data->nRow = 0 ;
  data->A_pc = NULL ;
  data->A_ia = NULL ;
  data->A_a  = NULL ;

  data->np = np ;
  data->ip = my_p ;
  data->recv_pc = NULL  ;
  data->send_pc = NULL ;
  data->send_id = NULL ;

  box_partition_rcb( np, my_p,
                     (const int (*)[2]) gbox, ghost,
                     & pbox ,
                     & map_local_ord ,
                     & data->recv_pc ,
                     & data->send_pc ,
                     & data->send_id );

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
    MATRIX_SCALAR  * const a  = (MATRIX_SCALAR *) malloc( sizeof(MATRIX_SCALAR) * nnz );

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

    data->nRow = irow ;
    data->A_pc = pc ;
    data->A_ia = ia ;
    data->A_a  = a ;
  }

  free( map_local_ord );
  free( pbox );
}

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
  int trials = 5 ;
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
    fprintf(print_file,"\"PROC\" , \"THREAD\" , \"EQUATION\" , \"NON-ZERO\" , \"MXV\"    , \"AXPBY\"  , \"DOT\" , \"Xerror\" , \"Iter\"\n");
    fprintf(print_file,"\"COUNT\" , \"COUNT\"  , \"COUNT\"    , \"COUNT\"    , \"Mflops\" , \"Mflops\" , \"Mflops\" , \"L2norm\" , \"COUNT\"\n");
  }

  for ( ntest = 0 ; ! ntest || ( ntest < max_cube && ncube[ntest] ) ; ++ntest ) {
    struct cgsolve_data cgdata ;

    if ( ncube[ntest] ) {
      gbox[0][1] = gbox[1][1] = gbox[2][1] = ncube[ntest] ;
    }

    hpccg_alloc_and_fill( np, my_p, (const int (*)[2]) gbox, ghost, &cgdata);

    cgdata.max_iter   = max_iter ;
    cgdata.print_iter = print_iter ;
    cgdata.tolerance  = tolerance ;

    {
      double dt_mxv[2] = { 0 , 0 };
      double dt_axpby[2] = { 0 , 0 };
      double dt_dot[2] = { 0 , 0 };
      VECTOR_SCALAR norm_resid = 0.0 ;
      int iter_count = 0 ;
      int iter_total = 0 ;
      int k ;

      VECTOR_SCALAR * const b      = (VECTOR_SCALAR *) malloc( sizeof(VECTOR_SCALAR) * cgdata.nRow );
      VECTOR_SCALAR * const x      = (VECTOR_SCALAR *) malloc( sizeof(VECTOR_SCALAR) * cgdata.nRow );
      VECTOR_SCALAR * const xexact = (VECTOR_SCALAR *) malloc( sizeof(VECTOR_SCALAR) * cgdata.nRow );

      {
        const VECTOR_SCALAR value = 1.0 /* 1.0 / 3.0 */ ;
        int i ;
        for ( i = 0 ; i < cgdata.nRow ; ++i ) xexact[i] = value ;
      }

      for ( k = 0 ; k < trials ; ++k ) {
        int i ;

        for ( i = 0 ; i < cgdata.nRow ; ++i ) { x[i] = 0.0 ; }

        cgsolve_set_lhs( & cgdata , xexact , b );

        cgsolve( & cgdata, b, x,
                 & iter_count, & norm_resid,
                 dt_mxv , dt_axpby , dt_dot );

        iter_total += iter_count ;
      }

      {
        int nnzGlobal = cgdata.A_pc[ cgdata.nRow ];
        double error[2] = { 0 , 0 };

        for ( k = 0 ; k < cgdata.nRow ; ++k ) {
          error[0] += ( x[k] - xexact[k] ) * ( x[k] - xexact[k] );
          error[1] += xexact[k] * xexact[k] ;
        }

#ifdef HAVE_MPI
        {
          double error_global[2] = { 0.0 , 0.0 };
          int nnz = nnzGlobal ;

          MPI_Allreduce( & nnz , & nnzGlobal , 1 , MPI_INT , MPI_SUM ,
                         MPI_COMM_WORLD );

          MPI_Allreduce( error , error_global , 2 , MPI_DOUBLE , MPI_SUM ,
                         MPI_COMM_WORLD );

          error[0] = error_global[0];
          error[1] = error_global[1];
        }
#endif

        error[0] = sqrt( error[0] );
        error[1] = sqrt( error[1] );

        if ( ! my_p ) {
          const int nRowGlobal = ( gbox[0][1] - gbox[0][0] ) *
                                 ( gbox[1][1] - gbox[1][0] ) *
                                 ( gbox[2][1] - gbox[2][0] );

          const double mflop_mxv =
             1.0e-6 * ( iter_total ) * 2 * nnzGlobal / dt_mxv[0] ;

          const double mflop_axpby =
             1.0e-6 * ( iter_total * 3 ) * 3 * nRowGlobal / dt_axpby[0] ;

          const double mflop_dot =
             1.0e-6 * ( iter_total * 2 ) * 2 * nRowGlobal / dt_dot[0] ;

          fprintf(print_file,"%8d , %8d , %8d , %8d , %10g , %10g , %10g , %g , %d\n",
                  np , nt , nRowGlobal , nnzGlobal ,
                  mflop_mxv , mflop_axpby , mflop_dot ,
                  error[0] / error[1] , iter_total );
          fflush(print_file);
        }
      }

      free( xexact );
      free( x );
      free( b );
    }
    free( cgdata.A_a );
    free( cgdata.A_ia );
    free( cgdata.A_pc );
    free( cgdata.recv_pc );
    free( cgdata.send_pc );
    free( cgdata.send_id );
  }

  if ( nt ) { TPI_Finalize(); }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return 0 ;
}

