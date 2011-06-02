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
#include <tpi_vector.h>
#include <CGSolver.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*--------------------------------------------------------------------*/

#ifdef HAVE_MPI

#define TIMER( DT , F )	\
  { double tb , te , tbg , teg , dt ; \
    tb = TPI_Walltime(); \
    F ; \
    te = TPI_Walltime(); \
    MPI_Allreduce(&tb, &tbg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); \
    MPI_Allreduce(&te, &teg, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); \
    DT[0] += dt = teg - tbg ; \
    DT[1] += dt * dt ; }

#else

#define TIMER( DT , F )	\
  { const double tb = TPI_Walltime(); double dt ; \
    F ; \
    DT[0] += dt = TPI_Walltime() - tb ; \
    DT[1] += dt * dt ; }

#endif

/*--------------------------------------------------------------------*/

static
VECTOR_SCALAR comm_sum( VECTOR_SCALAR v )
{
#ifdef HAVE_MPI
  VECTOR_SCALAR result = 0 ;
  if ( sizeof(VECTOR_SCALAR) == sizeof(double) ) {
    MPI_Allreduce( & v , & result , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
  }
  else {
    MPI_Allreduce( & v , & result , 1 , MPI_FLOAT , MPI_SUM , MPI_COMM_WORLD );
  }
  return result ;
#else
  return v ;
#endif
}

#ifdef HAVE_MPI
static
void comm_rhs_vector( const struct cgsolve_data * const data ,
                      VECTOR_SCALAR * const vec )
{
  const int np = data->np ;
  const int my_p = data->ip ;
  const int * const recv_pc = data->recv_pc ;
  const int * const send_pc = data->send_pc ;
  const int * const send_id = data->send_id ;
  int i , irecv ;

  for ( irecv = 0 , i = 1 ; i < np ; ++i ) {
    if ( recv_pc[i] < recv_pc[i+1] ) ++irecv ;
  }

#ifdef DEBUG_PRINT
  fflush(stdout);
  MPI_Barrier( MPI_COMM_WORLD );
  fflush(stdout);
#endif

  {
    VECTOR_SCALAR * const send_buf =
      (VECTOR_SCALAR *) malloc( sizeof(VECTOR_SCALAR) * send_pc[np] );

    MPI_Request * const recv_request =
      (MPI_Request *) malloc( sizeof(MPI_Request) * irecv );

    MPI_Status * const recv_status =
      (MPI_Status *) malloc( sizeof(MPI_Status) * irecv );

    for ( irecv = 0 , i = 1 ; i < np ; ++i ) {
      const int ip = ( i + my_p ) % np ;
      const int recv_beg    = recv_pc[i];
      const int recv_length = recv_pc[i+1] - recv_beg ;
      if ( recv_length ) {
#ifdef DEBUG_PRINT
        fprintf(stdout,"  comm_rhs_vector P%d Irecv P%d : %d\n",
                       my_p, ip, recv_length );
        fflush(stdout);
#endif
        MPI_Irecv( vec + recv_beg ,
                   recv_length * sizeof(VECTOR_SCALAR), MPI_BYTE ,
                   ip , 0 , MPI_COMM_WORLD , recv_request + irecv );
        ++irecv ;
      }
    }

    /* Gather components into send buffer */

    for ( i = 0 ; i < send_pc[np] ; ++i ) {
      send_buf[i] = vec[ send_id[i] ];
    }

    MPI_Barrier( MPI_COMM_WORLD );

    for ( i = 1 ; i < np ; ++i ) {
      const int ip = ( i + my_p ) % np ;
      const int send_beg    = send_pc[i];
      const int send_length = send_pc[i+1] - send_beg ;
      if ( send_length ) { /* Send to 'i' */
#ifdef DEBUG_PRINT
        fprintf(stdout,"  comm_rhs_vector P%d Rsend P%d : %d\n",
                       my_p, ip, send_length );
        fflush(stdout);
#endif
        MPI_Rsend( send_buf + send_beg ,
                   send_length * sizeof(VECTOR_SCALAR), MPI_BYTE ,
                   ip , 0 , MPI_COMM_WORLD );
      }
    }

    MPI_Waitall( irecv , recv_request , recv_status );

    free( recv_status );
    free( recv_request );
    free( send_buf );
  }
}
#else
#define comm_rhs_vector( D , V ) /* */
#endif

/*--------------------------------------------------------------------*/

void cgsolve_set_lhs( const struct cgsolve_data * const data ,
                      const VECTOR_SCALAR * const x ,
                            VECTOR_SCALAR * const b )
{
  const int nRow = data->nRow ;
  const int nVec = data->recv_pc[ data->np ] ;
  const int   * const A_pc = data->A_pc ;
  const int   * const A_ia = data->A_ia ;
  const MATRIX_SCALAR * const A_a  = data->A_a ;

  VECTOR_SCALAR * const p = (VECTOR_SCALAR *) malloc( nVec * sizeof(VECTOR_SCALAR) );

  tpi_copy( nRow , x , p );

  comm_rhs_vector( data , p );

  tpi_crs_matrix_apply( nRow, A_pc, A_ia, A_a, p, b );

  free( p );
}

/*--------------------------------------------------------------------*/

void cgsolve( const struct cgsolve_data * const data ,
              const VECTOR_SCALAR * const b ,
                    VECTOR_SCALAR * const x ,
                    int    * const iter_count ,
                    VECTOR_SCALAR * const norm_resid ,
                    double * const dt_mxv ,  
                    double * const dt_axpby ,
                    double * const dt_dot )
{
  const int nRow = data->nRow ;
  const int nVec = data->recv_pc[ data->np ] ;
  const int max_iter = data->max_iter ;
  const int print_iter = data->print_iter ;
  const int   * const A_pc = data->A_pc ;
  const int   * const A_ia = data->A_ia ;
  const MATRIX_SCALAR * const A_a  = data->A_a ;
  const VECTOR_SCALAR tolerance = data->tolerance ;

  const VECTOR_SCALAR tol_2 = tolerance * tolerance ;

  VECTOR_SCALAR * const r  = (VECTOR_SCALAR *) malloc( nRow * sizeof(VECTOR_SCALAR) );
  VECTOR_SCALAR * const p  = (VECTOR_SCALAR *) malloc( nVec * sizeof(VECTOR_SCALAR) );
  VECTOR_SCALAR * const Ap = (VECTOR_SCALAR *) malloc( nRow * sizeof(VECTOR_SCALAR) );

  VECTOR_SCALAR rtrans = 0.0 ;

  int k ;

  tpi_copy( nRow , b , r );
  tpi_copy( nRow , x , p );

  comm_rhs_vector( data , p ); tpi_crs_matrix_apply( nRow, A_pc, A_ia, A_a, p, Ap );

  tpi_axpby( nRow , -1.0, Ap, 1.0 , r );

  /* Include timing dot product for 2 * #iter dot products */
  TIMER( dt_dot , rtrans = comm_sum( tpi_dot( nRow , r , r ) ) );

  for ( k = 0 ; k < max_iter && tol_2 < rtrans ; ++k ) {
    VECTOR_SCALAR alpha ;
    VECTOR_SCALAR beta = 0.0 ;
    VECTOR_SCALAR pAp = 0.0 ;

    if ( k ) {
      const VECTOR_SCALAR oldrtrans = rtrans ;
      TIMER( dt_dot , rtrans = comm_sum( tpi_dot( nRow , r , r ) ) );
      beta = rtrans / oldrtrans ;
    }

    TIMER( dt_axpby , tpi_axpby( nRow, 1.0, r, beta, p ) );

    TIMER( dt_mxv , comm_rhs_vector( data , p ); tpi_crs_matrix_apply( nRow, A_pc, A_ia, A_a, p, Ap ) );

    TIMER( dt_dot , pAp = comm_sum( tpi_dot( nRow , p , Ap ) ) );

    if ( 0 < fabs( pAp ) ) {
      alpha = rtrans / pAp ;
    }
    else {
      alpha = rtrans = 0.0 ; /* Orthogonal, cannot continue */
    }

    if ( ! ( ( k + 1 ) % print_iter ) ) {
      fprintf(stdout,"  cgsolve | r(%d) | = %g\n",k,sqrt(rtrans));
      fflush(stdout);
    }
  
    TIMER( dt_axpby , tpi_axpby( nRow , alpha,  p,  1.0, x) );
    TIMER( dt_axpby , tpi_axpby( nRow , -alpha, Ap, 1.0, r) );
  }

  *norm_resid = sqrt( rtrans );
  *iter_count = k ;

  free( Ap );
  free( p );
  free( r );
}

/*--------------------------------------------------------------------*/

