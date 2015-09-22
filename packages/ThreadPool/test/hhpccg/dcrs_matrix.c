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


#include <stdlib.h>
#include <math.h>

#include <ThreadPool_config.h>
#include <TPI.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <dcrs_matrix.h>

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#if ! defined( HAVE_MPI )

static
double comm_sum( double v ) { return v ; }

#define get_off_process_entries( M , V )  /* */

/*--------------------------------------------------------------------*/
#else /* defined( HAVE_MPI ) */
/*--------------------------------------------------------------------*/

static
double comm_sum( double v )
{
  double result = 0 ;
  MPI_Allreduce( & v , & result , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
  return result ;
}

static
void get_off_process_entries(
  const struct distributed_crs_matrix * const matrix ,
  VECTOR_SCALAR * const vec )
{
  const int np   = matrix->p_size ;
  const int my_p = matrix->p_rank ;
  const int * const recv_pc = matrix->p_recv_pc ;
  const int * const send_pc = matrix->p_send_pc ;
  const int * const send_id = matrix->p_send_id ;
  int i , irecv ;

  for ( irecv = 0 , i = 1 ; i < np ; ++i ) {
    if ( recv_pc[i] < recv_pc[i+1] ) ++irecv ;
  }

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

#endif

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static void dcrs_apply_and_dot_span(
  const struct distributed_crs_matrix * const matrix ,
  const int span_begin ,
  const int span_end ,
  const VECTOR_SCALAR * const x ,
        VECTOR_SCALAR * const y ,
        double        * const result )
{
  const int           * const A_pc  = matrix->A_pc ;
  const int           * const A_ia  = matrix->A_ia ;
  const MATRIX_SCALAR * const A_a   = matrix->A_a ;

  double dot_x_y = *result ;

  int row = span_begin ;

  for ( ; row < span_end ; ++row ) {
    const int pcBeg = A_pc[ row ];
    const int pcEnd = A_pc[ row + 1 ];

    const int           *       ia    = A_ia + pcBeg ;
    const MATRIX_SCALAR *       a     = A_a  + pcBeg ;
    const MATRIX_SCALAR * const a_end = A_a  + pcEnd ;

    VECTOR_SCALAR y_tmp = 0 ;
    for ( ; a != a_end ; ++a , ++ia ) {
      y_tmp += *a * x[ *ia ];
    }
    dot_x_y += x[ row ] * y_tmp ;
    y[ row ] = y_tmp ;
  }

  *result = dot_x_y ;
}

static void dcrs_apply_span(
  const struct distributed_crs_matrix * const matrix ,
  const int span_begin ,
  const int span_end ,
  const VECTOR_SCALAR * const x ,
        VECTOR_SCALAR * const y )
{
  const int           * const A_pc  = matrix->A_pc ;
  const int           * const A_ia  = matrix->A_ia ;
  const MATRIX_SCALAR * const A_a   = matrix->A_a ;

  int row = span_begin ;

  for ( ; row < span_end ; ++row ) {
    const int pcBeg = A_pc[ row ];
    const int pcEnd = A_pc[ row + 1 ];

    const int           *       ia    = A_ia + pcBeg ;
    const MATRIX_SCALAR *       a     = A_a  + pcBeg ;
    const MATRIX_SCALAR * const a_end = A_a  + pcEnd ;

    VECTOR_SCALAR y_tmp = 0 ;
    for ( ; a != a_end ; ++a , ++ia ) {
      y_tmp += *a * x[ *ia ];
    }
    y[ row ] = y_tmp ;
  }
}

static void work_span( const int count , const int rank ,
                       int * jBeg , int * jEnd )
{
  const int length = *jEnd - *jBeg ;
  const int chunk  = ( length + count - 1 ) / count ;
  const int begin  = chunk * rank ;
        int end    = begin + chunk ;

  if ( length < end ) { end = length ; }

  *jEnd  = *jBeg + end ;
  *jBeg += begin ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static void tpi_work_dot_join( TPI_Work * work , const void * src  )
{ *((double *) ( work->reduce) ) += *((const double *) src); }

static void tpi_work_dot_init( TPI_Work * work )
{ *((double *) ( work->reduce) ) = 0 ; }

/*--------------------------------------------------------------------*/

struct work_dcrs {
  const struct distributed_crs_matrix * matrix ;
  const VECTOR_SCALAR * x ;
        VECTOR_SCALAR * y ;
  int   jBeg ;
  int   jEnd ;
};

/*--------------------------------------------------------------------*/

static void tpi_work_dcrs_apply_and_dot( TPI_Work * work )
{
  const struct work_dcrs * const info = (const struct work_dcrs *) work->info ;

  int local_begin = info->jBeg ;
  int local_end   = info->jEnd ;

  work_span( work->count , work->rank , & local_begin , & local_end );

  dcrs_apply_and_dot_span( info->matrix , local_begin , local_end ,
                           info->x , info->y , (double *) work->reduce );
}

double dcrs_apply_and_dot(
  const struct distributed_crs_matrix * matrix ,
  VECTOR_SCALAR * x ,
  VECTOR_SCALAR * y ,
  const int overlap_communication )
{
  struct work_dcrs info ;

  double result = 0.0 ;

  info.matrix = matrix ;
  info.x      = x ;
  info.y      = y ;

  if ( overlap_communication &&
       matrix->n_internal_row < matrix->n_local_row ) {

    double remote_result = 0 ;

    /* Start the internal matrix-vector multiply */
    /* result += dot( output = A * input , input ); */

    info.jBeg = 0 ;
    info.jEnd = matrix->n_internal_row ;

    /*  Divide internal work evenly among worker threads.
     *  This leave the primary thread completely out of the computation.
     */
    TPI_Start_threads_reduce( tpi_work_dcrs_apply_and_dot , & info , 
                              tpi_work_dot_join ,
                              tpi_work_dot_init ,
                              sizeof(result) , & result );

    get_off_process_entries( matrix , x );

    TPI_Wait(); /* Wait for internal result */

    info.jBeg = matrix->n_internal_row ;
    info.jEnd = matrix->n_local_row ;

    TPI_Run_threads_reduce( tpi_work_dcrs_apply_and_dot , & info , 
                            tpi_work_dot_join ,
                            tpi_work_dot_init ,
                            sizeof(remote_result) , & remote_result );

    result += remote_result ;
  }
  else {
    info.jBeg = 0 ;
    info.jEnd = matrix->n_local_row ;

    get_off_process_entries( matrix , x );

    TPI_Run_threads_reduce( tpi_work_dcrs_apply_and_dot , & info , 
                            tpi_work_dot_join ,
                            tpi_work_dot_init ,
                            sizeof(result) , & result );
  }

  result = comm_sum( result );

  return result ;
}

/*--------------------------------------------------------------------*/

static void tpi_work_dcrs_apply( TPI_Work * work )
{
  const struct work_dcrs * const info = (const struct work_dcrs *) work->info ;

  int local_begin = info->jBeg ;
  int local_end   = info->jEnd ;

  work_span( work->count , work->rank , & local_begin , & local_end );

  dcrs_apply_span( info->matrix , local_begin , local_end ,
                   info->x , info->y );
}

void dcrs_apply(
  const struct distributed_crs_matrix * matrix ,
  VECTOR_SCALAR * x ,
  VECTOR_SCALAR * y )
{
  struct work_dcrs info ;

  info.matrix = matrix ;
  info.x      = x ;
  info.y      = y ;
  info.jBeg   = 0 ;
  info.jEnd   = matrix->n_local_row ;

  get_off_process_entries( matrix , x );

  TPI_Run_threads( tpi_work_dcrs_apply , & info , 0 );
}

