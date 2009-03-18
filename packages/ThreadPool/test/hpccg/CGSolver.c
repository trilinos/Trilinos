
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <TPI.h>
#include <tpi_vector.h>
#include <CGSolver.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*--------------------------------------------------------------------*/

#ifdef HAVE_MPI

#define TIMER( DT , F )	\
  { double tb , te , tbg , teg ; \
    tb = TPI_Walltime(); \
    F ; \
    te = TPI_Walltime(); \
    MPI_Allreduce(&tb, &tbg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD); \
    MPI_Allreduce(&te, &teg, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); \
    DT += teg - tbg ; }

#else

#define TIMER( DT , F )	\
  { const double t = TPI_Walltime(); F ; DT += TPI_Walltime() - t ; }

#endif

/*--------------------------------------------------------------------*/

static
double comm_sum( double v )
{
#ifdef HAVE_MPI
  double result = 0 ;
  MPI_Allreduce( & v , & result , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );
  return result ;
#else
  return v ;
#endif
}

#ifdef HAVE_MPI
static
void comm_rhs_vector( const struct cgsolve_data * const data ,
                      double * const vec )
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
    double * const send_buf =
      (double *) malloc( sizeof(double) * send_pc[np] );

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
        MPI_Irecv( vec + recv_beg , recv_length , MPI_DOUBLE ,
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
        MPI_Rsend( send_buf + send_beg , send_length , MPI_DOUBLE ,
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
                      const double * const x ,
                            double * const b )
{
  const int nRow = data->nRow ;
  const int nVec = data->recv_pc[ data->np ] ;
  const int   * const A_pc = data->A_pc ;
  const int   * const A_ia = data->A_ia ;
  const float * const A_a  = data->A_a ;

  double * const p = (double *) malloc( nVec * sizeof(double) );

  tpi_copy( nRow , x , p );

  comm_rhs_vector( data , p );

  tpi_crs_matrix_apply( nRow, A_pc, A_ia, A_a, p, b );
}

/*--------------------------------------------------------------------*/

void cgsolve( const struct cgsolve_data * const data ,
              const double * const b ,
                    double * const x ,
                    int    * const iter_count ,
                    double * const norm_resid ,
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
  const float * const A_a  = data->A_a ;
  const double tolerance = data->tolerance ;

  const double tol_2 = tolerance * tolerance ;

  double * const r  = (double *) malloc( nRow * sizeof(double) );
  double * const p  = (double *) malloc( nVec * sizeof(double) );
  double * const Ap = (double *) malloc( nRow * sizeof(double) );

  double rtrans = 0.0 ;

  int k ;

  tpi_copy( nRow , b , r );
  tpi_copy( nRow , x , p );

  TIMER( *dt_mxv , comm_rhs_vector( data , p ); tpi_crs_matrix_apply( nRow, A_pc, A_ia, A_a, p, Ap ) );

  TIMER( *dt_axpby , tpi_axpby( nRow , -1.0, Ap, 1.0 , r ) );

  TIMER( *dt_dot , rtrans = comm_sum( tpi_dot( nRow , r , r ) ) );

  for ( k = 0 ; k < max_iter && tol_2 < rtrans ; ++k ) {
    double alpha ;
    double beta = 0.0 ;
    double pAp = 0.0 ;

    if ( k ) {
      const double oldrtrans = rtrans ;
      TIMER( *dt_dot , rtrans = comm_sum( tpi_dot( nRow , r , r ) ) );
      beta = 0 < rtrans ? rtrans / oldrtrans : 0.0 ;
    }

    TIMER( *dt_axpby , tpi_axpby( nRow, 1.0, r, beta, p ) );

    TIMER( *dt_mxv , comm_rhs_vector( data , p ); tpi_crs_matrix_apply( nRow, A_pc, A_ia, A_a, p, Ap ) );

    TIMER( *dt_dot , pAp = comm_sum( tpi_dot( nRow , p , Ap ) ) );

    alpha = 0 < rtrans ? rtrans / pAp : 0.0 ;

    if ( ! ( ( k + 1 ) % print_iter ) ) {
      fprintf(stdout,"  cgsolve | r(%d) | = %g\n",k,sqrt(rtrans));
      fflush(stdout);
    }
  
    TIMER( *dt_axpby , tpi_axpby( nRow , alpha,  p,  1.0, x) );
    TIMER( *dt_axpby , tpi_axpby( nRow , -alpha, Ap, 1.0, r) );
  }

  *norm_resid = sqrt( rtrans );
  *iter_count = k ;

  free( Ap );
  free( p );
  free( r );
}

/*--------------------------------------------------------------------*/

