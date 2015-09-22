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

void cgsolve_set_lhs( const struct distributed_crs_matrix * const matrix ,
                      const VECTOR_SCALAR * const x ,
                            VECTOR_SCALAR * const b )
{
  const int nRow = matrix->n_local_row ;
  const int nVec = matrix->p_recv_pc[ matrix->p_size ] ;

  VECTOR_SCALAR * const p =
    (VECTOR_SCALAR *) malloc( nVec * sizeof(VECTOR_SCALAR) );

  tpi_copy( nRow , x , p );

  dcrs_apply( matrix , p , b );

  free( p );
}

/*--------------------------------------------------------------------*/

/*  x += alpha * p ;
 *  r -= alpha * Ap ;
 *  return dot( r , r );
 */
static
double cgsolver_update( const int length ,
                        const VECTOR_SCALAR alpha ,
                        const VECTOR_SCALAR * p ,
                        const VECTOR_SCALAR * Ap ,
                              VECTOR_SCALAR * x ,
                              VECTOR_SCALAR * r );

/*--------------------------------------------------------------------*/

void cgsolve_blas( const struct distributed_crs_matrix * matrix ,
                   const VECTOR_SCALAR * const b ,
                         VECTOR_SCALAR * const x ,
                   const VECTOR_SCALAR tolerance ,
                   const int max_iter ,
                   const int print_iter ,
                         int    * const iter_count ,
                         VECTOR_SCALAR * const norm_resid ,
                         double * const solve_dt )
{
  const int nRow = matrix->n_local_row ;
  const int nVec = matrix->p_recv_pc[ matrix->p_size ] ;

  const VECTOR_SCALAR tol_2 = tolerance * tolerance ;

  VECTOR_SCALAR * const r  =
    (VECTOR_SCALAR *) malloc( nRow * sizeof(VECTOR_SCALAR) );
  VECTOR_SCALAR * const p  =
    (VECTOR_SCALAR *) malloc( nVec * sizeof(VECTOR_SCALAR) );
  VECTOR_SCALAR * const Ap =
    (VECTOR_SCALAR *) malloc( nRow * sizeof(VECTOR_SCALAR) );

  VECTOR_SCALAR rtrans = 0.0 ;
  VECTOR_SCALAR beta = 0.0 ;
  VECTOR_SCALAR pAp = 0.0 ;
  VECTOR_SCALAR alpha ;
  double time_begin , time_end ;

  int k ;

  tpi_copy( nRow , b , r );
  tpi_copy( nRow , x , p );

  /*  Ap = matrix * p ; */
  dcrs_apply( matrix , p , Ap );

  /*  r -= Ap ; */
  tpi_axpy( nRow , -1.0 , Ap , r );

  rtrans = tpi_dot( nRow , r , r );

  time_begin = TPI_Walltime();

  for ( k = 0 ; k < max_iter && tol_2 < rtrans ; ++k ) {

    /*  p = r + beta * p ; */
    tpi_xpby( nRow, r, beta, p ); /* parallel */

    dcrs_apply( matrix , p , Ap );

    pAp = tpi_dot( nRow , p , Ap );

    /* If orthogonal then cannot update */
    alpha = 0 < fabs( pAp ) ? rtrans / pAp : 0.0 ;

    /*  x += alpha * p ;
     *  r -= alpha * Ap ;
     *  return dot( r , r );
     */
    beta = rtrans ;

    tpi_axpy( nRow ,  alpha , p , x );
    tpi_axpy( nRow , -alpha , Ap , r );
    rtrans = tpi_dot( nRow , r , r );
    beta = rtrans / beta ;
  }

  time_end = TPI_Walltime();

#ifdef HAVE_MPI
  {
    double tb = time_begin ;
    double te = time_end ;
    MPI_Allreduce(&tb, &time_begin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&te, &time_end,   1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
#endif

  *solve_dt += time_end - time_begin ;

  *norm_resid = sqrt( rtrans );
  *iter_count = k ;

  free( Ap );
  free( p );
  free( r );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

void cgsolve( const struct distributed_crs_matrix * matrix ,
              const VECTOR_SCALAR * const b ,
                    VECTOR_SCALAR * const x ,
              const int overlap_comm ,
              const VECTOR_SCALAR tolerance ,
              const int max_iter ,
              const int print_iter ,
                    int    * const iter_count ,
                    VECTOR_SCALAR * const norm_resid ,
                    double * const solve_dt )
{
  const int nRow = matrix->n_local_row ;
  const int nVec = matrix->p_recv_pc[ matrix->p_size ] ;

  const VECTOR_SCALAR tol_2 = tolerance * tolerance ;

  VECTOR_SCALAR * const r  =
    (VECTOR_SCALAR *) malloc( nRow * sizeof(VECTOR_SCALAR) );
  VECTOR_SCALAR * const p  =
    (VECTOR_SCALAR *) malloc( nVec * sizeof(VECTOR_SCALAR) );
  VECTOR_SCALAR * const Ap =
    (VECTOR_SCALAR *) malloc( nRow * sizeof(VECTOR_SCALAR) );

  VECTOR_SCALAR rtrans = 0.0 ;
  VECTOR_SCALAR beta = 0.0 ;
  VECTOR_SCALAR pAp = 0.0 ;
  VECTOR_SCALAR alpha ;
  double time_begin , time_end ;

  int k ;

  tpi_copy( nRow , b , r );
  tpi_copy( nRow , x , p );

  /*  gather off-processor components of 'p'.
   *  Ap = matrix * p ;
   *  return dot( Ap , p );
   */
  pAp = dcrs_apply_and_dot( matrix , p , Ap , overlap_comm );

  /*  r -= 1.0 * Ap ;
   *  return dot( r , r );
   */
  alpha = 1.0 ;
  rtrans = cgsolver_update( nRow, alpha, NULL, Ap, NULL, r ); /* parallel */

  time_begin = TPI_Walltime();

  for ( k = 0 ; k < max_iter && tol_2 < rtrans ; ++k ) {

    /*  p = r + beta * p ; */
    tpi_xpby( nRow, r, beta, p ); /* parallel */

    /*  gather off-processor components of 'p'.
     *  Ap = matrix * p ;
     *  return dot( Ap , p );
     */
    pAp = dcrs_apply_and_dot( matrix , p , Ap , overlap_comm ); /* parallel */

    /* If orthogonal then cannot update */
    alpha = 0 < fabs( pAp ) ? rtrans / pAp : 0.0 ;

    /*  x += alpha * p ;
     *  r -= alpha * Ap ;
     *  return dot( r , r );
     */
    beta = rtrans ;
    rtrans = cgsolver_update( nRow , alpha , p , Ap , x , r ); /* parallel */
    beta = rtrans / beta ;
  }

  time_end = TPI_Walltime();

#ifdef HAVE_MPI
  {
    double tb = time_begin ;
    double te = time_end ;
    MPI_Allreduce(&tb, &time_begin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&te, &time_end,   1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  }
#endif

  *solve_dt += time_end - time_begin ;

  *norm_resid = sqrt( rtrans );
  *iter_count = k ;

  free( Ap );
  free( p );
  free( r );
}

/*--------------------------------------------------------------------*/

struct tpi_work_cgsolve {
  const VECTOR_SCALAR * p ;
  const VECTOR_SCALAR * Ap ;
        VECTOR_SCALAR * x ;
        VECTOR_SCALAR * r ;
        VECTOR_SCALAR alpha ;
  int length ;
};

static void tpi_work_dot_join( TPI_Work * work , const void * src  )
{ *((double *) work->reduce ) += *((const double *) src); }
 
static void tpi_work_dot_init( TPI_Work * work )
{ *((double *) work->reduce ) = 0 ; }

static void tpi_work_update( TPI_Work * work )
{
  const struct tpi_work_cgsolve * const cg_work = 
    (const struct tpi_work_cgsolve *) work->info ;

  const int           length = cg_work->length ;
  const VECTOR_SCALAR alpha  = cg_work->alpha ;
  const VECTOR_SCALAR * const p  = cg_work->p ;
  const VECTOR_SCALAR * const Ap = cg_work->Ap ;
        VECTOR_SCALAR * const x  = cg_work->x ;
        VECTOR_SCALAR * const r  = cg_work->r ;

  double mag = 0 ;
  int iBeg , iEnd , i ;

  tpi_work_span( work , length , & iBeg , & iEnd );

  if ( x ) { for ( i = iBeg ; i < iEnd ; ++i ) { x[i] += alpha * p[i]; } }

  for ( i = iBeg ; i < iEnd ; ++i ) {
    const VECTOR_SCALAR val = ( r[i] -= alpha * Ap[i] );
    mag += val * val ;
  }

  *((double*) work->reduce ) = mag ;
}

double cgsolver_update( const int length ,
                        const VECTOR_SCALAR alpha ,
                        const VECTOR_SCALAR * p ,
                        const VECTOR_SCALAR * Ap ,
                              VECTOR_SCALAR * x ,
                              VECTOR_SCALAR * r )
{
  struct tpi_work_cgsolve work ;

  double result = 0.0 ;

  work.length = length ;
  work.alpha  = alpha ;
  work.p  = p ;
  work.Ap = Ap ;
  work.x  = x ;
  work.r  = r ;

  TPI_Run_threads_reduce( tpi_work_update , & work ,
                          tpi_work_dot_join , tpi_work_dot_init ,
                          sizeof(result) , & result );

#ifdef HAVE_MPI
  {
    double local = result ;
    MPI_Allreduce( & local, & result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  }
#endif

  return result ;
}

/*--------------------------------------------------------------------*/

