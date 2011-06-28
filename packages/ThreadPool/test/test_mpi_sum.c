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

/**
 * @author H. Carter Edwards
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <TPI.h>
#include <ThreadPool_config.h>

int rand_r( unsigned int * );

/*--------------------------------------------------------------------*/

#if defined(HAVE_MPI)

#include <mpi.h>

typedef MPI_Comm COMM ;

#else

typedef int COMM ;

#endif

static int comm_size( COMM );
static int comm_rank( COMM );
static void comm_reduce_dmax( COMM , double * );
static void comm_reduce_dsum( COMM , double * );
static void comm_reduce_d4_sum( COMM , double * );

/*--------------------------------------------------------------------*/

static void my_span( const unsigned count , const unsigned rank ,
                     const unsigned size ,
                     unsigned * begin , unsigned * length )
{
  const unsigned int max = ( size + count - 1 ) / count ;
  const unsigned int end = size - max * ( count - ( rank + 1 ) );
  if ( rank ) {
    *begin  = end - max ;
    *length = max ;
  }
  else {
    *begin  = 0 ;
    *length = end ;
  }
}

/*--------------------------------------------------------------------*/

#define LESS_ABS( X , Y )	( ( X < 0 ? -X : X ) < ( Y < 0 ? -Y : Y ) )

static void d2_add_d( double v[] , const double a )
{
  const int AltV = a < 0 ? ( - a < ( v[0] < 0 ? - v[0] : v[0] ) )
                         : (   a < ( v[0] < 0 ? - v[0] : v[0] ) );

  const double VpA = v[0] + a ;

  v[1] += AltV ? ( a - ( VpA - v[0] ) ) : ( v[0] - ( VpA - a ) );
  v[0]  = VpA + v[1] ;
  v[1] += VpA - v[0] ;
}

void d4_dot( double v[] , unsigned n , const double * x , const double * y )
{
  double * pos = v ;
  double * neg = v + 2 ;
  const double * const x_end = x + n ;
  for ( ; x < x_end ; ++x , ++y ) {
    const double a = *x * *y ;
    if ( a < 0 ) { d2_add_d( neg , a ); }
    else         { d2_add_d( pos , a ); }
  }
}

double ddot( unsigned n , const double * x , const double * y )
{
  double val = 0 ;
  const double * const x_end = x + n ;
  for ( ; x < x_end ; ++x , ++y ) { val += *x * *y ; }
  return val ;
}

/*--------------------------------------------------------------------*/

struct TaskXY {
  unsigned int   nreduce ;
  unsigned int   n ;
  const double * x ;
  const double * y ;
};

static
void reduce_init( TPI_Work * work )
{
  struct TaskXY * const info = (struct TaskXY *) work->info ;
  double        * const dst  = (double *)        work->reduce ;

  if ( info->nreduce == 4 ) {
    dst[0] = 0 ;
    dst[1] = 0 ;
    dst[2] = 0 ;
    dst[3] = 0 ;
  }
  else if ( info->nreduce == 1 ) {
    dst[0] = 0 ;
  }
}

static
void reduce_join( TPI_Work * work , const void * arg_src )
{
  struct TaskXY * const info = (struct TaskXY *) work->info ;
  double        * const dst  = (double *)        work->reduce ;
  const double  * const src  = (const double *)  arg_src ;

  if ( info->nreduce == 4 ) {
    d2_add_d( dst ,     src[0] );
    d2_add_d( dst ,     src[1] );
    d2_add_d( dst + 2 , src[2] );
    d2_add_d( dst + 2 , src[3] );
  }
  else if ( info->nreduce == 1 ) {
    dst[0] += src[0] ;
  }
}

/*--------------------------------------------------------------------*/

static
void work_d4_dot_tp( TPI_Work * work )
{
  struct TaskXY * const info = (struct TaskXY *) work->info ;
  double        * const dst  = (double *)        work->reduce ;

  unsigned int begin , length ;

  my_span( work->count , work->rank , info->n , & begin , & length );

  d4_dot( dst , length , info->x + begin , info->y + begin );
}

double d4_dot_tp( COMM comm, unsigned nwork, unsigned n,
                  const double * x, const double * y )
{
  struct TaskXY info = { 4 , 0 , NULL , NULL };
  double result[4] = { 0 , 0 , 0 , 0 };
  info.n = n ;
  info.x = x ;
  info.y = y ;

  if ( nwork ) {
    TPI_Run_reduce( work_d4_dot_tp , & info , nwork ,
                    reduce_join, reduce_init, sizeof(result) , result );
  }
  else {
    TPI_Run_threads_reduce( work_d4_dot_tp , & info ,
                            reduce_join, reduce_init, sizeof(result), result);
  }

  comm_reduce_d4_sum( comm , result );

  d2_add_d( result , result[2] );
  d2_add_d( result , result[3] );

  return result[0] ;
}

static
void task_ddot_tp( TPI_Work * work )
{
  struct TaskXY * const info = (struct TaskXY *) work->info ;
  double        * const dst  = (double *) work->reduce ;
  unsigned int begin , length ;

  my_span( work->count , work->rank , info->n , & begin , & length );

  *dst += ddot( length , info->x + begin , info->y + begin );

  return ;
}

double ddot_tp( COMM comm, unsigned nwork, unsigned n,
                const double * x, const double * y )
{
  struct TaskXY info = { 1 , 0 , NULL , NULL };
  double result = 0 ;
  info.n = n ;
  info.x = x ;
  info.y = y ;

  if ( nwork ) {
    TPI_Run_reduce( task_ddot_tp , & info , nwork ,
                    reduce_join, reduce_init, sizeof(result), & result);
  }
  else {
    TPI_Run_threads_reduce( task_ddot_tp , & info ,
                            reduce_join, reduce_init, sizeof(result), & result);
  }

  comm_reduce_dsum( comm , & result );

  return result ;
}

/*--------------------------------------------------------------------*/

void dfill_rand( unsigned seed , unsigned n , double * x , double mag )
{
  const double scale = 2.0 * mag / (double) RAND_MAX ;
  double * const xe = x + n ;
  for ( ; xe != x ; ++x , ++seed ) {
    unsigned s = seed ;
    *x = scale * ((double) rand_r( & s )) - mag ;
  }
}

struct FillWork {
  double   mag ;
  double * beg ;
  unsigned length ;
  unsigned seed ;
};

static void task_dfill_rand( TPI_Work * work )
{
  struct FillWork * const w = (struct FillWork *) work->info ;

  unsigned int begin , length ;

  my_span( work->count, work->rank, w->length, & begin , & length );

  dfill_rand( w->seed + begin , length , w->beg + begin , w->mag );
}

void dfill_rand_tp( unsigned nblock , unsigned seed ,
                    unsigned n , double * x , double mag )
{
  struct FillWork data ;
  data.mag    = mag ;
  data.beg    = x ;
  data.length = n ;
  data.seed   = seed ;
  if ( nblock ) {
    const int nwork = ( n + nblock - 1 ) / nblock ;
    TPI_Run( & task_dfill_rand , & data , nwork , 0 );
  }
  else {
    TPI_Run_threads( & task_dfill_rand , & data , 0 );
  }
}

/*--------------------------------------------------------------------*/

static
void test_ddot_performance(
  COMM comm ,
  const int nthreads ,
  const int nblock ,
  const unsigned int num_trials ,
  const unsigned int num_tests ,
  const unsigned int length_array[]  /* Global array length for each test */ ,
  const double   mag )
{
  const unsigned int ddot_flop   = 2 ;  /* 1 mult, 1 sum */
  const unsigned int d4_dot_flop = 12 ; /* 1 mult, 7 sum, 4 compare */

  const unsigned int p_rank = comm_rank( comm );
  const unsigned int p_size = comm_size( comm );

  const unsigned int max_array = length_array[ num_tests - 1 ];

  unsigned int local_max_size = 0 ;
  unsigned int i_test ;

  TPI_Init( nthreads );

  if ( 0 == p_rank ) {
    fprintf(stdout,"\n\"DDOT and D4DOT Performance testing\"\n");
    fprintf(stdout,"\"MPI size = %u , TPI size = %d , BlockSize = %d , #Trials = %u\"\n",p_size,nthreads,nblock,num_trials);
    fprintf(stdout,"\"TEST\" , \"LENGTH\" , \"#CYCLE\" , \"DT-MEAN\" , \"DT-STDDEV\" , \"MFLOP-MEAN\" , \"MFLOP-STDDEV\"\n");
  }

  for ( i_test = 0 ; i_test < num_tests ; ++i_test ) {
    const unsigned length = length_array[ i_test ]; /* Global */
    const unsigned ncycle = 2 * max_array / length ;
    const unsigned local_max = ncycle * ( ( length + p_size - 1 ) / p_size );
    if ( local_max_size < local_max ) { local_max_size = local_max ; }
  }

  {
    double * const x = (double*) malloc(local_max_size * 2 * sizeof(double));
    double * const y = x + local_max_size ;

    unsigned int i , j ;

    dfill_rand_tp( nblock, 0,              local_max_size, x, mag );
    dfill_rand_tp( nblock, local_max_size, local_max_size, y, mag );

    for ( i_test = 0 ; i_test < num_tests ; ++i_test ) {
      const unsigned length = length_array[ i_test ]; /* Global */
      const unsigned ncycle = 2 * max_array / length ;

      unsigned int local_begin , local_length , local_nwork ;

      double dt_sum = 0.0 ;
      double dt_sum_2 = 0.0 ;

      my_span( p_size, p_rank, length, & local_begin , & local_length );

      local_nwork = nblock ? ( local_length + nblock - 1 ) / nblock : 0 ;

      /*--------------------------------------------------------------*/

      for ( i = 0 ; i < num_trials ; ++i ) {
        double dt = TPI_Walltime();
        for ( j = 0 ; j < ncycle ; ++j ) {
            ddot_tp( comm, local_nwork, local_length,
                     x + j * local_length ,
                     y + j * local_length );
        }
        dt = TPI_Walltime() - dt ;
        comm_reduce_dmax( comm , & dt );
        dt_sum   += dt ;
        dt_sum_2 += dt * dt ;
      }

      if ( 0 == p_rank ) {
        const double mflop = ((double)( ddot_flop * length * ncycle ) ) / ((double) 1e6 );

        const double dt_mean = dt_sum / num_trials ;
        const double dt_sdev = sqrt( ( num_trials * dt_sum_2 - dt_sum * dt_sum ) /
                                     ( num_trials * ( num_trials - 1 ) ) );
        const double mflop_mean = mflop / dt_mean ;
        const double mflop_sdev = mflop_mean * dt_sdev / ( dt_mean + dt_sdev );

        fprintf(stdout,"\"DDOT\"  , %8u , %8u , %9.5g , %9.5g , %9.5g , %9.5g\n",
                length, ncycle, dt_mean, dt_sdev, mflop_mean, mflop_sdev );
        fflush(stdout);
      }
    }

    for ( i_test = 0 ; i_test < num_tests ; ++i_test ) {
      const unsigned length = length_array[ i_test ]; /* Global */
      const unsigned ncycle = 2 * max_array / length ;

      unsigned int local_begin , local_length , local_nwork ;

      double dt_sum = 0 ;
      double dt_sum_2 = 0 ;

      my_span( p_size, p_rank, length, & local_begin , & local_length );

      local_nwork = nblock ? ( local_length + nblock - 1 ) / nblock : 0 ;

      /*--------------------------------------------------------------*/

      for ( i = 0 ; i < num_trials ; ++i ) {
        double dt = TPI_Walltime();
        for ( j = 0 ; j < ncycle ; ++j ) {
            d4_dot_tp( comm, local_nwork, local_length,
                       x + j * local_length ,
                       y + j * local_length );
        }
        dt = TPI_Walltime() - dt ;
        comm_reduce_dmax( comm , & dt );
        dt_sum   += dt ;
        dt_sum_2 += dt * dt ;
      }

      if ( 0 == p_rank ) {
        const double mflop = ((double)( d4_dot_flop * length * ncycle ) ) / ((double) 1e6 );

        const double dt_mean = dt_sum / num_trials ;
        const double dt_sdev = sqrt( ( num_trials * dt_sum_2 - dt_sum * dt_sum ) /
                                     ( num_trials * ( num_trials - 1 ) ) );
        const double mflop_mean = mflop / dt_mean ;
        const double mflop_sdev = mflop_mean * dt_sdev / ( dt_mean + dt_sdev );

        fprintf(stdout,"\"D4DOT\" , %8u , %8u , %9.5g , %9.5g , %9.5g , %9.5g\n",
                length, ncycle, dt_mean, dt_sdev, mflop_mean, mflop_sdev );
        fflush(stdout);
      }
    }

    /*--------------------------------------------------------------*/

    free( x );
  }

  TPI_Finalize();

  return ;
}

/*--------------------------------------------------------------------*/

static
void test_ddot_accuracy(
  COMM comm ,
  const int nthreads ,
  const int nblock ,
  const unsigned int num_tests ,
  const unsigned int length_array[]  /* Global array length for each test */ ,
  const double   mag )
{
  const unsigned int p_rank = comm_rank( comm );
  const unsigned int p_size = comm_size( comm );

  const unsigned int max_array = length_array[ num_tests - 1 ];
  const unsigned int local_max_size = ( max_array + p_size - 1 ) / p_size ;

  unsigned int i_test ;

  TPI_Init( nthreads );

  if ( 0 == p_rank ) {
    fprintf(stdout,"\n\"DDOT and D4DOT Accuracy testing\"\n");
    fprintf(stdout,"\"MPI size = %u , TPI size = %d , BlockSize = %d\"\n",p_size,nthreads,nblock);
    fprintf(stdout,"\"TEST\" , \"LENGTH\" , \"VALUE\"\n");
  }

  {
    double * const x = (double*) malloc(local_max_size * 2 * sizeof(double));
    double * const y = x + local_max_size ;

    for ( i_test = 0 ; i_test < num_tests ; ++i_test ) {
      const unsigned length      = length_array[ i_test ]; /* Global */
      const unsigned length_half = length / 2 ;

      unsigned local_begin , local_length , local_nwork ;

      double val_ddot ;

      my_span( p_size, p_rank, length, & local_begin , & local_length );

      local_nwork = nblock ? ( local_length + nblock - 1 ) / nblock : 0 ;

      /*--------------------------------------------------------------*/

      if ( local_begin < length_half ) {
        const unsigned len = local_length < length_half - local_begin
                           ? local_length : length_half - local_begin ;

        dfill_rand_tp( nblock,          local_begin, len, x, mag );
        dfill_rand_tp( nblock, length + local_begin, len, y, mag );
      }

      if ( length_half < local_begin + local_length ) {
        const unsigned beg = length_half > local_begin
                           ? length_half : local_begin ;
        const unsigned off = beg - local_begin ;
        const unsigned len = local_length - off ;

        dfill_rand_tp( nblock,          beg - length_half, len, x + off, mag );
        dfill_rand_tp( nblock, length + beg - length_half, len, y + off, - mag );
      }

      /*--------------------------------------------------------------*/

      val_ddot = ddot_tp( comm, local_nwork, local_length, x, y );

      if ( 0 == p_rank ) {
        fprintf(stdout,"\"DDOT\"  , %8u , %9.3g\n", length , val_ddot );
        fflush(stdout);
      }
    }

    for ( i_test = 0 ; i_test < num_tests ; ++i_test ) {
      const unsigned length      = length_array[ i_test ]; /* Global */
      const unsigned length_half = length / 2 ;

      unsigned local_begin , local_length , local_nwork ;

      double val_d4_dot ;

      my_span( p_size, p_rank, length, & local_begin , & local_length );

      local_nwork = nblock ? ( local_length + nblock - 1 ) / nblock : 0 ;

      /*--------------------------------------------------------------*/

      if ( local_begin < length_half ) {
        const unsigned len = local_length < length_half - local_begin
                           ? local_length : length_half - local_begin ;

        dfill_rand_tp( nblock,          local_begin, len, x, mag );
        dfill_rand_tp( nblock, length + local_begin, len, y, mag );
      }

      if ( length_half < local_begin + local_length ) {
        const unsigned beg = length_half > local_begin
                           ? length_half : local_begin ;
        const unsigned off = beg - local_begin ;
        const unsigned len = local_length - off ;

        dfill_rand_tp( nblock,          beg - length_half, len, x + off, mag );
        dfill_rand_tp( nblock, length + beg - length_half, len, y + off, - mag );
      }

      /*--------------------------------------------------------------*/

      val_d4_dot = d4_dot_tp( comm, local_nwork, local_length, x , y );

      if ( 0 == p_rank ) {
        fprintf(stdout,"\"DDOT\"  , %8u , %9.3g\n", length , val_d4_dot );
        fflush(stdout);
      }
    }

    /*--------------------------------------------------------------*/

    free( x );
  }

  TPI_Finalize();

  return ;
}

/*--------------------------------------------------------------------*/

const unsigned test_lengths[] = 
  { 1e4 , 2e4 , 5e4 ,
    1e5 , 2e5 , 5e5 ,
    1e6 , 2e6 , 5e6 , 1e7 };

const unsigned test_count = sizeof(test_lengths) / sizeof(unsigned);
const unsigned nblock = 2500 ;

const double test_mag = 1e4 ;

static void test_performance(
  COMM comm , const int test_thread_count , const int test_thread[] )
{
  const unsigned num_trials = 11 ;

  int i ;

  for ( i = 0 ; i < test_thread_count ; ++i ) {

    test_ddot_performance( comm , test_thread[i] , nblock,
                           num_trials , test_count , test_lengths , test_mag );

    test_ddot_performance( comm , test_thread[i] , 0,
                           num_trials , test_count , test_lengths , test_mag );
  }
}

static void test_accuracy(
  COMM comm , const int test_thread_count , const int test_thread[] ,
              unsigned test_do )
{
  int i ;

  if ( test_count < test_do ) { test_do = test_count ; }

  for ( i = 0 ; i < test_thread_count ; ++i ) {

    test_ddot_accuracy( comm, test_thread[i], nblock,
                        test_do, test_lengths, test_mag );

    test_ddot_accuracy( comm, test_thread[i], 0,
                        test_do, test_lengths, test_mag );
  }
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#define TEST_THREAD_MAX 128

#if defined(HAVE_MPI)

int main( int argc , char **argv )
{
  int nthread[ TEST_THREAD_MAX ];
  int i ;

  MPI_Init( & argc , & argv );

  for ( i = 0 ; i < TEST_THREAD_MAX ; ++i ) { nthread[i] = 0 ; }

  if ( 0 == comm_rank( MPI_COMM_WORLD ) ) {
    if ( 1 < argc && argc < TEST_THREAD_MAX ) {
      nthread[0] = 1 ;
      nthread[1] = argc - 1 ;
      for ( i = 1 ; i < argc ; ++i ) { nthread[i+1] = atoi( argv[i] ); }
    }
    else {
      nthread[0] = 0 ;
      nthread[1] = 1 ;
      nthread[2] = 1 ;
    }
  }

  MPI_Bcast( nthread , TEST_THREAD_MAX , MPI_INT , 0 , MPI_COMM_WORLD );

  if ( nthread[0] ) {
    test_accuracy(    MPI_COMM_WORLD , nthread[1] , nthread + 2 , test_count );
    test_performance( MPI_COMM_WORLD , nthread[1] , nthread + 2 );
  }
  else {
    test_accuracy(    MPI_COMM_WORLD , nthread[1] , nthread + 2 , 3 );
  }

  MPI_Finalize();

  return 0 ;
}

static int comm_size( COMM comm )
{
  int size = 0 ;
  MPI_Comm_size( comm , & size );
  return size ;
}

static int comm_rank( COMM comm )
{
  int rank = 0 ;
  MPI_Comm_rank( comm , & rank );
  return rank ;
}

static void comm_reduce_dmax( COMM comm , double * val )
{
  double tmp ;
  if ( MPI_SUCCESS ==
       MPI_Allreduce( val , & tmp , 1 , MPI_DOUBLE , MPI_MAX , comm ) ) {
    *val = tmp ;
  }
  else {
    *val = 0 ;
  }
}

static void comm_reduce_dsum( COMM comm , double * val )
{
  double tmp ;
  if ( MPI_SUCCESS ==
       MPI_Allreduce( val , & tmp , 1 , MPI_DOUBLE , MPI_SUM , comm ) ) {
    *val = tmp ;
  }
  else {
    *val = 0 ;
  }
}

static void comm_reduce_d4_op( void * argin ,
                               void * argout ,
                               int * n ,
                               MPI_Datatype * d )
{
  if ( d && n && *n == 4 ) {
    double * const in  = (double*) argin ;
    double * const out = (double*) argout ;
    d2_add_d( out ,     in[0] );
    d2_add_d( out ,     in[1] );
    d2_add_d( out + 2 , in[2] );
    d2_add_d( out + 2 , in[3] );
  }
  return ; 
}

static void comm_reduce_d4_sum( COMM comm , double * val )
{
  double tmp[4] ;
  MPI_Op mpi_op = MPI_OP_NULL ;

  /* Use Reduce->Bcast instead of Allreduce due to a bug with the SUN MPI. */

  MPI_Op_create( comm_reduce_d4_op , 0 , & mpi_op );
  MPI_Reduce( val , tmp , 4 , MPI_DOUBLE , mpi_op , 0 , comm );
  MPI_Bcast(        tmp , 4 , MPI_DOUBLE ,          0 , comm );
  MPI_Op_free( & mpi_op );

  val[0] = tmp[0] ;
  val[1] = tmp[1] ;
  val[2] = tmp[2] ;
  val[3] = tmp[3] ;
}

#else

int main( int argc , char **argv )
{
  int nthread[ TEST_THREAD_MAX ];
  int i ;

  for ( i = 0 ; i < TEST_THREAD_MAX ; ++i ) { nthread[i] = 0 ; }

  if ( 1 < argc && argc < TEST_THREAD_MAX ) {
    nthread[0] = 1 ;
    nthread[1] = argc - 1 ;
    for ( i = 1 ; i < argc ; ++i ) { nthread[i+1] = atoi( argv[i] ); }
  }
  else {
    nthread[0] = 0 ;
    nthread[1] = 4 ;
    nthread[2] = 1 ;
    nthread[3] = 2 ;
    nthread[4] = 4 ;
    nthread[5] = 8 ;
  }

  if ( nthread[0] ) {
    test_accuracy(    0 , nthread[1] , nthread + 2 , test_count );
    test_performance( 0 , nthread[1] , nthread + 2 );
  }
  else {
    test_accuracy(    0 , nthread[1] , nthread + 2 , 3 );
  }

  return 0 ;
}

static int comm_size( COMM comm ) { return comm ? -1 : 1 ; }
static int comm_rank( COMM comm ) { return comm ? -1 : 0 ; }
static void comm_reduce_dmax( COMM comm , double * val )
{
  if ( comm ) { *val = 0 ; }
  return ;
}
static void comm_reduce_dsum( COMM comm , double * val )
{
  if ( comm ) { *val = 0 ; }
  return ;
}
static void comm_reduce_d4_sum( COMM comm , double * val )
{
  if ( comm ) { val[0] = val[1] = val[2] = val[3] = 0 ; }
  return ;
}

#endif

/*--------------------------------------------------------------------*/

