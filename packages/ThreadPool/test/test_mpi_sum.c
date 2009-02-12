/*------------------------------------------------------------------------*/
/*                    TPI: Thread Pool Interface                          */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
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
  const unsigned int len = size - ( *begin = max * rank );
  *length = max < len ? max : len ;
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
  double pos[2] = { 0 , 0 };
  double neg[2] = { 0 , 0 };
  const double * const x_end = x + n ;
  for ( ; x < x_end ; ++x , ++y ) {
    const double a = *x * *y ;
    if ( a < 0 ) { d2_add_d( neg , a ); }
    else         { d2_add_d( pos , a ); }
  }
  v[0] = pos[0] ;
  v[1] = pos[1] ;
  v[2] = neg[0] ;
  v[3] = neg[1] ;
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
  unsigned n ;
  const double * x ;
  const double * y ;
        double   v[4] ;
};

static
void task_d4_dot_tp( TPI_Work * work )
{
  struct TaskXY * const data = (struct TaskXY *) work->shared ;

  double tmp[4] = { 0 , 0 , 0 , 0 };

  unsigned int begin , length ;

  my_span( work->work_count , work->work_rank , data->n , & begin , & length );

  d4_dot( tmp , length , data->x + begin , data->y + begin );

  TPI_Lock(0);
  d2_add_d( data->v ,     tmp[0] );
  d2_add_d( data->v ,     tmp[1] );
  d2_add_d( data->v + 2 , tmp[2] );
  d2_add_d( data->v + 2 , tmp[3] );
  TPI_Unlock(0);
}

double d4_dot_tp( COMM comm, unsigned n, unsigned nblock ,
                  const double * x, const double * y )
{
  struct TaskXY data = { 0 , NULL , NULL , { 0 , 0 , 0 , 0 } };
  data.n = n ;
  data.x = x ;
  data.y = y ;

  if ( nblock ) {
    const int nwork = ( n + nblock - 1 ) / nblock ;
    TPI_Run( task_d4_dot_tp , & data , nwork , 1 );
  }
  else {
    TPI_Run_threads( task_d4_dot_tp , & data , 1 );
  }

  comm_reduce_d4_sum( comm , data.v );

  d2_add_d( data.v , data.v[2] );
  d2_add_d( data.v , data.v[3] );

  return data.v[0] ;
}

static
void task_ddot_tp( TPI_Work * work )
{
  struct TaskXY * const data = (struct TaskXY *) work->shared ;
  double tmp ;
  unsigned int begin , length ;

  my_span( work->work_count , work->work_rank , data->n , & begin , & length );

  tmp = ddot( length , data->x + begin , data->y + begin );

  TPI_Lock(0);
  data->v[0] += tmp ;
  TPI_Unlock(0);
  return ;
}

double ddot_tp( COMM comm, unsigned n, unsigned nblock,
                const double * x, const double * y )
{
  struct TaskXY data = { 0 , NULL , NULL , { 0 , 0 , 0 , 0 } };
  data.n = n ;
  data.x = x ;
  data.y = y ;

  if ( nblock ) {
    const int nwork = ( n + nblock - 1 ) / nblock ;
    TPI_Run( task_ddot_tp , & data , nwork , 1 );
  }
  else {
    TPI_Run_threads( task_ddot_tp , & data , 1 );
  }

  comm_reduce_dsum( comm , data.v );

  return data.v[0] ;
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
  struct FillWork * const w = (struct FillWork *) work->shared ;

  unsigned int begin , length ;

  my_span( work->work_count, work->work_rank, w->length, & begin , & length );

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

      unsigned int local_begin , local_length ;

      double dt_sum = 0.0 ;
      double dt_sum_2 = 0.0 ;

      my_span( p_size, p_rank, length, & local_begin , & local_length );

      /*--------------------------------------------------------------*/

      for ( i = 0 ; i < num_trials ; ++i ) {
        double dt = TPI_Walltime();
        for ( j = 0 ; j < ncycle ; ++j ) {
            ddot_tp( comm, local_length, nblock,
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

      unsigned int local_begin , local_length ;

      double dt_sum = 0 ;
      double dt_sum_2 = 0 ;

      my_span( p_size, p_rank, length, & local_begin , & local_length );

      /*--------------------------------------------------------------*/

      for ( i = 0 ; i < num_trials ; ++i ) {
        double dt = TPI_Walltime();
        for ( j = 0 ; j < ncycle ; ++j ) {
            d4_dot_tp( comm, local_length, nblock,
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

        fprintf(stdout,"\"DDOT\"  , %8u , %8u , %9.5g , %9.5g , %9.5g , %9.5g\n",
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
      const int      length_half = length / 2 ;

      const unsigned local_begin = ( length * p_rank ) / p_size ;
      const unsigned local_end   = ( length * ( p_rank + 1 ) ) / p_size ;
      const unsigned local_size  = local_end - local_begin ;

      double val_ddot ;

      /*--------------------------------------------------------------*/

      if ( local_begin < length_half ) {
        const unsigned len = local_size < length_half - local_begin
                           ? local_size : length_half - local_begin ;

        dfill_rand_tp( nblock,          local_begin, len, x, mag );
        dfill_rand_tp( nblock, length + local_begin, len, y, mag );
      }

      if ( length_half < local_begin + local_size ) {
        const unsigned beg = length_half > local_begin
                           ? length_half : local_begin ;
        const unsigned off = beg - local_begin ;
        const unsigned len = local_size - off ;

        dfill_rand_tp( nblock,          beg - length_half, len, x + off, mag );
        dfill_rand_tp( nblock, length + beg - length_half, len, y + off, - mag );
      }

      /*--------------------------------------------------------------*/

      val_ddot = ddot_tp( comm, local_size, nblock, x, y );

      if ( 0 == p_rank ) {
        fprintf(stdout,"\"DDOT\"  , %8u , %9.3g\n", length , val_ddot );
        fflush(stdout);
      }
    }

    for ( i_test = 0 ; i_test < num_tests ; ++i_test ) {
      const unsigned length      = length_array[ i_test ]; /* Global */
      const int      length_half = length / 2 ;

      const unsigned local_begin = ( length * p_rank ) / p_size ;
      const unsigned local_end   = ( length * ( p_rank + 1 ) ) / p_size ;
      const unsigned local_size  = local_end - local_begin ;

      double val_d4_dot ;

      /*--------------------------------------------------------------*/

      if ( local_begin < length_half ) {
        const unsigned len = local_size < length_half - local_begin
                           ? local_size : length_half - local_begin ;

        dfill_rand_tp( nblock,          local_begin, len, x, mag );
        dfill_rand_tp( nblock, length + local_begin, len, y, mag );
      }

      if ( length_half < local_begin + local_size ) {
        const unsigned beg = length_half > local_begin
                           ? length_half : local_begin ;
        const unsigned off = beg - local_begin ;
        const unsigned len = local_size - off ;

        dfill_rand_tp( nblock,          beg - length_half, len, x + off, mag );
        dfill_rand_tp( nblock, length + beg - length_half, len, y + off, - mag );
      }

      /*--------------------------------------------------------------*/

      val_d4_dot = d4_dot_tp( comm, local_size, nblock, x , y );

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

static void test_main( COMM comm , int nthreads )
{
  const unsigned lengths[] = { 1e4 , 2e4 , 5e4 ,
                               1e5 , 2e5 , 5e5 ,
                               1e6 , 2e6 , 5e6 , 1e7 };

  const unsigned nblock = 1000 ;
  const unsigned num_tests = sizeof(lengths) / sizeof(unsigned);
  const unsigned num_trials = 11 ;
  const double mag = 1e4 ;

  test_ddot_accuracy( comm , nthreads , nblock, num_tests, lengths, mag );
  test_ddot_accuracy( comm , nthreads , 0, num_tests, lengths, mag );

  test_ddot_performance( comm , nthreads , nblock,
                         num_trials , num_tests , lengths , mag );

  test_ddot_performance( comm , nthreads , 0,
                         num_trials , num_tests , lengths , mag );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#if defined(HAVE_MPI)

int main( int argc , char **argv )
{
  MPI_Init( & argc , & argv );
  {
    int i ;
    for ( i = 1 ; i < argc ; ++i ) {
      int nthreads = atoi( argv[i] );
      MPI_Bcast( & nthreads , 1 , MPI_INT , 0 , MPI_COMM_WORLD );
      test_main( MPI_COMM_WORLD , nthreads );
    }
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
  if ( 1 < argc ) {
    int i ;
    for ( i = 1 ; i < argc ; ++i ) {
      int nthreads = atoi( argv[i] );
      test_main( 0 , nthreads );
    }
  }
  else {
    int nthreads ;
    for ( nthreads = 1 ; nthreads <= 16 ; nthreads *= 2 ) {
      test_main( 0 , nthreads );
    }
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

