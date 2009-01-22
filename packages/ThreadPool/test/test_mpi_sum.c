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

#include <stdlib.h>
#include <stdio.h>
#include <TPI.h>

int rand_r( unsigned int * );

/*--------------------------------------------------------------------*/

#if defined(TEST_WITH_MPI)

#include <mpi.h>

typedef MPI_Comm COMM ;

#else

typedef int COMM ;

#endif

static int comm_size( COMM );
static int comm_rank( COMM );
static void comm_reduce_dsum( COMM , double * );
static void comm_reduce_d4_sum( COMM , double * );

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
  unsigned nb ;
  const double * x ;
  const double * y ;
        double * v ;
};

static
void task_d4_dot_tp( TPI_Work * work )
{
  struct TaskXY * const data = (struct TaskXY *) work->shared ;

  const int begin  = data->nb * work->work_rank ;
  const int length = data->n - begin < data->nb ?
                     data->n - begin : data->nb ;

  d4_dot( data->v + 4 * work->work_rank ,
          length , data->x + begin , data->y + begin );
}

double d4_dot_tp( COMM comm, unsigned n, unsigned nblock ,
                  const double * x, const double * y )
{
  const int nwork = ( n + nblock - 1 ) / nblock ;
  const int ntmp = 4 * nwork ;
  double tmp[ ntmp ];

  struct TaskXY data = { n , nblock , x , y , tmp };

  int i ;

  for ( i = 0 ; i < ntmp ; ++i ) { tmp[i] = 0 ; }

  TPI_Run( task_d4_dot_tp , & data , nwork , 0 );

  for ( i = 1 ; i < nwork ; ++i ) {
    d2_add_d( tmp ,   tmp[  i*4] ); 
    d2_add_d( tmp ,   tmp[1+i*4] ); 
    d2_add_d( tmp+2 , tmp[2+i*4] ); 
    d2_add_d( tmp+2 , tmp[3+i*4] ); 
  }

  comm_reduce_d4_sum( comm , tmp );

  d2_add_d( tmp , tmp[2] );
  d2_add_d( tmp , tmp[3] );

  return tmp[0] ;
}

static
void task_ddot_tp( TPI_Work * work )
{
  struct TaskXY * const data = (struct TaskXY *) work->shared ;

  const int begin  = data->nb * work->work_rank ;
  const int length = data->n - begin < data->nb ?
                     data->n - begin : data->nb ;

  data->v[ work->work_rank ] =
    ddot( length , data->x + begin , data->y + begin );
  return ;
}

double ddot_tp( COMM comm, unsigned n, unsigned nblock,
                const double * x, const double * y )
{
  const int nwork = ( n + nblock - 1 ) / nblock ;
  double tmp[ nwork ];

  struct TaskXY data = { n , nblock , x , y , tmp };

  int i ;

  for ( i = 0 ; i < nwork ; ++i ) { tmp[i] = 0 ; }

  TPI_Run( task_ddot_tp , & data , nwork , 0 );

  for ( i = 1 ; i < nwork ; ++i ) { tmp[0] += tmp[i] ; }

  comm_reduce_dsum( comm , tmp );

  return tmp[0] ;
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
  unsigned nblock ;
  unsigned seed ;
};

static void task_dfill_rand( TPI_Work * work )
{
  struct FillWork * const w = (struct FillWork *) work->shared ;

  const int begin  = w->nblock * work->work_rank ;
  const int length = w->length - begin < w->nblock ?
                     w->length - begin : w->nblock ;

  dfill_rand( w->seed + begin , length , w->beg + begin , w->mag );
}

void dfill_rand_tp( unsigned seed , unsigned n , unsigned nblock ,
                    double * x , double mag )
{
  const int nwork = ( n + nblock - 1 ) / nblock ;
  struct FillWork data = { mag , x , n , nblock , seed };
  TPI_Run( & task_dfill_rand , & data , nwork , 0 );
}

/*--------------------------------------------------------------------*/

static
void test_tpi_ddot_driver(
  COMM comm ,
  const int nthreads ,
  const int nblock ,
  const unsigned Mflop_target  /* Per concurrent thread */ ,
  const unsigned num_trials ,
  const unsigned num_tests ,
  const unsigned length_array[]  /* Global array length for each test */ ,
  const double   mag )
{
  const unsigned num_times = num_trials * num_tests ;
  double dt_ddot[ num_times ];
  double dt_d4_dot[ num_times ];
  double val_d4_dot[ num_tests ] ;
  double val_ddot[ num_tests ] ;

  const unsigned ddot_flop   = 2 ;  /* 1 mult, 1 sum */
  const unsigned d4_dot_flop = 12 ; /* 1 mult, 7 sum, 4 compare */

  const unsigned p_rank = comm_rank( comm );
  const unsigned p_size = comm_size( comm );
  const unsigned np = p_size * nthreads ;

  const unsigned max_array = length_array[ num_tests - 1 ];

  const unsigned local_max_begin = ( max_array * p_rank ) / p_size ;
  const unsigned local_max_end   = ( max_array * ( p_rank + 1 ) ) / p_size ;
  const unsigned local_max_size  = local_max_end - local_max_begin ;

  TPI_Init( nthreads );

  {
    double * const x = (double*) malloc(local_max_size * 2 * sizeof(double));
    double * const y = x + local_max_size ;

    unsigned i_test , i , j ;

    for ( i_test = 0 ; i_test < num_tests ; ++i_test ) {
      const unsigned length      = length_array[ i_test ]; /* Global */
      const int      length_half = length / 2 ;
      const unsigned num_array   = max_array / length ;

      const unsigned ddot_ncycle =
        1 + (unsigned)( ( Mflop_target * np * 1e6 ) / ( ddot_flop * length ) );

      const unsigned d4_dot_ncycle =
        1 + (unsigned)( ( Mflop_target * np * 1e6 ) / ( d4_dot_flop * length));

      const unsigned local_begin = ( length * p_rank ) / p_size ;
      const unsigned local_end   = ( length * ( p_rank + 1 ) ) / p_size ;
      const unsigned local_size  = local_end - local_begin ;

      /*--------------------------------------------------------------*/

      if ( local_begin < length_half ) {
        const unsigned len = local_size < length_half - local_begin
                           ? local_size : length_half - local_begin ;

        for ( i = 0 ; i < num_array ; ++i ) {
          const int n = i * local_size ;
          dfill_rand_tp(          local_begin, len, nblock, x + n, mag );
          dfill_rand_tp( length + local_begin, len, nblock, y + n, mag );
        }
      }

      if ( length_half < local_begin + local_size ) {
        const unsigned beg = length_half > local_begin
                           ? length_half : local_begin ;
        const unsigned off = beg - local_begin ;
        const unsigned len = local_size - off ;

        for ( i = 0 ; i < num_array ; ++i ) {
          const int n = i * local_size + off ;
          dfill_rand_tp(          beg - length_half, len, nblock, x + n, mag );
          dfill_rand_tp( length + beg - length_half, len, nblock, y + n, - mag );
        }
      }

      /*--------------------------------------------------------------*/

      for ( i = 0 ; i < num_trials ; ++i ) {
        {
          const unsigned n = i * d4_dot_ncycle ;
          const double   t = TPI_Walltime();
          for ( j = 0 ; j < d4_dot_ncycle ; ++j ) {
            const unsigned k = ( ( j + n ) % num_array ) * local_size ;
            val_d4_dot[i_test] =
              d4_dot_tp( comm, local_size, nblock, x + k, y + k);
          }
          dt_d4_dot[ i + i_test * num_trials ] = TPI_Walltime() - t ;
        }
      }

      /*--------------------------------------------------------------*/

      for ( i = 0 ; i < num_trials ; ++i ) {
        {
          const unsigned n = i * ddot_ncycle ;
          const double   t = TPI_Walltime();
          for ( j = 0 ; j < ddot_ncycle ; ++j ) {
            const unsigned k = ( ( j + n ) % num_array ) * local_size ;
            val_ddot[i_test] = ddot_tp( comm, local_size, nblock, x + k, y + k);
          }
          dt_ddot[ i + i_test * num_trials ] = TPI_Walltime() - t ;
        }
      }
    }

    free( x );
  }

  TPI_Finalize();

  /*------------------------------------------------------------------*/

  if ( 0 == p_rank ) {

    unsigned i_test , i ;

    fprintf(stdout,"\n\"DDOT and D4DOT Performance testing\"\n");
    fprintf(stdout,"\"MPI size = %u , TPI size = %d , Bounds = %g , #Trials = %u\"\n",p_size,nthreads,mag,num_trials);
    fprintf(stdout,"\"TEST\" , \"LENGTH\" , \"VALUE\" , \"DT-MAX\" , \"DT-AVG\" , \"DT-MIN\" , \"MFLOP-MIN\" , \"MFLOP-AVG\" , \"MFLOP-MAX\"\n");

    for ( i_test = 0 ; i_test < num_tests ; ++i_test ) {
      const unsigned length = length_array[ i_test ];
      const unsigned ncycle =
        1 + (unsigned)( ( Mflop_target * np * 1e6 ) / ( ddot_flop * length ));
      const double mflop = ((double)( ddot_flop * length ) ) / ((double) 1e6 );

      const double * const dt = dt_ddot + i_test * num_trials ;
      double dt_max = dt[0] ;
      double dt_min = dt[0] ;
      double dt_avg = dt[0]  ;

      for ( i = 1 ; i < num_trials ; ++i ) {
        if ( dt_max < dt[i] ) { dt_max = dt[i] ; }
        if ( dt_min > dt[i] ) { dt_min = dt[i] ; }
        dt_avg += dt[i] ;
      }
      dt_avg /= (double) num_trials ;

      dt_max /= (double) ncycle ;
      dt_avg /= (double) ncycle ;
      dt_min /= (double) ncycle ;

      fprintf(stdout,"\"DDOT\"  , %8u , %9.3g , %9.5g , %9.5g , %9.5g , %9.5g , %9.5g , %9.5g\n",
              length , val_ddot[i_test] ,
              dt_max , dt_avg , dt_min,
              ( mflop / dt_max ) , ( mflop / dt_avg ) , ( mflop / dt_min ) );
    }

    for ( i_test = 0 ; i_test < num_tests ; ++i_test ) {
      const unsigned length = length_array[ i_test ];
      const unsigned ncycle =
        1 + (unsigned)(( Mflop_target * np * 1e6 ) / ( ddot_flop * length ));
      const double mflop = ((double)( d4_dot_flop * length ) ) /
                           ((double) 1e6 );

      const double * const dt = dt_d4_dot + i_test * num_trials ;
      double dt_max = dt[0] ;
      double dt_min = dt[0] ;
      double dt_avg = dt[0] ;

      for ( i = 1 ; i < num_trials ; ++i ) {
        if ( dt_max < dt[i] ) { dt_max = dt[i] ; }
        if ( dt_min > dt[i] ) { dt_min = dt[i] ; }
        dt_avg += dt[i] ;
      }
      dt_avg /= (double) num_trials ;

      dt_max /= (double) ncycle ;
      dt_avg /= (double) ncycle ;
      dt_min /= (double) ncycle ;

      fprintf(stdout,"\"D4DOT\" , %8u , %9.3g , %9.5g , %9.5g , %9.5g , %9.5g , %9.5g , %9.5g\n",
              length , val_d4_dot[i_test] ,
              dt_max , dt_avg , dt_min ,
              ( mflop / dt_max ) , ( mflop / dt_avg ) , ( mflop / dt_min ) );
    }
  }

  fflush(stdout);

  return ;
}

/*--------------------------------------------------------------------*/

static void test_main( COMM comm , int nthreads )
{
  const unsigned lengths[] = { 1e3 , 2e3 , 5e3 ,
                               1e4 , 2e4 , 5e4 ,
                               1e5 , 2e5 , 5e5 ,
                               1e6 , 2e6 , 5e6 , 1e7 };

  const unsigned nblock = 1000 ;
  const unsigned num_tests = sizeof(lengths) / sizeof(unsigned);
  const unsigned num_trials = 5 ;
  const unsigned mflop_target = 100 ;
  const double mag = 1e4 ;

  test_tpi_ddot_driver( comm , nthreads , nblock,
                        mflop_target , num_trials ,
                        num_tests , lengths , mag );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#if defined(TEST_WITH_MPI)

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

