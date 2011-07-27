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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <TPI.h>

#if defined( HAVE_MPI )
#include <mpi.h>
#endif

/*--------------------------------------------------------------------*/

static void test_work( TPI_Work * );
static void test_reduce_work( TPI_Work * );
static void test_reduce_init( TPI_Work * );
static void test_reduce_join( TPI_Work * , const void * );
static void test_reduce_via_lock( TPI_Work * );
static void test_reduce_via_nolock( TPI_Work * );

void test_tpi_init(   const int ntest, const int nthread[], const int ntrial);
void test_tpi_block(  const int ntest, const int nthread[], const int ntrial);
void test_tpi_reduce( const int ntest, const int nthread[], const int ntrial);
void test_tpi_work(   const int ntest, const int nthread[],
                      const int nwork , const int ntrial );
void test_tpi_work_async(
  const int ntest , const int nthread[] , const int nwork , const int ntrial );

int main( int argc , char ** argv )
{
  int num_thread[] = { 1 , 2 , 4 , 6 , 8 , 12 , 16 };
  int num_test = sizeof(num_thread) / sizeof(int);

#if defined( HAVE_MPI )
  int rank ;

  MPI_Init( & argc , & argv );
  MPI_Comm_rank( MPI_COMM_WORLD , & rank );
  if ( 0 == rank ) {
#endif
 
  const int ntrial = 1 < argc ? atoi( argv[1] ) : 5 ;
  const int nwork  = 2 < argc ? atoi( argv[2] ) : 100 ;
 
  /* Get the configuration print message out. */
  fprintf( stdout , "\"%s\"\n" , TPI_Version() );
  fprintf( stdout , "\"Unit Testing: ntrial = %d , nwork = %d\"\n" , ntrial , nwork );
 
  test_tpi_init(   num_test , num_thread , ntrial );
  test_tpi_block(  num_test , num_thread , ntrial );
  test_tpi_reduce( num_test , num_thread , ntrial );
  test_tpi_work(   num_test , num_thread , nwork , ntrial );
  test_tpi_work_async( num_test , num_thread , nwork , ntrial );
 
#if defined( HAVE_MPI )
  }
  MPI_Finalize();
#endif

  return 0 ;
}

/*--------------------------------------------------------------------*/

void test_tpi_init( const int ntest , const int nthread[] , const int ntrial )
{
  int j ;

  fprintf( stdout , "\n\"TEST TPI_Init / TPI_Finalize\"\n" );
  fprintf( stdout , "\"#Thread\" , \"#Trial\" , \"TPI_Init(avg-msec)\" , \"TPI_Init(stddev-msec)\" , \"TPI_Finalize(avg-msec)\" , \"TPI_Finalize(stddev-msec)\"\n");

  for ( j = 0 ; j < ntest ; ++j ) {
    const int nth = nthread[j];
    double dt_init_total   = 0.0 ;
    double dt_init_total_2 = 0.0 ;
    double dt_fin_total    = 0.0 ;
    double dt_fin_total_2  = 0.0 ;
    int i ;
    int result ;

    for ( i = 0 ; i < ntrial ; ++i ) {
      double t , dt ;

      t = TPI_Walltime();
      result = TPI_Init( nth );
      dt = TPI_Walltime() - t ;
      dt_init_total += dt ;
      dt_init_total_2 += dt * dt ;

      if ( result != nth ) {
        fprintf(stderr,"%d != TPI_Init(%d) : FAILED at trial %d\n",
                result , nth , i );
        abort();
      }

      t = TPI_Walltime();
      TPI_Finalize();
      dt = TPI_Walltime() - t ;
      dt_fin_total += dt ;
      dt_fin_total_2 += dt * dt ;
    }

    if ( 1 < ntrial ) {
      const double init_mean = 1.0e6 * dt_init_total / ntrial ;
      const double init_sdev = 1.0e6 * sqrt( ( ntrial * dt_init_total_2 -
                                       dt_init_total * dt_init_total ) /
                                     ( ntrial * ( ntrial - 1 ) ) );

      const double fin_mean = 1.0e6 * dt_fin_total / ntrial ;
      const double fin_sdev = 1.0e6 * sqrt( ( ntrial * dt_fin_total_2 -
                                      dt_fin_total * dt_fin_total ) /
                                    ( ntrial * ( ntrial - 1 ) ) );
      
      fprintf(stdout,"%d , %d , %10g , %10g , %10g , %10g\n",
              nth , ntrial , init_mean , init_sdev , fin_mean , fin_sdev );
    }
  }
}

/*--------------------------------------------------------------------*/

void test_tpi_block( const int ntest , const int nthread[] , const int ntrial )
{
  int i, j ;

  fprintf( stdout , "\n\"TEST TPI_Block / TPI_Unblock\"\n" );
  fprintf( stdout , "\"#Thread\" , \"#Trial\" , \"TPI_Block(avg-msec)\" , \"TPI_Block(stddev-msec)\" , \"TPI_Unblock(avg-msec)\" , \"TPI_Unblock(stddev-msec)\"\n");

  for ( j = 0 ; j < ntest ; ++j ) {
    const int nth = nthread[j];

    double dt_block_total   = 0.0 ;
    double dt_block_total_2 = 0.0 ;
    double dt_unblock_total    = 0.0 ;
    double dt_unblock_total_2  = 0.0 ;

    int result = TPI_Init( nth );

    if ( result != nth ) {
      fprintf(stderr,"%d != TPI_Init(%d) : FAILED\n", result , nth );
      abort();
    }

    for ( i = 0 ; i < ntrial ; ++i ) {
      double t , dt ;

      t = TPI_Walltime();
      TPI_Block();
      dt = TPI_Walltime() - t ;
      dt_block_total += dt ;
      dt_block_total_2 += dt * dt ;


      t = TPI_Walltime();
      TPI_Unblock();
      dt = TPI_Walltime() - t ;
      dt_unblock_total += dt ;
      dt_unblock_total_2 += dt * dt ;
    }

    TPI_Finalize();

    if ( 1 < ntrial ) {
      const double block_mean = 1.0e6 * dt_block_total / ntrial ;
      const double block_sdev = 1.0e6 * sqrt( ( ntrial * dt_block_total_2 -
                                        dt_block_total * dt_block_total ) /
                                      ( ntrial * ( ntrial - 1 ) ) );

      const double unblock_mean = 1.0e6 * dt_unblock_total / ntrial ;
      const double unblock_sdev = 1.0e6 * sqrt( ( ntrial * dt_unblock_total_2 -
                                          dt_unblock_total * dt_unblock_total) /
                                        ( ntrial * ( ntrial - 1 ) ) );
      
      fprintf(stdout,"%d , %d , %10g , %10g , %10g , %10g\n",
              nth , ntrial , block_mean , block_sdev , unblock_mean , unblock_sdev );
    }
  }
}

/*--------------------------------------------------------------------*/

void test_tpi_reduce( const int ntest , const int nthread[] , const int ntrial )
{
  int j ;

  fprintf( stdout , "\n\"TEST TPI_Run_threads(reduce) / TPI_Run_threads_reduce\"\n" );
  fprintf( stdout , "\"#Thread\" , \"#Trial\" , \"TPI_Run_threads(avg-msec)\" , \"TPI_Run_threads(stddev-msec)\" , \"TPI_Run_threads_reduce(avg-msec)\" , \"TPI_Run_threads_reduce(stddev-msec)\"\n");

  for ( j = 0 ; j < ntest ; ++j ) {
    const int nth = nthread[j];

    double dt_lock_total   = 0.0 ;
    double dt_lock_total_2 = 0.0 ;
    double dt_reduce_total    = 0.0 ;
    double dt_reduce_total_2  = 0.0 ;
    int i ;

    int result = TPI_Init( nth );

    if ( result != nth ) {
      fprintf(stderr,"%d != TPI_Init(%d) : FAILED\n", result , nth );
    }

    for ( i = 0 ; i < ntrial ; ++i ) {
      double t , dt ;
      int value = 0 ;
      int * const ptr = & value ;

      t = TPI_Walltime();
      TPI_Run_threads( test_reduce_via_lock , & ptr , 1 );
      dt = TPI_Walltime() - t ;
      dt_lock_total += dt ;
      dt_lock_total_2 += dt * dt ;

      if ( value != nth ) {
        fprintf(stderr,
                "TPI_Run_threads(reduce,...) : FAILED at trial %d\n",
                i );
        abort();
      }

      value = 0 ;

      t = TPI_Walltime();
      TPI_Run_threads_reduce( test_reduce_via_nolock , NULL ,
                              test_reduce_join , test_reduce_init ,
                              sizeof(value) , & value );
  
      dt = TPI_Walltime() - t ;
      dt_reduce_total += dt ;
      dt_reduce_total_2 += dt * dt ;

      if ( value != nth ) {
        fprintf(stderr,
                "TPI_Run_threads_reduce(...) : FAILED at trial %d\n",
                i );
        abort();
      }
    }

    TPI_Finalize();

    if ( 1 < ntrial ) {
      const double lock_mean = 1.0e6 * dt_lock_total / ntrial ;
      const double lock_sdev = 1.0e6 * sqrt( ( ntrial * dt_lock_total_2 -
                                       dt_lock_total * dt_lock_total ) /
                                     ( ntrial * ( ntrial - 1 ) ) );

      const double reduce_mean = 1.0e6 * dt_reduce_total / ntrial ;
      const double reduce_sdev = 1.0e6 * sqrt( ( ntrial * dt_reduce_total_2 -
                                         dt_reduce_total * dt_reduce_total) /
                                       ( ntrial * ( ntrial - 1 ) ) );
      
      fprintf(stdout,"%d , %d , %10g , %10g , %10g , %10g\n",
              nth, ntrial, lock_mean, lock_sdev, reduce_mean, reduce_sdev);
    }
  }
}

/*--------------------------------------------------------------------*/

void test_tpi_work( const int ntest , const int nthread[] , const int nwork ,
                    const int ntrial )
{
  int * const flags = (int *) malloc( sizeof(int) * nwork );
  int j ;

  fprintf( stdout , "\n\"TEST TPI_Run / TPI_Run_reduce\"\n" );
  fprintf( stdout , "\"#Thread\" , \"#Work\" , \"#Trial\" , \"TPI_Run(avg-msec)\" , \"TPI_Run(stddev-msec)\" , \"TPI_Run_reduce(avg-msec)\" , \"TPI_Run_reduce(stddev-msec)\"\n");

  for ( j = 0 ; j < ntest ; ++j ) {
    const int nth = nthread[j];

    double dt_work_total   = 0.0 ;
    double dt_work_total_2 = 0.0 ;
    double dt_reduce_total    = 0.0 ;
    double dt_reduce_total_2  = 0.0 ;
    int i , k ;

    int result = TPI_Init( nth );

    if ( result != nth ) {
      fprintf(stderr,"%d != TPI_Init(%d) : FAILED\n", result , nth );
    }

    for ( i = 0 ; i < ntrial ; ++i ) {
      double t , dt ;
      int value = 0 ;

      for ( k = 0 ; k < nwork ; ++k ) { flags[k] = 0 ; }

      t = TPI_Walltime();
      TPI_Run( test_work , & flags , nwork , 0 );
      dt = TPI_Walltime() - t ;
      dt_work_total += dt ;
      dt_work_total_2 += dt * dt ;

      for ( k = 0 ; k < nwork && flags[k] ; ++k );

      if ( k < nwork ) {
        fprintf(stderr, "TPI_Run(...) : FAILED at trial %d\n", i );
        abort();
      }

      for ( k = 0 ; k < nwork ; ++k ) { flags[k] = 0 ; }

      t = TPI_Walltime();
      TPI_Run_reduce( test_reduce_work , & flags , nwork ,
                      test_reduce_join , test_reduce_init ,
                      sizeof(value) , & value );
  
      dt = TPI_Walltime() - t ;
      dt_reduce_total += dt ;
      dt_reduce_total_2 += dt * dt ;

      for ( k = 0 ; k < nwork && flags[k] ; ++k );

      if ( value != nwork || k < nwork ) {
        fprintf(stderr, "TPI_Run_reduce(...) : FAILED at trial %d\n", i );
        abort();
      }
    }

    TPI_Finalize();

    if ( 1 < ntrial ) {
      const double work_mean = 1.0e6 * dt_work_total / ntrial ;
      const double work_sdev = 1.0e6 * sqrt( ( ntrial * dt_work_total_2 -
                                       dt_work_total * dt_work_total ) /
                                     ( ntrial * ( ntrial - 1 ) ) );

      const double reduce_mean = 1.0e6 * dt_reduce_total / ntrial ;
      const double reduce_sdev = 1.0e6 * sqrt( ( ntrial * dt_reduce_total_2 -
                                         dt_reduce_total * dt_reduce_total) /
                                       ( ntrial * ( ntrial - 1 ) ) );
      
      fprintf(stdout,"%d , %d , %d , %10g , %10g , %10g , %10g\n",
              nth, ntrial, nwork, work_mean, work_sdev, reduce_mean, reduce_sdev);
    }
  }

  free( flags );
}

/*--------------------------------------------------------------------*/

void test_tpi_work_async(
  const int ntest , const int nthread[] , const int nwork , const int ntrial )
{
  int * const flags = (int *) malloc( sizeof(int) * nwork );
  int j ;

  fprintf( stdout , "\n\"TEST TPI_Start / TPI_Start_reduce\"\n" );
  fprintf( stdout , "\"#Thread\" , \"#Work\" , \"#Trial\" , \"TPI_Start(avg-msec)\" , \"TPI_Start(stddev-msec)\" , \"TPI_Start_reduce(avg-msec)\" , \"TPI_Start_reduce(stddev-msec)\"\n");

  for ( j = 0 ; j < ntest ; ++j ) {
    const int nth = nthread[j];

    double dt_work_total   = 0.0 ;
    double dt_work_total_2 = 0.0 ;
    double dt_reduce_total    = 0.0 ;
    double dt_reduce_total_2  = 0.0 ;
    int i , k ;

    int result = TPI_Init( nth );

    if ( result != nth ) {
      fprintf(stderr,"%d != TPI_Init(%d) : FAILED\n", result , nth );
    }

    for ( i = 0 ; i < ntrial ; ++i ) {
      double t , dt ;
      int value = 0 ;

      for ( k = 0 ; k < nwork ; ++k ) { flags[k] = 0 ; }

      t = TPI_Walltime();
      TPI_Start( test_work , & flags , nwork , 0 );
      TPI_Wait();
      dt = TPI_Walltime() - t ;
      dt_work_total += dt ;
      dt_work_total_2 += dt * dt ;

      for ( k = 0 ; k < nwork && flags[k] ; ++k );

      if ( k < nwork ) {
        fprintf(stderr, "TPI_Run(...) : FAILED at trial %d\n", i );
        abort();
      }

      for ( k = 0 ; k < nwork ; ++k ) { flags[k] = 0 ; }

      t = TPI_Walltime();

      TPI_Start_reduce( test_reduce_work , & flags , nwork ,
                        test_reduce_join , test_reduce_init ,
                        sizeof(value) , & value );
      TPI_Wait();
  
      dt = TPI_Walltime() - t ;
      dt_reduce_total += dt ;
      dt_reduce_total_2 += dt * dt ;

      for ( k = 0 ; k < nwork && flags[k] ; ++k );

      if ( value != nwork || k < nwork ) {
        fprintf(stderr, "TPI_Run_reduce(...) : FAILED at trial %d\n", i );
        abort();
      }
    }

    TPI_Finalize();

    if ( 1 < ntrial ) {
      const double work_mean = 1.0e6 * dt_work_total / ntrial ;
      const double work_sdev = 1.0e6 * sqrt( ( ntrial * dt_work_total_2 -
                                       dt_work_total * dt_work_total ) /
                                     ( ntrial * ( ntrial - 1 ) ) );

      const double reduce_mean = 1.0e6 * dt_reduce_total / ntrial ;
      const double reduce_sdev = 1.0e6 * sqrt( ( ntrial * dt_reduce_total_2 -
                                         dt_reduce_total * dt_reduce_total) /
                                       ( ntrial * ( ntrial - 1 ) ) );
      
      fprintf(stdout,"%d , %d , %d , %10g , %10g , %10g , %10g\n",
              nth, ntrial, nwork, work_mean, work_sdev, reduce_mean, reduce_sdev);
    }
  }

  free( flags );
}

/*--------------------------------------------------------------------*/

static void test_work( TPI_Work * work )
{
  int * const flags = * (int *const*) work->info ;
  flags[ work->rank ] = 1 ;
}

static void test_reduce_work( TPI_Work * work )
{
  int * const flags = * (int *const*) work->info ;
  flags[ work->rank ] = 1 ;

  *((int *) work->reduce) += 1 ;
}

static void test_reduce_init( TPI_Work * work )
{
  *((int *) work->reduce) = 0 ;
}

static void test_reduce_join( TPI_Work * work , const void * src )
{
  *((int *) work->reduce) += *( (const int *) src );
}

static void test_reduce_via_lock( TPI_Work * work )
{
  int * const value = * ((int *const*) work->info );
  int result ;
  if ( ( result = TPI_Lock(0) ) ) {
    fprintf(stderr,"TPI_Lock(0) = %d : FAILED\n", result);
    abort();
  }
  *value += 1 ;
  if ( ( result = TPI_Unlock(0) ) ) {
    fprintf(stderr,"TPI_Unlock(0) = %d : FAILED\n", result);
    abort();
  }
}

static void test_reduce_via_nolock( TPI_Work * work )
{
  int * const value = (int *) work->reduce ;
  *value += 1 ;
}

/*--------------------------------------------------------------------*/

