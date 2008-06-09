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

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <TPI.h>

/*------------------------------------------------------------------------*/
/* Test various ways of controling worker threads */

typedef struct TestPthreads_struct {
  pthread_mutex_t  m_lock ;
  pthread_mutex_t  m_lock_run ;
  pthread_cond_t   m_cond ;
  pthread_cond_t   m_cond_run ;
  int              m_blocking ;
  int              m_thread_count ;
  int              m_work_size ;
} TestPthreads ;

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

#if 1

static void test_work( TestPthreads * const p , const int am_thread )
{
  const int blocking = am_thread && p->m_blocking ;

  int working = 1 ;

  while ( working ) {
    const int work_size = p->m_work_size ; /* may change */

    if ( ( 0 >  work_size ) ||
         ( 0 == work_size && ! am_thread ) ) {
      working = 0 ;
    }
    else if ( blocking || 0 < work_size ) {
      pthread_mutex_lock( & p->m_lock_run );

      if ( blocking ) {
        while ( ! p->m_work_size ) {
          pthread_cond_wait( & p->m_cond_run , & p->m_lock_run );
        }
      }

      if ( 0 < p->m_work_size ) { --( p->m_work_size ); }

      pthread_mutex_unlock( & p->m_lock_run );
    }
  }
}

#else

static void test_work( TestPthreads * const p , const int am_thread )
{
  const int blocking = am_thread && p->m_blocking ;

  int working = 1 ;

  while ( working ) {
    int work_size = p->m_work_size ; /* may change */

    if ( ( 0 >  work_size ) ||
         ( 0 == work_size && ! am_thread ) ) {
      working = 0 ;
    }
    else if ( ! blocking && 0 == work_size ) {
      /* No work to do and not blocking, just spin */ ;
    }
    else {
      pthread_mutex_lock( & p->m_lock_run );

      /* I have the lock */

      if ( 0 < p->m_work_size ) { --( p->m_work_size ); }

      work_size = p->m_work_size ;

      if ( blocking && 0 == work_size ) {
        pthread_cond_wait( & p->m_cond_run , & p->m_lock_run );
      }

      pthread_mutex_unlock( & p->m_lock_run );
    }
  }
}

#endif

static void * test_driver( void * arg )
{
  TestPthreads * const shared = (TestPthreads*) arg ;

  /*------------------------------*/
  /* Initializing */

  pthread_mutex_lock( & shared->m_lock );
  ++( shared->m_thread_count );
  pthread_mutex_unlock( & shared->m_lock );

  pthread_cond_signal(  & shared->m_cond );

  /*------------------------------*/

  test_work( shared , 1 );

  /*------------------------------*/
  /* Terminating */

  pthread_mutex_lock( & shared->m_lock );
  if ( 0 == --( shared->m_thread_count ) ) {
    pthread_cond_signal( & shared->m_cond );
  }
  pthread_mutex_unlock( & shared->m_lock );

  return NULL ;
}


static void test_run( pthread_attr_t * thread_attr ,
                      int blocking ,
                      int number_threads ,
                      int number_trials ,
                      int number_loops ,
                      int number_work ,
                      double * dt_start_stop ,
                      double * dt_loop )
{
  TestPthreads data = { PTHREAD_MUTEX_INITIALIZER ,
                        PTHREAD_MUTEX_INITIALIZER ,
                        PTHREAD_COND_INITIALIZER ,
                        PTHREAD_COND_INITIALIZER ,
                        0 , 0 , 0 };
  double dt_total ;
  double dt_run = 0 ;
  int j ;

  dt_total = TPI_Walltime();

  for ( j = 0 ; j < number_trials ; ++j ) {
    int i , k ;

    data.m_thread_count = 1 ;
    data.m_blocking = blocking ;

    pthread_mutex_lock( & data.m_lock );

    for ( i = 0 ; i < number_threads ; ++i ) {
      pthread_t pt ;
      pthread_create( & pt, thread_attr, & test_driver , & data );
      pthread_cond_wait( & data.m_cond , & data.m_lock );
    }

    /* Running */

    {
      double dt = TPI_Walltime();

      for ( k = 0 ; k < number_loops ; ++k ) {

        pthread_mutex_lock( & data.m_lock_run );
        data.m_work_size = number_work ;
        pthread_mutex_unlock( & data.m_lock_run );

        if ( data.m_blocking ) {
          pthread_cond_broadcast( & data.m_cond_run );
        }

        test_work( & data , 0 );
      }

      dt_run += TPI_Walltime() - dt ;
    }

    /* Termination */

    --( data.m_thread_count );

    if ( data.m_thread_count ) {
      pthread_mutex_lock( & data.m_lock_run );
      data.m_work_size = -1 ;
      if ( data.m_blocking ) {
        pthread_cond_broadcast( & data.m_cond_run );
      }
      pthread_mutex_unlock( & data.m_lock_run );
      pthread_cond_wait( & data.m_cond , & data.m_lock );
    }

    pthread_mutex_unlock( & data.m_lock );

    /* All threads have terminated */

    data.m_work_size = 0 ;
  }

  dt_total = TPI_Walltime() - dt_total ;

  *dt_loop       = 1.0e6 * dt_run / (double) ( number_trials * number_loops );
  *dt_start_stop = 1.0e6 * ( dt_total - dt_run ) / (double) number_trials ;

  pthread_mutex_destroy( & data.m_lock );
  pthread_mutex_destroy( & data.m_lock_run );
  pthread_cond_destroy(  & data.m_cond );
  pthread_cond_destroy(  & data.m_cond_run );
}

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

static double test_mutex_init_destroy( const int number )
{
  pthread_mutex_t mutex ;
  double dt ;
  int i ;
  dt = TPI_Walltime();
  for ( i = 0 ; i < number ; ++i ) {
    pthread_mutex_init( & mutex , NULL );
    pthread_mutex_destroy( & mutex );
  }
  dt = ( TPI_Walltime() - dt ) / (double) number ;
  return dt ;
}

static double test_mutex_lock_unlock( const int number )
{
  pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER ;
  double dt ;
  int i ;

  dt = TPI_Walltime();
  for ( i = 0 ; i < number ; ++i ) {
    pthread_mutex_lock( & mutex );
    pthread_mutex_unlock( & mutex );
  }
  dt = ( TPI_Walltime() - dt ) / (double) number ;

  pthread_mutex_destroy( & mutex );
  return dt ;
}

/*------------------------------------------------------------------------*/

void test_pthreads_performance( int n_test , int * n_concurrent )
{
  const int n_mutex = 1e4 /* 1e8 */ ;
  const int n_trial = 1e3 /* 1e4 */ ;
  const int n_loop  = 1e3 /* 1e4 */ ;
  const int n_work  = 1e2 ;

  {
    const double dt = 1e6 * test_mutex_init_destroy( n_mutex );
    fprintf(stdout,"\n\"test pthreads mutex init/destroy (microsec)\" , %g\n",dt);
  }

  {
    const double dt = 1e6 * test_mutex_lock_unlock( n_mutex );
    fprintf(stdout,"\n\"test pthreads mutex lock/unlock (microsec)\" , %g\n",dt);
  }

  /*------------------------------------------------------------------*/

  {
    int i ;

    fprintf(stdout,"\n\"test pthreads SCOPE_SYSTEM run-blocking\"\n");
    fprintf(stdout,"\"#Threads\" , \"#Spawned\" \"#WorkLoop\" \"Spawn (microsec)\" , \"Loop (microsec)\"\n");

    pthread_attr_t thread_attr ;

    pthread_attr_init( & thread_attr );
    pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
    pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

    for ( i = 0 ; i < n_test ; ++i ) {
      const int blocking = 1 ;
      const int nthread = n_concurrent[i] - 1 ;
      double dt_start_stop , dt_loop ;

      test_run( & thread_attr, blocking, nthread, n_trial, n_loop, n_work,
                & dt_start_stop , & dt_loop );

      fprintf( stdout, "%d , %d , %d , %g , %g\n",
               n_concurrent[i] , nthread , n_work , dt_start_stop , dt_loop );
      fflush( stdout );
    }

    pthread_attr_destroy( & thread_attr );
  }

  /*------------------------------------------------------------------*/

  {
    int i ;

    fprintf(stdout,"\n\"test pthreads SCOPE_PROCESS run-blocking\"\n");
    fprintf(stdout,"\"#Threads\" , \"#Spawned\" \"#WorkLoop\" \"Spawn (microsec)\" , \"Loop (microsec)\"\n");

    pthread_attr_t thread_attr ;

    pthread_attr_init( & thread_attr );
    pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_PROCESS );
    pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

    for ( i = 0 ; i < n_test ; ++i ) {
      const int blocking = 1 ;
      const int nthread = n_concurrent[i] - 1 ;
      double dt_start_stop , dt_loop ;

      test_run( & thread_attr, blocking, nthread, n_trial, n_loop, n_work,
                & dt_start_stop , & dt_loop );

      fprintf( stdout, "%d , %d , %d , %g , %g\n",
               n_concurrent[i] , nthread , n_work , dt_start_stop , dt_loop );
      fflush( stdout );
    }

    pthread_attr_destroy( & thread_attr );
  }

  /*------------------------------------------------------------------*/

  {
    int i ;

    fprintf(stdout,"\n\"test pthreads SCOPE_SYSTEM run-spinning\"\n");
    fprintf(stdout,"\"#Threads\" , \"#Spawned\" \"#WorkLoop\" \"Spawn (microsec)\" , \"Loop (microsec)\"\n");

    pthread_attr_t thread_attr ;

    pthread_attr_init( & thread_attr );
    pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
    pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

    for ( i = 0 ; i < n_test ; ++i ) {
      const int blocking = 0 ;
      const int nthread = n_concurrent[i] - 1 ;
      double dt_start_stop , dt_loop ;

      test_run( & thread_attr, blocking, nthread, n_trial, n_loop, n_work,
                & dt_start_stop , & dt_loop );

      fprintf( stdout, "%d , %d , %d , %g , %g\n",
               n_concurrent[i] , nthread , n_work , dt_start_stop , dt_loop );
      fflush( stdout );
    }

    pthread_attr_destroy( & thread_attr );
  }

  /*------------------------------------------------------------------*/

  {
    int i ;

    fprintf(stdout,"\n\"test pthreads SCOPE_PROCESS run-spinning\"\n");
    fprintf(stdout,"\"#Threads\" , \"#Spawned\" \"#WorkLoop\" \"Spawn (microsec)\" , \"Loop (microsec)\"\n");

    pthread_attr_t thread_attr ;

    pthread_attr_init( & thread_attr );
    pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_PROCESS );
    pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

    for ( i = 0 ; i < n_test ; ++i ) {
      const int blocking = 0 ;
      const int nthread = n_concurrent[i] - 1 ;
      double dt_start_stop , dt_loop ;

      test_run( & thread_attr, blocking, nthread, n_trial, n_loop, n_work,
                & dt_start_stop , & dt_loop );

      fprintf( stdout, "%d , %d , %d , %g , %g\n",
               n_concurrent[i] , nthread , n_work , dt_start_stop , dt_loop );
      fflush( stdout );
    }

    pthread_attr_destroy( & thread_attr );
  }

  /*------------------------------------------------------------------*/

  fflush( stdout );
}

/*------------------------------------------------------------------------*/


