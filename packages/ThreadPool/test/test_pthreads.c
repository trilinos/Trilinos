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
  pthread_mutex_t              m_lock ;
  pthread_mutex_t              m_lock_b ;
  pthread_cond_t               m_cond ;
  pthread_cond_t               m_cond_b ;
  struct TestPthreads_struct * m_next ;
  int                          m_work_flag ;
} TestPthreads ;

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

static void * test_spawn_driver_spin( void * arg )
{
  TestPthreads local = { PTHREAD_MUTEX_INITIALIZER ,
                         PTHREAD_MUTEX_INITIALIZER ,
                         PTHREAD_COND_INITIALIZER ,
                         PTHREAD_COND_INITIALIZER ,
                         NULL , 0 };

  TestPthreads * const shared = (TestPthreads*) arg ;

  /*------------------------------*/
  /* Initializing */

  /* Insert myself into the pool */

  pthread_mutex_lock( & shared->m_lock );

  local.m_next = shared->m_next ;

  shared->m_next = & local ;

  /* Signal that I've initialized */
  pthread_cond_signal(  & shared->m_cond );
  pthread_mutex_unlock( & shared->m_lock );

  /* Block until all peers are released to work */
  pthread_mutex_lock(   & shared->m_lock_b );
  pthread_mutex_unlock( & shared->m_lock_b );

  /*------------------------------*/

  {
    volatile int * const work_flag = & local.m_work_flag ;
    int flag ;

    while ( 0 <= ( flag = *work_flag ) ) {

      if ( 0 < flag ) {
        pthread_mutex_lock( & local.m_lock );

        /* My own parallel work goes here. */


        /* I completed my work */

        *work_flag = 0 ;

        pthread_mutex_unlock( & local.m_lock );
      }
    }
  }

  /*------------------------------*/
  /* Terminating */

  /* Remove myself from the pool */

  pthread_mutex_lock( & shared->m_lock );

  if ( shared->m_next == & local ) {
    shared->m_next = local.m_next ;
  }
  else {
    TestPthreads * p = shared->m_next ;
    while ( NULL != p && p->m_next != & local ) { p = p->m_next ; }
    if ( NULL == p ) {
      fprintf(stderr,"test_spawn_driver_cond TERMINATION FAILED\n");
      abort();
    }
    p->m_next = local.m_next ;
  }

  {
    const int finished = NULL == shared->m_next ;

    pthread_mutex_unlock( & shared->m_lock );

    if ( finished ) { pthread_cond_signal( & shared->m_cond ); }
  }

  /*------------------------------*/

  /* Clean up */

  pthread_mutex_destroy( & local.m_lock );
  pthread_mutex_destroy( & local.m_lock_b );
  pthread_cond_destroy(  & local.m_cond );
  pthread_cond_destroy(  & local.m_cond_b );

  return NULL ;
}


static double test_spawn_spin( int number_threads ,
                               int number_loops ,
                               int number_trials )
{
  TestPthreads data = { PTHREAD_MUTEX_INITIALIZER ,
                        PTHREAD_MUTEX_INITIALIZER ,
                        PTHREAD_COND_INITIALIZER ,
                        PTHREAD_COND_INITIALIZER ,
                        NULL , 0 };
  pthread_attr_t thread_attr ;
  double dt ;
  int j ;

  pthread_attr_init( & thread_attr );

  pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
  pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

  dt = TPI_Walltime();

  for ( j = 0 ; j < number_trials ; ++j ) {
    TestPthreads * p ;
    int i , k ;

    /* Initialization */

    pthread_mutex_lock( & data.m_lock_b );

    pthread_mutex_lock( & data.m_lock );
    for ( i = 0 ; i < number_threads ; ++i ) {
      pthread_t pt ;
      pthread_create( & pt, & thread_attr, & test_spawn_driver_spin , & data );
      pthread_cond_wait( & data.m_cond , & data.m_lock );
    }
    pthread_mutex_unlock( & data.m_lock );

    pthread_mutex_unlock( & data.m_lock_b ); /* Release */

    /* Running */

    for ( k = 0 ; k < number_loops ; ++k ) {

      for ( p = data.m_next ; p ; p = p->m_next ) {
        p->m_work_flag = 1 ; /* Can start running */
      }

      /* Parallel work happens here */


      /* Work is done, make sure all threads are stopped */

      for ( p = data.m_next ; p ; p = p->m_next ) {
        volatile int * const flag = & p->m_work_flag ;
        if ( *flag ) {
          pthread_mutex_lock( & p->m_lock );
          *flag = 0 ;
          pthread_mutex_unlock( & p->m_lock );
        }
      }
    }

    /* Termination */

    if ( NULL != data.m_next ) {
      pthread_mutex_lock( & data.m_lock );

      for ( p = data.m_next ; p ; p = p->m_next ) {
        p->m_work_flag = -1 ; /* Can start running to terminate */
      }

      pthread_cond_wait( & data.m_cond , & data.m_lock );

      pthread_mutex_unlock( & data.m_lock );
    }
  }

  dt = TPI_Walltime() - dt ;

  pthread_mutex_destroy( & data.m_lock );
  pthread_mutex_destroy( & data.m_lock_b );
  pthread_cond_destroy(  & data.m_cond );
  pthread_cond_destroy(  & data.m_cond_b );

  return dt ;
}

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

static void * test_spawn_driver_cond( void * arg )
{
  TestPthreads local = { PTHREAD_MUTEX_INITIALIZER ,
                         PTHREAD_MUTEX_INITIALIZER ,
                         PTHREAD_COND_INITIALIZER ,
                         PTHREAD_COND_INITIALIZER ,
                         NULL , 0 };

  TestPthreads * const shared = (TestPthreads*) arg ;

  /*------------------------------*/
  /* Initializing */

  /* Insert myself into the pool */

  local.m_next = shared->m_next ;

  shared->m_next = & local ;

  /*------------------------------*/
  /* Running */

  pthread_mutex_lock( & local.m_lock );

  /* Signal that I have started and am ready to run */

  pthread_mutex_lock(   & shared->m_lock );
  pthread_cond_signal(  & shared->m_cond );
  pthread_mutex_unlock( & shared->m_lock );

  do {

    pthread_cond_wait( & local.m_cond , & local.m_lock );

    if ( local.m_work_flag ) {
      /* Pass the non-zero work signal along to the next thread */

      if ( NULL != local.m_next ) {
        pthread_cond_signal( & local.m_next->m_cond );
      }

      if ( 0 < local.m_work_flag ) {

        /* My own parallel work goes here. */


        /* I completed my work */

        local.m_work_flag = 0 ;
      }
    }
  } while ( 0 <= local.m_work_flag );

  pthread_mutex_unlock( & local.m_lock );

  /*------------------------------*/
  /* Terminating */

  /* Remove myself from the pool */

  pthread_mutex_lock( & shared->m_lock );

  if ( shared->m_next == & local ) {
    shared->m_next = local.m_next ;
  }
  else {
    TestPthreads * p = shared->m_next ;
    while ( NULL != p && p->m_next != & local ) { p = p->m_next ; }
    if ( NULL == p ) {
      fprintf(stderr,"test_spawn_driver_cond TERMINATION FAILED\n");
      abort();
    }
    p->m_next = local.m_next ;
  }

  {
    const int finished = NULL == shared->m_next ;

    pthread_mutex_unlock( & shared->m_lock );

    if ( finished ) { pthread_cond_signal( & shared->m_cond ); }
  }

  /*------------------------------*/

  /* Clean up */

  pthread_mutex_destroy( & local.m_lock );
  pthread_mutex_destroy( & local.m_lock_b );
  pthread_cond_destroy(  & local.m_cond );
  pthread_cond_destroy(  & local.m_cond_b );

  return NULL ;
}


static double test_spawn_cond( int number_threads ,
                               int number_loops ,
                               int number_trials )
{
  TestPthreads data = { PTHREAD_MUTEX_INITIALIZER ,
                        PTHREAD_MUTEX_INITIALIZER ,
                        PTHREAD_COND_INITIALIZER ,
                        PTHREAD_COND_INITIALIZER ,
                        NULL , 0 };
  pthread_attr_t thread_attr ;
  double dt ;
  int j ;

  pthread_attr_init( & thread_attr );

  pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
  pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

  dt = TPI_Walltime();

  for ( j = 0 ; j < number_trials ; ++j ) {
    TestPthreads * p ;
    int i , k ;

    /* Initialization */

    pthread_mutex_lock( & data.m_lock );
    for ( i = 0 ; i < number_threads ; ++i ) {
      pthread_t pt ;
      pthread_create( & pt, & thread_attr, & test_spawn_driver_cond , & data );
      pthread_cond_wait( & data.m_cond , & data.m_lock );
    }
    pthread_mutex_unlock( & data.m_lock );

    /* Running */

    for ( k = 0 ; k < number_loops ; ++k ) {

      for ( p = data.m_next ; p ; p = p->m_next ) {
        pthread_mutex_lock( & p->m_lock );
        p->m_work_flag = 1 ;
        pthread_mutex_unlock( & p->m_lock );
      }

      if ( NULL != data.m_next ) {
        pthread_cond_signal( & data.m_next->m_cond );
      }

      /* Parallel work happens here */


      /* Work is done, make sure all threads are stopped */

      for ( p = data.m_next ; p ; p = p->m_next ) {
        if ( p->m_work_flag ) {
          pthread_mutex_lock( & p->m_lock );
          p->m_work_flag = 0 ;
          pthread_mutex_unlock( & p->m_lock );
        }
      }
    }

    /* Termination */

    for ( p = data.m_next ; p ; p = p->m_next ) {
      pthread_mutex_lock( & p->m_lock );
      p->m_work_flag = -1 ;
      pthread_mutex_unlock( & p->m_lock );
    }

    if ( NULL != data.m_next ) {
      pthread_mutex_lock( & data.m_lock );

      pthread_cond_signal( & data.m_next->m_cond );

      pthread_cond_wait( & data.m_cond , & data.m_lock );

      pthread_mutex_unlock( & data.m_lock );
    }
  }

  dt = TPI_Walltime() - dt ;

  pthread_mutex_destroy( & data.m_lock );
  pthread_mutex_destroy( & data.m_lock_b );
  pthread_cond_destroy(  & data.m_cond );
  pthread_cond_destroy(  & data.m_cond_b );

  return dt ;
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

typedef struct TestCondSignal_Struct {
  pthread_mutex_t mutex ;
  pthread_cond_t  cond  ;
} TestCondSignal ;

static void * test_cond_signal_driver( void * arg )
{
  TestCondSignal * const data = (TestCondSignal *) arg ;

  pthread_cond_signal(  & data->cond ); /* Start */

  pthread_mutex_lock(   & data->mutex );
  pthread_cond_wait(    & data->cond , & data->mutex );
  pthread_mutex_unlock( & data->mutex );

  pthread_cond_signal(  & data->cond ); /* End */

  return NULL ;
}

static double test_cond_signal( const int number )
{
  TestCondSignal data = { PTHREAD_MUTEX_INITIALIZER ,
                          PTHREAD_COND_INITIALIZER };
  pthread_attr_t thread_attr ;
  pthread_t pt ;

  double dt ;
  int i ;

  pthread_attr_init( & thread_attr );

  pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
  pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

  pthread_mutex_lock( & data.mutex );

  pthread_create( & pt, & thread_attr, & test_cond_signal_driver , & data );

  pthread_cond_wait( & data.cond , & data.mutex );

  dt = TPI_Walltime();
  for ( i = 0 ; i < number ; ++i ) {
    /* Signal but do not allow it to run */
    pthread_cond_signal( & data.cond );
  }
  dt = ( TPI_Walltime() - dt ) / (double) number ;

  pthread_mutex_unlock( & data.mutex );

  pthread_cond_wait( & data.cond , & data.mutex );

  pthread_mutex_destroy( & data.mutex );
  pthread_cond_destroy(  & data.cond );
  return dt ;
}

/*------------------------------------------------------------------------*/

void test_pthreads_performance( int n_test , int * n_concurrent )
{
  const int n_mutex = 1e8 ;
  const int n_trial = 1e4 ;
  const int n_loop  = 1e4 ;

  {
    const double dt = 1e6 * test_mutex_init_destroy( n_mutex );
    fprintf(stdout,"\n\"test pthreads mutex init/destroy (microsec)\" , %g\n",dt);
  }

  {
    const double dt = 1e6 * test_mutex_lock_unlock( n_mutex );
    fprintf(stdout,"\n\"test pthreads mutex lock/unlock (microsec)\" , %g\n",dt);
  }

  {
    const double dt = 1e6 * test_cond_signal( n_mutex );
    fprintf(stdout,"\n\"test pthreads cond_signal (microsec)\" , %g\n",dt);
  }

  {
    int i ;

    fprintf(stdout,"\n\"test pthreads spawn and signal-loop\"\n");
    fprintf(stdout,"\"#Threads\" , \"#Spawned\" \"Spawn (microsec)\" , \"Loop (microsec)\"\n");

    for ( i = 0 ; i < n_test ; ++i ) {
      const int nthread = n_concurrent[i] - 1 ;
      double dt_spawn = test_spawn_cond( nthread , 0 , n_trial );
      double dt_loop  = test_spawn_cond( nthread , n_loop , 1 );

      dt_spawn /= (double) n_trial ;
      dt_loop  -= dt_spawn ;

      dt_spawn *= 1.0e6 ;
      dt_loop  *= 1.0e6 / (double) n_loop ;

      fprintf(stdout,"%d , %d , %g , %g\n", n_concurrent[i] , nthread , dt_spawn , dt_loop );
      fflush( stdout );
    }
  }

  {
    int i ;

    fprintf(stdout,"\n\"test pthreads spawn and spin-loop\"\n");
    fprintf(stdout,"\"#Threads\" , \"#Spawned\" \"Spawn (microsec)\" , \"Loop (microsec)\"\n");

    for ( i = 0 ; i < n_test ; ++i ) {
      const int nthread = n_concurrent[i] - 1 ;
      double dt_spawn = test_spawn_spin( nthread , 0 , n_trial );
      double dt_loop  = test_spawn_spin( nthread , n_loop , 1 );

      dt_spawn /= (double) n_trial ;
      dt_loop  -= dt_spawn ;

      dt_spawn *= 1.0e6 ;
      dt_loop  *= 1.0e6 / (double) n_loop ;

      fprintf(stdout,"%d , %d , %g , %g\n", n_concurrent[i] , nthread , dt_spawn , dt_loop );
      fflush( stdout );
    }
  }

  fflush( stdout );
}

/*------------------------------------------------------------------------*/


