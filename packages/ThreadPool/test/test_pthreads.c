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
  pthread_cond_t   m_cond ;
  int              m_thread_rank ;
  int              m_thread_count ;
} TestPthreads ;

/*------------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/

static void * test_driver( void * arg )
{
  TestPthreads * const data = (TestPthreads*) arg ;
  TestPthreads * const root = data - data->m_thread_rank ;

  /*------------------------------*/
  /* Initializing */

  pthread_mutex_lock(   & data->m_lock );

  pthread_mutex_lock(   & root->m_lock );
  pthread_cond_signal(  & root->m_cond );
  pthread_mutex_unlock( & root->m_lock );

  /*------------------------------*/

  while ( data->m_thread_rank ) { 
    pthread_cond_wait( & data->m_cond , & data->m_lock );
  } 
  pthread_mutex_unlock( & data->m_lock );

  /*------------------------------*/
  /* Terminating */

  pthread_mutex_lock( & root->m_lock );
  if ( 0 == --( root->m_thread_count ) ) {
    pthread_cond_signal( & root->m_cond );
  }
  pthread_mutex_unlock( & root->m_lock );

  return NULL ;
}


static void test_run( pthread_attr_t * const thread_attr ,
                      const int number_threads ,
                      const int number_trials ,
                      const int number_loops ,
                      double * const dt_start_stop ,
                      double * const dt_loop )
{
  TestPthreads data[ number_threads ];
  double dt_total ;
  double dt_run = 0 ;
  int j ;

  dt_total = TPI_Walltime();

  for ( j = 0 ; j < number_trials ; ++j ) {
    int i ;

    for ( i = 0 ; i < number_threads ; ++i ) {
      pthread_cond_init( & data[i].m_cond , NULL );
      pthread_mutex_init( & data[i].m_lock , NULL );
      data[i].m_thread_rank = i ;
      data[i].m_thread_count = number_threads ;
    }

    pthread_mutex_lock( & data->m_lock );

    for ( i = 1 ; i < number_threads ; ++i ) {
      pthread_t pt ;
      pthread_create( & pt, thread_attr, & test_driver , data + i );
      pthread_cond_wait( & data->m_cond , & data->m_lock );
      pthread_mutex_lock( & data[i].m_lock );
    }

    /* Running */

    {
      double dt = TPI_Walltime();
      int k ;

      for ( k = 1 ; k < number_loops ; ++k ) {
        for ( i = 1 ; i < number_threads ; ++i ) {
          pthread_cond_signal(  & data[i].m_cond );
          pthread_mutex_unlock( & data[i].m_lock );
        }

        /* Work goes here */

        for ( i = 1 ; i < number_threads ; ++i ) {
          pthread_mutex_lock( & data[i].m_lock );
        }
      }

      dt_run += TPI_Walltime() - dt ;
    }

    /* Termination */

    --( data->m_thread_count );

    if ( data->m_thread_count ) {
      for ( i = 1 ; i < number_threads ; ++i ) {
        data[i].m_thread_rank = 0 ;
        pthread_cond_signal(  & data[i].m_cond );
        pthread_mutex_unlock( & data[i].m_lock );
      }

      pthread_cond_wait( & data->m_cond , & data->m_lock );
    }

    pthread_mutex_unlock( & data->m_lock );

    for ( i = 0 ; i < number_threads ; ++i ) {
      pthread_cond_destroy( & data[i].m_cond );
      pthread_mutex_destroy( & data[i].m_lock );
    }
  }

  dt_total = TPI_Walltime() - dt_total ;

  *dt_loop       = 1.0e6 * dt_run / (double) ( number_trials * number_loops );
  *dt_start_stop = 1.0e6 * ( dt_total - dt_run ) / (double) number_trials ;
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
  const int n_trial = 1e2 /* 1e4 */ ;
  const int n_loop  = 1e3 /* 1e4 */ ;

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
    fprintf(stdout,"\"#Threads\" , \"#Spawned\" \"Spawn (microsec)\" , \"Loop (microsec)\"\n");

    pthread_attr_t thread_attr ;

    pthread_attr_init( & thread_attr );
    pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
    pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

    for ( i = 0 ; i < n_test ; ++i ) {
      const int nthread = n_concurrent[i] ;
      double dt_start_stop , dt_loop ;

      test_run( & thread_attr, nthread, n_trial, n_loop,
                & dt_start_stop , & dt_loop );

      fprintf( stdout, "%d , %d , %g , %g\n",
               nthread , nthread - 1 , dt_start_stop , dt_loop );
      fflush( stdout );
    }

    pthread_attr_destroy( & thread_attr );
  }

  /*------------------------------------------------------------------*/

  {
    int i ;

    fprintf(stdout,"\n\"test pthreads SCOPE_PROCESS run-blocking\"\n");
    fprintf(stdout,"\"#Threads\" , \"#Spawned\" \"Spawn (microsec)\" , \"Loop (microsec)\"\n");

    pthread_attr_t thread_attr ;

    pthread_attr_init( & thread_attr );
    pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_PROCESS );
    pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

    for ( i = 0 ; i < n_test ; ++i ) {
      const int nthread = n_concurrent[i] ;
      double dt_start_stop , dt_loop ;

      test_run( & thread_attr, nthread, n_trial, n_loop,
                & dt_start_stop , & dt_loop );

      fprintf( stdout, "%d , %d , %g , %g\n",
               nthread , nthread - 1 , dt_start_stop , dt_loop );
      fflush( stdout );
    }

    pthread_attr_destroy( & thread_attr );
  }

  /*------------------------------------------------------------------*/

  fflush( stdout );
}

/*------------------------------------------------------------------------*/


