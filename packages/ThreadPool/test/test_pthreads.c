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
#include <pthread.h>
#include <TPI.h>

/*------------------------------------------------------------------------*/

typedef struct TestPthreads_struct {
  pthread_mutex_t              m_lock ;
  pthread_cond_t               m_cond ;
  struct TestPthreads_struct * m_next ;
} TestPthreads ;

static void * test_spawn_driver( void * arg )
{
  TestPthreads local = { PTHREAD_MUTEX_INITIALIZER ,
                         PTHREAD_COND_INITIALIZER ,
                         NULL };

  TestPthreads * const shared = (TestPthreads*) arg ;

  pthread_mutex_lock( & local.m_lock );

  local.m_next = shared->m_next ;

  shared->m_next = & local ;

  /* Signal that I have started */

  pthread_mutex_lock(   & shared->m_lock );
  pthread_cond_signal(  & shared->m_cond );
  pthread_mutex_unlock( & shared->m_lock );

  /* Wait for my signal to terminate */

  do {
    pthread_cond_wait( & local.m_cond , & local.m_lock );
  } while ( shared->m_next );

  /* Signal that I am terminating */

  pthread_mutex_lock(   & shared->m_lock );
  pthread_cond_signal(  & shared->m_cond );
  pthread_mutex_unlock( & shared->m_lock );

  /* Clean up */

  pthread_mutex_unlock(  & local.m_lock );
  pthread_mutex_destroy( & local.m_lock );
  pthread_cond_destroy(  & local.m_cond );

  return NULL ;
}

static double test_spawn( int number_threads ,
                          int number_loops ,
                          int number_trials )
{
  TestPthreads data = { PTHREAD_MUTEX_INITIALIZER ,
                        PTHREAD_COND_INITIALIZER ,
                        NULL };
  pthread_attr_t thread_attr ;
  double dt ;
  int j ;

  pthread_attr_init( & thread_attr );

  pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
  pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

  pthread_mutex_lock( & data.m_lock );

  dt = TPI_Walltime();

  for ( j = 0 ; j < number_trials ; ++j ) {
    int i , k ;
    for ( i = 0 ; i < number_threads ; ++i ) {
      pthread_t pt ;
      pthread_create( & pt, & thread_attr, & test_spawn_driver , & data );
      pthread_cond_wait( & data.m_cond , & data.m_lock );
      pthread_mutex_lock( & data.m_next->m_lock );
    }

    for ( k = 0 ; k < number_loops ; ++k ) {
      TestPthreads * p ;
      for ( p = data.m_next ; p ; p = p->m_next ) {
        pthread_cond_signal(  & p->m_cond );
        pthread_mutex_unlock( & p->m_lock );
      }

      /* Parallel work happens here */

      for ( p = data.m_next ; p ; p = p->m_next ) {
        pthread_mutex_lock( & p->m_lock );
      }
    }

    {
      TestPthreads * p_next = data.m_next ; data.m_next = NULL ;
      while ( p_next ) {
        TestPthreads * const p = p_next ;
        p_next = p->m_next ;
        pthread_cond_signal(  & p->m_cond );
        pthread_mutex_unlock( & p->m_lock );
        pthread_cond_wait( & data.m_cond , & data.m_lock );
      }
    }
  }

  dt = TPI_Walltime() - dt ;

  pthread_mutex_unlock(  & data.m_lock );
  pthread_mutex_destroy( & data.m_lock );
  pthread_cond_destroy(  & data.m_cond );

  return dt ;
}

/*------------------------------------------------------------------------*/

static double test_mutex_init_destroy( int number_mutex )
{
  pthread_mutex_t mutex ;
  double dt ;
  int i ;
  dt = TPI_Walltime();
  for ( i = 0 ; i < number_mutex ; ++i ) {
    pthread_mutex_init( & mutex , NULL );
    pthread_mutex_destroy( & mutex );
  }
  dt = TPI_Walltime() - dt ;
  return dt ;
}

/*------------------------------------------------------------------------*/

void test_pthreads_performance( int n_test , int * n_concurrent )
{
  const int n_mutex = 1e7 ;
  const int n_spawn = 1e5 ;
  const int n_loop  = 1e2 ;

  {
    double dt_mutex = 1e6 * test_mutex_init_destroy( n_mutex ) / n_mutex ;
    fprintf(stdout,"\n\"test pthreads mutex (microsec)\" , %g\n",dt_mutex);
  }

  {
    int i ;

    fprintf(stdout,"\n\"test pthreads spawn and loop\"\n");
    fprintf(stdout,"\"#Threads\" , \"#Spawned\" \"Spawn (microsec)\" , \"Loop (microsec)\"\n");

    for ( i = 0 ; i < n_test ; ++i ) {
      const int nthread = n_concurrent[i] - 1 ;
      const int ntrial  = n_spawn / nthread ;
      double dt_spawn = test_spawn( nthread , 0 , ntrial );
      double dt_loop  = test_spawn( nthread , n_loop , ntrial ) - dt_spawn ;

      dt_spawn *= 1.0e6 / ntrial ;
      dt_loop  *= 1.0e6 / ( ntrial * n_loop );

      fprintf(stdout,"%d , %d , %g , %g\n", n_concurrent[i] , nthread , dt_spawn , dt_loop );
      fflush( stdout );
    }
  }
  fflush( stdout );
}

/*------------------------------------------------------------------------*/


