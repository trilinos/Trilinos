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

/* #define DEBUG_PRINT */

#include <stdio.h>
#include <stdlib.h>

#include <TPI.h>
#include <ThreadPool_config.h>

enum { THREAD_COUNT_MAX = 1024 };
enum { LOCK_COUNT_MAX   = 256 };

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#ifdef HAVE_PTHREAD

#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include <sys/types.h>

struct ThreadPool_Data ;

typedef void (*Thread_Work)( struct ThreadPool_Data * const , unsigned int );

typedef struct Thread_Data {
  pthread_cond_t  m_cond ;
  pthread_mutex_t m_lock ;
  unsigned int    m_pending ;
} Thread ;

typedef struct ThreadPool_Data {
  pthread_mutex_t       m_work_lock ;
  pthread_cond_t        m_work_cond ;
  TPI_work_subprogram   m_work_routine ;
  void                * m_work_argument ;
  int                   m_work_count_total ;
  int                   m_work_count_claim ;
  int                   m_work_count_pending ;

  Thread_Work           m_thread_driver ;
  int                   m_thread_count ;
  int                   m_thread_upper ;
  int                   m_lock_init ;
  int                   m_lock_count ;

  Thread                m_thread[ THREAD_COUNT_MAX ];
  pthread_mutex_t       m_lock[ LOCK_COUNT_MAX ];
} ThreadPool ;

/*--------------------------------------------------------------------*/

static ThreadPool * local_thread_pool()
{
  static ThreadPool thread_pool = {
    /* m_work_lock           */  PTHREAD_MUTEX_INITIALIZER ,
    /* m_work_cond           */  PTHREAD_COND_INITIALIZER ,
    /* m_work_routine        */  NULL ,
    /* m_work_argument       */  NULL ,
    /* m_work_count_total    */  0 ,
    /* m_work_count_claim    */  0 ,
    /* m_work_count_pending  */  0 ,
    /* m_thread_driver       */  NULL ,
    /* m_thread_count        */  0 ,
    /* m_thread_upper        */  0 ,
    /* m_lock_init           */  0 ,
    /* m_lock_count          */  0 };

  return & thread_pool ;
}

/*--------------------------------------------------------------------*/

int TPI_Lock( int i )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result = i < 0 || thread_pool->m_lock_count <= i ? TPI_ERROR_SIZE : 0 ;

  if ( ! result && pthread_mutex_lock( thread_pool->m_lock + i ) ) {
    result = TPI_ERROR_LOCK ;
  }

  return result ;
}

int TPI_Trylock( int i )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result = i < 0 || thread_pool->m_lock_count <= i ? TPI_ERROR_SIZE : 0 ;

  if ( ! result ) {
    result = pthread_mutex_trylock( thread_pool->m_lock + i );

    if ( EBUSY == result ) { result = TPI_LOCK_BUSY ; }
    else if ( result )     { result = TPI_ERROR_LOCK ; }
  }
  return result ;
}

int TPI_Unlock( int i )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result = i < 0 || thread_pool->m_lock_count <= i ? TPI_ERROR_SIZE : 0 ;

  if ( ! result && pthread_mutex_unlock( thread_pool->m_lock + i ) ) {
    result = TPI_ERROR_LOCK ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*  Run the work queue until it is empty.
 *  The input 'thread_number' is deliberately ignored.
 */
static void local_run_work( ThreadPool * const thread_pool ,
                            unsigned int work_item )
{
  struct TPI_Work_Struct work ;

  work.shared     = thread_pool->m_work_argument ;
  work.lock_count = thread_pool->m_lock_count ;
  work.work_count = thread_pool->m_work_count_total ;

  work_item = 0 ;

  do {

    if ( work_item ) { /* Have a work item to do */
      work.work_rank = work.work_count - work_item ;

      (* thread_pool->m_work_routine)( & work );
    }

    /* Finished my work item, update the work queue */

    pthread_mutex_lock( & thread_pool->m_work_lock );

    /* Record work just performed, if any */
    if ( work_item ) { --( thread_pool->m_work_count_pending ); }

    /* Claim more work, if any */
    work_item = thread_pool->m_work_count_claim ;

    if ( work_item ) { thread_pool->m_work_count_claim = work_item - 1 ; }

    pthread_mutex_unlock( & thread_pool->m_work_lock );

  } while ( work_item );
}

/*--------------------------------------------------------------------*/

static
void local_broadcast( ThreadPool * const control ,
                      const unsigned int thread_rank ,
                      unsigned int thread_cycle )
{
  Thread * const     thread_data  = control->m_thread ;
  const unsigned int thread_count = control->m_thread_count ;
 
  for ( ; thread_cycle < thread_count ; thread_cycle <<= 1 ) {
    const unsigned int next_rank = thread_rank + thread_cycle ;
    if ( next_rank < thread_count ) {
      Thread * const next_thread = thread_data + next_rank ;
 
#ifdef DEBUG_PRINT
  fprintf(stdout,"  broadcast run %u to %u\n",thread_rank,next_rank);
  fflush(stdout);
#endif

      next_thread->m_pending = thread_data->m_pending ;
 
      pthread_mutex_unlock( & next_thread->m_lock );
      pthread_cond_signal(  & next_thread->m_cond );
    }
  }
}
 
static
void local_barrier( ThreadPool * const control ,
                    const unsigned int thread_rank )
{
  Thread * const     thread_data  = control->m_thread ;
  const unsigned int thread_count = control->m_thread_count ;
        unsigned int thread_cycle = control->m_thread_upper ;
 
  for ( ; thread_rank < thread_cycle ; thread_cycle >>= 1 ) {
    const unsigned next_rank = thread_rank + thread_cycle ;
    if ( next_rank < thread_count ) {
      Thread * const next_thread = thread_data + next_rank ;
 
#ifdef DEBUG_PRINT
  fprintf(stdout,"  barrier halt %u from %u lock\n",thread_rank,next_rank);
  fflush(stdout);
#endif
 
      pthread_mutex_lock( & next_thread->m_lock );

      if ( next_thread->m_pending ) {
 
#ifdef DEBUG_PRINT
  fprintf(stdout,"  barrier halt %u from %u wait\n",thread_rank,next_rank);
  fflush(stdout);
#endif

        pthread_cond_wait( & next_thread->m_cond , & next_thread->m_lock );
      }

#ifdef DEBUG_PRINT
  fprintf(stdout,"  barrier halt %u from %u done\n",thread_rank,next_rank);
  fflush(stdout);
#endif
 
    }
  }

  if ( thread_data->m_pending ) {
    Thread * const my_thread = thread_data + thread_rank ;
    my_thread->m_pending = 0 ;
    pthread_cond_signal( & my_thread->m_cond );
  }
}
 
static unsigned int local_fan_cycle( const unsigned int i )
{
  unsigned j ;
  for ( j = 1 ; j <= i ; j <<= 1 );
  return j ;
}

/*--------------------------------------------------------------------*/
/*  The driver given to 'pthread_create'.
 *  Run work until told to terminate.
 */
static void * local_driver( void * arg )
{
  ThreadPool * const thread_pool = (ThreadPool *) arg ;
  Thread     * const root_thread = thread_pool->m_thread ;

  const unsigned my_rank   = --thread_pool->m_work_count_pending ;
  const unsigned my_cycle  = local_fan_cycle( my_rank );
  Thread * const my_thread = root_thread + my_rank ;

  /* Aquire my lock */
  pthread_mutex_lock( & my_thread->m_lock );

  /* Signal master thread that I have started */

  pthread_mutex_lock( & root_thread->m_lock );
  pthread_cond_signal(  & root_thread->m_cond );
  pthread_mutex_unlock( & root_thread->m_lock );

  local_barrier( thread_pool , my_rank );

  /*------------------------------*/

  do {
    /* Wait for the run signal, giving up my lock */
    pthread_cond_wait( & my_thread->m_cond , & my_thread->m_lock );

    /* Continue broadcasting the run signal */
    local_broadcast( thread_pool , my_rank , my_cycle );

    if ( thread_pool->m_thread_driver ) {

#ifdef DEBUG_PRINT
  fprintf(stdout,"  thread %u running\n",my_rank);
  fflush(stdout);
#endif

      /* Participate in the work */
      (*thread_pool->m_thread_driver)( thread_pool , my_rank );

      /* Resume waiting */
      local_barrier( thread_pool , my_rank );
    }
  } while ( thread_pool->m_thread_driver );

  /*------------------------------*/
  /* Termination, unlock my thread and count down */

  pthread_mutex_unlock( & my_thread->m_lock );

  pthread_mutex_lock(  & root_thread->m_lock );
  if ( ! --thread_pool->m_work_count_pending ) {
    pthread_cond_signal( & root_thread->m_cond );
  }
  pthread_mutex_unlock(  & root_thread->m_lock );

#ifdef DEBUG_PRINT
  fprintf( stdout , "  thread %u exiting\n",my_rank);
  fflush( stdout );
#endif

  return NULL ;
}

/*--------------------------------------------------------------------*/
/*  Run the work subprogram exactly once for each thread.
 *  The last thread signals the root thread.
 */
static void local_run_thread( ThreadPool * const thread_pool ,
                              unsigned int thread_rank )
{
  struct TPI_Work_Struct work ;

  work.shared     = thread_pool->m_work_argument ;
  work.lock_count = thread_pool->m_lock_count ;
  work.work_count = thread_pool->m_thread_count ;
  work.work_rank  = thread_rank ;

  (* thread_pool->m_work_routine)( & work );
}

/*--------------------------------------------------------------------*/

static int set_lock_count( ThreadPool * const thread_pool ,
                           const int lock_count )
{
  int result = lock_count < 0 || LOCK_COUNT_MAX < lock_count
             ? TPI_ERROR_SIZE : 0 ;

  while ( ! result && thread_pool->m_lock_init < lock_count ) {

    pthread_mutex_t * const lock = thread_pool->m_lock +
                                   thread_pool->m_lock_init ;

    if ( pthread_mutex_init( lock , NULL ) ) {
      result = TPI_ERROR_INTERNAL ;
    }
    else {
      ++( thread_pool->m_lock_init );
    }
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Run( TPI_work_subprogram subprogram  ,
             void * const        shared_data ,
             const int           work_count  ,
             const int           lock_count  )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result =
    NULL != thread_pool->m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == subprogram                  ? TPI_ERROR_NULL : (
    0 ) );

  if ( ! result ) {

    result = set_lock_count( thread_pool , lock_count );

    if ( ! result ) {

      thread_pool->m_lock_count         = lock_count ;
      thread_pool->m_work_argument      = shared_data ;
      thread_pool->m_work_routine       = subprogram ;
      thread_pool->m_work_count_total   = work_count ;
      thread_pool->m_work_count_pending = work_count ;
      thread_pool->m_work_count_claim   = work_count ;
      thread_pool->m_thread_driver      = & local_run_work ;

      if ( 1 < thread_pool->m_thread_count ) {
        /* Activate the blocked worker threads */
        local_broadcast( thread_pool , 0 , 1 );

        /* Participate in work queue */
        local_run_work( thread_pool , 0 );

        /* Block the worker threads */
        local_barrier( thread_pool , 0 );
      }
      else {
        /* No worker threads, can bypass work queue locking */
        struct TPI_Work_Struct work ;

        work.shared     = shared_data ;
        work.lock_count = lock_count ;
        work.work_count = work_count ;
        work.work_rank  = 0 ;

        for ( ; work.work_rank < work.work_count ; ++( work.work_rank ) ) {
          (* subprogram)( & work );
        }
      }

      thread_pool->m_thread_driver      = NULL ;
      thread_pool->m_work_count_claim   = 0 ;
      thread_pool->m_work_count_pending = 0 ;
      thread_pool->m_work_count_total   = 0 ;
      thread_pool->m_work_routine       = NULL ;
      thread_pool->m_work_argument      = NULL ;
      thread_pool->m_lock_count         = 0 ;
    }
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Run_threads( TPI_work_subprogram subprogram  ,
                     void * const        shared_data ,
                     const int           lock_count  )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result =
    NULL != thread_pool->m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == subprogram                  ? TPI_ERROR_NULL : (
    0 ) );

  if ( ! result ) {

    result = set_lock_count( thread_pool , lock_count );

    if ( ! result ) {

      thread_pool->m_lock_count    = lock_count ;
      thread_pool->m_work_argument = shared_data ;
      thread_pool->m_work_routine  = subprogram ;
      thread_pool->m_thread_driver = & local_run_thread ;

      if ( 1 < thread_pool->m_thread_count ) {
        thread_pool->m_thread->m_pending = 1 ;

        /* Activate the blocked worker threads */
        local_broadcast( thread_pool , 0 , 1 );

        /* Participate in work */
        local_run_thread( thread_pool , 0 );

        /* Block the worker threads */
        local_barrier( thread_pool , 0 );

        thread_pool->m_thread->m_pending = 1 ;
      }
      else {
        /* No worker threads */
        struct TPI_Work_Struct work ;

        work.shared     = shared_data ;
        work.lock_count = lock_count ;
        work.work_count = 1 ;
        work.work_rank  = 0 ;

        (* subprogram)( & work );
      }

      thread_pool->m_thread_driver = NULL ;
      thread_pool->m_work_routine  = NULL ;
      thread_pool->m_work_argument = NULL ;
      thread_pool->m_lock_count    = 0 ;
    }
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Init( int n )
{
  ThreadPool * const thread_pool = local_thread_pool();
  pthread_attr_t attr ;
  int result , i ;

  result = thread_pool->m_thread_count ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && ( n < 1 || THREAD_COUNT_MAX <= n ) ) {
    result = TPI_ERROR_SIZE ;
  }

  for ( i = 0 ; ! result && i < n ; ++i ) {
    thread_pool->m_thread[i].m_pending = 0 ;
    if ( pthread_mutex_init( & thread_pool->m_thread[i].m_lock , NULL ) ||
         pthread_cond_init(  & thread_pool->m_thread[i].m_cond , NULL ) ) {
      result = TPI_ERROR_INTERNAL ;
    }
  }

  if ( ! result &&
       ( pthread_attr_init( & attr ) ||
         pthread_attr_setscope(       & attr, PTHREAD_SCOPE_SYSTEM ) ||
         pthread_attr_setdetachstate( & attr, PTHREAD_CREATE_DETACHED ) ) ) {
    result = TPI_ERROR_INTERNAL ;
  }

  if ( ! result ) {

    pthread_mutex_lock( & thread_pool->m_thread->m_lock );

    for ( i = 1 ; i <= n ; i <<= 1 );

    thread_pool->m_thread_upper = i ;
    thread_pool->m_thread_count = n ;
    thread_pool->m_work_count_pending = n ;

    for ( i = 1 ; ! result && i < n ; ++i ) {
      pthread_t pt ;

      if ( pthread_create( & pt, & attr, & local_driver, thread_pool ) ) {
        result = TPI_ERROR_INTERNAL ;
      }
      else {
        /* Wait for worker thread to start */
        pthread_cond_wait( & thread_pool->m_thread->m_cond ,
                           & thread_pool->m_thread->m_lock );
      }
    }

    thread_pool->m_work_count_pending = 0 ;

    pthread_attr_destroy( & attr );

    if ( ! result ) {
      local_barrier( thread_pool , 0 );
    }
    else {
      TPI_Finalize();
    }
  }

  if ( ! result ) { result = n ; }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  ThreadPool * const thread_pool = local_thread_pool();
  int result , i ;

  result = NULL != thread_pool->m_thread_driver ? TPI_ERROR_ACTIVE : 0 ;

#ifdef DEBUG_PRINT
  fprintf(stdout,"TPI_Finalize\n");
  fflush(stdout);
#endif

  if ( ! result ) {

    thread_pool->m_work_count_pending = thread_pool->m_thread_count - 1 ;

    local_broadcast( thread_pool , 0 , 1 );

    if ( 0 < thread_pool->m_work_count_pending ) {
      /* Wait for signal from last worker thread */
      pthread_cond_wait( & thread_pool->m_thread->m_cond ,
                         & thread_pool->m_thread->m_lock );
    }

    pthread_mutex_unlock( & thread_pool->m_thread->m_lock );
 
    for ( i = 0 ; i < thread_pool->m_thread_count ; ++i ) {
      pthread_mutex_destroy( & thread_pool->m_thread[i].m_lock );
      pthread_cond_destroy(  & thread_pool->m_thread[i].m_cond );
    }
 
    thread_pool->m_thread_count = 0 ;
    thread_pool->m_thread_upper = 0 ;

    while ( thread_pool->m_lock_init ) {
      pthread_mutex_destroy( thread_pool->m_lock + thread_pool->m_lock_init );
      --( thread_pool->m_lock_init );
    }
  }

#ifdef DEBUG_PRINT
  fprintf(stdout,"TPI_Finalize DONE\n");
  fflush(stdout);
#endif

  return result ;
}

/* #ifdef HAVE_PTHREAD */
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#else

typedef struct LocalWorkData {
  int m_active ;
  int m_lock_count ;
  int m_lock[ LOCK_COUNT_MAX ];
} LocalWork ;

static LocalWork * local_work()
{
  static struct LocalWorkData data = { 0 , 0 };
  return & data ;
}

int TPI_Lock( int i )
{
  LocalWork * const work = local_work();

  const int result = 
    ( i < 0 || work->m_lock_count <= i ? TPI_ERROR_SIZE :
    ( work->m_lock[i] )                ? TPI_ERROR_LOCK : 0 );
  if ( ! result ) { work->m_lock[i] = 1 ; }
  return result ;
}

int TPI_Trylock( int i )
{
  LocalWork * const work = local_work();

  const int result = 
    ( i < 0 || work->m_lock_count <= i ? TPI_ERROR_SIZE :
    ( work->m_lock[i] )                ? TPI_LOCK_BUSY : 0 );
  if ( ! result ) { work->m_lock[i] = 1 ; }
  return result ;
}

int TPI_Unlock( int i )
{
  LocalWork * const work = local_work();

  const int result = 
    ( i < 0 || work->m_lock_count <= i ? TPI_ERROR_SIZE :
    ( ! work->m_lock[i] )              ? TPI_ERROR_LOCK : 0 );
  if ( ! result ) { work->m_lock[i] = 0 ; }

  return result ;
}

static int set_lock_count( LocalWork * const work , const int lock_count )
{
  const int result = lock_count < 0 || LOCK_COUNT_MAX < lock_count
                   ? TPI_ERROR_SIZE : 0 ;

  if ( ! result ) {
    int i ;
    for ( i = 0 ; i < lock_count ; ++i ) { work->m_lock[i] = 0 ; }
  }

  return result ; 
}

int TPI_Run( TPI_work_subprogram subprogram  ,
             void * const        shared_data ,
             const int           work_count  ,
             const int           lock_count  )
{
  LocalWork * const work = local_work();

  int result = work->m_active ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {

    result = set_lock_count( work , lock_count );

    if ( ! result ) {
      struct TPI_Work_Struct w ;

      work->m_active     = 1 ;
      work->m_lock_count = lock_count ;

      w.shared     = shared_data ;
      w.lock_count = lock_count ;
      w.work_count = work_count ;
      w.work_rank  = 0 ;

      for ( w.work_rank = 0 ; w.work_rank < work_count ; ++(w.work_rank) ) {
        (*subprogram)( & w );
      }

      work->m_active     = 0 ;
      work->m_lock_count = 0 ;
    }
  }

  return result ;
}

int TPI_Run_threads( TPI_work_subprogram subprogram  ,
                     void * const        shared_data ,
                     const int           lock_count  )
{
  return TPI_Run( subprogram , shared_data , 1 , lock_count );
}

int TPI_Init( int )
{
  LocalWork * const work = local_work();

  int result = work->m_active ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) { result = 1 ; }

  return result ;
}

int TPI_Finalize()
{
  LocalWork * const work = local_work();
  const int result = work->m_active ? TPI_ERROR_ACTIVE : 0 ;
  return result ;
}

#endif

