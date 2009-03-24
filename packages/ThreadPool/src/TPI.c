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
#include <string.h>

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

/*--------------------------------------------------------------------*/
/*  Performance is heavily impacted by an
 *  atomic decrement of the work counter.
 *  Optimize this if at all possible.
 */

#if defined(__GNUC__) && defined(__linux__) && ( 4 <= __GNUC__ )

#define atomic_fetch_and_decrement( VALUE_PTR )	\
	__sync_fetch_and_sub(VALUE_PTR,1)

#elif defined( HAVE_PTHREAD_SPIN )

static pthread_once_t     atomic_once = PTHREAD_ONCE_INIT ;
static pthread_spinlock_t atomic_lock ;

static void atomic_init()
{
  pthread_spin_init( & atomic_lock , PTHREAD_PROCESS_PRIVATE );
}

static
int atomic_fetch_and_decrement( volatile int * value )
{
  pthread_once( & atomic_once , once_init );
  int result ;
  pthread_spin_lock( & lock );
  result = ( *value )-- ;
  pthread_spin_unlock( & lock );
  return result ;
}

#else

static pthread_mutex_t atomic_lock = PTHREAD_MUTEX_INITIALIZER ;

static
int atomic_fetch_and_decrement( volatile int * value )
{
  int result ;
  while ( pthread_mutex_trylock( & atomic_lock ) );
  result = ( *value )-- ;
  pthread_mutex_unlock( & atomic_lock );
  return result ;
}

#endif

/*--------------------------------------------------------------------*/

struct ThreadPool_Data ;

typedef void (*Thread_Work)( struct ThreadPool_Data * const , unsigned int );

typedef struct Thread_Data {
  pthread_cond_t         m_cond ;
  pthread_mutex_t        m_lock ;
  int                    m_must_run ;
  struct TPI_Work_Struct m_work ;
} Thread ;

typedef struct ThreadPool_Data {
  TPI_work_subprogram   m_work_routine ;
  int                   m_work_count_claim ;

  TPI_reduce_subprogram m_reduce_routine ;
  unsigned char       * m_reduce_alloc ;
  int                   m_reduce_alloc_size ;
  int                   m_reduce_size ;
  int                   m_reduce_grain ;

  Thread_Work           m_thread_driver ;
  int                   m_thread_count ;
  int                   m_thread_upper ;
  int                   m_lock_init ;

  Thread                m_thread[ THREAD_COUNT_MAX ];
  pthread_mutex_t       m_lock[ LOCK_COUNT_MAX ];
} ThreadPool ;

/*--------------------------------------------------------------------*/

static ThreadPool * local_thread_pool()
{
  static ThreadPool thread_pool = {
    /* m_work_routine        */  NULL ,
    /* m_work_count_claim    */  0 ,

    /* m_reduce_routine      */  NULL ,
    /* m_reduce_alloc        */  NULL ,
    /* m_reduce_alloc_size   */  0 ,
    /* m_reduce_size         */  0 ,
    /* m_reduce_grain        */  0 ,

    /* m_thread_driver       */  NULL ,
    /* m_thread_count        */  0 ,
    /* m_thread_upper        */  0 ,
    /* m_lock_init           */  0 };

  return & thread_pool ;
}

/*--------------------------------------------------------------------*/

int TPI_Lock( int i )
{
  ThreadPool * const thread_pool = local_thread_pool();
  const int lock_count = thread_pool->m_thread->m_work.lock_count ;

  int result = i < 0 || lock_count <= i ? TPI_ERROR_SIZE : 0 ;

  if ( ! result && pthread_mutex_lock( thread_pool->m_lock + i ) ) {
    result = TPI_ERROR_LOCK ;
  }

  return result ;
}

int TPI_Trylock( int i )
{
  ThreadPool * const thread_pool = local_thread_pool();
  const int lock_count = thread_pool->m_thread->m_work.lock_count ;

  int result = i < 0 || lock_count <= i ? TPI_ERROR_SIZE : 0 ;

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
  const int lock_count = thread_pool->m_thread->m_work.lock_count ;

  int result = i < 0 || lock_count <= i ? TPI_ERROR_SIZE : 0 ;

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
                            unsigned int thread_rank )
{
  struct TPI_Work_Struct * const work =
    &( thread_pool->m_thread[ thread_rank ].m_work );

  int * const claim = & thread_pool->m_work_count_claim ;

  while ( 0 < ( work->rank = atomic_fetch_and_decrement(claim))) {

    work->rank = work->count - work->rank ;

    (* thread_pool->m_work_routine)( work );
  }
  return ;
}

/*--------------------------------------------------------------------*/
/*  Run the work subprogram exactly once for each thread.
 *  The last thread signals the root thread.
 */
static void local_run_thread( ThreadPool * const thread_pool ,
                              unsigned int thread_rank )
{
  struct TPI_Work_Struct * const work =
    &( thread_pool->m_thread[ thread_rank ].m_work );

  work->rank = thread_rank ;

  (* thread_pool->m_work_routine)( work );
}

/*--------------------------------------------------------------------*/

static
void local_broadcast( ThreadPool * const control ,
                      const unsigned int thread_rank )
{
  Thread * const     thread_data  = control->m_thread ;
  const unsigned int thread_count = control->m_thread_count ;
        unsigned int thread_cycle = 1 ;
        unsigned int next_rank ;

  for ( ; thread_cycle <= thread_rank ; thread_cycle <<= 1 );

  for ( ; ( next_rank = thread_rank + thread_cycle ) < thread_count ;
          thread_cycle <<= 1 ) {
    Thread * const this_thread = thread_data + thread_rank ;
    Thread * const next_thread = thread_data + next_rank ;
 
#ifdef DEBUG_PRINT
  fprintf(stdout,"  broadcast run %u to %u\n",thread_rank,next_rank);
  fflush(stdout);
#endif

    next_thread->m_must_run        = this_thread->m_must_run ;
    next_thread->m_work.info       = this_thread->m_work.info ;
    next_thread->m_work.reduce     = NULL ;
    next_thread->m_work.rank       = -1 ;
    next_thread->m_work.count      = this_thread->m_work.count ;
    next_thread->m_work.lock_count = this_thread->m_work.lock_count ;

    if ( control->m_reduce_routine ) {
      next_thread->m_work.reduce = control->m_reduce_alloc +
                                   control->m_reduce_grain * next_rank ;

      memcpy( next_thread->m_work.reduce ,
              this_thread->m_work.reduce ,
              control->m_reduce_size );
    }

    pthread_cond_signal(  & next_thread->m_cond );
    pthread_mutex_unlock( & next_thread->m_lock );
  }
}
 
static
void local_barrier( ThreadPool * const control ,
                    const unsigned int thread_rank )
{
  Thread * const     thread_data  = control->m_thread ;
  Thread * const     this_thread  = thread_data + thread_rank ;
  const unsigned int thread_count = control->m_thread_count ;
        unsigned int thread_cycle = control->m_thread_upper ;

  for ( ; thread_count <= thread_rank + thread_cycle ; thread_cycle >>= 1 );

  for ( ; thread_rank < thread_cycle ; thread_cycle >>= 1 ) {
    const unsigned next_rank = thread_rank + thread_cycle ;
    Thread * const next_thread = thread_data + next_rank ;
 
#ifdef DEBUG_PRINT
  fprintf(stdout,"  barrier halt %u from %u lock\n",thread_rank,next_rank);
  fflush(stdout);
#endif
 
    pthread_mutex_lock( & next_thread->m_lock );

    if ( next_thread->m_must_run ) {
 
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

    /* Reduce next_thread's data into this thread's data */

    if ( control->m_reduce_routine ) {
      this_thread->m_work.count = thread_count ;
      this_thread->m_work.rank  = thread_rank ;

      (* control->m_reduce_routine)( & this_thread->m_work ,
                                       next_thread->m_work.reduce );
    }

    /* next_thread is done */

    next_thread->m_work.info       = NULL ;
    next_thread->m_work.reduce     = NULL ;
    next_thread->m_work.rank       = -1 ;
    next_thread->m_work.count      = 0 ;
    next_thread->m_work.lock_count = 0 ;
  }

  if ( thread_data->m_must_run ) { /* My parent is waiting for my signal */
    this_thread->m_must_run = 0 ;
    pthread_cond_signal( & this_thread->m_cond );
  }
}
 
/*--------------------------------------------------------------------*/
/*  The driver given to 'pthread_create'.
 *  Run work until told to terminate.
 */
static void * local_driver( void * arg )
{
  ThreadPool * const thread_pool = (ThreadPool *) arg ;
  Thread     * const root_thread = thread_pool->m_thread ;

  const unsigned my_rank   = --thread_pool->m_work_count_claim ;
  Thread * const my_thread = root_thread + my_rank ;

  /* Aquire my lock */
  pthread_mutex_lock( & my_thread->m_lock );

  /* Signal root thread that I have started */

  pthread_mutex_lock(   & root_thread->m_lock );
  pthread_cond_signal(  & root_thread->m_cond );
  pthread_mutex_unlock( & root_thread->m_lock );

  /*------------------------------*/

  do {
    /* Acquire locks for fan-in / fan-out */
    local_barrier( thread_pool , my_rank );

    /* Wait for the run signal, giving up my lock */
    pthread_cond_wait( & my_thread->m_cond , & my_thread->m_lock );

    /* Continue broadcasting the run signal */
    local_broadcast( thread_pool , my_rank );

    if ( thread_pool->m_thread_driver ) {

#ifdef DEBUG_PRINT
  fprintf(stdout,"  thread %u working\n",my_rank);
  fflush(stdout);
#endif

      /* Participate in the work */
      (*thread_pool->m_thread_driver)( thread_pool , my_rank );
    }
  } while ( thread_pool->m_thread_driver );

  /*------------------------------*/
  /* Termination, unlock my thread and count down */

  pthread_mutex_unlock( & my_thread->m_lock );

  pthread_mutex_lock(  & root_thread->m_lock );
  if ( ! --thread_pool->m_work_count_claim ) {
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

int TPI_Run( TPI_work_subprogram work_subprogram  ,
             const void *        work_info ,
             int                 work_count  ,
             int                 lock_count  )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result =
    NULL != thread_pool->m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram             ? TPI_ERROR_NULL : (
    work_count  < 0                     ? TPI_ERROR_SIZE : (
    lock_count  < 0                     ? TPI_ERROR_SIZE : (
    0 ) ) ) );

  if ( ! result && ! ( result = set_lock_count(thread_pool,lock_count) ) ) {
    Thread * const root_thread = thread_pool->m_thread ;

    root_thread->m_work.info       = work_info ;
    root_thread->m_work.reduce     = NULL ;
    root_thread->m_work.rank       = -1 ;
    root_thread->m_work.count      = work_count ;
    root_thread->m_work.lock_count = lock_count ;

    if ( 1 < thread_pool->m_thread_count ) {
      thread_pool->m_work_routine     = work_subprogram ;
      thread_pool->m_work_count_claim = work_count ;
      thread_pool->m_thread_driver    = & local_run_work ;

      /* Activate the blocked worker threads */
      local_broadcast( thread_pool , 0 );

      /* Participate in work queue */
      local_run_work( thread_pool , 0 );

      /* Block the worker threads */
      local_barrier( thread_pool , 0 );

      thread_pool->m_thread_driver    = NULL ;
      thread_pool->m_work_count_claim = 0 ;
      thread_pool->m_work_routine     = NULL ;
    }
    else {
      /* No worker threads, can bypass work queue locking */
      struct TPI_Work_Struct * const work = & root_thread->m_work ;

      for ( work->rank = 0 ; work->rank < work->count ; ++( work->rank ) ) {
        (* work_subprogram)( work );
      }
    }

    root_thread->m_work.info       = NULL ;
    root_thread->m_work.reduce     = NULL ;
    root_thread->m_work.rank       = -1 ;
    root_thread->m_work.count      = 0 ;
    root_thread->m_work.lock_count = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/

static void alloc_reduce( ThreadPool * const thread_pool ,
                          int                reduce_size )
{
  const int grain = 128 ; /* Byte grain size */
  const int grain_count  = ( reduce_size + grain - 1 ) / grain ;
  const int reduce_grain = grain * grain_count ; 
  const int alloc_size   = reduce_grain * thread_pool->m_thread_count ;

  if ( thread_pool->m_reduce_alloc_size < alloc_size ) {
    thread_pool->m_reduce_alloc_size = alloc_size ;

    if ( thread_pool->m_reduce_alloc ) {
      thread_pool->m_reduce_alloc =
        (unsigned char *) realloc( thread_pool->m_reduce_alloc , alloc_size );
    }
    else {
      thread_pool->m_reduce_alloc = (unsigned char *) malloc( alloc_size );
    }
  }

  thread_pool->m_reduce_size  = reduce_size ;
  thread_pool->m_reduce_grain = reduce_grain ;
}

int TPI_Run_reduce( TPI_work_subprogram   work_subprogram  ,
                    const void *          work_info ,
                    int                   work_count  ,
                    TPI_reduce_subprogram reduce_subprogram ,
                    void *                reduce_data ,
                    int                   reduce_size )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result =
    NULL != thread_pool->m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram             ? TPI_ERROR_NULL : (
    NULL == reduce_subprogram           ? TPI_ERROR_NULL : (
    NULL == reduce_data                 ? TPI_ERROR_NULL : (
    work_count  < 0                     ? TPI_ERROR_SIZE : (
    reduce_size < 0                     ? TPI_ERROR_SIZE : (
    0 ) ) ) ) ) );

  if ( ! result ) {
    Thread * const root_thread = thread_pool->m_thread ;

    root_thread->m_work.info       = work_info ;
    root_thread->m_work.reduce     = reduce_data ;
    root_thread->m_work.rank       = -1 ;
    root_thread->m_work.count      = work_count ;
    root_thread->m_work.lock_count = 0 ;

    if ( 1 < thread_pool->m_thread_count ) {
      alloc_reduce( thread_pool, reduce_size );

      thread_pool->m_thread_driver    = & local_run_work ;
      thread_pool->m_reduce_routine   = reduce_subprogram ;
      thread_pool->m_work_routine     = work_subprogram ;
      thread_pool->m_work_count_claim = work_count ;

      /* Activate the blocked worker threads */
      local_broadcast( thread_pool , 0 );

      /* Participate in work queue */
      local_run_work( thread_pool , 0 );

      /* Block the worker threads */
      local_barrier( thread_pool , 0 );

      thread_pool->m_thread_driver    = NULL ;
      thread_pool->m_work_count_claim = 0 ;
      thread_pool->m_work_routine     = NULL ;
      thread_pool->m_reduce_routine   = NULL ;
      thread_pool->m_reduce_grain     = 0 ;
      thread_pool->m_reduce_size      = 0 ;
    }
    else {
      /* No worker threads, can bypass work queue locking */
      struct TPI_Work_Struct * const work = & root_thread->m_work ;

      for ( work->rank = 0 ; work->rank < work->count ; ++( work->rank ) ) {
        (* work_subprogram)( work );
      }
    }

    root_thread->m_work.info       = NULL ;
    root_thread->m_work.reduce     = NULL ;
    root_thread->m_work.rank       = -1 ;
    root_thread->m_work.count      = 0 ;
    root_thread->m_work.lock_count = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Run_threads( TPI_work_subprogram work_subprogram  ,
                     const void *        work_info ,
                     int                 lock_count  )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result =
    NULL != thread_pool->m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram             ? TPI_ERROR_NULL : (
    0 ) );

  if ( ! result && ! ( result = set_lock_count(thread_pool,lock_count) ) ) {
    Thread * const root_thread = thread_pool->m_thread ;

    root_thread->m_work.info       = work_info ;
    root_thread->m_work.reduce     = NULL ;
    root_thread->m_work.rank       = 0 ;
    root_thread->m_work.lock_count = lock_count ;

    if ( 1 < thread_pool->m_thread_count ) {
      thread_pool->m_work_routine  = work_subprogram ;
      thread_pool->m_thread_driver = & local_run_thread ;

      root_thread->m_work.count = thread_pool->m_thread_count ;
      root_thread->m_must_run   = 1 ;

      /* Activate the blocked worker threads */
      local_broadcast( thread_pool , 0 );

      /* Participate in work */
      local_run_thread( thread_pool , 0 );

      /* Block the worker threads */
      local_barrier( thread_pool , 0 );

      root_thread->m_must_run = 1 ;

      thread_pool->m_thread_driver = NULL ;
      thread_pool->m_work_routine  = NULL ;
    }
    else {
      /* No worker threads */
      root_thread->m_work.count = 1 ;
      (* work_subprogram)( & root_thread->m_work );
    }

    root_thread->m_work.info       = NULL ;
    root_thread->m_work.reduce     = NULL ;
    root_thread->m_work.rank       = -1 ;
    root_thread->m_work.count      = 0 ;
    root_thread->m_work.lock_count = 0 ;
  }

  return result ;
}

int TPI_Run_threads_reduce( TPI_work_subprogram   work_subprogram  ,
                            const void *          work_info ,
                            TPI_reduce_subprogram reduce_subprogram ,
                            void *                reduce_data ,
                            int                   reduce_size )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result =
    NULL != thread_pool->m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram             ? TPI_ERROR_NULL : (
    NULL == reduce_subprogram           ? TPI_ERROR_NULL : (
    NULL == reduce_data                 ? TPI_ERROR_NULL : (
    reduce_size < 0                     ? TPI_ERROR_SIZE : (
    0 ) ) ) ) );

  if ( ! result ) {
    Thread * const root_thread = thread_pool->m_thread ;

    root_thread->m_work.info       = work_info ;
    root_thread->m_work.reduce     = reduce_data ;
    root_thread->m_work.rank       = 0 ;
    root_thread->m_work.lock_count = 0 ;

    if ( 1 < thread_pool->m_thread_count ) {

      alloc_reduce( thread_pool , reduce_size );

      thread_pool->m_thread_driver  = & local_run_thread ;
      thread_pool->m_reduce_routine = reduce_subprogram ;
      thread_pool->m_work_routine   = work_subprogram ;

      root_thread->m_work.count = thread_pool->m_thread_count ;
      root_thread->m_must_run   = 1 ;

      /* Activate the blocked worker threads */
      local_broadcast( thread_pool , 0 );

      /* Participate in work */
      local_run_thread( thread_pool , 0 );

      /* Block the worker threads */
      local_barrier( thread_pool , 0 );

      root_thread->m_must_run = 1 ;

      thread_pool->m_thread_driver  = NULL ;
      thread_pool->m_work_routine   = NULL ;
      thread_pool->m_reduce_routine = NULL ;
      thread_pool->m_reduce_grain   = 0 ;
      thread_pool->m_reduce_size    = 0 ;
    }
    else {
      /* No worker threads */
      root_thread->m_work.count = 1 ;
      (* work_subprogram)( & root_thread->m_work );
    }

    root_thread->m_work.info   = NULL ;
    root_thread->m_work.reduce = NULL ;
    root_thread->m_work.rank   = -1 ;
    root_thread->m_work.count  = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
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
    thread_pool->m_thread[i].m_must_run = 0 ;
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
    thread_pool->m_work_count_claim = n ;

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

    thread_pool->m_work_count_claim = 0 ;

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

    thread_pool->m_work_count_claim = thread_pool->m_thread_count - 1 ;

    local_broadcast( thread_pool , 0 );

    if ( 0 < thread_pool->m_work_count_claim ) {
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
    if ( thread_pool->m_reduce_alloc ) {
      free( thread_pool->m_reduce_alloc );
      thread_pool->m_reduce_alloc = NULL ;
      thread_pool->m_reduce_alloc_size = 0 ;
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

int TPI_Run( TPI_work_subprogram work_subprogram  ,
             const void *        work_info ,
             int                 work_count  ,
             int                 lock_count  )
{
  LocalWork * const work = local_work();

  int result = work->m_active ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {

    result = set_lock_count( work , lock_count );

    if ( ! result ) {
      struct TPI_Work_Struct w ;

      work->m_active     = 1 ;
      work->m_lock_count = lock_count ;

      w.info       = work_info ;
      w.reduce     = NULL ;
      w.count      = work_count ;
      w.rank       = 0 ;
      w.lock_count = lock_count ;

      for ( w.rank = 0 ; w.rank < work_count ; ++(w.rank) ) {
        (*work_subprogram)( & w );
      }

      work->m_active     = 0 ;
      work->m_lock_count = 0 ;
    }
  }

  return result ;
}

int TPI_Run_reduce( TPI_work_subprogram   work_subprogram  ,
                    const void *          work_info ,
                    int                   work_count  ,
                    TPI_reduce_subprogram reduce_subprogram ,
                    void *                reduce_data ,
                    int                   reduce_size )
{
  LocalWork * const work = local_work();

  int result =
    0    != work->m_active    ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram   ? TPI_ERROR_NULL : (
    NULL == reduce_subprogram ? TPI_ERROR_NULL : (
    NULL == reduce_data       ? TPI_ERROR_NULL : (
    work_count  < 0           ? TPI_ERROR_SIZE : (
    reduce_size < 0           ? TPI_ERROR_SIZE : (
    0 ) ) ) ) ) );

  if ( ! result ) {
    struct TPI_Work_Struct w ;

    work->m_active = 1 ;

    w.info       = work_info ;
    w.reduce     = reduce_data ;
    w.count      = work_count ;
    w.rank       = 0 ;
    w.lock_count = 0 ;

    for ( ; w.rank < w.count ; ++( w.rank ) ) {
      (* work_subprogram)( & w );
    }

    work->m_active = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Run_threads( TPI_work_subprogram work_subprogram  ,
                     const void *        work_info ,
                     int                 lock_count  )
{
  return TPI_Run( work_subprogram , work_info , 1 , lock_count );
}

int TPI_Run_threads_reduce( TPI_work_subprogram   work_subprogram  ,
                            const void *          work_info ,
                            TPI_reduce_subprogram reduce_subprogram ,
                            void *                reduce_data ,
                            int                   reduce_size )
{
  return TPI_Run_reduce( work_subprogram , work_info , 1 ,
                         reduce_subprogram , reduce_data , reduce_size );
}

int TPI_Init( int nthread )
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

