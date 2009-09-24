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

/*--------------------------------------------------------------------*/

#include <TPI.h>
#include <stdlib.h>
#include <ThreadPool_config.h>

enum { THREAD_COUNT_MAX = 1024 };
enum { LOCK_COUNT_MAX   = 256 };

/*--------------------------------------------------------------------*/

#if defined( HAVE_PTHREAD )

#if	defined( __linux__ ) && defined( __GNUC__ ) && ( 4 <= __GNUC__ )

#define HAVE_FUTEX

#if	! defined( __INTEL_COMPILER )

#define HAVE_ATOMIC_SYNC

#else

#define HAVE_PTHREAD_SPIN

#endif /* ! defined( __INTEL_COMPILER ) */
#endif /* ! ( defined(__linux__) && defined(__GNUC__) && ( 4 <= __GNUC__ ) ) */

/*--------------------------------------------------------------------*/
/* #define DEBUG_PRINT */

#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/time.h>

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*  Performance is heavily impacted by an
 *  atomic decrement of the work counter.
 *  Optimize this if at all possible.
 */

#if defined( HAVE_ATOMIC_SYNC )

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
  int result ;
  pthread_once( & atomic_once , atomic_init );
  pthread_spin_lock( & atomic_lock );
  result = ( *value )-- ;
  pthread_spin_unlock( & atomic_lock );
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
/*--------------------------------------------------------------------*/

struct ThreadPool_Data ;

static void * local_driver( void * );

#if defined( HAVE_FUTEX )

#include <sys/syscall.h>
#include <linux/futex.h>
 
/* FUTEX_WAIT = 0 */
/* FUTEX_WAKE = 1 */
 
typedef struct Thread_Data {
  struct TPI_Work_Struct   m_work ;
  struct ThreadPool_Data * m_pool ;
  long                     m_rank ;
  long                     m_cycle ;
  long                     m_active ;
} Thread ;

static void wake_thread( Thread * thread , int val )
{
  thread->m_active = val ;
  syscall( SYS_futex, & thread->m_active, FUTEX_WAKE, 1, NULL, NULL, 0 );
}

static void wait_thread( Thread * thread , int val )
{
  syscall( SYS_futex, & thread->m_active, FUTEX_WAIT, val, NULL, NULL, 0 );
}

static int create_thread( pthread_attr_t * attr , Thread * thread )
{
  int result = 0 ;
  pthread_t pt ;

  thread->m_work.info       = NULL ;
  thread->m_work.reduce     = NULL ;
  thread->m_work.count      = 0 ;
  thread->m_work.rank       = 0 ;
  thread->m_work.lock_count = 0 ;

  thread->m_active          = 1 ;

  if ( thread->m_rank != 0 ) {
    if ( pthread_create( & pt, attr, & local_driver, thread ) ) {
      thread->m_active = 0 ;
      result = TPI_ERROR_INTERNAL ;
    }
  }

  return result ;
}

static void destroy_thread( Thread * thread )
{
  thread->m_work.info       = NULL ;
  thread->m_work.reduce     = NULL ;
  thread->m_work.count      = 0 ;
  thread->m_work.rank       = 0 ;
  thread->m_work.lock_count = 0 ;
}

/*--------------------------------------------------------------------*/

#else /* ! defined( HAVE_FUTEX ) */

typedef struct Thread_Data {
  pthread_cond_t           m_cond ;
  pthread_mutex_t          m_lock ;
  struct TPI_Work_Struct   m_work ;
  struct ThreadPool_Data * m_pool ;
  long                     m_rank ;
  long                     m_cycle ;
  long                     m_active ;
} Thread ;

static void wake_thread( Thread * thread , int val )
{
  pthread_mutex_lock(   & thread->m_lock );
  thread->m_active = val ;
  pthread_mutex_unlock( & thread->m_lock );
  pthread_cond_signal(  & thread->m_cond );
}

static void wait_thread( Thread * thread , int val )
{
  pthread_mutex_lock( & thread->m_lock );
  if ( thread->m_active == val ) {
    pthread_cond_wait( & thread->m_cond , & thread->m_lock );
  }
  pthread_mutex_unlock( & thread->m_lock );
}

static int create_thread( pthread_attr_t * attr , Thread * thread )
{
  int result = 0 ;
  pthread_t pt ;

  thread->m_work.info       = NULL ;
  thread->m_work.reduce     = NULL ;
  thread->m_work.count      = 0 ;
  thread->m_work.rank       = 0 ;
  thread->m_work.lock_count = 0 ;

  thread->m_active          = 1 ;

  if ( thread->m_rank != 0 ) {
    if ( pthread_mutex_init( & thread->m_lock , NULL ) ||
         pthread_cond_init(  & thread->m_cond , NULL ) ||
         pthread_create( & pt, attr, & local_driver, thread ) ) {
      thread->m_active = 0 ;
      result = TPI_ERROR_INTERNAL ;
    }
  }

  return result ;
}

static void destroy_thread( Thread * thread )
{
  if ( 0 < thread->m_rank ) {
    pthread_mutex_destroy( & thread->m_lock );
    pthread_cond_destroy(  & thread->m_cond );
  }

  thread->m_work.info       = NULL ;
  thread->m_work.reduce     = NULL ;
  thread->m_work.count      = 0 ;
  thread->m_work.rank       = 0 ;
  thread->m_work.lock_count = 0 ;
}

#endif

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

typedef void (*Thread_Work)( Thread * const );

typedef struct ThreadPool_Data {
  TPI_work_subprogram   m_work_routine ;
  int                   m_work_count_claim ;

  TPI_reduce_join       m_reduce_join ;
  TPI_reduce_init       m_reduce_init ;
  unsigned char       * m_reduce_alloc ;
  int                   m_reduce_alloc_size ;
  int                   m_reduce_size ;
  int                   m_reduce_grain ;

  Thread_Work           m_thread_driver ;
  int                   m_thread_count ;
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

    /* m_reduce_join         */  NULL ,
    /* m_reduce_init         */  NULL ,
    /* m_reduce_alloc        */  NULL ,
    /* m_reduce_alloc_size   */  0 ,
    /* m_reduce_size         */  0 ,
    /* m_reduce_grain        */  0 ,

    /* m_thread_driver       */  NULL ,
    /* m_thread_count        */  0 ,
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
static void local_run_work( Thread * const thread )
{
  struct TPI_Work_Struct * const work = &( thread->m_work );

  ThreadPool * const thread_pool = thread->m_pool ;

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
static void local_run_thread( Thread * const thread )
{
  struct TPI_Work_Struct * const work = &( thread->m_work );

  ThreadPool * const thread_pool = thread->m_pool ;

  work->rank = thread->m_rank ;

  (* thread_pool->m_work_routine)( work );
}

/*--------------------------------------------------------------------*/

static
void local_broadcast( Thread * const this_thread )
{
  ThreadPool * const thread_pool  = this_thread->m_pool ;
  const unsigned int thread_rank  = this_thread->m_rank ;
  const unsigned int thread_count = thread_pool->m_thread_count ;
        unsigned int thread_cycle = 1 ;
        unsigned int next_rank ;

  for ( ; thread_cycle <= thread_rank ; thread_cycle <<= 1 );

  for ( ; ( next_rank = thread_rank + thread_cycle ) < thread_count ;
          thread_cycle <<= 1 ) {
    Thread * const next_thread = this_thread + thread_cycle ;
 
    next_thread->m_work.info   = this_thread->m_work.info ;
    next_thread->m_work.reduce = NULL ;

    if ( thread_pool->m_reduce_init ) {
      next_thread->m_work.reduce = thread_pool->m_reduce_alloc +
                                   thread_pool->m_reduce_grain * next_rank ;

      thread_pool->m_reduce_init( & next_thread->m_work );
    }

    next_thread->m_work.rank       = -1 ;
    next_thread->m_work.count      = this_thread->m_work.count ;
    next_thread->m_work.lock_count = this_thread->m_work.lock_count ;

    /*  Activate the next thread.  It is known to be blocked. */
    wake_thread( next_thread , 1 );
  }
}
 
static
void local_barrier( Thread * const this_thread )
{
  ThreadPool * const thread_pool  = this_thread->m_pool ;
  const unsigned int thread_rank  = this_thread->m_rank ;
        unsigned int thread_cycle = this_thread->m_cycle ;

  for ( ; thread_rank < thread_cycle ; thread_cycle >>= 1 ) {
    Thread * const next_thread = this_thread + thread_cycle ;
 
    /*  Block if next_thread is still active. */
    wait_thread( next_thread , 1 );

    /* next_thread is done */

    next_thread->m_work.rank       = -1 ;
    next_thread->m_work.count      = 0 ;
    next_thread->m_work.lock_count = 0 ;

    if ( thread_pool->m_reduce_join ) {
      /* Reduce next_thread's data into this thread's data */

      (* thread_pool->m_reduce_join)( & this_thread->m_work ,
                                      next_thread->m_work.reduce );
    }

    next_thread->m_work.reduce = NULL ;
    next_thread->m_work.info   = NULL ;
  }

  /* My parent should be waiting for me to finish */
  wake_thread( this_thread , 0 );
}
 
/*--------------------------------------------------------------------*/
/*  The driver given to 'pthread_create'.
 *  Run work until told to terminate.
 */
static void * local_driver( void * arg )
{
  Thread     * const this_thread = (Thread *) arg ;
  ThreadPool * const thread_pool = this_thread->m_pool ;

  do {
    local_barrier( this_thread );   /*  Wait for my tree threads to block */

    wait_thread( this_thread , 0 ); /*  Block until I am activated. */

    local_broadcast( this_thread ); /*  Broadcast the activation */

    if ( thread_pool->m_thread_driver ) { /* Participate in work */
      (*thread_pool->m_thread_driver)( this_thread );
    }
  } while ( thread_pool->m_thread_driver );

  local_barrier( this_thread ); /* Termination barrier */

  return NULL ;
}

/*--------------------------------------------------------------------*/
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
      local_broadcast( root_thread );

      /* Participate in work queue */
      local_run_work( root_thread );

      /* Block the worker threads */
      local_barrier( root_thread );

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
                    TPI_reduce_join       reduce_join ,
                    TPI_reduce_init       reduce_init ,
                    int                   reduce_size ,
                    void *                reduce_data )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result =
    NULL != thread_pool->m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram             ? TPI_ERROR_NULL : (
    NULL == reduce_join                 ? TPI_ERROR_NULL : (
    NULL == reduce_init                 ? TPI_ERROR_NULL : (
    NULL == reduce_data                 ? TPI_ERROR_NULL : (
    work_count  <= 0                    ? TPI_ERROR_SIZE : (
    reduce_size <= 0                    ? TPI_ERROR_SIZE : (
    0 ) ) ) ) ) ) );

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
      thread_pool->m_reduce_join      = reduce_join ;
      thread_pool->m_reduce_init      = reduce_init ;
      thread_pool->m_work_routine     = work_subprogram ;
      thread_pool->m_work_count_claim = work_count ;

      /* Activate the blocked worker threads */
      local_broadcast( root_thread );

      /* Participate in work queue */
      local_run_work( root_thread );

      /* Block the worker threads */
      local_barrier( root_thread );

      thread_pool->m_thread_driver    = NULL ;
      thread_pool->m_work_count_claim = 0 ;
      thread_pool->m_work_routine     = NULL ;
      thread_pool->m_reduce_join      = NULL ;
      thread_pool->m_reduce_init      = NULL ;
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

      /* Activate the blocked worker threads */
      local_broadcast( root_thread );

      /* Participate in work */
      local_run_thread( root_thread );

      /* Block the worker threads */
      local_barrier( root_thread );

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
                            TPI_reduce_join       reduce_join ,
                            TPI_reduce_init       reduce_init ,
                            int                   reduce_size ,
                            void *                reduce_data )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result =
    NULL != thread_pool->m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram             ? TPI_ERROR_NULL : (
    NULL == reduce_join                 ? TPI_ERROR_NULL : (
    NULL == reduce_init                 ? TPI_ERROR_NULL : (
    NULL == reduce_data                 ? TPI_ERROR_NULL : (
    reduce_size <= 0                    ? TPI_ERROR_SIZE : (
    0 ) ) ) ) ) );

  if ( ! result ) {
    Thread * const root_thread = thread_pool->m_thread ;

    root_thread->m_work.info       = work_info ;
    root_thread->m_work.reduce     = reduce_data ;
    root_thread->m_work.rank       = 0 ;
    root_thread->m_work.lock_count = 0 ;

    if ( 1 < thread_pool->m_thread_count ) {

      alloc_reduce( thread_pool , reduce_size );

      thread_pool->m_thread_driver  = & local_run_thread ;
      thread_pool->m_reduce_join    = reduce_join ;
      thread_pool->m_reduce_init    = reduce_init ;
      thread_pool->m_work_routine   = work_subprogram ;

      root_thread->m_work.count = thread_pool->m_thread_count ;

      /* Activate the blocked worker threads */
      local_broadcast( root_thread );

      /* Participate in work */
      local_run_thread( root_thread );

      /* Block the worker threads */
      local_barrier( root_thread );

      thread_pool->m_thread_driver  = NULL ;
      thread_pool->m_work_routine   = NULL ;
      thread_pool->m_reduce_join    = NULL ;
      thread_pool->m_reduce_init    = NULL ;
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
  int result = thread_pool->m_thread_count ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && ( n < 1 || THREAD_COUNT_MAX <= n ) ) {
    result = TPI_ERROR_SIZE ;
  }

  if ( ! result ) {
    pthread_attr_t attr ;

    if ( pthread_attr_init( & attr ) ||
         pthread_attr_setscope(       & attr, PTHREAD_SCOPE_SYSTEM ) ||
         pthread_attr_setdetachstate( & attr, PTHREAD_CREATE_DETACHED ) ) {
      result = TPI_ERROR_INTERNAL ;
    }

    if ( ! result ) {
      int thread_upper ;
      int i ;

      for ( thread_upper = 1 ; thread_upper <= n ; thread_upper <<= 1 );

      thread_pool->m_thread_count = n ;

      /* Create threads last-to-first for start up fan-in barrier */

      for ( i = n ; ! result && 0 < i ; ) {
        int cycle = thread_upper ;

        --i ;
        for ( ; n <= i + cycle ; cycle >>= 1 );

        thread_pool->m_thread[i].m_pool  = thread_pool ;
        thread_pool->m_thread[i].m_rank  = i ;
        thread_pool->m_thread[i].m_cycle = cycle ;

        result = create_thread( & attr , thread_pool->m_thread + i );
      }

      if ( result ) {
        while ( i < --( thread_pool->m_thread_count ) ) {
          Thread * thread = thread_pool->m_thread +
                            thread_pool->m_thread_count ;
          wait_thread( thread , 0 ); /* Wait for thread to block */
          wake_thread( thread , 1 ); /* Activate thread */
          wait_thread( thread , 0 ); /* Wait for thread to terminate */
          destroy_thread( thread );
        }
        thread_pool->m_thread_count = 0 ;
      }

      pthread_attr_destroy( & attr );
    }
  }

  if ( ! result ) {
    local_barrier( thread_pool->m_thread );
    result = n ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  ThreadPool * const thread_pool = local_thread_pool();
  int result , i ;

  result = NULL != thread_pool->m_thread_driver ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {

    local_broadcast( thread_pool->m_thread );
    local_barrier(   thread_pool->m_thread );

    while ( 0 <= --( thread_pool->m_thread_count ) ) {
      destroy_thread( thread_pool->m_thread + thread_pool->m_thread_count );
    }

    thread_pool->m_thread_count = 0 ;

    while ( thread_pool->m_lock_init ) {
      --( thread_pool->m_lock_init );
      pthread_mutex_destroy( thread_pool->m_lock + thread_pool->m_lock_init );
    }

    if ( thread_pool->m_reduce_alloc ) {
      free( thread_pool->m_reduce_alloc );
      thread_pool->m_reduce_alloc = NULL ;
      thread_pool->m_reduce_alloc_size = 0 ;
    }
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#else /* ! defined( HAVE_PTHREAD ) */

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
                    TPI_reduce_join       reduce_join ,
                    TPI_reduce_init       reduce_init ,
                    int                   reduce_size ,
                    void *                reduce_data )
{
  LocalWork * const work = local_work();

  int result =
    0    != work->m_active    ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram   ? TPI_ERROR_NULL : (
    NULL == reduce_join       ? TPI_ERROR_NULL : (
    NULL == reduce_init       ? TPI_ERROR_NULL : (
    NULL == reduce_data       ? TPI_ERROR_NULL : (
    work_count  <= 0          ? TPI_ERROR_SIZE : (
    reduce_size <= 0          ? TPI_ERROR_SIZE : (
    0 ) ) ) ) ) ) );

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
                            TPI_reduce_join       reduce_join ,
                            TPI_reduce_init       reduce_init ,
                            int                   reduce_size ,
                            void *                reduce_data )
{
  return TPI_Run_reduce( work_subprogram, work_info, 1,
                         reduce_join, reduce_init, reduce_size, reduce_data );
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


