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
#include <stdio.h>
#include <ThreadPool_config.h>

/*--------------------------------------------------------------------*/
/*----------- PTHREAD CONFIGURATION (BEGIN) --------------------------*/
/*--------------------------------------------------------------------*/

#if	defined( HAVE_PTHREAD )

#include <errno.h>
#include <pthread.h>
#include <sched.h>

/*--------------------------------------------------------------------*/
/*---------------- COMPILER SPECIFICS (BEGIN) ------------------------*/
/*--------------------------------------------------------------------*/

/*  Performance is heavily impacted by an
 *  atomic decrement of the work counter.
 *  Optimize this if at all possible.
 */

#if	defined( __INTEL_COMPILER )

#define THREADPOOL_CONFIG "PTHREAD SCHED_YIELD"

#elif	defined( __linux__ ) && \
	defined( __GNUC__ ) && ( 4 <= __GNUC__ )

#define THREADPOOL_CONFIG "PTHREAD SCHED_YIELD ATOMIC_SYNC"

#define atomic_fetch_and_decrement( VALUE_PTR )	\
	__sync_fetch_and_sub( VALUE_PTR , 1 )

#else

#define THREADPOOL_CONFIG "PTHREAD SCHED_YIELD"

#endif

#if ! defined( atomic_fetch_and_decrement )

static int atomic_fetch_and_decrement( volatile int * value )
{
  static pthread_mutex_t atomic_lock = PTHREAD_MUTEX_INITIALIZER ;
  int result ;
  while ( EBUSY == pthread_mutex_trylock( & atomic_lock ) );
  result = ( *value )-- ;
  pthread_mutex_unlock( & atomic_lock );
  return result ;
}

#endif

/*--------------------------------------------------------------------*/
/*---------------- COMPILER SPECIFICS (END) --------------------------*/
/*--------------------------------------------------------------------*/

typedef pthread_mutex_t  local_lock_type ;

#else /* ! defined( HAVE_PTHREAD ) */

#define THREADPOOL_CONFIG "NO THREADING"

typedef int  local_lock_type ;

#endif

/*--------------------------------------------------------------------*/
/*----------- PTHREAD CONFIGURATION (END) ----------------------------*/
/*--------------------------------------------------------------------*/

const char * TPI_Version()
{
  static const char version_string[] =
    "TPI Version 1.1 , November 2009 , Configuration = " THREADPOOL_CONFIG ;

  return version_string ;
}

/*--------------------------------------------------------------------*/

enum { THREAD_COUNT_MAX = 256 };
enum { LOCK_COUNT_MAX   = 32 };

struct ThreadPool_Data ;

typedef struct Thread_Data {
  struct Thread_Data * m_thread_fan ; /* Fan-in / fan-out begin */
  void               * m_reduce ;     /* Reduction memory */
  long                 m_rank ;
  long                 m_barrier_wait_max ;
  long                 m_barrier_wait_total ;
  long                 m_barrier_wait_count ;
  volatile long        m_control ;
} Thread ;

typedef struct ThreadPool_Data {
  TPI_work_subprogram   m_work_routine ;
  const void *          m_work_info ;
  TPI_reduce_join       m_reduce_join ;
  TPI_reduce_init       m_reduce_init ;
  unsigned char       * m_reduce_alloc ;
  int                   m_reduce_alloc_size ;
  int                   m_thread_count ;
  int                   m_lock_init ;
  int                   m_lock_count ;
  int                   m_work_thread_count ;
  int                   m_work_count ;
  int                   m_work_count_claim ;

  Thread                m_thread[ THREAD_COUNT_MAX ];
  local_lock_type       m_lock[ LOCK_COUNT_MAX ];
} ThreadPool ;


static ThreadPool thread_pool =
{
  /* m_work_routine        */  NULL ,
  /* m_work_info           */  NULL ,
  /* m_reduce_join         */  NULL ,
  /* m_reduce_init         */  NULL ,
  /* m_reduce_alloc        */  NULL ,
  /* m_reduce_alloc_size   */  0 ,
  /* m_thread_count        */  0 ,
  /* m_lock_init           */  0 ,
  /* m_lock_count          */  0 ,
  /* m_work_thread_count   */  0 ,
  /* m_work_count          */  0 ,
  /* m_work_count_claim    */  0
};

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#if defined( HAVE_PTHREAD )

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Lock( int i )
{
  int result = i < 0 || thread_pool.m_lock_count <= i ? TPI_ERROR_SIZE : 0 ;

  if ( ! result ) {
    pthread_mutex_t * const lock = thread_pool.m_lock + i ;

    while ( EBUSY == ( result = pthread_mutex_trylock( lock ) ) );

    if ( result ) { result = TPI_ERROR_LOCK ; }
  }
  return result ;
}

int TPI_Unlock( int i )
{
  int result = i < 0 || thread_pool.m_lock_count <= i ? TPI_ERROR_SIZE : 0 ;

  if ( ! result && pthread_mutex_unlock( thread_pool.m_lock + i ) ) {
    result = TPI_ERROR_LOCK ;
  }

  return result ;
}

static int local_set_lock_count( const int lock_count )
{
  int result = lock_count < 0 || LOCK_COUNT_MAX < lock_count
             ? TPI_ERROR_SIZE : 0 ;

  while ( ! result && thread_pool.m_lock_init < lock_count ) {

    pthread_mutex_t * const lock = thread_pool.m_lock +
                                   thread_pool.m_lock_init ;

    if ( pthread_mutex_init( lock , NULL ) ) {
      result = TPI_ERROR_INTERNAL ;
    }
    else {
      ++( thread_pool.m_lock_init );
    }
  }

  return result ;
}

static void local_destroy_locks()
{
  while ( thread_pool.m_lock_init ) {
    --( thread_pool.m_lock_init );
    pthread_mutex_destroy( thread_pool.m_lock + thread_pool.m_lock_init );
  }
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/* Run work if any, then wait for child threads to block. */

static void local_run( Thread * const this_thread , void * reduce )
{
  struct TPI_Work_Struct work ;

  work.info       = thread_pool.m_work_info ;
  work.reduce     = reduce ;
  work.count      = thread_pool.m_work_count ;
  work.lock_count = thread_pool.m_lock_count ;

  if ( work.count <= thread_pool.m_work_thread_count ) {

    work.rank = ( thread_pool.m_thread_count - 1 ) - this_thread->m_rank ;

    if ( work.rank < work.count ) {
      (*thread_pool.m_work_routine)( & work );
    }
  }
  else {

    int * const claim = & thread_pool.m_work_count_claim ;

    while ( 0 < ( work.rank = atomic_fetch_and_decrement( claim ))) {

      work.rank = work.count - work.rank ;

      (*thread_pool.m_work_routine)( & work );
    }
  }
}

static int wait_thread( volatile long * const control , const int val )
{
  int count = 0 ;
  while ( val == *control ) {
    sched_yield();
    ++count ;
  }
  return count ;
}

static void local_barrier_wait( Thread * const this_thread ,
                                Thread * const thread )
{
  const long count = wait_thread( & thread->m_control , 1 );

  ++( this_thread->m_barrier_wait_count );

  this_thread->m_barrier_wait_total += count ;

  if ( this_thread->m_barrier_wait_max < count ) {
    this_thread->m_barrier_wait_max = count ;
  }
}

static void local_barrier( Thread * const this_thread )
{
  Thread * const thread_beg = this_thread[0].m_thread_fan ;
  Thread *       thread     = this_thread[1].m_thread_fan ;

  if ( ! thread_pool.m_work_routine ) {
    while ( thread_beg < thread ) {
      --thread ; local_barrier_wait( this_thread , thread );
    }
  }
  else if ( ! thread_pool.m_reduce_join ) {

    local_run( this_thread , NULL );

    while ( thread_beg < thread ) {
      --thread ; local_barrier_wait( this_thread , thread );
    }
  }
  else {

    /* Work data for the reduction initialization and join */

    struct TPI_Work_Struct work ;

    work.info       = thread_pool.m_work_info ;
    work.reduce     = this_thread->m_reduce ;
    work.count      = -1 ;
    work.rank       = -1 ;
    work.lock_count = -1 ;

    /* Initialize reduction value for non-root thread */

    if ( this_thread->m_rank ) { (*thread_pool.m_reduce_init)( & work ); }

    /* Run the work routine with barrier blocking */

    local_run( this_thread , work.reduce );

    /* Reduction of thread's contributions */

    while ( thread_beg < thread ) {
      --thread ; local_barrier_wait( this_thread , thread );
      (*thread_pool.m_reduce_join)( & work , thread->m_reduce );
    }
  }
}
 
/*--------------------------------------------------------------------*/
/*  The driver given to 'pthread_create'.
 *  Run work until told to terminate.
 */
static void * local_driver( void * arg )
{
  Thread * const this_thread = (Thread *) arg ;

  do {
    /* Wait for my subtree of threads to complete */
    local_barrier( this_thread );

    this_thread->m_control = 0 ;

    /*  Spin until I am activated. */
    wait_thread( & this_thread->m_control , 0 );

  } while ( thread_pool.m_work_routine );

  local_barrier( this_thread ); /* Termination barrier */

  this_thread->m_control = 0 ;

  return NULL ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
 
static void alloc_reduce( int reduce_size )
{
  const int alloc_count = thread_pool.m_thread_count - 1 ;

  if ( thread_pool.m_reduce_alloc_size < alloc_count * reduce_size ) {

    const int grain_shift  = 8 ; /* grain_size = 0x80 */
    const int grain_size   = 1 << grain_shift ; /* Byte grain size */
    const int grain_count  = ( reduce_size + grain_size - 1 ) >> grain_shift ;
    const int reduce_grain = grain_size * grain_count ; 
    const int alloc_size   = alloc_count * reduce_grain ;

    int i ;

    if ( thread_pool.m_reduce_alloc ) {
      thread_pool.m_reduce_alloc =
        (unsigned char *) realloc( thread_pool.m_reduce_alloc , alloc_size );
    }
    else {
      thread_pool.m_reduce_alloc = (unsigned char *) malloc( alloc_size );
    }

    thread_pool.m_reduce_alloc_size = alloc_size ;

    for ( i = 0 ; i < alloc_count ; ++i ) {
      thread_pool.m_thread[i+1].m_reduce =
        thread_pool.m_reduce_alloc + reduce_grain * i ;
    }
  }
}

static int local_start(
  int                   work_thread_count ,
  TPI_work_subprogram   work_subprogram  ,
  const void *          work_info ,
  int                   work_count  ,
  int                   lock_count ,
  TPI_reduce_join       reduce_join ,
  TPI_reduce_init       reduce_init ,
  int                   reduce_size ,
  void *                reduce_data )
{
  const int result = lock_count ? local_set_lock_count( lock_count ) : 0 ;

  if ( ! result ) {

    thread_pool.m_work_routine     = work_subprogram ;
    thread_pool.m_work_info        = work_info ;
    thread_pool.m_work_count       = work_count ;
    thread_pool.m_lock_count       = lock_count ;
    thread_pool.m_thread->m_reduce = reduce_data ;

    if ( 1 < thread_pool.m_thread_count ) {

      if ( reduce_size ) { alloc_reduce( reduce_size ); }

      thread_pool.m_reduce_join       = reduce_join ;
      thread_pool.m_reduce_init       = reduce_init ;
      thread_pool.m_work_thread_count = work_thread_count ;
      thread_pool.m_work_count_claim  = work_count ;

      /* Activate the spinning worker threads */
      {
        Thread * const thread_beg = thread_pool.m_thread + 1 ;
        Thread *       thread     = thread_pool.m_thread +
                                    thread_pool.m_thread_count ;

        while ( thread_beg < thread ) { (--thread)->m_control = 1 ; }
      }
    }
  }

  return result ;
}

static void local_wait()
{
  if ( 1 < thread_pool.m_thread_count ) {

    local_barrier( thread_pool.m_thread );

    thread_pool.m_reduce_join       = NULL ;
    thread_pool.m_reduce_init       = NULL ;
    thread_pool.m_work_thread_count = 0 ;
    thread_pool.m_work_count_claim  = 0 ;
  }
  else {
    struct TPI_Work_Struct w = { NULL , NULL , 0 , 0 , 0 };

    w.info       = thread_pool.m_work_info ;
    w.count      = thread_pool.m_work_count ;
    w.lock_count = thread_pool.m_lock_count ;
    w.reduce     = thread_pool.m_thread->m_reduce ;

    for ( w.rank = 0 ; w.rank < w.count ; ++( w.rank ) ) {
      (* thread_pool.m_work_routine )( & w );
    }
  }

  thread_pool.m_work_routine     = NULL ;
  thread_pool.m_work_info        = NULL ;
  thread_pool.m_work_count       = 0 ;
  thread_pool.m_lock_count       = 0 ;
  thread_pool.m_thread->m_reduce = NULL ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Init( int n )
{
  int result = thread_pool.m_thread_count ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && ( n < 1 || THREAD_COUNT_MAX + 1 <= n ) ) {
    result = TPI_ERROR_SIZE ;
  }

  if ( ! result ) {
    pthread_attr_t attr ;

    if ( pthread_attr_init( & attr )
         || pthread_attr_setscope(       & attr, PTHREAD_SCOPE_SYSTEM )
         || pthread_attr_setdetachstate( & attr, PTHREAD_CREATE_DETACHED ) ) {
      result = TPI_ERROR_INTERNAL ;
    }

    if ( ! result ) {
      int thread_rank = 0 ;
      int count = 1 ;

      /* Initialize one lock for blocking and unblocking */

      local_set_lock_count( 1 );

      /* Initialize threads with fan-in / fan-out span of threads */

      for ( thread_rank = 0 ; thread_rank <= n ; ++thread_rank ) {
        Thread * const thread = thread_pool.m_thread + thread_rank ;

        thread->m_thread_fan         = thread_pool.m_thread + count ;
        thread->m_reduce             = NULL ;
        thread->m_rank               = thread_rank ;
        thread->m_barrier_wait_max   = 0 ;
        thread->m_barrier_wait_total = 0 ;
        thread->m_barrier_wait_count = 0 ;
        thread->m_control            = 1 ;

        {
          int up = 1 ;
          while ( up <= thread_rank )    { up <<= 1 ; }
          while ( thread_rank + up < n ) { up <<= 1 ; ++count ; }
        }
      }

      thread_pool.m_thread_count = n ;

      /* Create threads last-to-first for start up fan-in barrier */

      for ( thread_rank = n ; ! result && 1 < thread_rank ; ) {
        Thread * const thread = thread_pool.m_thread + --thread_rank ;

        pthread_t pt ;

        if ( pthread_create( & pt, & attr, & local_driver, thread ) ) {
          thread->m_control = 0 ;
          result = TPI_ERROR_INTERNAL ;
        }
      }

      /* If a thread-spawn failed, terminate the created threads */

      if ( result ) {
        while ( thread_rank < --( thread_pool.m_thread_count ) ) {
          Thread * thread = thread_pool.m_thread + thread_pool.m_thread_count ;
          wait_thread( & thread->m_control , 1 ); /* Wait for blocking */
          thread->m_control = 1 ; /* Reactivate thread */
          wait_thread( & thread->m_control , 1 ); /* Wait for termination */
        }
        thread_pool.m_thread_count = 0 ;
      }

      pthread_attr_destroy( & attr );
    }
  }

  if ( ! result ) {
    local_barrier( thread_pool.m_thread );
    result = n ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  static int print_statistics = 0 ;

  int result ;

  result = NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {

    /* Wake up threads then wait for them to terminate */
    local_start( 0 , NULL , NULL , 0 ,
                 0 , NULL , NULL , 0 , NULL );

    local_wait();

    if ( print_statistics ) {
      int i = 0 ;
      for ( ; i < thread_pool.m_thread_count ; ++i ) {
        if ( thread_pool.m_thread[i].m_barrier_wait_count ) {
          const long mean =
            (long) ( ( thread_pool.m_thread[i].m_barrier_wait_total + 0.5 ) /
                       thread_pool.m_thread[i].m_barrier_wait_count );
          fprintf(stdout,"Thread[%d] barrier_wait( max %ld , mean %ld )\n", i ,
                   thread_pool.m_thread[i].m_barrier_wait_max , mean );
        }
      }
    }

    thread_pool.m_thread_count = 0 ;

    local_destroy_locks();

    if ( thread_pool.m_reduce_alloc ) {
      free( thread_pool.m_reduce_alloc );
      thread_pool.m_reduce_alloc = NULL ;
      thread_pool.m_reduce_alloc_size = 0 ;
    }
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static void local_block( TPI_Work * work )
{
  if ( work->rank ) {
    pthread_mutex_lock(   thread_pool.m_lock );
    pthread_mutex_unlock( thread_pool.m_lock );
  }
}

int TPI_Block()
{
  const int result =
    NULL != thread_pool.m_work_routine       ? TPI_ERROR_ACTIVE : (
    pthread_mutex_lock( thread_pool.m_lock ) ? TPI_ERROR_INTERNAL :

    local_start( thread_pool.m_thread_count ,
                 local_block , NULL ,
                 thread_pool.m_thread_count ,
                 0 /* lock_count */ ,
                 NULL , NULL , 0 , NULL ) );

  return result ;
}

int TPI_Unblock()
{
  const int result =
    local_block != thread_pool.m_work_routine  ? TPI_ERROR_ACTIVE : (
    pthread_mutex_unlock( thread_pool.m_lock ) ? TPI_ERROR_INTERNAL : 0 );

  if ( ! result ) { local_wait(); }

  return result ;
}

int TPI_Isblocked()
{
  return local_block == thread_pool.m_work_routine ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#else /* ! defined( HAVE_PTHREAD ) */

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Lock( int i )
{
  int result = i < 0 || thread_pool.m_lock_count <= i ? TPI_ERROR_SIZE : 0 ;

  if ( ! result ) {
    if ( 0 != thread_pool.m_lock[i] ) {
      result = TPI_ERROR_LOCK ;
    }
    else {
      thread_pool.m_lock[i] = 1 ;
    }
  }
  return result ;
}

int TPI_Unlock( int i )
{
  int result = i < 0 || thread_pool.m_lock_count <= i ? TPI_ERROR_SIZE : 0 ;

  if ( ! result ) {
    if ( 0 == thread_pool.m_lock[i] ) {
      result = TPI_ERROR_LOCK ;
    }
    else {
      thread_pool.m_lock[i] = 0 ;
    }
  }
  return result ;
}

static int local_set_lock_count( const int lock_count )
{
  int result = lock_count < 0 || LOCK_COUNT_MAX < lock_count
             ? TPI_ERROR_SIZE : 0 ;

  while ( thread_pool.m_lock_init < lock_count ) {

    thread_pool.m_lock[ thread_pool.m_lock_init ] = 0 ;

    ++( thread_pool.m_lock_init );
  }

  return result ;
}

/*--------------------------------------------------------------------*/

static int local_start(
  int                   work_thread_count ,
  TPI_work_subprogram   work_subprogram  ,
  const void *          work_info ,
  int                   work_count  ,
  int                   lock_count ,
  TPI_reduce_join       reduce_join ,
  TPI_reduce_init       reduce_init ,
  int                   reduce_size ,
  void *                reduce_data )
{
  const int result = lock_count ? local_set_lock_count( lock_count ) : 0 ;

  if ( ! result ) {
    thread_pool.m_work_routine     = work_subprogram ;
    thread_pool.m_work_info        = work_info ;
    thread_pool.m_work_count       = work_count ;
    thread_pool.m_lock_count       = lock_count ;
    thread_pool.m_thread->m_reduce = reduce_data ;
  }

  return result ;
}

static void local_wait()
{
  struct TPI_Work_Struct w = { NULL , NULL , 0 , 0 , 0 };

  w.info       = thread_pool.m_work_info ;
  w.count      = thread_pool.m_work_count ;
  w.lock_count = thread_pool.m_lock_count ;
  w.reduce     = thread_pool.m_thread->m_reduce ;

  for ( w.rank = 0 ; w.rank < w.count ; ++( w.rank ) ) {
    (* thread_pool.m_work_routine )( & w );
  }

  thread_pool.m_work_routine     = NULL ;
  thread_pool.m_work_info        = NULL ;
  thread_pool.m_work_count       = 0 ;
  thread_pool.m_lock_count       = 0 ;
  thread_pool.m_thread->m_reduce = NULL ;
}

/*--------------------------------------------------------------------*/

static void local_block( TPI_Work * work ) {}

int TPI_Block()
{
  const int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE :

    local_start( thread_pool.m_thread_count ,
                 local_block , NULL ,
                 thread_pool.m_thread_count ,
                 0 /* lock_count */ ,
                 NULL , NULL , 0 , NULL ) ;

  return result ;
}

int TPI_Unblock()
{
  const int result =
    local_block != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) { local_wait(); }

  return result ;
}

int TPI_Isblocked()
{
  return local_block == thread_pool.m_work_routine ;
}

/*--------------------------------------------------------------------*/

int TPI_Init( int n )
{
  int result = thread_pool.m_thread_count ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && ( n < 1 || THREAD_COUNT_MAX + 1 <= n ) ) {
    result = TPI_ERROR_SIZE ;
  }
  else {
    Thread * const thread = thread_pool.m_thread ;

    thread->m_thread_fan         = NULL ;
    thread->m_reduce             = NULL ;
    thread->m_rank               = 0 ;
    thread->m_barrier_wait_max   = 0 ;
    thread->m_barrier_wait_total = 0 ;
    thread->m_barrier_wait_count = 0 ;
    thread->m_control            = 1 ;

    thread_pool.m_thread_count = result = n ;

    /* Initialize one lock for blocking and unblocking */

    local_set_lock_count( 1 );
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  int result =  NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {
    thread_pool.m_thread_count = 0 ;
    thread_pool.m_lock_init = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#endif

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Wait()
{
  const int result =
    ( NULL        == thread_pool.m_work_routine ||
      local_block == thread_pool.m_work_routine ) ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) { local_wait(); }

  return result ;
}

int TPI_Start( TPI_work_subprogram work_subprogram  ,
               const void *        work_info ,
               int                 work_count ,
               int                 lock_count )
{
  const int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram            ? TPI_ERROR_NULL : (
    work_count  < 0                    ? TPI_ERROR_SIZE :
    local_start( thread_pool.m_thread_count - 1 ,
                 work_subprogram , work_info , work_count , lock_count ,
                 NULL , NULL , 0 , NULL ) ) );

  return result ;
}

int TPI_Run( TPI_work_subprogram work_subprogram  ,
             const void *        work_info ,
             int                 work_count ,
             int                 lock_count )
{
  const int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram            ? TPI_ERROR_NULL : (
    work_count  < 0                    ? TPI_ERROR_SIZE :
    local_start( thread_pool.m_thread_count ,
                 work_subprogram , work_info , work_count , lock_count ,
                 NULL , NULL , 0 , NULL ) ) );

  if ( ! result ) { local_wait(); }

  return result ;
}

int TPI_Run_threads( TPI_work_subprogram work_subprogram  ,
                     const void *        work_info ,
                     int                 lock_count  )
{
  const int work_count = 0 < thread_pool.m_thread_count ?
                             thread_pool.m_thread_count : 1 ;

  const int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram            ? TPI_ERROR_NULL : (
    local_start( thread_pool.m_thread_count ,
                 work_subprogram , work_info , work_count , lock_count ,
                 NULL , NULL , 0 , NULL ) ) );

  if ( ! result ) { local_wait(); }

  return result ;
}

int TPI_Start_threads( TPI_work_subprogram work_subprogram  ,
                       const void *        work_info ,
                       int                 lock_count  )
{
  const int work_count = 1 < thread_pool.m_thread_count ?
                             thread_pool.m_thread_count - 1 : 1 ;

  const int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram            ? TPI_ERROR_NULL : (
    local_start( thread_pool.m_thread_count - 1 ,
                 work_subprogram , work_info , work_count , lock_count ,
                 NULL , NULL , 0 , NULL ) ) );

  if ( ! result ) { local_wait(); }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Run_reduce( TPI_work_subprogram   work_subprogram  ,
                    const void *          work_info ,
                    int                   work_count  ,
                    TPI_reduce_join       reduce_join ,
                    TPI_reduce_init       reduce_init ,
                    int                   reduce_size ,
                    void *                reduce_data )
{
  const int lock_count = 0 ;

  const int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram            ? TPI_ERROR_NULL : (
    NULL == reduce_join                ? TPI_ERROR_NULL : (
    NULL == reduce_init                ? TPI_ERROR_NULL : (
    NULL == reduce_data                ? TPI_ERROR_NULL : (
    work_count  <= 0                   ? TPI_ERROR_SIZE : (
    reduce_size <= 0                   ? TPI_ERROR_SIZE : 

    local_start( thread_pool.m_thread_count ,
                 work_subprogram, work_info, work_count, lock_count,
                 reduce_join, reduce_init, reduce_size, reduce_data )))))));

  if ( ! result ) { local_wait(); }

  return result ;
}

int TPI_Run_threads_reduce( TPI_work_subprogram   work_subprogram  ,
                            const void *          work_info ,
                            TPI_reduce_join       reduce_join ,
                            TPI_reduce_init       reduce_init ,
                            int                   reduce_size ,
                            void *                reduce_data )
{
  const int lock_count = 0 ;
  const int work_count = 0 < thread_pool.m_thread_count ?
                             thread_pool.m_thread_count : 1 ;

  const int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram            ? TPI_ERROR_NULL : (
    NULL == reduce_join                ? TPI_ERROR_NULL : (
    NULL == reduce_init                ? TPI_ERROR_NULL : (
    NULL == reduce_data                ? TPI_ERROR_NULL : (
    reduce_size <= 0                   ? TPI_ERROR_SIZE : 

    local_start( thread_pool.m_thread_count ,
                 work_subprogram , work_info , work_count , lock_count ,
                 reduce_join, reduce_init, reduce_size, reduce_data ))))));

  if ( ! result ) { local_wait(); }

  return result ;
}

int TPI_Start_reduce( TPI_work_subprogram   work_subprogram  ,
                      const void *          work_info ,
                      int                   work_count  ,
                      TPI_reduce_join       reduce_join ,
                      TPI_reduce_init       reduce_init ,
                      int                   reduce_size ,
                      void *                reduce_data )
{
  const int lock_count = 0 ;

  const int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram            ? TPI_ERROR_NULL : (
    NULL == reduce_join                ? TPI_ERROR_NULL : (
    NULL == reduce_init                ? TPI_ERROR_NULL : (
    NULL == reduce_data                ? TPI_ERROR_NULL : (
    work_count  <= 0                   ? TPI_ERROR_SIZE : (
    reduce_size <= 0                   ? TPI_ERROR_SIZE : 

    local_start( thread_pool.m_thread_count - 1 ,
                 work_subprogram , work_info , work_count , lock_count ,
                 reduce_join, reduce_init, reduce_size, reduce_data )))))));

  return result ;
}

int TPI_Start_threads_reduce( TPI_work_subprogram   work_subprogram  ,
                              const void *          work_info ,
                              TPI_reduce_join       reduce_join ,
                              TPI_reduce_init       reduce_init ,
                              int                   reduce_size ,
                              void *                reduce_data )
{
  const int lock_count = 0 ;
  const int work_count = 1 < thread_pool.m_thread_count ?
                             thread_pool.m_thread_count - 1 : 1 ;

  const int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram            ? TPI_ERROR_NULL : (
    NULL == reduce_join                ? TPI_ERROR_NULL : (
    NULL == reduce_init                ? TPI_ERROR_NULL : (
    NULL == reduce_data                ? TPI_ERROR_NULL : (
    reduce_size <= 0                   ? TPI_ERROR_SIZE : 

    local_start( thread_pool.m_thread_count - 1 ,
                 work_subprogram , work_info , work_count , lock_count ,
                 reduce_join, reduce_init, reduce_size, reduce_data ))))));

  return result ;
}


