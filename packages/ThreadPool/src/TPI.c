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

#if	defined( __linux__ ) && \
	defined( __GNUC__ ) && \
	( 4 <= __GNUC__ ) && \
	! defined( __INTEL_COMPILER )

#define HAVE_ATOMIC_SYNC

#define THREADPOOL_CONFIG "PTHREAD SCHED_YIELD ATOMIC_SYNC"

#else

#define THREADPOOL_CONFIG "PTHREAD SCHED_YIELD"

#endif

#else /* ! defined( HAVE_PTHREAD ) */

#define THREADPOOL_CONFIG "NO THREADING"

#endif

const char * TPI_Version()
{
  static const char version_string[] =
    "TPI Version 1.1 , October 2009 , Configuration = " THREADPOOL_CONFIG ;

  return version_string ;
}

/*--------------------------------------------------------------------*/
/*----------- PTHREAD CONFIGURATION (END) ----------------------------*/
/*--------------------------------------------------------------------*/

enum { THREAD_COUNT_MAX = 256 };
enum { LOCK_COUNT_MAX   = 32 };

#if defined( HAVE_PTHREAD )

#include <errno.h>
#include <pthread.h>
#include <sched.h>

/*  Performance is heavily impacted by an
 *  atomic decrement of the work counter.
 *  Optimize this if at all possible.
 */

#if defined( HAVE_ATOMIC_SYNC )

#define atomic_fetch_and_decrement( VALUE_PTR )	\
	__sync_fetch_and_sub( VALUE_PTR , 1 )

#else

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
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------*/

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
  int                   m_work_count ;
  int                   m_work_count_claim ;

  Thread                m_thread[ THREAD_COUNT_MAX ];
  pthread_mutex_t       m_lock[ LOCK_COUNT_MAX ];
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
  /* m_work_count          */  0 ,
  /* m_work_count_claim    */  0
};

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

  if ( work.count <= thread_pool.m_thread_count ) {

    work.rank = this_thread->m_rank ;

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
    work.count      = thread_pool.m_thread_count ;
    work.rank       = this_thread->m_rank ;
    work.lock_count = 0 ;

    /* Initialize reduction value for non-root thread */

    if ( work.rank ) { (*thread_pool.m_reduce_init)( & work ); }

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
/* Activate child threads.  They are known to be blocked. */

static void local_activate()
{
  Thread * const thread_beg = thread_pool.m_thread + 1 ;
  Thread *       thread     = thread_pool.m_thread +
                              thread_pool.m_thread_count ;

  while ( thread_beg < thread ) { (--thread)->m_control = 1 ; }
}
 
/*--------------------------------------------------------------------*/
/*  The driver given to 'pthread_create'.
 *  Run work until told to terminate.
 */
static void * local_driver( void * arg )
{
  Thread * const this_thread = (Thread *) arg ;

  do {
    local_barrier( this_thread );   /*  Wait for my tree threads to block */

    this_thread->m_control = 0 ;

    /*  Block until I am activated. */
    wait_thread( & this_thread->m_control , 0 );

  } while ( thread_pool.m_work_routine );

  local_barrier( this_thread ); /* Termination barrier */

  this_thread->m_control = 0 ;

  return NULL ;
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
  int result = NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && 1 < thread_pool.m_thread_count ) {
    local_set_lock_count( 1 );

    thread_pool.m_work_routine = local_block ;
    thread_pool.m_work_count   = thread_pool.m_thread_count ;

    pthread_mutex_lock( thread_pool.m_lock );

    local_activate();
  }
  return result ;
}

int TPI_Unblock()
{
  int result = 0 ;

  if ( thread_pool.m_work_routine ) {

    result = local_block != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : 0 ;

    if ( ! result && 1 < thread_pool.m_thread_count ) {

      pthread_mutex_unlock( thread_pool.m_lock );

      local_barrier( thread_pool.m_thread );

      thread_pool.m_work_routine = NULL ;
      thread_pool.m_work_count   = 0 ;
    }
  }
  return result ;
}

int TPI_Isblocked()
{
  return local_block == thread_pool.m_work_routine ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Run( TPI_work_subprogram work_subprogram  ,
             const void *        work_info ,
             int                 work_count ,
             int                 lock_count )
{
  int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram            ? TPI_ERROR_NULL : (
    work_count  < 0                    ? TPI_ERROR_SIZE : 0 ) );

  if ( ! result && lock_count ) {
    result = local_set_lock_count( lock_count );
  }

  if ( ! result ) {

    if ( 1 < thread_pool.m_thread_count ) {

      thread_pool.m_work_routine     = work_subprogram ;
      thread_pool.m_work_info        = work_info ;
      thread_pool.m_work_count       = work_count ;
      thread_pool.m_work_count_claim = work_count ;
      thread_pool.m_lock_count       = lock_count ;

      /* Activate the blocked worker threads */
      local_activate();

      /* Work and then block the worker threads */
      local_barrier( thread_pool.m_thread );

      thread_pool.m_work_routine     = NULL ;
      thread_pool.m_work_info        = NULL ;
      thread_pool.m_work_count       = 0 ;
      thread_pool.m_work_count_claim = 0 ;
      thread_pool.m_lock_count       = 0 ;
    }
    else {
      struct TPI_Work_Struct w = { NULL , NULL , 0 , 0 , 0 };

      w.info       = work_info ;
      w.count      = work_count ;
      w.lock_count = lock_count ;

      for ( w.rank = 0 ; w.rank < w.count ; ++( w.rank ) ) {
        (* work_subprogram)( & w );
      }
    }
  }

  return result ;
}

int TPI_Run_threads( TPI_work_subprogram work_subprogram  ,
                     const void *        work_info ,
                     int                 lock_count  )
{
  const int count = thread_pool.m_thread_count ?
                    thread_pool.m_thread_count : 1 ;

  return TPI_Run( work_subprogram , work_info , count , lock_count );
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
    thread_pool.m_thread[0].m_reduce = NULL ;

    for ( i = 0 ; i < alloc_count ; ++i ) {
      thread_pool.m_thread[i+1].m_reduce =
        thread_pool.m_reduce_alloc + reduce_grain * i ;
    }
  }
}

int TPI_Run_reduce( TPI_work_subprogram   work_subprogram  ,
                    const void *          work_info ,
                    int                   work_count  ,
                    TPI_reduce_join       reduce_join ,
                    TPI_reduce_init       reduce_init ,
                    int                   reduce_size ,
                    void *                reduce_data )
{
  int result =
    NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : (
    NULL == work_subprogram            ? TPI_ERROR_NULL : (
    NULL == reduce_join                ? TPI_ERROR_NULL : (
    NULL == reduce_init                ? TPI_ERROR_NULL : (
    NULL == reduce_data                ? TPI_ERROR_NULL : (
    work_count  <= 0                   ? TPI_ERROR_SIZE : (
    reduce_size <= 0                   ? TPI_ERROR_SIZE : 0 ) ) ) ) ) );

  if ( ! result ) {

    if ( 1 < thread_pool.m_thread_count ) {

      alloc_reduce( reduce_size );

      thread_pool.m_reduce_join      = reduce_join ;
      thread_pool.m_reduce_init      = reduce_init ;
      thread_pool.m_work_routine     = work_subprogram ;
      thread_pool.m_work_info        = work_info ;
      thread_pool.m_work_count       = work_count ;
      thread_pool.m_work_count_claim = work_count ;

      {
        Thread * const root_thread = thread_pool.m_thread ;

        root_thread->m_reduce = reduce_data ;

        /* Activate the blocked worker threads */
        local_activate();

        /* Work and then block the worker threads */
        local_barrier( root_thread );

        root_thread->m_reduce = NULL ;
      }

      thread_pool.m_reduce_join      = NULL ;
      thread_pool.m_reduce_init      = NULL ;
      thread_pool.m_work_routine     = NULL ;
      thread_pool.m_work_info        = NULL ;
      thread_pool.m_work_count       = 0 ;
      thread_pool.m_work_count_claim = 0 ;
    }
    else {
      struct TPI_Work_Struct w = { NULL , NULL , 0 , 0 , 0 };

      w.info   = work_info ;
      w.reduce = reduce_data ;
      w.count  = work_count ;

      for ( w.rank = 0 ; w.rank < w.count ; ++( w.rank ) ) {
        (* work_subprogram)( & w );
      }
    }
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
  const int count = thread_pool.m_thread_count ?
                    thread_pool.m_thread_count : 1 ;

  return TPI_Run_reduce( work_subprogram , work_info , count ,
                         reduce_join , reduce_init ,
                         reduce_size , reduce_data );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Init( int n )
{
  static int first_pass = 1 ;
  int result = thread_pool.m_thread_count ? TPI_ERROR_ACTIVE : 0 ;

  if ( first_pass ) {
    fprintf( stdout , TPI_Version() );
    first_pass = 0 ;
  }

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

      /* Initialize threads with fan-in / fan-out span of threads */

      int count = 1 ;
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
/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  static int print_statistics = 0 ;

  int result ;

  result = NULL != thread_pool.m_work_routine ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {

    /* Wake up threads then wait for them to terminate */
    local_activate();
    local_barrier( thread_pool.m_thread );

    if ( print_statistics ) {
      int i = 0 ;
      for ( ; i < thread_pool.m_thread_count ; ++i ) {
        if ( thread_pool.m_thread[i].m_barrier_wait_count ) {
          long mean = ( thread_pool.m_thread[i].m_barrier_wait_total + 0.5 ) /
                        thread_pool.m_thread[i].m_barrier_wait_count ;
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

int TPI_Block()     { work->m_active ? TPI_ERROR_ACTIVE : 0 ; }
int TPI_Unblock()   { work->m_active ? TPI_ERROR_ACTIVE : 0 ; }
int TPI_Isblocked() { work->m_active ? TPI_ERROR_ACTIVE : 0 ; }

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
      struct TPI_Work_Struct w = { NULL , NULL , 0 , 0 , 0 };

      work->m_active     = 1 ;
      work->m_lock_count = lock_count ;

      w.info       = work_info ;
      w.count      = work_count ;
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
    struct TPI_Work_Struct w = { NULL , NULL , 0 , 0 , 0 };

    work->m_active = 1 ;

    w.info   = work_info ;
    w.reduce = reduce_data ;
    w.count  = work_count ;

    for ( w.rank = 0 ; w.rank < w.count ; ++( w.rank ) ) {
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


