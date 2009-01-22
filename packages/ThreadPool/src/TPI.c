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

#ifndef HAVE_PTHREAD
#define HAVE_PTHREAD 1
#endif

#include <TPI.h>

#include <stdlib.h>

enum { MAXIMUM_LOCK_COUNT = 256 };

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#if HAVE_PTHREAD

#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include <sys/types.h>

struct ThreadPool_Data ;

typedef struct TPI_ThreadData_Struct {
  pthread_cond_t           m_run_cond ;
  pthread_mutex_t          m_run_lock ;
  pthread_mutex_t          m_work_lock ;
  struct ThreadPool_Data * m_thread_pool ;
  int                      m_thread_rank ;
  int                      m_work_begin ;
  int                      m_work_count ;
  int                      m_work_actual ;
} ThreadData ;

typedef struct ThreadPool_Data {
  ThreadData * const    m_thread_data ;
  int                   m_thread_count ;
  TPI_work_subprogram   m_work_routine ;
  void                * m_work_argument ;
  int                   m_work_count_total ;
  pthread_mutex_t     * m_lock ;
  int                   m_lock_count ;
} ThreadPool ;

/*--------------------------------------------------------------------*/

enum { MAXIMUM_THREAD_COUNT = 256 };

static ThreadPool * local_thread_pool()
{
  static ThreadData data_pool[ MAXIMUM_THREAD_COUNT ];

  static ThreadPool thread_pool = {
    /* m_thread_data      */  data_pool ,
    /* m_thread_count     */  0 ,
    /* m_work_routine     */  NULL ,
    /* m_work_argument    */  NULL ,
    /* m_work_count_total */  0 ,
    /* m_lock             */  NULL ,
    /* m_lock_count       */  0 };

  return & thread_pool ;
}

/*--------------------------------------------------------------------*/

int TPI_Lock( int i )
{
  ThreadPool * const thread_pool = local_thread_pool();
  const int result = 
    ( NULL == thread_pool                           ? TPI_ERROR_NULL :
    ( i < 0 || thread_pool->m_lock_count <= i       ? TPI_ERROR_SIZE :
    ( pthread_mutex_lock( thread_pool->m_lock + i ) ? TPI_ERROR_LOCK : 0 ) ) );
  return result ;
}

int TPI_Trylock( int i )
{
  ThreadPool * const thread_pool = local_thread_pool();
  int p ;
  const int result = 
    ( NULL == thread_pool                     ? TPI_ERROR_NULL :
    ( i < 0 || thread_pool->m_lock_count <= i ? TPI_ERROR_SIZE :
    ( EBUSY == ( p = pthread_mutex_trylock(thread_pool->m_lock+i) )
               ? TPI_LOCK_BUSY : ( p ? TPI_ERROR_LOCK : 0 ) ) ) );
  return result ;
}

int TPI_Unlock( int i )
{
  ThreadPool * const thread_pool = local_thread_pool();
  const int result = 
    ( NULL == thread_pool                           ? TPI_ERROR_NULL :
    ( i < 0 || thread_pool->m_lock_count <= i       ? TPI_ERROR_SIZE :
    ( pthread_mutex_unlock(thread_pool->m_lock + i) ? TPI_ERROR_LOCK : 0 ) ) );
  return result ;
}

/*--------------------------------------------------------------------*/
/*  Steal work to replenish my queue. */

static void steal_work( ThreadData * my_pool )
{
  ThreadPool * const thread_pool = my_pool->m_thread_pool ;
  ThreadData * const all_pool    = thread_pool->m_thread_data ;
  const int my_thread_rank       = my_pool->m_thread_rank ;
  const int thread_count =
    thread_pool->m_thread_count < thread_pool->m_work_count_total ?
    thread_pool->m_thread_count : thread_pool->m_work_count_total ;

  int guess_max_work ;
  int guess_total_work ;
  int my_work_count ;

  do {
    ThreadData * steal_pool = NULL ;
    int p ;

    guess_max_work = 0 ;
    guess_total_work = 0 ;
    my_work_count = 0 ;

    /*  Quickly, but not guaranteed accurate due to not locking,
     *  get a count of remaining total work and who has the maximum.
     */
    for ( p = 1 ; p < thread_count ; ++p ) {
      ThreadData * const try_pool =
        all_pool + ( my_thread_rank + p ) % thread_count ;

      /* The following query is not thread safe, but only need a quick guess. */

      const int thread_work = try_pool->m_work_count ;

      /* Keep up with guessed total work to determine how much work to steal */

      guess_total_work += thread_work ;

      if ( guess_max_work < thread_work ) {
        guess_max_work = thread_work ;
        steal_pool     = try_pool ;
      }
    }

    /*  If a task pool (might) have work remaining steal some. */

    if ( guess_max_work ) {
      /*  There is work available to steal.
       *  How much work to steal ?
       *  At least my fraction of the total remaining work, rounded up.
       */
      my_work_count = ( guess_total_work + thread_count - 1 ) / thread_count ;

      pthread_mutex_lock( & steal_pool->m_work_lock );

      {
        const int max_steal = ( steal_pool->m_work_count + 1 ) / 2 ;

        if ( max_steal < my_work_count ) { my_work_count = max_steal ; }
      }

      if ( my_work_count ) {

        my_pool   ->m_work_count =  my_work_count ;
        steal_pool->m_work_count -= my_work_count ;

        my_pool->m_work_begin = steal_pool->m_work_begin +
                                steal_pool->m_work_count ;
      }
      pthread_mutex_unlock( & steal_pool->m_work_lock );
    }

    /*  If ( guess_max_work < guess_total_work ) then
     *    There might be more work to steal.
     *  If ( 0 == my_work_count ) then
     *    The steal failed because another thread 
     *    completed or stole the work first.
     */
  } while ( ! my_work_count && guess_total_work );
}

/*--------------------------------------------------------------------*/
/*  Run the task queue.
 *  A worker runs until commanded to terminate.
 *  The control runs until the current queue is empty.
 */
static void local_run_task( ThreadData * const my_pool )
{
  pthread_mutex_t * const my_work_lock = & my_pool->m_work_lock ;
  ThreadPool      * const thread_pool  =   my_pool->m_thread_pool ;

  struct TPI_Work_Struct work =
    { thread_pool->m_work_argument ,
      thread_pool->m_lock_count ,
      thread_pool->m_work_count_total ,
      0 };

  int work_count ;

  do {
    pthread_mutex_lock( my_work_lock );

    if ( ! my_pool->m_work_count ) { steal_work( my_pool ); }

    work_count     = my_pool->m_work_count ;
    work.work_rank = my_pool->m_work_begin + work_count - 1 ;

    if ( work_count ) {
      --( my_pool->m_work_count );
      ++( my_pool->m_work_actual );
    }

    pthread_mutex_unlock( my_work_lock );

    if ( work_count ) {
      (* thread_pool->m_work_routine)( & work );
    }
  } while ( work_count );

  return ;
}

/*--------------------------------------------------------------------*/
/*  The driver given to 'pthread_create'.
 *  Run work until told to terminate.
 */

static void * local_thread_pool_driver( void * arg )
{
  ThreadData * const my_pool     = (ThreadData*) arg ;
  ThreadPool * const thread_pool = my_pool->m_thread_pool ;
  ThreadData * const root        = thread_pool->m_thread_data ;

  /*------------------------------*/
  /* Acquire my run lock and signal main thread that I started */

  pthread_mutex_lock( & my_pool->m_run_lock );

  pthread_mutex_lock(   & root->m_run_lock );
  pthread_cond_signal(  & root->m_run_cond );
  pthread_mutex_unlock( & root->m_run_lock );

  /*------------------------------*/
  /*  While I'm active wait for run signal.
   *  I give up my run lock while I am waiting.
   */

  while ( my_pool->m_thread_rank ) {
    pthread_cond_wait( & my_pool->m_run_cond , & my_pool->m_run_lock );
    if ( my_pool->m_thread_rank ) { local_run_task( my_pool ); }
  }
  pthread_mutex_unlock( & my_pool->m_run_lock );

  /*------------------------------*/
  /* Termination, the main thread is waiting for the last thread to signal */

  pthread_mutex_lock( & root->m_run_lock );
  if ( ! --( thread_pool->m_thread_count ) ) {
    pthread_cond_signal( & root->m_run_cond );
  }
  pthread_mutex_unlock( & root->m_run_lock );

  return NULL ;
}

/*--------------------------------------------------------------------*/

int TPI_Run( TPI_work_subprogram subprogram  ,
             void * const        shared_data ,
             const int           work_count  ,
             const int           lock_count  )
{
  ThreadPool * const thread_pool = local_thread_pool();
  ThreadData * const all_pool    = thread_pool->m_thread_data ;

  int result = 0 ;

  if ( ! result && ( NULL == subprogram || NULL == all_pool->m_thread_pool ) ) {
    result = TPI_ERROR_NULL ;
  }

  if ( ! result && ( lock_count < 0 || MAXIMUM_LOCK_COUNT < lock_count ) ) {
    result = TPI_ERROR_SIZE ;
  }

  if ( ! result && NULL == thread_pool ) {
    result = TPI_ERROR_ACTIVE ;
  }

  if ( ! result ) {
    pthread_mutex_t locks[ lock_count ];

    /* Don't bother activating more threads than their is work */
    const int np = thread_pool->m_thread_count < work_count ?
                   thread_pool->m_thread_count : work_count ;
    int p ;

    int lock_count_init = 0 ;

    while ( ! result && lock_count_init < lock_count ) {
      if ( pthread_mutex_init( locks + lock_count_init , NULL ) ) {
        result = TPI_ERROR_INTERNAL ;
      }
      else {
        ++lock_count_init ;
      }
    }
    if ( ! result && lock_count_init ) {
      thread_pool->m_lock = locks ;
      thread_pool->m_lock_count = lock_count_init ;
    }

    /*  I have all of the threads' run locks, safe to
     *  set the initial work for each thread.
     */
    {
      int work_offset = 0 ;

      for ( p = 0 ; p < np ; ) {
        ThreadData * const pool = all_pool + p ; ++p ;

        pool->m_work_begin = work_offset ;

        work_offset = ( work_count * p ) / np ;

        pool->m_work_count = work_offset - pool->m_work_begin ;
        pool->m_work_actual = 0 ;
      }
    }

    thread_pool->m_work_argument    = shared_data ;
    thread_pool->m_work_routine     = subprogram ;
    thread_pool->m_work_count_total = work_count ;

    /* Signal and unlock all of the threads' run locks */
    for ( p = 1 ; p < np ; ++p ) {
      ThreadData * const pool = all_pool + p ;
      pthread_cond_signal(  & pool->m_run_cond );
      pthread_mutex_unlock( & pool->m_run_lock );
    }

    local_run_task( all_pool ); /* Participate in the work */

    /* All work has been claimed, wait to reclaim threads' run locks */

    for ( p = 1 ; p < np ; ++p ) {
      pthread_mutex_lock( & all_pool[p].m_run_lock );
    }

    thread_pool->m_work_count_total    = 0 ;
    thread_pool->m_work_routine  = NULL ;
    thread_pool->m_work_argument = NULL ;

    while ( lock_count_init ) {
      pthread_mutex_destroy( locks + lock_count_init );
      --lock_count_init ;
    }

    /* Compute level of parallelism that has occured */
    if ( ! result ) {
      int max = 0 ;
      for ( p = 0 ; p < np ; ++p ) {
        ThreadData * const pool = all_pool + p ;
        if ( max < pool->m_work_actual ) { max = pool->m_work_actual ; }
      }

      result = ( work_count + max - 1 ) / max ;
    }
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Init( int n )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result = ! thread_pool ||
                 thread_pool->m_work_count_total ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && ( n <= 0 || MAXIMUM_THREAD_COUNT < n ) ) {
    result = TPI_ERROR_SIZE ;
  }

  if ( ! result ) {
    pthread_attr_t thread_attr ;

    if ( pthread_attr_init( & thread_attr ) ) {
      result = TPI_ERROR_INTERNAL ;
    }
    else {
      ThreadData * const task_pool = thread_pool->m_thread_data ;

      int p ;
      for ( p = 0 ; ! result && p < n ; ++p ) {
        ThreadData * task = task_pool + p ;
        task->m_thread_pool = thread_pool ;
        task->m_thread_rank = p ;
        task->m_work_begin = 0 ;
        task->m_work_count = 0 ;

        if ( pthread_cond_init(  & task->m_run_cond , NULL ) ||
             pthread_mutex_init( & task->m_run_lock , NULL ) ||
             pthread_mutex_init( & task->m_work_lock , NULL ) ) {
          result = TPI_ERROR_INTERNAL ;
        }
      }

      pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
      pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

      pthread_mutex_lock( & task_pool->m_run_lock );

      thread_pool->m_thread_count = 1 ;

      for ( p = 1 ; p < n && ! result ; ) {
        ThreadData * task = thread_pool->m_thread_data + p ;
        pthread_t pt ;

        if ( pthread_create( & pt, & thread_attr,
                             & local_thread_pool_driver, task ) ) {
          result = TPI_ERROR_INTERNAL ;
        }
        else {
          /* Wait for start and claim its run lock. */
          thread_pool->m_thread_count = ++p ;
          pthread_cond_wait( & task_pool->m_run_cond ,
                             & task_pool->m_run_lock );
          pthread_mutex_lock( & task->m_run_lock );
        }
      }

      pthread_attr_destroy( & thread_attr );

      pthread_mutex_unlock( & task_pool->m_run_lock );
    }

    if ( result ) { TPI_Finalize(); }
  }

  if ( ! result ) { result = thread_pool->m_thread_count ; }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result = ! thread_pool ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {
    ThreadData * const task_pool = thread_pool->m_thread_data ;

    const int np = thread_pool->m_thread_count ;
    int p ;

    pthread_mutex_lock( & task_pool->m_run_lock );

    --( thread_pool->m_thread_count );

    /* Clear active flag and activate threads */

    for ( p = 1 ; p < np ; ++p ) {
      ThreadData * task = task_pool + p ;
      task->m_thread_rank = 0 ;
      pthread_cond_signal(  & task->m_run_cond );
      pthread_mutex_unlock( & task->m_run_lock );
    }

    if ( 0 < thread_pool->m_thread_count ) {
      /* Wait for last thread to signal termination */
      pthread_cond_wait( & task_pool->m_run_cond ,
                         & task_pool->m_run_lock );
    }

    pthread_mutex_unlock( & task_pool->m_run_lock );

    /* Destroy the task pool mutexes and conditions */

    for ( p = 0 ; p < np ; ++p ) {
      ThreadData * task = task_pool + p ;
      task->m_thread_pool = NULL ;
      task->m_thread_rank = 0 ;
      task->m_work_begin = 0 ;
      task->m_work_count = 0 ;

      pthread_mutex_destroy( & task->m_work_lock );
      pthread_mutex_destroy( & task->m_run_lock );
      pthread_cond_destroy(  & task->m_run_cond );
    }
  }

  return result ;
}

/* #if HAVE_PTHREAD */
/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#else

typedef struct LocalWorkData {
  int * m_lock ;
  int   m_lock_count ;
  int   m_active ;
} LocalWork ;

static LocalWork * local_work()
{
  static struct LocalWorkData data = { NULL , 0 , 0 };
  return & data ;
}

int TPI_Lock( int i )
{
  LocalWork * const work = local_work();

  const int result = 
    ( NULL == work                     ? TPI_ERROR_NULL :
    ( i < 0 || work->m_lock_count <= i ? TPI_ERROR_SIZE :
    ( work->m_lock[i] )                ? TPI_ERROR_LOCK : 0 ) );
  if ( ! result ) { work->m_lock[i] = 1 ; }
  return result ;
}

int TPI_Trylock( int i )
{
  LocalWork * const work = local_work();

  const int result = 
    ( NULL == work                     ? TPI_ERROR_NULL :
    ( i < 0 || work->m_lock_count <= i ? TPI_ERROR_SIZE :
    ( work->m_lock[i] )                ? TPI_LOCK_BUSY : 0 ) );
  if ( ! result ) { work->m_lock[i] = 1 ; }
  return result ;
}

int TPI_Unlock( int i )
{
  LocalWork * const work = local_work();

  const int result = 
    ( NULL == work                     ? TPI_ERROR_NULL :
    ( i < 0 || work->m_lock_count <= i ? TPI_ERROR_SIZE :
    ( ! work->m_lock[i] )              ? TPI_ERROR_LOCK : 0 ) );
  if ( ! result ) { work->m_lock[i] = 0 ; }
}

int TPI_Run( TPI_work_subprogram subprogram  ,
             void * const        shared_data ,
             const int           work_count  ,
             const int           lock_count  )
{
  LocalWork * const work = local_work();
  int i ;
  int result = 0 ;

  if ( ! work || work->m_active ) { result = TPI_ERROR_ACTIVE ; }

  if ( ! result && MAXIMUM_LOCK_COUNT <= lock_count ) {
    result = TPI_ERROR_SIZE ;
  }

  if ( ! result ) {
    int locks[ lock_count ];

    work->m_active     = 1 ;
    work->m_lock_count = lock_count ;
    work->m_locks      = locks ;

    for ( i = 0 ; i < lock_count ; ++i ) { locks[i] = 0 ; }

    for ( i = 0 ; i < work_count ; ++i ) {
      TPI_Work work = { shared_data , lock_count , work_count , i };
      (* thread_pool->m_work_routine)( & work );
    }

    work->m_active     = 0 ;
    work->m_lock_count = 0 ;
    work->m_locks      = NULL ;

    result = 1 ;
  }

  return result ;
}

int TPI_Init( int )
{
  LocalWork * const work = local_work();
  int result = 0 ;

  if ( ! work || work->m_active ) { result = TPI_ERROR_ACTIVE ; }

  if ( ! result ) { result = 1 ; }

  return result ;
}

int TPI_Finalize()
{
  LocalWork * const work = local_work();
  int result = 0 ;

  if ( ! work || work->m_active ) { result = TPI_ERROR_ACTIVE ; }

  return result ;
}

}

#endif

