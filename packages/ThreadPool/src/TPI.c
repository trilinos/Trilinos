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

#include <TPI.h>
#include <ThreadPool_config.h>

enum { MAXIMUM_LOCK_COUNT = 256 };

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#ifdef HAVE_PTHREAD

#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include <sys/types.h>

struct ThreadPool_Data ;

typedef void (*Thread_Work)( struct ThreadPool_Data * const , unsigned int );

struct Thread_Data {
  pthread_cond_t       m_run_cond ;
  pthread_mutex_t      m_run_lock ;
  struct Thread_Data * m_next ;
  unsigned int         m_number ;
};

typedef struct ThreadPool_Data {
  pthread_mutex_t       m_work_lock ;
  pthread_cond_t        m_work_cond ;
  pthread_mutex_t       m_master_lock ;
  pthread_cond_t        m_master_cond ;
  Thread_Work           m_thread_work ;
  struct Thread_Data  * m_thread_next ;
  int                   m_thread_count ;
  pthread_mutex_t     * m_lock ;
  int                   m_lock_init ;
  int                   m_lock_count ;
  TPI_work_subprogram   m_work_routine ;
  void                * m_work_argument ;
  int                   m_work_count_total ;
  int                   m_work_count_claim ;
  int                   m_work_count_pending ;
} ThreadPool ;

/*--------------------------------------------------------------------*/

static ThreadPool * local_thread_pool()
{
  static pthread_mutex_t lock_pool[ MAXIMUM_LOCK_COUNT ];

  static ThreadPool thread_pool = {
    /* m_work_lock           */  PTHREAD_MUTEX_INITIALIZER ,
    /* m_work_cond           */  PTHREAD_COND_INITIALIZER ,
    /* m_master_lock         */  PTHREAD_MUTEX_INITIALIZER ,
    /* m_master_cond         */  PTHREAD_COND_INITIALIZER ,
    /* m_thread_work         */  NULL ,
    /* m_thread_next         */  NULL ,
    /* m_thread_count        */  0 ,
    /* m_lock                */  lock_pool ,
    /* m_lock_init           */  0 ,
    /* m_lock_count          */  0 ,
    /* m_work_routine        */  NULL ,
    /* m_work_argument       */  NULL ,
    /* m_work_count_total    */  0 ,
    /* m_work_count_claim    */  0 ,
    /* m_work_count_pending  */  0 };

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
/*  Run the work queue until it is empty.
 *  The input 'thread_number' is deliberately ignored.
 */
static void local_run_work( ThreadPool * const thread_pool ,
                            unsigned int work_item )
{
  work_item = 0 ;

  do {

    /* Have a work item to do */

    if ( work_item ) {
      struct TPI_Work_Struct work ;

      work.shared     = thread_pool->m_work_argument ;
      work.lock_count = thread_pool->m_lock_count ;
      work.work_count = thread_pool->m_work_count_total ;
      work.work_rank  = thread_pool->m_work_count_total - work_item ;

      (* thread_pool->m_work_routine)( & work );
    }

    /* Finished my work item, update the work queue */

    pthread_mutex_lock( & thread_pool->m_work_lock );

    if ( work_item ) {
      /* Record work just performed */
      if ( ! --( thread_pool->m_work_count_pending ) ) {
        /* Signal master thread I completed the last work item */
        pthread_cond_signal( & thread_pool->m_work_cond );
      }
    }

    /* Claim more work, if any */
    work_item = thread_pool->m_work_count_claim ;

    if ( work_item ) {
      thread_pool->m_work_count_claim = work_item - 1 ;
    }

    pthread_mutex_unlock( & thread_pool->m_work_lock );

  } while ( work_item );
}

/*--------------------------------------------------------------------*/
/*  Run the work subprogram exactly once for each thread.
 */
static void local_run_thread( ThreadPool * const thread_pool ,
                              unsigned int thread_number )
{
  struct TPI_Work_Struct work ;

  work.shared     = thread_pool->m_work_argument ;
  work.lock_count = thread_pool->m_lock_count ;
  work.work_count = thread_pool->m_thread_count ;
  work.work_rank  = thread_pool->m_thread_count - thread_number ;

  (* thread_pool->m_work_routine)( & work );

  /* Finished my work item, update the work queue */

  pthread_mutex_lock( & thread_pool->m_work_lock );

  /* Record work just performed */
  if ( ! --( thread_pool->m_work_count_pending ) ) {
    /* Signal master thread I completed the last work item */
    pthread_cond_signal( & thread_pool->m_work_cond );
  }

  pthread_mutex_unlock( & thread_pool->m_work_lock );
}

/*--------------------------------------------------------------------*/
/*  The driver given to 'pthread_create'.
 *  Run work until told to terminate.
 */

static void * local_thread_pool_driver( void * arg )
{
  ThreadPool * const thread_pool = (ThreadPool *) arg ;

  struct Thread_Data my_data = {
    /* m_run_cond */  PTHREAD_COND_INITIALIZER ,
    /* m_run_lock */  PTHREAD_MUTEX_INITIALIZER ,
    /* m_next     */  NULL ,
    /* m_number   */  0 };

  /*------------------------------*/
  /* Signal master thread that I have started */

  pthread_mutex_lock( & my_data.m_run_lock );
  pthread_mutex_lock( & thread_pool->m_master_lock );

  my_data.m_next   = thread_pool->m_thread_next ;
  my_data.m_number = my_data.m_next ? my_data.m_next->m_number - 1
                                    : thread_pool->m_thread_count - 1 ;

  thread_pool->m_thread_next = & my_data ;

  pthread_cond_signal(  & thread_pool->m_master_cond );
  pthread_mutex_unlock( & thread_pool->m_master_lock );

  /*------------------------------*/
  /*  While I'm active wait for run signal.
   *  I give up my run lock while I am waiting.
   */
  do {
    pthread_cond_wait( & my_data.m_run_cond , & my_data.m_run_lock );

    if ( my_data.m_next &&
        ! pthread_mutex_trylock( & my_data.m_next->m_run_lock ) ) {
      pthread_mutex_unlock( & my_data.m_next->m_run_lock );
      pthread_cond_signal(  & my_data.m_next->m_run_cond );
    }

    if ( thread_pool->m_thread_work ) {
      (*thread_pool->m_thread_work)( thread_pool , my_data.m_number );
    }
  } while ( thread_pool->m_thread_work );

  /*------------------------------*/
  /* Termination, the main thread is waiting for the last thread to signal */

  pthread_mutex_unlock(  & my_data.m_run_lock );
  pthread_mutex_destroy( & my_data.m_run_lock );
  pthread_cond_destroy(  & my_data.m_run_cond );

  pthread_mutex_lock(   & thread_pool->m_master_lock );
  pthread_cond_signal(  & thread_pool->m_master_cond );
  pthread_mutex_unlock( & thread_pool->m_master_lock );

  return NULL ;
}

/*--------------------------------------------------------------------*/

static int set_lock_count( ThreadPool * const thread_pool ,
                           const int lock_count )
{
  int result = 0 ;

  if ( lock_count < 0 || MAXIMUM_LOCK_COUNT < lock_count ) {
    result = TPI_ERROR_SIZE ;
  }

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
    struct Thread_Data * const worker = thread_pool->m_thread_next ;

    if ( ! worker ) {
      /* No worker threads */
      result = set_lock_count( thread_pool , lock_count );

      if ( ! result ) {
        struct TPI_Work_Struct work ;

        work.shared     = shared_data ;
        work.lock_count = lock_count ;
        work.work_count = work_count ;
        work.work_rank  = 0 ;

        for ( ; work.work_rank < work.work_count ; ++( work.work_rank ) ) {
          (* subprogram)( & work );
        }
      }
    }
    else {
      /* Have worker threads */

      pthread_mutex_lock( & thread_pool->m_work_lock );

      result = set_lock_count( thread_pool , lock_count );

      if ( ! result ) {

        thread_pool->m_lock_count         = lock_count ;
        thread_pool->m_work_argument      = shared_data ;
        thread_pool->m_work_routine       = subprogram ;
        thread_pool->m_work_count_total   = work_count ;
        thread_pool->m_work_count_pending = work_count ;
        thread_pool->m_work_count_claim   = work_count ;

        pthread_mutex_unlock( & thread_pool->m_work_lock );

        /* If first worker waiting for a signal then signal it.
         * This signal will be passed along thread-to-thread.
         */
        if ( ! pthread_mutex_trylock( & worker->m_run_lock ) ) {
          pthread_mutex_unlock( & worker->m_run_lock );
          pthread_cond_signal(  & worker->m_run_cond );
        }

        local_run_work( thread_pool , thread_pool->m_thread_count );

        pthread_mutex_lock( & thread_pool->m_work_lock );

        if ( thread_pool->m_work_count_pending ) {
          /* Wait for worker threads to complete last claimed work item */
          pthread_cond_wait( & thread_pool->m_work_cond ,
                             & thread_pool->m_work_lock );
        }

        thread_pool->m_work_count_claim   = 0 ;
        thread_pool->m_work_count_pending = 0 ;
        thread_pool->m_work_count_total   = 0 ;
        thread_pool->m_work_routine       = NULL ;
        thread_pool->m_work_argument      = NULL ;
        thread_pool->m_lock_count         = 0 ;
      }

      pthread_mutex_unlock( & thread_pool->m_work_lock );
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

    struct Thread_Data * w ;

    /* Lock all worker threads */

    for ( w = thread_pool->m_thread_next ; w ; w = w->m_next ) {
      pthread_mutex_lock( & w->m_run_lock );
    }

    result = set_lock_count( thread_pool , lock_count );

    if ( ! result ) {
      /* Use the run_threads driver */

      thread_pool->m_thread_work        = & local_run_thread ;
      thread_pool->m_lock_count         = lock_count ;
      thread_pool->m_work_argument      = shared_data ;
      thread_pool->m_work_routine       = subprogram ;
      thread_pool->m_work_count_pending = thread_pool->m_thread_count ;

      /* Unlock all threads */
      for ( w = thread_pool->m_thread_next ; w ; w = w->m_next ) {
        pthread_mutex_unlock( & w->m_run_lock );
      }

      /* Signal the first thread to go.
       * This signal will be passed along thread-to-thread.
       */
      if ( thread_pool->m_thread_next ) {
        pthread_cond_signal( & thread_pool->m_thread_next->m_run_cond );
      }

      local_run_thread( thread_pool , thread_pool->m_thread_count );

      pthread_mutex_lock( & thread_pool->m_work_lock );

      if ( thread_pool->m_work_count_pending ) {
        /* Wait for worker threads to complete last claimed work item */
        pthread_cond_wait( & thread_pool->m_work_cond ,
                           & thread_pool->m_work_lock );
      }

      pthread_mutex_unlock( & thread_pool->m_work_lock );

      /* Lock all threads, they must complete to be locked. */

      for ( w = thread_pool->m_thread_next ; w ; w = w->m_next ) {
        pthread_mutex_lock( & w->m_run_lock );
      }

      /* Change the driver back to run-work */

      thread_pool->m_work_routine     = NULL ;
      thread_pool->m_work_argument    = NULL ;
      thread_pool->m_lock_count       = 0 ;
      thread_pool->m_thread_work      = & local_run_work ;
    }

    /* Unlock all worker threads */
    for ( w = thread_pool->m_thread_next ; w ; w = w->m_next ) {
      pthread_mutex_unlock( & w->m_run_lock );
    }
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Init( int n )
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result = thread_pool->m_thread_count ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {
    pthread_attr_t thread_attr ;

    if ( pthread_attr_init( & thread_attr ) ) {
      result = TPI_ERROR_INTERNAL ;
    }
    else {
      int p ;

      pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
      pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

      thread_pool->m_thread_count = n ;
      thread_pool->m_thread_work  = & local_run_work ;

      pthread_mutex_lock( & thread_pool->m_master_lock );

      for ( p = 1 ; p < n && ! result ; ++p ) {
        pthread_t pt ;

        if ( pthread_create( & pt, & thread_attr,
                             & local_thread_pool_driver, thread_pool ) ) {
          result = TPI_ERROR_INTERNAL ;
        }
        else {
          /* Wait for start */
          pthread_cond_wait( & thread_pool->m_master_cond ,
                             & thread_pool->m_master_lock );
        }
      }

      pthread_attr_destroy( & thread_attr );
    }

    if ( result ) { TPI_Finalize(); }
  }

  if ( ! result ) { result = n ; }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  ThreadPool * const thread_pool = local_thread_pool();

  int result = NULL != thread_pool->m_work_routine ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {

    struct Thread_Data * t ;

    for ( t = thread_pool->m_thread_next ; NULL != t ; t = t->m_next ) {
      pthread_mutex_lock( & t->m_run_lock );
    }

    thread_pool->m_thread_work = NULL ;

    for ( t = thread_pool->m_thread_next ; NULL != t ; ) {
      struct Thread_Data * const t_now = t ; t = t->m_next ;
      pthread_mutex_unlock( & t_now->m_run_lock );
      pthread_cond_signal(  & t_now->m_run_cond );

      /* Wait for termination */
      pthread_cond_wait( & thread_pool->m_master_cond ,
                         & thread_pool->m_master_lock );
    }

    thread_pool->m_thread_next  = NULL ;
    thread_pool->m_thread_count = 0 ;

    while ( thread_pool->m_lock_init ) {
      pthread_mutex_destroy( thread_pool->m_lock + thread_pool->m_lock_init );
      --( thread_pool->m_lock_init );
    }

    pthread_mutex_unlock( & thread_pool->m_master_lock );
  }

  return result ;
}

/* #ifdef HAVE_PTHREAD */
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
  static int locks[ MAXIMUM_LOCK_COUNT ];
  static struct LocalWorkData data = { locks , 0 , 0 };
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

  return result ;
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
    work->m_active     = 1 ;
    work->m_lock_count = lock_count ;

    for ( i = 0 ; i < lock_count ; ++i ) { work->m_lock[i] = 0 ; }

    for ( i = 0 ; i < work_count ; ++i ) {
      TPI_Work w = { shared_data , lock_count , work_count , i };
      (*subprogram)( & w );
    }

    work->m_active     = 0 ;
    work->m_lock_count = 0 ;
  }

  return result ;
}

int TPI_Run_threads( TPI_work_subprogram subprogram  ,
                     void * const        shared_data ,
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

    work->m_active     = 1 ;
    work->m_lock_count = lock_count ;

    for ( i = 0 ; i < lock_count ; ++i ) { work->m_lock[i] = 0 ; }

    {
      TPI_Work w = { shared_data , lock_count , 1 , 0 };
      (*subprogram)( & w );
    }

    work->m_active     = 0 ;
    work->m_lock_count = 0 ;
  }

  return result ;
}

int TPI_Init( int nthread )
{
  LocalWork * const work = local_work();
  int result = 0 ;

  nthread = 1 ;

  if ( ! work || work->m_active ) { result = TPI_ERROR_ACTIVE ; }

  if ( ! result ) { result = nthread ; }

  return result ;
}

int TPI_Finalize()
{
  LocalWork * const work = local_work();
  int result = 0 ;

  if ( ! work || work->m_active ) { result = TPI_ERROR_ACTIVE ; }

  return result ;
}

#endif

