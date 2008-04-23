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
#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include <sys/types.h>

#include <TPI.h>

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

typedef struct TPI_ThreadPool_Private {
  pthread_mutex_t * m_lock ;
  int               m_lock_size ;
  int               m_group_size ;
  int               m_group_rank ;
  int               m_size ;
  int               m_rank ;
} ThreadWork ;

typedef struct ThreadPool_Data {
  pthread_mutex_t           m_pool_lock ;
  pthread_mutex_t           m_pool_lock_run ;
  pthread_cond_t            m_pool_cond ;
  pthread_cond_t            m_pool_cond_run ;
  TPI_parallel_subprogram * m_work_group_routine ;
  void                   ** m_work_group_argument ;
  const int               * m_work_group_number ;
  pthread_mutex_t         * m_work_lock ;
  int                       m_work_lock_size ;
  int                       m_work_group_size ;
  int                       m_number_threads ;
  int                       m_work_size ;
  int                       m_work_begin_count ;
  int                       m_work_end_count ;
  int                       m_thread_work ;
} ThreadPool ;

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Group_rank( TPI_ThreadPool local , int * rank , int * size )
{
  const int result = NULL == local ? TPI_ERROR_NULL : 0 ;

  if ( ! result ) {
    if ( rank ) { *rank = local->m_group_rank ; }
    if ( size ) { *size = local->m_group_size ; }
  }
  return result ;
}

int TPI_Rank( TPI_ThreadPool local , int * rank , int * size )
{
  const int result = NULL == local ? TPI_ERROR_NULL : 0 ;

  if ( ! result ) {
    if ( rank ) { *rank = local->m_rank ; }
    if ( size ) { *size = local->m_size ; }
  }
  return result ;
}

int TPI_Partition( int Thread_Rank ,
                   int Thread_Size ,
                   int Number ,
                   int * Local_Begin ,
                   int * Local_Number )
{
  const int result =
    Local_Begin && Local_Number
      ? ( 0 <= Thread_Rank && Thread_Rank < Thread_Size ? 0 : TPI_ERROR_SIZE )
      :  TPI_ERROR_NULL ;

  if ( ! result ) {
    const int rem  = Number % Thread_Size ;
    const int base = Number / Thread_Size ;
    const int add  = Thread_Rank < rem ;
    *Local_Begin  = base * Thread_Rank + ( add ? Thread_Rank : rem );
    *Local_Number = base +               ( add ? 1 : 0 );
  }
  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Lock( TPI_ThreadPool local , int i )
{
  const int result = 
    ( NULL == local                           ? TPI_ERROR_NULL :
    ( i < 0 || local->m_lock_size <= i        ? TPI_ERROR_SIZE :
    ( pthread_mutex_lock( local->m_lock + i ) ? TPI_ERROR_LOCK : 0 ) ) );
  return result ;
}

int TPI_Trylock( TPI_ThreadPool local , int i )
{
  int p ;
  const int result = 
    ( NULL == local                    ? TPI_ERROR_NULL :
    ( i < 0 || local->m_lock_size <= i ? TPI_ERROR_SIZE :
    ( EBUSY == ( p = pthread_mutex_trylock(local->m_lock+i) ) ? TPI_LOCK_BUSY :
    ( p ? TPI_ERROR_LOCK : 0 ) ) ) );
  return result ;
}

int TPI_Unlock( TPI_ThreadPool local , int i )
{
  const int result = 
    ( NULL == local                             ? TPI_ERROR_NULL :
    ( i < 0 || local->m_lock_size <= i          ? TPI_ERROR_SIZE :
    ( pthread_mutex_unlock( local->m_lock + i ) ? TPI_ERROR_LOCK : 0 ) ) );
  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*  Run the work queue.
 *  A worker runs until commanded to terminate.
 *  The control runs until the current queue is empty.
 */
static void local_thread_pool_run_work( ThreadPool * const pool ,
                                        const int am_worker )
{
  pthread_mutex_t * const lock = & pool->m_pool_lock_run ;

  int working = 1 ;
  int work_n = 0 ;

  while ( working ) {
    if ( ! pthread_mutex_trylock( lock ) ) {

      if ( work_n ) {
        /* In the previous iteration I finished some work. */
        --( pool->m_work_end_count );
      }

      if ( am_worker && 0 <= pool->m_work_end_count &&
                        0 == pool->m_work_begin_count ) {
        /*  Worker with no work to do, unlock and block */
        work_n = 0 ;
        pthread_cond_wait( & pool->m_pool_cond_run , lock );
        /*  Now running and have the lock; however, 'pool' has changed */
      }

      {
        const int work_end_count = pool->m_work_end_count ;

        if ( 0 < work_end_count && 0 < pool->m_work_begin_count ) {
          /* There is work for me to do. */
          if ( ! work_n ) {
            /* First time I have obtained work from this queue */
            ++( pool->m_thread_work );
          }
          work_n = pool->m_work_begin_count ;
          --( pool->m_work_begin_count );
        }
        else {
          /* There is no work for me to do. */
          work_n = 0 ;
          if ( ( 0 >  work_end_count ) ||
               ( 0 == work_end_count && ! am_worker ) ) {
            /* All work is done */
            working = 0 ;
          }
        }
      }

      pthread_mutex_unlock( lock );

      /*----------------------*/

      if ( work_n ) { /* Have work to do */
        const int         work_group_size   = pool->m_work_group_size ;
        const int * const work_group_number = pool->m_work_group_number ;
        int i = pool->m_work_size - work_n ;
        int j = 0 ;
        while ( j < work_group_size && work_group_number[j] <= i ) {
          i -= work_group_number[j] ;
          ++j ;
        }

        {
          const TPI_parallel_subprogram work_routine =
            pool->m_work_group_routine[j] ;
          void * const work_argument =
            pool->m_work_group_argument[j] ;

          ThreadWork work = {
            /*  m_lock       */  pool->m_work_lock ,
            /*  m_lock_size  */  pool->m_work_lock_size ,
            /*  m_group_size */  pool->m_work_group_size ,
            /*  m_group_rank */  j ,
            /*  m_size       */  pool->m_work_group_number[j] ,
            /*  m_rank       */  i };

          (*work_routine)( work_argument , & work );
        }
      }
    }
  }

  return ;
}

/*--------------------------------------------------------------------*/
/*  The driver given to 'pthread_create'.
 *  Run work until told to terminate.
 */

static void * local_thread_pool_driver( void * arg )
{
  ThreadPool * const pool = (ThreadPool*) arg ;

  /*------------------------------*/
  /* Start up */

  pthread_mutex_lock(   & pool->m_pool_lock );
  ++( pool->m_number_threads );
  pthread_cond_signal(  & pool->m_pool_cond );
  pthread_mutex_unlock( & pool->m_pool_lock );

  /*------------------------------*/

  local_thread_pool_run_work( pool , 1 );

  /*------------------------------*/
  /* Termination:
   *   Remove myself from the pool.
   *   Signal the pool when it is empty.
   */

  pthread_mutex_lock( & pool->m_pool_lock );

  if ( ! --( pool->m_number_threads ) ) {
    pthread_cond_signal( & pool->m_pool_cond );
  }
  pthread_mutex_unlock( & pool->m_pool_lock );

  return NULL ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/* The work queue shared by all threads */

static ThreadPool * local_thread_pool()
{
  static ThreadPool pool = {
    /* m_pool_lock           */  PTHREAD_MUTEX_INITIALIZER ,
    /* m_pool_lock_run       */  PTHREAD_MUTEX_INITIALIZER ,
    /* m_pool_cond           */  PTHREAD_COND_INITIALIZER ,
    /* m_pool_cond_run       */  PTHREAD_COND_INITIALIZER ,
    /* m_work_group_routine  */  NULL ,
    /* m_work_group_argument */  NULL ,
    /* m_work_group_number   */  NULL ,
    /* m_work_lock           */  NULL ,
    /* m_work_lock_size      */  0 ,
    /* m_work_group_size     */  0 ,
    /* m_number_threads      */  0 ,
    /* m_work_size           */  0 ,
    /* m_work_begin_count    */  0 ,
    /* m_work_end_count      */  0 ,
    /* m_thread_work         */  0 };

  /* Guard against recursive call */

  return pool.m_work_group_routine ? NULL : & pool ;
}

/*--------------------------------------------------------------------*/

int TPI_Run_many( const int number_routine ,
                  TPI_parallel_subprogram routine[] ,
                  void * routine_data[] ,
                  const int number[] )
{
  ThreadPool * const pool = local_thread_pool();

  int i , nwork ;

  int result =
    ( NULL == routine ||
      NULL == routine_data ||
      NULL == number ? TPI_ERROR_NULL :
    ( NULL == pool ? TPI_ERROR_ACTIVE : 0 ) );  

  for ( nwork = i = 0 ; ! result && i < number_routine ; ++i ) {
    if ( ! routine[i] ) {
      result = TPI_ERROR_NULL ;
    }
    else if ( number[i] <= 0 ) {
      result = TPI_ERROR_SIZE ;
    }
    else {
      nwork += number[i] ;
    }
  }

  if ( ! result ) {

    const int number_total = nwork ;
    const int number_locks = pool->m_work_lock_size ;

    pthread_mutex_t lock[ number_locks ];

    int nlocks = 0 ;

    while ( nlocks < number_locks && ! result ) {
      if ( pthread_mutex_init( lock + nlocks , NULL ) ) {
        result = TPI_ERROR_INTERNAL ;
      }
      else {
        ++nlocks ;
      }
    }

    if ( ! result ) {

      pool->m_work_group_routine  = routine ;
      pool->m_work_group_argument = routine_data ;
      pool->m_work_lock           = lock ;
      pool->m_work_lock_size      = number_locks ;
      pool->m_work_group_size     = number_routine ;
      pool->m_work_group_number   = number ;
      pool->m_work_size           = number_total ;
      pool->m_work_begin_count    = number_total ;
      pool->m_work_end_count      = number_total ; /* Trigger to start */

      pthread_cond_broadcast( & pool->m_pool_cond_run ); /* Unblock workers */

      local_thread_pool_run_work(pool,0); /* Participate in the work */

      pool->m_work_size           = 0 ;
      pool->m_work_group_routine  = NULL ;
      pool->m_work_group_argument = NULL ;
      pool->m_work_lock           = NULL ;
      pool->m_work_lock_size      = 0 ;
      pool->m_work_group_size     = 0 ;
      pool->m_work_group_number   = 0 ;
    }

    while ( nlocks-- ) { pthread_mutex_destroy( lock + nlocks ); }
  }

  return result ;
}

/* Run one routine on all allocated threads. */

int TPI_Run( TPI_parallel_subprogram routine , void * routine_data )
{
  ThreadPool * const pool = local_thread_pool();

  int result = ! pool ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {
    const int number = pool->m_number_threads ;
    result = TPI_Run_many( 1 , & routine , & routine_data , & number );
  }

  return result ;
}

int TPI_Set_lock_size( int number )
{
  ThreadPool * const pool = local_thread_pool();

  int result = ! pool ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && number < 0 ) { result = TPI_ERROR_SIZE ; }

  if ( ! result ) { pool->m_work_lock_size = number ; }

  return result ;
}

int TPI_Run_count( int * count )
{
  ThreadPool * const pool = local_thread_pool();

  int result = ! pool ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && ! count ) { result = TPI_ERROR_NULL ; }

  if ( ! result ) {
    *count = pool->m_thread_work ;
    pool->m_thread_work = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Init( int n )
{
  ThreadPool * const pool = local_thread_pool();

  int result = ! pool || pool->m_number_threads ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && n <= 0 ) { result = TPI_ERROR_SIZE ; }

  if ( ! result ) {
    pthread_attr_t thread_attr ;

    if ( pthread_attr_init( & thread_attr ) ) {
      result = TPI_ERROR_INTERNAL ;
    }
    else {
      int n_thread = 1 ; /* Count myself among the threads */

      pool->m_number_threads = 1 ;

      pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
      pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

      pthread_mutex_lock( & pool->m_pool_lock );

      for ( ; n_thread < n && ! result ; ++n_thread ) {
        pthread_t pt ;

        if ( pthread_create( & pt, & thread_attr,
                             & local_thread_pool_driver, pool ) ) {
          result = TPI_ERROR_INTERNAL ;
        }
        else {
          /* Wait for start */
          pthread_cond_wait( & pool->m_pool_cond , & pool->m_pool_lock );
        }
      }

      pthread_attr_destroy( & thread_attr );

      pthread_mutex_unlock( & pool->m_pool_lock );
    }

    if ( result ) { TPI_Finalize(); }
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  ThreadPool * const pool = local_thread_pool();

  int result = ! pool ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) {
    --( pool->m_number_threads );

    if ( 0 < pool->m_number_threads ) {

      pthread_mutex_lock( & pool->m_pool_lock );

      pool->m_work_end_count = -1 ; /* Trigger to terminate */

      /* Make sure workers are waiting for signal before signaling */
      pthread_mutex_lock(     & pool->m_pool_lock_run );
      pthread_cond_broadcast( & pool->m_pool_cond_run );
      pthread_mutex_unlock(   & pool->m_pool_lock_run );

      pthread_cond_wait(    & pool->m_pool_cond , & pool->m_pool_lock );
      pthread_mutex_unlock( & pool->m_pool_lock );
    }

    pool->m_work_size  = 0 ;
    pool->m_work_begin_count = 0 ;
    pool->m_work_end_count = 0 ;
    pool->m_number_threads = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Size( int * number_allocated )
{
  int result = number_allocated ? 0 : TPI_ERROR_NULL ;

  if ( ! result ) {
    ThreadPool * const pool = local_thread_pool();

    if ( pool ) {
      *number_allocated = pool->m_number_threads ;
    }
    else {
      result = TPI_ERROR_ACTIVE ;
    }
  }

  return result ;
}

