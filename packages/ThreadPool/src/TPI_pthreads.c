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

#include <unistd.h>
#include <errno.h>
#include <pthread.h>
#include <sys/types.h>

#include <TPI.h>

#define TPI_NO_SCHED 0
#define TPI_USE_SPINLOCK 0

/*  IDEA: try <pthread.h> advanced realtime 'pthread_spinlock_t'. */

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/* Choose between a hard lock or spin lock */

static int local_thread_pool_lock( pthread_mutex_t * T )
{
  int result ;

#if TPI_USE_SPINLOCK

  while ( EBUSY == ( result = pthread_mutex_trylock( T ) ) ); 

#else

  result = pthread_mutex_lock( T );

#endif

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

typedef struct TPI_ThreadPool_Private {
  pthread_mutex_t          m_work_lock ;
  TPI_parallel_subprogram  m_routine ;
  void                   * m_argument ;
  pthread_mutex_t        * m_lock ;
  int                      m_lock_size ;
  int                      m_size ;
  int                      m_rank ;
  int                      m_group_size ;
  int                      m_group_rank ;
} ThreadWork ;

typedef struct ThreadControl_Data {
  pthread_mutex_t             m_thread_lock ;
  pthread_cond_t              m_thread_cond ;
  struct ThreadControl_Data * m_thread_next ;
  unsigned                    m_thread_rank ;
  int                         m_work_count ;
} ThreadControl ;
  
typedef struct ThreadPool_Data {
  pthread_mutex_t m_pool_lock ;
  pthread_cond_t  m_pool_wait ;
  ThreadControl * m_thread_first ;
  ThreadWork    * m_work_begin ;
  int             m_work_size ;
  int             m_number_threads ;
  int             m_number_locks ;
  int             m_work_count ;    /* Root thread */
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
      ? ( 0 <= Thread_Rank && Thread_Rank < Thread_Size ? 0 : TPI_ERROR_SIZE )       :  TPI_ERROR_NULL ;

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
  int result = NULL == local ? TPI_ERROR_NULL : 0 ;

  if ( ! result ) {
    if ( i < 0 || local->m_lock_size <= i ) {
      result = TPI_ERROR_SIZE ;
    }
    else if ( local_thread_pool_lock( local->m_lock + i ) ) {
      result = TPI_ERROR_LOCK ;
    }
  }
  return result ;
}

int TPI_Trylock( TPI_ThreadPool local , int i )
{
  int result = NULL == local ? TPI_ERROR_NULL : 0 ;

  if ( ! result ) {
    if ( i < 0 || local->m_lock_size <= i ) {
      result = TPI_ERROR_SIZE ;
    }
    else {
      switch( pthread_mutex_trylock( local->m_lock + i ) ) {
      case 0 : break ;
      case EBUSY : result = TPI_LOCK_BUSY ; break ;
      default    : result = TPI_ERROR_LOCK ; break ;
      }
    }
  }
  return result ;
}

int TPI_Unlock( TPI_ThreadPool local , int i )
{
  int result = NULL == local ? TPI_ERROR_NULL : 0 ;

  if ( ! result ) {
    if ( i < 0 || local->m_lock_size <= i ) {
      result = TPI_ERROR_SIZE ;
    }
    else if ( pthread_mutex_unlock( local->m_lock + i ) ) {
      result = TPI_ERROR_LOCK ;
    }
  }
  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/*  Run the work queue as long as the queue has work.
 *  Return how many work tasks were actually run.
 */
static int local_thread_pool_run_work( ThreadWork * const work ,
                                       const unsigned     begin ,
                                       const unsigned     number )
{
  int counter = 0 ;
  unsigned i = 0 ;

  for ( ; i < number ; ++i ) {
     ThreadWork * const w = work + ( i + begin ) % number ;

     if ( w->m_routine && ! pthread_mutex_trylock( & w->m_work_lock ) ) {

      /*  Cannot be sure of 'work->m_routine' until I have it locked.
       *  Claim the work by setting the 'work->m_routine' to NULL.
       *  Do this quickly so my sibling threads know it has been claimed.
       */

      const TPI_parallel_subprogram routine = w->m_routine ;

      w->m_routine = NULL ;

      pthread_mutex_unlock( & w->m_work_lock );

      if ( routine ) {
        (*routine)( w->m_argument , w );
        ++counter ;
      }
    }
  }
  return counter ;
}

/*--------------------------------------------------------------------*/
/*  The driver given to 'pthread_create'.
 *  Increment number of threads on start up.
 *  Run work until told to terminate.
 *  Decrement number of threads on shut down.
 */

static void * local_thread_pool_driver( void * arg )
{
  ThreadControl control = {
    /* m_thread_lock */  PTHREAD_MUTEX_INITIALIZER ,
    /* m_thread_cond */  PTHREAD_COND_INITIALIZER ,
    /* m_thread_next */  NULL ,
    /* m_thread_rank */  0 ,
    /* m_work_count  */  0 };

  ThreadPool * pool = (ThreadPool*) arg ;

  pthread_mutex_lock( & control.m_thread_lock );

  control.m_thread_next = pool->m_thread_first ;

  pool->m_thread_first = & control ;

  /* Signal I've started */
  pthread_mutex_lock(   & pool->m_pool_lock );
  pthread_cond_signal(  & pool->m_pool_wait );
  pthread_mutex_unlock( & pool->m_pool_lock );

  while ( pool && ! pthread_cond_wait( & control.m_thread_cond ,
                                       & control.m_thread_lock ) ) {
    if ( pool->m_work_begin ) {
      control.m_work_count +=
        local_thread_pool_run_work( pool->m_work_begin ,
                                    control.m_thread_rank ,
                                    pool->m_work_size );
    }
    else {
      /* Signal I'm terminating */
      pool->m_thread_first = control.m_thread_next ;
      pthread_mutex_lock(   & pool->m_pool_lock );
      pthread_cond_signal(  & pool->m_pool_wait );
      pthread_mutex_unlock( & pool->m_pool_lock );
      pool = NULL ;
    }
  }

  pthread_mutex_unlock(  & control.m_thread_lock );
  pthread_mutex_destroy( & control.m_thread_lock );
  pthread_cond_destroy(  & control.m_thread_cond );

  return NULL ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/* The work queue shared by all threads */

static ThreadPool * local_thread_pool()
{
  static ThreadPool pool = {
    /* m_pool_lock      */  PTHREAD_MUTEX_INITIALIZER ,
    /* m_pool_wait      */  PTHREAD_COND_INITIALIZER ,
    /* m_thread_first   */  NULL ,
    /* m_work_begin     */  NULL ,
    /* m_work_size      */  0 ,
    /* m_number_threads */  0 ,
    /* m_number_locks   */  0 ,
    /* m_work_count     */  0 };

  /* Guard against recursive call */

  return pool.m_work_begin ? NULL : & pool ;
}

/*--------------------------------------------------------------------*/

int TPI_Run_many( const int number_routine ,
                  TPI_parallel_subprogram routine[] ,
                  void * routine_data[] ,
                  const int number[] )
{
  int i , nwork ;

  int result =
    ! number || ! routine || ! routine_data || ! number ? TPI_ERROR_NULL : 0 ;

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

    ThreadPool * const pool = local_thread_pool();

    if ( ! pool ) { result = TPI_ERROR_ACTIVE ; }

    if ( ! result ) {

      const int number_total = nwork ;
      const int number_locks = pool->m_number_locks ;

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

        ThreadWork work[ number_total ];

        { /* Fill the work queue */
          ThreadWork * w = work ;

          for ( i = 0 ; i < number_routine ; ++i ) {
            int k = 0 ;
            for ( ; k < number[i] ; ++k , ++w ) {
              pthread_mutex_init( & w->m_work_lock , NULL );
              w->m_routine   = routine[i] ;
              w->m_argument  = routine_data[i] ;
              w->m_lock      = lock ;
              w->m_lock_size = number_locks ;
              w->m_size      = number[i] ;
              w->m_rank      = k ;
              w->m_group_rank = i ;
              w->m_group_size = number_routine ;
            }
          }
        }

        pool->m_work_begin = work ;
        pool->m_work_size  = number_total ;

        { /* Signal all thread to get to work */
          ThreadControl * c = pool->m_thread_first ;
          for ( ; c ; c = c->m_thread_next ) {
            pthread_cond_signal(  & c->m_thread_cond );
            pthread_mutex_unlock( & c->m_thread_lock );
          }
        }

        /* Participate in the work */
        pool->m_work_count += local_thread_pool_run_work(work,0,number_total);

        { /* Reacquire all thread locks - to be sure they are blocked */
          ThreadControl * c = pool->m_thread_first ;
          for ( ; c ; c = c->m_thread_next ) {
            pthread_mutex_lock( & c->m_thread_lock );
          }
        }

        pool->m_work_begin = NULL ;
        pool->m_work_size  = 0 ;
        pool->m_number_locks = 0 ;

        { /* Destroy the work queue locks */
          ThreadWork * const work_end = work + number_total ;
          ThreadWork * w = work ;
          for ( ; w < work_end ; ++w ) {
            pthread_mutex_destroy( & w->m_work_lock );
          }
        }
      }

      while ( nlocks-- ) { pthread_mutex_destroy( lock + nlocks ); }
    }
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

  if ( ! result ) { pool->m_number_locks = number ; }

  return result ;
}

int TPI_Run_count( int number , int * count )
{
  ThreadPool * const pool = local_thread_pool();

  int result = ! pool ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && ! count ) { result = TPI_ERROR_NULL ; }

  if ( ! result ) {
    ThreadControl * c ;

    if ( number ) {
      *count = pool->m_work_count ;
      for ( c = pool->m_thread_first ; c && --number ; c = c->m_thread_next ) {
        *++count = c->m_work_count ;
      }
      while ( --number ) { *++count = 0 ; }
    }

    pool->m_work_count = 0 ;
    for ( c = pool->m_thread_first ; c ; c = c->m_thread_next ) {
      c->m_work_count = 0 ;
    }
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
          ThreadControl * volatile * const control = & pool->m_thread_first ;
          pthread_cond_wait( & pool->m_pool_wait , & pool->m_pool_lock );
          pthread_mutex_lock( & (*control)->m_thread_lock );
          (*control)->m_thread_rank = n_thread ;
        }
      }

      pthread_mutex_unlock( & pool->m_pool_lock );

      pthread_attr_destroy( & thread_attr );

      pool->m_number_threads = n_thread ;
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
    ThreadControl * volatile * const control = & pool->m_thread_first ;

    pthread_mutex_lock( & pool->m_pool_lock );

    while ( *control ) {
      pthread_cond_signal(  & (*control)->m_thread_cond );
      pthread_mutex_unlock( & (*control)->m_thread_lock );
      pthread_cond_wait( & pool->m_pool_wait , & pool->m_pool_lock );
    }

    pthread_mutex_unlock( & pool->m_pool_lock );

    pool->m_number_threads = 0 ;
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Size( int * number_allocated , int * number_concurrent )
{
  int result = 0 ;
  int count = 0 ;

  if ( number_concurrent ) {

#if ! TPI_NO_SCHED

    enum { NTMP = 8 };

    extern int sched_getaffinity( pid_t, unsigned int , unsigned long * );

    unsigned long tmp[ NTMP ] = { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 };

    if ( ! sched_getaffinity( 0 , sizeof(tmp) , tmp ) ) {

      int i ;

      for ( i = 0 ; i < NTMP ; ++i ) {
        for ( ; tmp[i] ; tmp[i] >>= 1 ) {
          if ( tmp[i] & 01 ) { ++count ; }
        }
      }
    }
    else {
      result = TPI_ERROR_INTERNAL ;
    }
#endif

    *number_concurrent = count ;
  }

  if ( number_allocated ) {
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

