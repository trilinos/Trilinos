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

#include <TPI.h>

#define NO_SCHED
#define USE_SPIN_LOCK 0
#define USE_MUTEX_TRYLOCK 1

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#if ! defined(NO_PTHREADS)

#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <pthread.h>

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/* Interact with thread scheduler to query and control thread affinity */

#if ! defined(NO_SCHED)

#include <sys/types.h>

static int set_affinity( int i )
{
  extern int sched_getaffinity( pid_t, unsigned int , unsigned long * );
  extern int sched_setaffinity( pid_t, unsigned int , unsigned long * );

  unsigned long tmp = 0 ;

  int result = sched_getaffinity( 0 , sizeof(unsigned long) , & tmp );

  if ( ! result ) {
    result = ( tmp &= 1 << i ) ? 0 : -1 ;
  }
/* For now simply verify that the rank is within the CPU count */
/*
  if ( ! result ) {
    result = sched_setaffinity( 0 , sizeof(unsigned long) , & tmp );
  }
  if ( ! result ) {
    result = sched_getaffinity( 0 , sizeof(unsigned long) , & tmp );
  }
  if ( ! result ) {
    result = ( tmp == 1 << i ) ? 0 : -1 ;
  }
*/
  return result ;
}

#else

static int set_affinity( int i ) { return i ? 0 : 0  ; }

#endif

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#if USE_SPIN_LOCK

/*  try <pthread.h> advanced realtime threads option. */
typedef pthread_spinlock_t fastlock_type ;

#define FASTLOCK_LOCK( T )     pthread_spin_lock( T )
#define FASTLOCK_TRYLOCK( T )  pthread_spin_trylock( T )
#define FASTLOCK_UNLOCK( T )   pthread_spin_unlock( T )
#define FASTLOCK_INIT( T )     pthread_spin_init( T , PTHREAD_PROCESS_PRIVATE )
#define FASTLOCK_DESTROY( T )  pthread_spin_destroy( T )

#elif USE_MUTEX_TRYLOCK

typedef pthread_mutex_t   fastlock_type ;

static int FASTLOCK_LOCK( fastlock_type * T )
{
  int result ;
  while ( EBUSY == ( result = pthread_mutex_trylock( T ) ) ); 
  return result ;
}

#define FASTLOCK_TRYLOCK( T )  pthread_mutex_trylock( T )
#define FASTLOCK_UNLOCK( T )   pthread_mutex_unlock( T )
#define FASTLOCK_INIT( T )     pthread_mutex_init( T , NULL )
#define FASTLOCK_DESTROY( T )  pthread_mutex_destroy( T )

#else

typedef pthread_mutext_t   fastlock_type ;

#define FASTLOCK_LOCK( T )     pthread_mutex_lock( T )
#define FASTLOCK_TRYLOCK( T )  pthread_mutex_trylock( T )
#define FASTLOCK_UNLOCK( T )   pthread_mutex_unlock( T )
#define FASTLOCK_INIT( T )     pthread_mutex_init( T , NULL )
#define FASTLOCK_DESTROY( T )  pthread_mutex_destroy( T )

#endif

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

struct ThreadPool_Data ;

typedef struct TPI_ThreadPool_Private {
  fastlock_type            m_lock ;
  int                   (* m_check )( TPI_ThreadPool );
  volatile int             m_status ;
  struct ThreadPool_Data * m_pool ;
  void                   * m_buffer ;
} ThreadData ;

typedef struct ThreadPool_Data {
  pthread_mutex_t          m_active ; /* Contention indicates an error */
  struct ThreadPool_Data * m_parent ;
  ThreadData             * m_thread ;
  int                      m_thread_size ;

  fastlock_type          * m_mutex ;
  int                      m_mutex_size ;
  TPI_parallel_subprogram  m_routine ;
  void                   * m_argument ;
  int                      m_buf_size ;
} ThreadPool ;

/*--------------------------------------------------------------------*/

static int local_thread_pool_check( TPI_ThreadPool local )
{
  return NULL == local ? TPI_ERROR_NULL :
       ( & local_thread_pool_check != local->m_check ? TPI_ERROR_POOL : 0 );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Rank( TPI_ThreadPool local , int * rank , int * size )
{
  const int result = local_thread_pool_check( local );
  if ( ! result ) {
    if ( NULL != rank ) { *rank = local - local->m_pool->m_thread ; }
    if ( NULL != size ) { *size = local->m_pool->m_thread_size ; }
  }
  else {
    if ( NULL != rank ) { *rank = 0 ; }
    if ( NULL != size ) { *size = 0 ; }
  }
  return result ;
}

int TPI_Partition( TPI_ThreadPool local , int N ,
                   int * I_local , int * N_local )
{
  const int result = local_thread_pool_check( local );
  if ( ! result ) {
    const int rank = local - local->m_pool->m_thread ;
    const int size = local->m_pool->m_thread_size ;
    const int rem  = N % size ;
    const int base = N / size ;
    const int add  = rank < rem ;
    if ( NULL != I_local ) { *I_local = base * rank + ( add ? rank : rem ); }
    if ( NULL != N_local ) { *N_local = base +        ( add ? 1 : 0 ); }
  }
  else {
    if ( NULL != I_local ) { *I_local = 0 ; }
    if ( NULL != N_local ) { *N_local = 0 ; }
  }
  return result ;
}

int TPI_Buffer( TPI_ThreadPool local , void ** buffer , int * size )
{
  const int result = local_thread_pool_check( local );
  if ( ! result ) {
    if ( NULL != buffer ) { *buffer = local->m_buffer ; }
    if ( NULL != size   ) { *size   = local->m_pool->m_buf_size ; }
  }
  else {
    if ( NULL != buffer ) { *buffer = NULL ; }
    if ( NULL != size   ) { *size   = 0 ; }
  }
  return result ;
}

int TPI_Lock_size( TPI_ThreadPool local , int * size )
{
  const int result = local_thread_pool_check( local );
  if ( ! result && NULL != size ) { *size = local->m_pool->m_mutex_size ; }
  return result ;
}

int TPI_Split_lock_size( TPI_ThreadPool local , int * size )
{
  const int result = local_thread_pool_check( local );
  if ( ! result && NULL != size ) {
    *size = local->m_pool->m_parent
          ? local->m_pool->m_parent->m_mutex_size : 0 ;
  }
  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static int local_thread_pool_lock( ThreadPool * pool , int i )
{
  int result = ( NULL != pool && 0 <= i && i < pool->m_mutex_size )
               ? 0 : TPI_ERROR_SIZE ;
  if ( ! result ) {
    result = FASTLOCK_LOCK( pool->m_mutex + i ) ? TPI_ERROR_LOCK : 0 ;
  }
  return result ;
}

static int local_thread_pool_trylock( ThreadPool * pool , int i )
{
  int result = ( NULL != pool && 0 <= i && i < pool->m_mutex_size )
               ? 0 : TPI_ERROR_SIZE ;
  if ( ! result ) {
    switch( FASTLOCK_TRYLOCK( pool->m_mutex + i ) ) {
    case 0 : break ;
    case EBUSY : result = TPI_LOCK_BUSY ; break ;
    default    : result = TPI_ERROR_LOCK ; break ;
    }
  }
  return result ;
}

static int local_thread_pool_unlock( ThreadPool * pool , int i )
{
  int result = ( NULL != pool && 0 <= i && i < pool->m_mutex_size )
               ? 0 : TPI_ERROR_SIZE ;
  if ( ! result ) {
    result = FASTLOCK_UNLOCK( pool->m_mutex + i ) ? TPI_ERROR_LOCK : 0 ;
  }
  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Lock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) { result = local_thread_pool_lock( local->m_pool , i ); }
  return result ;
}

int TPI_Split_lock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) {
    result = local_thread_pool_lock( local->m_pool->m_parent , i );
  }
  return result ;
}

int TPI_Trylock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) {
    result = local_thread_pool_trylock( local->m_pool , i );
  }
  return result ;
}

int TPI_Split_trylock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) {
    result = local_thread_pool_trylock( local->m_pool->m_parent , i );
  }
  return result ;
}

int TPI_Unlock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) { result = local_thread_pool_unlock( local->m_pool , i ); }
  return result ;
}

int TPI_Split_unlock( TPI_ThreadPool local , int i )
{
  int result = local_thread_pool_check( local );
  if ( ! result ) {
    result = local_thread_pool_unlock( local->m_pool->m_parent , i );
  }
  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static void local_thread_pool_die( void * arg , TPI_ThreadPool pool )
{
  if ( NULL == pool || arg != pool->m_pool ) {
    fprintf(stderr,"TPI internal corruption\n");
    fflush(stderr);
  }
}

static int local_thread_pool_run_routine( ThreadData * const td )
{
  int working = 1 ;

  if ( ! FASTLOCK_TRYLOCK( & td->m_lock ) ) {

    const ThreadPool * const p = td->m_pool ;
    const unsigned size = p->m_buf_size ;

    if ( td->m_status ) {
      char buffer[ size ]; /* Requested buffer for this run */
      void * const old = td->m_buffer ; /* Save the old buffer if nesting */
      td->m_buffer = size ? buffer : NULL ;
      (* p->m_routine)( p->m_argument, td );
      td->m_buffer = old ;

      td->m_status = 0 ; /* Have completed */
    }

    working = local_thread_pool_die != td->m_pool->m_routine ;

    FASTLOCK_UNLOCK( & td->m_lock );
  }

  return working ;
}

/*--------------------------------------------------------------------*/
/* Driver: call the routine, check for thread termination */

static void * local_thread_pool_driver( void * arg )
{
  ThreadData * const my_td = (ThreadData*)( arg );
  ThreadPool * const root_pool = my_td->m_pool ; /* started on root */
  const unsigned root_size = root_pool->m_thread_size ;
  const unsigned root_rank = my_td - root_pool->m_thread ;

  int working = 1 ;

  while ( working ) {

    /* Always do my work first, I may terminate */

    if ( my_td->m_status ) {

      working = local_thread_pool_run_routine( my_td );

      if ( working ) {
        /* Finished my work, try to steal some unstarted work */
        unsigned i = 1 ;
        for ( ; i < root_size ; ++i ) {
          const unsigned p = ( i + root_rank ) % root_size ;
          if ( p ) {
            ThreadData * const td = root_pool->m_thread + p ;
            if ( td->m_status ) {
              local_thread_pool_run_routine( td );
            }
          }
        }
      }
    }
  }

  return NULL ;
}

/*--------------------------------------------------------------------*/

static int local_thread_pool_run( ThreadPool * const pool )
{
  const unsigned number_thread = pool->m_thread_size ;
  const unsigned number_locks  = pool->m_mutex_size ;

  int result = 0 ;

  {
    fastlock_type locks[ number_locks ];

    unsigned nlocks = 0 ;

    while ( nlocks < number_locks && ! result ) {
      if ( FASTLOCK_INIT( locks + nlocks ) ) {
        result = TPI_ERROR_INTERNAL ;
      }
      else {
        ++nlocks ;
      }
    }

    if ( ! result ) {

      pool->m_mutex = locks ;

      {
        unsigned i = 0 ;
        for ( i = 0 ; i < number_thread ; ++i ) {
          ThreadData * const td = pool->m_thread + i ;
          td->m_status = NULL != td->m_pool ; /* Start running now */
        }
      }

      local_thread_pool_run_routine( pool->m_thread );

      {
        int waiting = 1 ;
        while ( waiting ) {
          waiting = 0 ;
          for ( unsigned i = 1 ; i < number_thread ; ++i ) {
            ThreadData * const td = pool->m_thread + i ;
            if ( td->m_status ) { /* Tried to start this, maybe it ran ? */
              waiting = 1 ;
              local_thread_pool_run_routine( td ); /* Make sure it ran */
            }
          }
        }
      }
    }

    while ( nlocks-- ) { FASTLOCK_DESTROY( locks + nlocks ); }
  }

  pool->m_mutex      = NULL ;
  pool->m_mutex_size = 0 ;
  pool->m_routine    = NULL ;
  pool->m_argument   = NULL ;
  pool->m_buf_size   = 0 ;

  return result ;
}

/*--------------------------------------------------------------------*/
/*  Acquire the active lock for the pool.
 *  This thread must be the first thread in the pool.
 */

static int local_thread_pool_acquire_active_lock( TPI_ThreadPool local )
{
  int result = local_thread_pool_check( local );

  if ( ! result ) {
    ThreadPool * const pool = local->m_pool ;
    if ( local != pool->m_thread || pthread_mutex_trylock(& pool->m_active) ) {
      result = TPI_ERROR_ACTIVE ;
    }
  }
  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Set_lock_size( TPI_ThreadPool local , int size )
{
  int result = local_thread_pool_acquire_active_lock( local );
  if ( ! result ) {
    local->m_pool->m_mutex_size = size ;
    pthread_mutex_unlock( & local->m_pool->m_active );
  }
  return result ;
}

int TPI_Set_buffer_size( TPI_ThreadPool local , int size )
{
  int result = local_thread_pool_acquire_active_lock( local );
  if ( ! result ) {
    local->m_pool->m_buf_size = size ;
    pthread_mutex_unlock( & local->m_pool->m_active );
  }
  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Run( TPI_ThreadPool local ,
             TPI_parallel_subprogram routine ,
             void * routine_data )
{
  int result = local_thread_pool_acquire_active_lock( local );

  if ( ! result ) { /* Now have the active lock for the whole pool */

    ThreadPool * const pool = local->m_pool ;

    pool->m_routine  = routine ;
    pool->m_argument = routine_data ;

    result = local_thread_pool_run( pool );

    pthread_mutex_unlock( & pool->m_active );
  }

  return result ;
}

/*--------------------------------------------------------------------*/

int TPI_Split( TPI_ThreadPool local ,
               const int number ,
               const int sizes[] ,
               TPI_parallel_subprogram routine[] ,
               void * routine_data[] )
{
  int result = local_thread_pool_acquire_active_lock( local );

  if ( ! result ) {

    ThreadPool * const parent = local->m_pool ;

    {
      int child_size = 0 ;
      int i = 0 ;
      for ( ; i < number && ! result ; ++i ) {
        if      ( NULL == routine[i] ) { result = TPI_ERROR_NULL ; }
        else if ( 0    >= sizes[i] )   { result = TPI_ERROR_SIZE ; }
        child_size += sizes[i] ;
      }

      if ( ! result && child_size != parent->m_thread_size ) {
        result = TPI_ERROR_SIZE ;
      }
    }

    if ( ! result ) {

      ThreadPool children[ number ];

      int child = 0 ;
      int i = 0 ;
      
      for ( i = 0 ; i < parent->m_thread_size ; ++i ) {
        parent->m_thread[i].m_pool = NULL ;
      }

      for ( i = 0 ; child < number && ! result ; ) {

        ThreadData * const td = parent->m_thread + i ;

        i += sizes[child] ;

        ThreadPool * const pool = children + child ;

        td->m_pool = pool ;

        pool->m_parent      = parent ;
        pool->m_thread      = td ;
        pool->m_thread_size = sizes[child] ;
        pool->m_mutex       = NULL ;
        pool->m_mutex_size  = 0 ;
        pool->m_routine     = routine[child] ;
        pool->m_argument    = routine_data[child] ;
        pool->m_buf_size    = 0 ;

        if ( pthread_mutex_init( & pool->m_active , NULL ) ) {
          result = TPI_ERROR_INTERNAL ;
        }
        else {
          ++child ;
        }
      }

      if ( ! result ) { result = local_thread_pool_run( parent ); }

      while ( child-- ) {
        pthread_mutex_destroy( & children[child].m_active );
      }

      for ( i = 0 ; i < parent->m_thread_size ; ++i ) {
        parent->m_thread[i].m_pool = parent ;
      }
    }

    pthread_mutex_unlock( & parent->m_active );
  }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

static void local_thread_pool_start_task( void * arg , TPI_ThreadPool local )
{
  const int p_rank = * ((int*) arg );

  if ( NULL  != arg &&
       NULL  != local &&
       NULL  != local->m_pool &&
       local == local->m_pool->m_thread + p_rank &&
       ! set_affinity( p_rank ) ) {
    local->m_check = & local_thread_pool_check ;
  }
}

static void local_thread_pool_terminate_threads( ThreadPool * pool )
{
  const int n = pool->m_thread_size ;

  pool->m_thread_size = 0 ;
  pool->m_mutex_size = 0 ;
  pool->m_routine  = & local_thread_pool_die ;
  pool->m_argument = pool ;
  pool->m_buf_size = 0 ;

  pool->m_thread[0].m_status = 0 ;

  {
    int i = 1 ;
    for ( ; i < n ; ++i ) {
      ThreadData * const td = pool->m_thread + i ;
      td->m_pool = pool ;
      td->m_status = 1 ; /* Start running */
    }
  }

  {
    int waiting = 1 ;
    while ( waiting ) {
      int i = 0 ;
      waiting = 0 ;
      for ( ; i < n ; ++i ) {
        ThreadData * const td = pool->m_thread + i ;
        if ( td->m_status ) {
          waiting = 1 ; /* Wait for thread to die */
        }
        else {
          FASTLOCK_LOCK(    & td->m_lock );
          FASTLOCK_UNLOCK(  & td->m_lock );
          FASTLOCK_DESTROY( & td->m_lock );
        }
      }
    }
  }
}

static int local_thread_pool_create_threads( int n , ThreadPool * pool )
{
  int result = pool->m_thread_size ? TPI_ERROR_ACTIVE : 0 ;

  int n_lock = 0 ;
  int n_thread = 0 ;

  /*------------------------------------------------------------------*/
  /* Initialize thread data and create all thread run-control locks */
  if ( ! result ) {

    for ( n_lock = 0 ; n_lock < n && ! result ; ) {
      ThreadData * const td = pool->m_thread + n_lock ;

      td->m_check  = NULL ;
      td->m_pool   = pool ;
      td->m_buffer = NULL ;
      td->m_status = 1 ;

      if ( FASTLOCK_INIT( & td->m_lock ) ) {
        result = TPI_ERROR_INTERNAL ;
      }
      else {
        ++n_lock ;
      }
    }
  }
  /*------------------------------------------------------------------*/
  /* Start up all threads and verify they started correctly */
  if ( ! result ) {
    pthread_attr_t thread_attr ;

    if ( pthread_attr_init( & thread_attr ) ) {
      result = TPI_ERROR_INTERNAL ;
    }
    else {

      /* pthread_setconcurrency( n ); */

      pthread_attr_setscope(       & thread_attr, PTHREAD_SCOPE_SYSTEM );
      pthread_attr_setdetachstate( & thread_attr, PTHREAD_CREATE_DETACHED );

      pool->m_routine  = & local_thread_pool_start_task ;
      pool->m_argument = & n_thread ;
      pool->m_thread->m_check = & local_thread_pool_check ;

      for ( n_thread = 1 ; n_thread < n && ! result ; ) {
        pthread_t pt ;

        ThreadData * const td = pool->m_thread + n_thread ;

        if ( pthread_create( & pt ,
                             & thread_attr ,
                             & local_thread_pool_driver ,
                             (void *) td ) ) {
          result = TPI_ERROR_INTERNAL ;
        }
        else {
          int waiting = 1 ;
          while ( waiting ) {
            if ( ! FASTLOCK_LOCK( & td->m_lock ) ) {
              if ( ! td->m_status ) { waiting = 0 ; /* Has run */ }
              FASTLOCK_UNLOCK( & td->m_lock );
            }
          }
          if ( & local_thread_pool_check != td->m_check ) {
            result = TPI_ERROR_INTERNAL ;
          }
          else {
            ++n_thread ;
          }
        }
      }
      pthread_attr_destroy( & thread_attr );
    }

    pool->m_thread_size = n_thread ;
    pool->m_routine  = NULL ;
    pool->m_argument = NULL ;
  }

  /*------------------------------------------------------------------*/
  if ( result ) {
    local_thread_pool_terminate_threads( pool );

    {
      int i = n_thread ;
      for ( ; i < n_lock ; ++i ) {
        FASTLOCK_DESTROY( & pool->m_thread[i].m_lock );
      }
    }
  }

  return result ;
}

static int local_thread_pool_root( int n , TPI_ThreadPool * local )
{
  enum { MAX_THREADS = 128 };

  static ThreadData threads[ MAX_THREADS ];

  static ThreadPool root = {
      /* m_active      */  PTHREAD_MUTEX_INITIALIZER ,
      /* m_parent      */  NULL ,
      /* m_thread      */  threads ,
      /* m_thread_size */  0 ,
      /* m_mutex       */  NULL ,
      /* m_mutex_size  */  0 ,
      /* m_routine     */  NULL ,
      /* m_argument    */  NULL ,
      /* m_buf_size    */  0 };

  int result = 0 < n && n < MAX_THREADS ? 0 : TPI_ERROR_SIZE ;

  if ( ! result && pthread_mutex_trylock( & root.m_active ) ) {
    result = TPI_ERROR_ACTIVE ;
  }

  if ( ! result ) {
    if ( local ) {
      result = local_thread_pool_create_threads( n , & root );
      *local = result ? NULL : threads ;
    }
    else {
      local_thread_pool_terminate_threads( & root );
    }
    pthread_mutex_unlock( & root.m_active );
  }

  return result ;
}

int TPI_Init( int n , TPI_ThreadPool * local )
{
  int result = NULL != local ? 0 : TPI_ERROR_NULL ;
  if ( ! result && set_affinity( 0 ) ) { result = TPI_ERROR_INTERNAL ; }
  if ( ! result ) { result = local_thread_pool_root( n , local ); }
  if ( result ) {
    fprintf(stderr,"TPI_Init FAILED: threads not activated\n");
  }
  return result ;
}

int TPI_Finalize()
{
  return local_thread_pool_root( 1 , NULL );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#else

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

#include <stddef.h>

struct TPI_ThreadPool_Private {
  TPI_ThreadPool m_parent ;
  int    m_lock_size ;
  int    m_buf_size ;
  void * m_buffer ;
};

int TPI_Init( int size , TPI_ThreadPool * local )
{
  static struct TPI_ThreadPool_Private root = { NULL , 0 , 0 , NULL };
  const int result = 1 == size ? 0 : TPI_ERROR_SIZE ;
  if ( ! result && local ) { *local = & root ; }
  return result ;
}

int TPI_Finalize() { return 0 ; }

int TPI_Rank( TPI_ThreadPool local , int * rank , int * size )
{
  const int result = NULL != local ? 0 : TPI_ERROR_NULL ;
  if ( rank ) { *rank = 0 ; }
  if ( size ) { *size = 1 ; }
  return result ;
}

int TPI_Partition( TPI_ThreadPool local , int N ,
                   int * I_local , int * N_local )
{
  const int result = NULL != local ? 0 : TPI_ERROR_NULL ;
  if ( I_local ) { *I_local = 0 ; }
  if ( N_local ) { *N_local = N ; }
  return result ;
}

int TPI_Buffer( TPI_ThreadPool local , void ** buffer , int * size )
{
  const int result = NULL != local ? 0 : TPI_ERROR_NULL ;
  if ( ! result ) {
    if ( buffer ) { *buffer = local->m_buffer ; }
    if ( size )   { *size   = local->m_buf_size ; }
  }
  return result ;
}

int TPI_Lock_size( TPI_ThreadPool local , int * size )
{
  const int result = NULL != local ? 0 : TPI_ERROR_NULL ;
  if ( ! result ) {
    if ( size ) { *size = local->m_lock_size ; }
  }
  return result ;
}

int TPI_Set_lock_size( TPI_ThreadPool local , int size )
{
  const int result = local->m_buffer ? TPI_ERROR_ACTIVE : 0 ;
  if ( ! result ) { local->m_lock_size = size ; }
  return result ;
}

int TPI_Set_buffer_size( TPI_ThreadPool local , int size )
{
  const int result = local->m_buffer ? TPI_ERROR_ACTIVE : 0 ;
  if ( ! result ) { local->m_buf_size = size ; }
  return result ;
}

int TPI_Run( TPI_ThreadPool local ,
             TPI_parallel_subprogram prog ,
             void * shared )
{
  const int result = local->m_buffer ? TPI_ERROR_ACTIVE : 0 ;
  const int buf_size = local->m_buf_size ;
  if ( ! result ) {
    char buffer[ buf_size ];
    local->m_buffer = & buffer[0] ;
    (*prog)( shared , local );
    local->m_buffer = NULL ;
  }
  return result ;
}

int TPI_Lock( TPI_ThreadPool local , int n )
{
  return 0 <= n && n < local->m_lock_size ? 0 : TPI_ERROR_SIZE ;
}

int TPI_Trylock( TPI_ThreadPool local , int n )
{
  return 0 <= n && n < local->m_lock_size ? 0 : TPI_ERROR_SIZE ;
}

int TPI_Unlock( TPI_ThreadPool local , int n )
{
  return 0 <= n && n < local->m_lock_size ? 0 : TPI_ERROR_SIZE ;
}

int TPI_Split( TPI_ThreadPool local ,
               const int nchild           /* Number of child pools     */ ,
               const int sizes[]          /* array of child pool sizes */ ,
               TPI_parallel_subprogram prog[] /* array of main functions  */ ,
               void * shared []           /* array of main functions' data */)
{
  int result = local->m_buffer ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && ( 1 != nchild || 1 != sizes[0] ) ) {
    result = TPI_ERROR_SIZE ;
  }

  if ( ! result ) {
    struct TPI_ThreadPool_Private child = { local , 0 , 0 , NULL };

    (*prog[0])( shared[0] , & child );
  }

  return result ;
}

int TPI_Split_lock_size( TPI_ThreadPool local , int * size )
{
  if ( size ) { *size = local->m_parent->m_lock_size ; }
  return 0 ;
}

int TPI_Split_lock(    TPI_ThreadPool local , int n )
{
  return 0 <= n && n < local->m_parent->m_lock_size ? 0 : TPI_ERROR_SIZE ;
}

int TPI_Split_trylock( TPI_ThreadPool local , int n )
{
  return 0 <= n && n < local->m_parent->m_lock_size ? 0 : TPI_ERROR_SIZE ;
}

int TPI_Split_unlock(  TPI_ThreadPool local , int n )
{
  return 0 <= n && n < local->m_parent->m_lock_size ? 0 : TPI_ERROR_SIZE ;
}

#endif

