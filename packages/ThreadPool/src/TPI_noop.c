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

#include <stddef.h>
#include <TPI.h>

/*--------------------------------------------------------------------*/

typedef struct TPI_ThreadPool_Private {
  TPI_parallel_subprogram  m_routine ;
  void                   * m_argument ;
  int                    * m_lock ;
  int                      m_lock_size ;
  int                      m_size ;
  int                      m_rank ;
  int                      m_group_size ;
  int                      m_group_rank ;
} ThreadWork ;

typedef struct ThreadPool_Data {
  ThreadWork    * m_work_begin ;
  int             m_work_size ;
  int             m_number_threads ;
  int             m_number_locks ;
  int             m_work_count ;
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
  int result = NULL == local ? TPI_ERROR_NULL : 0 ;

  if ( ! result ) {
    if ( i < 0 || local->m_lock_size <= i ) {
      result = TPI_ERROR_SIZE ;
    }
    else if ( local->m_lock[i] ) {
      result = TPI_ERROR_LOCK ;
    }
    else {
      local->m_lock[i] = 1 ;
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
    else if ( local->m_lock[i] ) {
      result = TPI_LOCK_BUSY ;
    }
    else {
      local->m_lock[i] = 1 ;
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
    else if ( ! local->m_lock[i] ) {
      result = TPI_ERROR_LOCK ;
    }
    else {
      local->m_lock[i] = 0 ;
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
                                       const unsigned     number )
{
  int counter = 0 ;
  unsigned i = 0 ;

  for ( ; i < number ; ++i ) {
    ThreadWork * const w = work + i ;
    (*w->m_routine)( w->m_argument , w );
  }
  return number ;
}

/*--------------------------------------------------------------------*/
/* The work queue shared by all threads */

static ThreadPool * local_thread_pool()
{
  static ThreadPool pool = {
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

      int lock[ number_locks ];

      int nlocks = 0 ;

      while ( nlocks < number_locks && ! result ) {
        lock[ nlocks ] = 0 ;
        ++nlocks ;
      }

      if ( ! result ) {

        ThreadWork work[ number_total ];

        { /* Fill the work queue */
          ThreadWork * w = work ;

          for ( i = 0 ; i < number_routine ; ++i ) {
            int k = 0 ;
            for ( ; k < number[i] ; ++k , ++w ) {
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

        /* Participate in the work */
        pool->m_work_count += local_thread_pool_run_work(work,number_total);

        pool->m_work_begin = NULL ;
        pool->m_work_size  = 0 ;
        pool->m_number_locks = 0 ;
      }
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

    if ( number ) {
      *count = pool->m_work_count ;
      while ( --number ) { *++count = 0 ; }
    }

    pool->m_work_count = 0 ;
  }
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Init( int n )
{
  ThreadPool * const pool = local_thread_pool();

  int result = ! pool || pool->m_number_threads ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && n <= 0 ) { result = TPI_ERROR_SIZE ; }

  if ( ! result ) { pool->m_number_threads = n ; }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  ThreadPool * const pool = local_thread_pool();

  int result = ! pool ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result ) { pool->m_number_threads = 0 ; }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Size( int * number_allocated , int * number_concurrent )
{
  int result = 0 ;

  if ( number_concurrent ) { *number_concurrent = 0 ; }

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

