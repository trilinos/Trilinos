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
  int * m_lock ;
  int   m_lock_size ;
  int   m_size ;
  int   m_rank ;
} ThreadWork ;

typedef struct ThreadPool_Data {
  int  m_work_size ;
  int  m_number_locks ;
} ThreadPool ;

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

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

int TPI_Rank( TPI_ThreadPool local , int * rank , int * size )
{
  const int result = NULL == local ? TPI_ERROR_NULL : 0 ;

  if ( ! result ) {
    if ( rank ) { *rank = local->m_rank ; }
    if ( size ) { *size = local->m_size ; }
  }
  return result ;
}

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
/* The work queue shared by all threads */

static ThreadPool * local_thread_pool()
{
  static ThreadPool pool = {
    /* m_work_size      */  0 ,
    /* m_number_locks   */  0 };

  /* Guard against recursive call */

  return pool.m_work_size ? NULL : & pool ;
}

/*--------------------------------------------------------------------*/

int TPI_Run( TPI_parallel_subprogram routine ,
             void * routine_data ,
             int number )
{
  int i ;

  int result = ! routine || ! routine_data ? TPI_ERROR_NULL : 0 ;

  if ( ! result ) {

    ThreadPool * const pool = local_thread_pool();

    if ( ! pool ) { result = TPI_ERROR_ACTIVE ; }

    if ( ! result ) {

      const int number_locks = pool->m_number_locks ;

      int lock[ number_locks ];

      int nlocks = 0 ;

      while ( nlocks < number_locks && ! result ) {
        lock[ nlocks ] = 0 ;
        ++nlocks ;
      }

      if ( ! result ) {
        if ( number <= 0 ) { number = 1 ; }

        ThreadWork work = { lock , number_locks , number , 0 };

        pool->m_work_size = number ;

        for ( i = 0 ; i < number ; ++i ) {
          work.m_rank = i ;
          (*routine)( routine_data , & work );
        }
        pool->m_work_size = 0 ;
      }
    }
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

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Init( int n )
{
  ThreadPool * const pool = local_thread_pool();

  int result = ! pool ? TPI_ERROR_ACTIVE : 0 ;

  if ( ! result && n <= 0 ) { result = TPI_ERROR_SIZE ; }

  return result ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

int TPI_Finalize()
{
  ThreadPool * const pool = local_thread_pool();

  int result = ! pool ? TPI_ERROR_ACTIVE : 0 ;

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
      *number_allocated = 1 ;
    }
    else {
      result = TPI_ERROR_ACTIVE ;
    }
  }

  return result ;
}

