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
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 *
 *  Thread Pool Interface (TPI).
 *
 *  A simple and miminalistic interface for executing subprograms
 *  from N parallel threads.
 *
 *  States: the underlying thread pool has three states.
 *    1) Uninitialized: no extra threads exist, this is the initial state.
 *    2) Paused: N-1 extra threads exists and are blocked.
 *    3) Active: N-1 extra threads are calling the subprograms.
 *
 *  Threads are created and blocked ready-for-use so that the thread
 *  creation cost is paid only once.  This cost is significantly larger
 *  than unblocking an already created thread.
 *
 *  The extra threads are blocked, as opposed to spinning, in the paused
 *  state so that they do not compete for compute cycles with other threads
 *  created external to the TPI interface.  For example, threads created
 *  by OpenMP.
 */

#ifndef ThreadPoolInterface_h
#define ThreadPoolInterface_h

#if defined( __cplusplus )
extern "C" {
#endif

/*--------------------------------------------------------------------*/
/** A utility to measure wall-clock time, which is frequently
 *  needed when performance testing HPC algorithms.
 */
double TPI_Walltime();

/*--------------------------------------------------------------------*/
/* All functions return zero for success. */

#define TPI_LOCK_BUSY      ((int)  1  /* trylock or unlock failed */ )
#define TPI_ERROR_NULL     ((int) -1  /* NULL input */ )
#define TPI_ERROR_SIZE     ((int) -2  /* BAD input: size or index */ )
#define TPI_ERROR_LOCK     ((int) -3  /* BAD lock or unlock */ )
#define TPI_ERROR_ACTIVE   ((int) -4  /* BAD input: the pool is active  */ )
#define TPI_ERROR_INTERNAL ((int) -5  /* internal resource error */ )

/*--------------------------------------------------------------------*/

struct TPI_ThreadPool_Private ; 

typedef struct TPI_ThreadPool_Private * TPI_ThreadPool ;

/*--------------------------------------------------------------------*/
/**  The interface for a thread pool parallel subprogram.
 *   When these subprograms are called the thread pool is
 *   in the 'active' state.
 */
typedef void (*TPI_parallel_subprogram)( void * shared_data ,
                                         TPI_ThreadPool pool );

/** Run a thread pool parallel subprogram on all TPI-created threads.
 *  Each thread in the pool will call the subprogram as:
 *
 *    (*subprogram)( shared_data , pool )
 *
 *  The thread pool must be in the 'paused' state when this
 *  function is called.  Thus recursive calls to TPI_Run or
 *  TPI_Run_many are illegal.
 */
int TPI_Run( TPI_parallel_subprogram /* subprogram  */ ,
             void *                  /* shared data */ );

/** Query the subprogram's execution "rank within size."
 *  The suprogram is called 'size' number of times where the
 *  current call is uniquely identified by 'rank', where
 *  0 <= rank < size.  
 *
 *  The thread pool must be in the active state.
 *  Thus this method may only be called from within
 *  a TPI_parallel_subprogram.
 */
int TPI_Rank( TPI_ThreadPool , int * /* rank */ , int * /* size */ );

/*--------------------------------------------------------------------*/
/** Run many possibly different thread pool parallel subprograms.
 *  Each subprogram is called within its own pool of the given size.
 *  If the sum of pool sizes is larger than the number of created
 *  threads then the threads use work-stealing to call 
 *  'subprogram[i]' 'size[i]' times.
 *
 *  The thread pool must be in the 'paused' state when this
 *  function is called.  Thus recursive calls to TPI_Run or
 *  TPI_Run_many are illegal.
 */
int TPI_Run_many( const int                  /* Number of groups */ ,
                  TPI_parallel_subprogram [] /* group of subprograms  */ ,
                  void * []                  /* group of subprograms' data */ ,
                  const int []               /* group of pool sizes   */ );

/** Query the current thread's "rank within size of group of pools."
 *  The 'TPI_Run_many' executes multiple subprograms within defined
 *  groups of a given size.  This method queries in which group the
 *  current call is taking place.
 *
 *  The thread pool must be in the active state.
 *  Thus this method may only be called from within
 *  a TPI_parallel_subprogram.
 */
int TPI_Group_rank( TPI_ThreadPool , int * /* rank */ , int * /* size */ );

/*--------------------------------------------------------------------*/
/** Query accumulated run counts for each created thread and then
 *  reset the counts to zero.  The run count is the number of
 *  times that a thread performed work since the thread was created
 *  or since the previous call to 'TPI_Run_count'.
 *
 *  This provides a measure of parallelism in that the work
 *  requested by 'TPI_Run' or 'TPI_Run_many' may be completed 
 *  by fewer threads than have been created.  This will occur
 *  if a thread fails to unblock before another created thread
 *  'steals' the blocked thread's work.
 */
int TPI_Run_count( int * /* count */ );

/*--------------------------------------------------------------------*/
/** Set the lock size for the next call to TPI_Run or TPI_Run_many.
 *  The thread pool must be in the 'paused' state.
 */
int TPI_Set_lock_size( int );

/** Blocks until lock # is obtained.
 *  The thread pool must be in the 'active' state.
 */
int TPI_Lock( TPI_ThreadPool , int );

/** Tries to lock #, returns 0 if successful. 
 *  The thread pool must be in the 'active' state.
 */
int TPI_Trylock( TPI_ThreadPool , int );

/** Unlocks lock #.
 *  The thread pool must be in the 'active' state.
 */
int TPI_Unlock( TPI_ThreadPool , int );

/*--------------------------------------------------------------------*/
/** Initialize the root thread pool to the specified size.
 *  The thread pool is transitioned from 'uninitialized' to 'paused' state.
 */
int TPI_Init( int /* size */ );

/** Query the number of threads created.
 *  The thread pool must be in the 'paused' state.
 */
int TPI_Size( int * );

/** Finalize (shut down) all threads and thread pools.
 *  The thread pool is transitioned from 'paused' to 'uninitialized' state.
 */
int TPI_Finalize();

/*--------------------------------------------------------------------*/
/** Query the number of threads known to be concurrently schedulable.
 *  If the concurrency is not detectable then return 0.
 */
int TPI_Concurrency();

/*--------------------------------------------------------------------*/
/** A frequently occuring parallel operation is to partition a
 *  collection of Number of items among Thread_Size threads.
 *  Each thread has an approximately equal span of Number
 *  defined by [ Local_Begin , Local_Begin + Local_Number - 1 ].
 */
int TPI_Partition( int   /* Thread_Rank */ ,
                   int   /* Thread_Size */ ,
                   int   /* Number */ ,
                   int * /* Local_Begin */ ,
                   int * /* Local_Number */ );

/*--------------------------------------------------------------------*/

#if defined( __cplusplus )
}
#endif

#endif

