/*------------------------------------------------------------------------*/
/*                    TPI: Thread Pool Interface                          */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/* Redistribution and use in source and binary forms, with or without     */
/* modification, are permitted provided that the following conditions are */
/* met:                                                                   */
/*                                                                        */
/* 1. Redistributions of source code must retain the above copyright      */
/* notice, this list of conditions and the following disclaimer.          */
/*                                                                        */
/* 2. Redistributions in binary form must reproduce the above copyright   */
/* notice, this list of conditions and the following disclaimer in the    */
/* documentation and/or other materials provided with the distribution.   */
/*                                                                        */
/* 3. Neither the name of the Corporation nor the names of the            */
/* contributors may be used to endorse or promote products derived from   */
/* this software without specific prior written permission.               */
/*                                                                        */
/* THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY        */
/* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE      */
/* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR     */
/* PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE    */
/* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,    */
/* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR     */
/* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF */
/* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING   */
/* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS     */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.           */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 *
 *  Thread Pool Interface (TPI).
 *
 *  A simple and miminalistic interface for executing subprograms
 *  in a thread parallel, shared memory mode.
 *
 *  States: the underlying thread pool has four states.
 *    1) Uninitialized: no extra threads exist, this is the initial state.
 *    2) Ready:  extra threads exist and are ready to run a subprogram.
 *    3) Active: extra threads are calling the subprogram.
 *    4) Blocked: extra threads blocked.
 *
 *  Threads are created on initialization and placed in the 'Ready' state.
 *  While in the 'Ready' state the threads are spin-waiting to minimize
 *  the cost of activating blocked threads.
 *  Threads can be blocked so that they do not compete for computatational
 *  resources with other threads created external to the TPI interface.
 *  For example, threads created by OpenMP or TBB.
 */

#ifndef ThreadPoolInterface_h
#define ThreadPoolInterface_h

#if defined( __cplusplus )
extern "C" {
#endif

/*--------------------------------------------------------------------*/
/** \brief  Version string. */
const char * TPI_Version();

/** Start up the requested number of threads, less the calling thread.
 *  Return the actual number of threads, including the calling thread,
 *  otherwise return an error.
 */
int TPI_Init( int thread_count );

/** Shut down all started threads. */
int TPI_Finalize();

/*--------------------------------------------------------------------*/
/** \brief  A utility to measure wall-clock time, which is frequently
 *          needed when performance testing HPC algorithms.
 */
double TPI_Walltime();

/*--------------------------------------------------------------------*/
/* All functions return zero for success. */

#define TPI_ERROR_NULL     ((int) -1)  /**<  NULL input */
#define TPI_ERROR_SIZE     ((int) -2)  /**<  BAD input: size or index */
#define TPI_ERROR_LOCK     ((int) -3)  /**<  BAD lock or unlock */
#define TPI_ERROR_ACTIVE   ((int) -4)  /**<  BAD input: the pool is active  */
#define TPI_ERROR_INTERNAL ((int) -5)  /**< internal resource error */

/*--------------------------------------------------------------------*/
/** \brief  Work information passed to a work subprogram. */
struct TPI_Work_Struct {
  const void * info ;       /**<  Shared info input to TPI_Run */
  void       * reduce ;     /**<  Data for reduce operation, if any */
  int          count ;      /**<  Count of work requested via TPI_Run */
  int          rank ;       /**<  Rank  of work for the current call */
  int          lock_count ; /**<  Count of locks requested via TPI_Run */
};

/** \brief  Typedef for work subprogram argument */
typedef const struct TPI_Work_Struct TPI_Work ;

/**  The interface for a parallel task */
typedef void (*TPI_work_subprogram)( TPI_Work * );

/**  The interface for a parallel reduction operation.
 *   Initialize  work->reduce value.
 */
typedef
void (*TPI_reduce_init)( TPI_Work * work );

/**  The interface for a parallel reduction operation.
 *   Perform reduction operation  work->reduce OP= reduce.
 *   Every initialized reduce value will appear exactly
 *   once as the 'reduce' argument of a call to the join function.
 */
typedef
void (*TPI_reduce_join)( TPI_Work * work , const void * reduce );

/*--------------------------------------------------------------------*/
/** \brief Run a work subprogram in thread parallel.
 *
 *  The thread pool must be in the 'paused' state when this
 *  function is called.  Thus a recursive call to TPI_Run is illegal.
 */
int TPI_Run( TPI_work_subprogram work_subprogram  ,
             const void *        work_info ,
             int                 work_count  ,
             int                 lock_count );

/** \brief Run a work and reduction subprograms in thread parallel.
 *
 *  Each call to the work_subprogram has exclusive (thread safe)
 *  access to its work->reduce data.
 *  The reduce_init and reduce_join subprograms have
 *  exclusive access to their arguments.
 */
int TPI_Run_reduce( TPI_work_subprogram   work_subprogram  ,
                    const void *          work_info ,
                    int                   work_count  ,
                    TPI_reduce_join       reduce_join ,
                    TPI_reduce_init       reduce_init ,
                    int                   reduce_size ,
                    void *                reduce_data );

/** \brief  Run a work subprogram exactly once on each thread.
 *
 *  The thread pool must be in the 'paused' state when this
 *  function is called.  Thus a recursive call to TPI_Run is illegal.
 */
int TPI_Run_threads( TPI_work_subprogram work_subprogram ,
                     const void *        work_info ,
                     int                 lock_count );

/** \brief Run a work and reduction subprograms in thread parallel.
 *
 *  Each call to the work_subprogram has exclusive (thread safe)
 *  access to its work->reduce data.
 *  The reduce_init and reduce_join subprograms have
 *  exclusive access to their arguments.
 */
int TPI_Run_threads_reduce( TPI_work_subprogram   work_subprogram ,
                            const void *          work_info ,
                            TPI_reduce_join       reduce_join ,
                            TPI_reduce_init       reduce_init ,
                            int                   reduce_size ,
                            void *                reduce_data );

/*--------------------------------------------------------------------*/
/** \brief  Start a work subprogram in thread parallel
 *          running on all but the 'main' calling thread;
 *          the 'main' calling thread returns immediately.
 *
 *  The thread pool must be in the 'paused' state when this
 *  function is called.  Thus a recursive call to TPI_Start is illegal.
 */
int TPI_Start( TPI_work_subprogram work_subprogram  ,
               const void *        work_info ,
               int                 work_count  ,
               int                 lock_count );

/** \brief  Start a work and reduction subprograms in thread parallel
 *          running on all but the 'main' calling thread;
 *          the 'main' calling thread returns immediately.
 *
 *  Each call to the work_subprogram has exclusive (thread safe)
 *  access to its work->reduce data.
 *  The reduce_init and reduce_join subprograms have
 *  exclusive access to their arguments.
 */
int TPI_Start_reduce( TPI_work_subprogram   work_subprogram  ,
                      const void *          work_info ,
                      int                   work_count  ,
                      TPI_reduce_join       reduce_join ,
                      TPI_reduce_init       reduce_init ,
                      int                   reduce_size ,
                      void *                reduce_data );

/** \brief  Run a work subprogram on each thread
 *          that is not the 'main' calling thread.
 *          The 'main' calling thread returns immediately.
 *
 *  The thread pool must be in the 'paused' state when this
 *  function is called.  Thus a recursive call to TPI_Start_threads is illegal.
 */
int TPI_Start_threads( TPI_work_subprogram work_subprogram ,
                       const void *        work_info ,
                       int                 lock_count );

/** \brief  Start a work / reduction subprogram 
 *          on each thread that is not the 'main' calling thread.
 *          The 'main' calling thread returns immediately.
 *
 *  Each call to the work_subprogram has exclusive (thread safe)
 *  access to its work->reduce data.
 *  The reduce_init and reduce_join subprograms have
 *  exclusive access to their arguments.
 */
int TPI_Start_threads_reduce( TPI_work_subprogram   work_subprogram ,
                              const void *          work_info ,
                              TPI_reduce_join       reduce_join ,
                              TPI_reduce_init       reduce_init ,
                              int                   reduce_size ,
                              void *                reduce_data );

/** \brief  Wait for a started work subprogram to complete. */
int TPI_Wait();

/*--------------------------------------------------------------------*/
/** \brief  Block threads within the operating system.
 *
 *  Normally the worker threads are unblocked and spinning for
 *  minimal start up overhead when running work subprograms.
 *  If no TPI work is to be performed for a long period of time
 *  then an application can block the worker threads.
 */
int TPI_Block();

/** \brief  Unblock blocked threads within the operating system */
int TPI_Unblock();

/** \brief  Query if threads are blocked */
int TPI_Isblocked();

/*--------------------------------------------------------------------*/
/** \brief  Blocks until lock lock_rank is obtained.
 *          The thread pool must be in the 'active' state.
 */
int TPI_Lock( int lock_rank );

/** \brief  Unlocks lock lock_rank.
 *          The thread pool must be in the 'active' state.
 */
int TPI_Unlock( int lock_rank );

/*--------------------------------------------------------------------*/

#if defined( __cplusplus )
}
#endif

#endif

