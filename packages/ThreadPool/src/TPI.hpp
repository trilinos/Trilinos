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
 */

#ifndef util_ThreadPool_hpp
#define util_ThreadPool_hpp

#include <TPI.h>

namespace TPI {

typedef TPI_Work Work ;

//----------------------------------------------------------------------
/** Run  worker.*method(work)  on all threads.
 */
template<class Worker>
int Run( Worker & worker , void (Worker::*method)(Work &) ,
         int work_count , int lock_count = 0 );

//----------------------------------------------------------------------

inline int Lock( int n )    { return TPI_Lock( n ); }
inline int Unlock( int n )  { return TPI_Unlock( n ); }

/** Lock guard to insure that a lock is released
 *  when control exists a block.
 *    {
 *      TPI::LockGuard local_lock( i );
 *    }
 */
class LockGuard {
private:
  LockGuard();
  LockGuard( const LockGuard & );
  LockGuard & operator = ( const LockGuard & );
  const int m_value ;
  const int m_result ;
public:
  operator int() const { return m_result ; }

  explicit LockGuard( unsigned i_lock )
    : m_value( i_lock ), m_result( TPI_Lock(i_lock) ) {}

  ~LockGuard() { TPI_Unlock( m_value ); }
};

//----------------------------------------------------------------------

inline
int Init( int n ) { return TPI_Init( n ); }

inline
int Finalize() { return TPI_Finalize(); }

inline
double Walltime() { return TPI_Walltime(); }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

template<class Worker>
class WorkerMethodHelper {
private:
  WorkerMethodHelper();
  WorkerMethodHelper( const WorkerMethodHelper & );
  WorkerMethodHelper & operator = ( const WorkerMethodHelper & );

public:

  typedef void (Worker::*Method)( Work & );

  Worker & worker ;
  Method   method ;

  WorkerMethodHelper( Worker & w , Method m ) : worker(w), method(m) {}

  static void run( TPI_Work * work )
    {
      try {
        const WorkerMethodHelper & wm =
          * reinterpret_cast<const WorkerMethodHelper*>(work->info);
        (wm.worker.*wm.method)(*work);
      } catch(...){}
    }
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<class Worker>
inline
int Run( Worker & worker, void (Worker::*method)(Work &) ,
         int work_count , int lock_count )
{
  typedef WorkerMethodHelper<Worker> WM ;

  WM tmp( worker , method );

  return TPI_Run( reinterpret_cast<TPI_work_subprogram>(& WM::run),&tmp,work_count,lock_count);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}

#endif

