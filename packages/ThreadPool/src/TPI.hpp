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
inline int Trylock( int n ) { return TPI_Trylock( n ); }
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

