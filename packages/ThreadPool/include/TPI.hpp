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

typedef TPI_ThreadPool ThreadPool ;

//----------------------------------------------------------------------
/** Run  (*func)(arg,pool)  on all threads.
 */
int Run( void (*func)(void*,ThreadPool), void * arg , int = 0 );

/** Run  worker.*method(pool)  on all threads.
 */
template<class Worker>
int Run( Worker & worker , void (Worker::*method)(ThreadPool) , int = 0 );

//----------------------------------------------------------------------

inline
int Set_lock_size( int n ) { return TPI_Set_lock_size( n ); }

inline
int Lock( ThreadPool pool , int n ) { return TPI_Lock( pool , n ); }

inline
int Trylock( ThreadPool pool , int n ) { return TPI_Trylock( pool , n ); }

inline
int Unlock( ThreadPool pool , int n ) { return TPI_Unlock( pool , n ); }

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
  const ThreadPool m_pool ;
  const int        m_value ;
  const int        m_result ;
public:
  operator int() const { return m_result ; }

  explicit LockGuard( ThreadPool pool , unsigned i_lock )
    : m_pool( pool ), m_value( i_lock ), m_result( TPI_Lock(pool,i_lock) ) {}

  ~LockGuard() { TPI_Unlock( m_pool , m_value ); }
};

//----------------------------------------------------------------------

inline
int Rank( ThreadPool pool , int & rank , int & size )
  { return TPI_Rank( pool , & rank , & size ); }

inline
int Partition( int Rank , int Size , int N , int & I_local , int & N_local )
  { return TPI_Partition( Rank , Size , N , & I_local , & N_local ); }

//----------------------------------------------------------------------

inline
int Init( int n ) { return TPI_Init( n ); }

inline
int Finalize() { return TPI_Finalize(); }

inline
int Size( int & number_allocated ) { return TPI_Size( & number_allocated ); }

inline
int Concurrency() { return TPI_Concurrency(); }

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

  typedef void (Worker::*Method)( ThreadPool );

  Worker & worker ;
  Method   method ;

  WorkerMethodHelper( Worker & w , Method m ) : worker(w), method(m) {}

  static void run( void * arg , ThreadPool pool )
    {
      try {
        WorkerMethodHelper & wm = * reinterpret_cast<WorkerMethodHelper*>(arg);
        (wm.worker.*wm.method)(pool);
      } catch(...){}
    }
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

inline
int Run( void (*func)( void * , ThreadPool ) , void * arg , int n )
{
  return TPI_Run( reinterpret_cast< TPI_parallel_subprogram >(func), arg , n );
}

template<class Worker>
inline
int Run( Worker & worker, void (Worker::*method)(ThreadPool) , int n )
{
  typedef WorkerMethodHelper<Worker> WM ;

  WM tmp( worker , method );

  return TPI_Run( reinterpret_cast<TPI_parallel_subprogram>(& WM::run),&tmp,n);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}

#endif

