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
int Run( void (*func)(void*,ThreadPool), void * arg );

/** Run  worker.*method(pool)  on all threads.
 */
template<class Worker>
int Run( Worker & worker , void (Worker::*method)(ThreadPool) );

//----------------------------------------------------------------------
/* Run  worker_#.*method_#(pool)  on  num_#  threads. */

template<class W0>
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 );

template<class W0, class W1>
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 );

template<class W0, class W1, class W2>
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 ,
         W2 & worker_2 , void (W2::*method_2)(ThreadPool) , int num_2 );

template<class W0, class W1, class W2, class W3>
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 ,
         W2 & worker_2 , void (W2::*method_2)(ThreadPool) , int num_2 ,
         W3 & worker_3 , void (W3::*method_3)(ThreadPool) , int num_3 );

template<class W0, class W1, class W2, class W3, class W4>
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 ,
         W2 & worker_2 , void (W2::*method_2)(ThreadPool) , int num_2 ,
         W3 & worker_3 , void (W3::*method_3)(ThreadPool) , int num_3 ,
         W4 & worker_4 , void (W3::*method_4)(ThreadPool) , int num_4 );

template<class W0, class W1, class W2, class W3, class W4, class W5>
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 ,
         W2 & worker_2 , void (W2::*method_2)(ThreadPool) , int num_2 ,
         W3 & worker_3 , void (W3::*method_3)(ThreadPool) , int num_3 ,
         W4 & worker_4 , void (W3::*method_4)(ThreadPool) , int num_4 ,
         W5 & worker_5 , void (W5::*method_5)(ThreadPool) , int num_5 );

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
int Run( void (*func)( void * , ThreadPool ) , void * arg )
{
  return TPI_Run( reinterpret_cast< TPI_parallel_subprogram >(func), arg );
}

template<class Worker>
inline
int Run( Worker & worker, void (Worker::*method)(ThreadPool) )
{
  typedef WorkerMethodHelper<Worker> WM ;

  WM tmp( worker , method );

  return TPI_Run( reinterpret_cast<TPI_parallel_subprogram>(& WM::run), &tmp);
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<class W0>
inline
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 )
{
  enum { N = 1 };
  typedef WorkerMethodHelper<W0> WM_0 ;

  WM_0 tmp_0( worker_0 , method_0 );

  TPI_parallel_subprogram sub[N] = {
    reinterpret_cast<TPI_parallel_subprogram>(& WM_0::run) };

  void * arg[N] = { & tmp_0 };
  int    num[N] = {   num_0 };

  return TPI_Run_many( N , sub , arg , num );
}

//----------------------------------------------------------------------

template<class W0, class W1>
inline
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 )
{
  enum { N = 2 };
  typedef WorkerMethodHelper<W0> WM_0 ;
  typedef WorkerMethodHelper<W1> WM_1 ;

  WM_0 tmp_0( worker_0 , method_0 );
  WM_1 tmp_1( worker_1 , method_1 );

  TPI_parallel_subprogram sub[N] = {
    reinterpret_cast<TPI_parallel_subprogram>(& WM_0::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_1::run) };

  void * arg[N] = { & tmp_0 , & tmp_1 };
  int    num[N] = {   num_0 ,   num_1 };

  return TPI_Run_many( N , sub , arg , num );
}

//----------------------------------------------------------------------

template<class W0, class W1, class W2>
inline
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 ,
         W2 & worker_2 , void (W2::*method_2)(ThreadPool) , int num_2 )
{
  enum { N = 3 };
  typedef WorkerMethodHelper<W0> WM_0 ;
  typedef WorkerMethodHelper<W1> WM_1 ;
  typedef WorkerMethodHelper<W2> WM_2 ;

  WM_0 tmp_0( worker_0 , method_0 );
  WM_1 tmp_1( worker_1 , method_1 );
  WM_2 tmp_2( worker_2 , method_2 );

  TPI_parallel_subprogram sub[N] = {
    reinterpret_cast<TPI_parallel_subprogram>(& WM_0::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_1::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_2::run) };

  void * arg[N] = { & tmp_0 , & tmp_1 , & tmp_2 };
  int    num[N] = {   num_0 ,   num_1 ,   num_2 };

  return TPI_Run_many( N , sub , arg , num );
}

//----------------------------------------------------------------------

template<class W0, class W1, class W2, class W3>
inline
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 ,
         W2 & worker_2 , void (W2::*method_2)(ThreadPool) , int num_2 ,
         W3 & worker_3 , void (W3::*method_3)(ThreadPool) , int num_3 )
{
  enum { N = 4 };
  typedef WorkerMethodHelper<W0> WM_0 ;
  typedef WorkerMethodHelper<W1> WM_1 ;
  typedef WorkerMethodHelper<W2> WM_2 ;
  typedef WorkerMethodHelper<W3> WM_3 ;

  WM_0 tmp_0( worker_0 , method_0 );
  WM_1 tmp_1( worker_1 , method_1 );
  WM_2 tmp_2( worker_2 , method_2 );
  WM_3 tmp_3( worker_3 , method_3 );

  TPI_parallel_subprogram sub[N] = {
    reinterpret_cast<TPI_parallel_subprogram>(& WM_0::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_1::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_2::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_3::run) };

  void * arg[N] = { & tmp_0 , & tmp_1 , & tmp_2 , & tmp_3 };
  int    num[N] = {   num_0 ,   num_1 ,   num_2 ,   num_3 };

  return TPI_Run_many( N , sub , arg , num );
}

//----------------------------------------------------------------------

template<class W0, class W1, class W2, class W3, class W4>
inline
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 ,
         W2 & worker_2 , void (W2::*method_2)(ThreadPool) , int num_2 ,
         W3 & worker_3 , void (W3::*method_3)(ThreadPool) , int num_3 ,
         W4 & worker_4 , void (W3::*method_4)(ThreadPool) , int num_4 )
{
  enum { N = 5 };
  typedef WorkerMethodHelper<W0> WM_0 ;
  typedef WorkerMethodHelper<W1> WM_1 ;
  typedef WorkerMethodHelper<W2> WM_2 ;
  typedef WorkerMethodHelper<W3> WM_3 ;
  typedef WorkerMethodHelper<W4> WM_4 ;

  WM_0 tmp_0( worker_0 , method_0 );
  WM_1 tmp_1( worker_1 , method_1 );
  WM_2 tmp_2( worker_2 , method_2 );
  WM_3 tmp_3( worker_3 , method_3 );
  WM_4 tmp_4( worker_4 , method_4 );

  TPI_parallel_subprogram sub[N] = {
    reinterpret_cast<TPI_parallel_subprogram>(& WM_0::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_1::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_2::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_3::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_4::run) };

  void * arg[N] = { & tmp_0 , & tmp_1 , & tmp_2 , & tmp_3 , & tmp_4 };
  int    num[N] = {   num_0 ,   num_1 ,   num_2 ,   num_3 ,   num_4 };

  return TPI_Run_many( N , sub , arg , num );
}

//----------------------------------------------------------------------

template<class W0, class W1, class W2, class W3, class W4, class W5>
inline
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 ,
         W2 & worker_2 , void (W2::*method_2)(ThreadPool) , int num_2 ,
         W3 & worker_3 , void (W3::*method_3)(ThreadPool) , int num_3 ,
         W4 & worker_4 , void (W3::*method_4)(ThreadPool) , int num_4 ,
         W5 & worker_5 , void (W5::*method_5)(ThreadPool) , int num_5 )
{
  enum { N = 6 };
  typedef WorkerMethodHelper<W0> WM_0 ;
  typedef WorkerMethodHelper<W1> WM_1 ;
  typedef WorkerMethodHelper<W2> WM_2 ;
  typedef WorkerMethodHelper<W3> WM_3 ;
  typedef WorkerMethodHelper<W4> WM_4 ;
  typedef WorkerMethodHelper<W4> WM_5 ;

  WM_0 tmp_0( worker_0 , method_0 );
  WM_1 tmp_1( worker_1 , method_1 );
  WM_2 tmp_2( worker_2 , method_2 );
  WM_3 tmp_3( worker_3 , method_3 );
  WM_4 tmp_4( worker_4 , method_4 );
  WM_5 tmp_5( worker_5 , method_5 );

  TPI_parallel_subprogram sub[N] = {
    reinterpret_cast<TPI_parallel_subprogram>(& WM_0::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_1::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_2::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_3::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_4::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_5::run) };

  void * arg[N] = { & tmp_0, & tmp_1, & tmp_2, & tmp_3, & tmp_4, & tmp_5 };
  int    num[N] = {   num_0,   num_1,   num_2,   num_3,   num_4,   num_5 };

  return TPI_Run_many( N , sub , arg , num );
}

//----------------------------------------------------------------------

}

#endif

