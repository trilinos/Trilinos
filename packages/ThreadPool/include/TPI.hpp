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

#include <string>
#include <stdexcept>

#include <util/TPI.h>

namespace TPI {

typedef TPI_ThreadPool ThreadPool ;

inline
int Run( void (*)(void*,ThreadPool), void * );

/** Run 'void Worker::method( TPI::ThreadPool )' */
template<class Worker>
inline
int Run( Worker & , void (Worker::*)(ThreadPool) );

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
int Partition( ThreadPool pool , int N , int & I_local , int & N_local )
  { return TPI_Partition( pool , N , & I_local , & N_local ); }

//----------------------------------------------------------------------

int Init( int n ) { return TPI_Init( n ); }

int Finalize() { return TPI_Finalize(); }

int Size( int & number_allocated , int & number_concurrent )
  { return TPI_Size( & number_allocated , & number_concurrent ); }

double Walltime() { return TPI_Walltime(); }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

template<class Worker>
class WorkerMethod {
private:
  WorkerMethod();
  WorkerMethod( const WorkerMethod & );
  WorkerMethod & operator = ( const WorkerMethod & );

public:

  typedef void (Worker::*Method)( ThreadPool );

  Worker & worker ;
  Method   method ;

  WorkerMethod( Worker & w , Method m ) : worker(w), method(m) {}

  static void run( void * arg , ThreadPool pool )
    {
      try {
        WorkerMethod & wm = * reinterpret_cast<WorkerMethod*>(arg);
        (wm.worker.*wm.method)(pool);
      } catch(...){}
    }
};

}

inline
void Run( ThreadPool pool ,
          void (*func)( void * , ThreadPool ) ,
          void * arg )
{
  TPI_Run( pool, reinterpret_cast< TPI_parallel_subprogram >(func), arg );
}

template<class Worker>
inline
int Run( Worker & worker, void (Worker::*method)(ThreadPool))
{
  typedef WorkerMethod<Worker> WM ;

  WM tmp( worker , method );

  return TPI_Run( reinterpret_cast<TPI_parallel_subprogram>(& WM::run), &tmp);
}

//----------------------------------------------------------------------

template<class W0>
inline
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 )
{
  typedef WorkerMethod<W0> WM_0 ;

  WM_0 tmp1( worker_0 , method_0 );

  TPI_parallel_subprogram sub[1] = {
    reinterpret_cast<TPI_parallel_subprogram>(& WM_0::run) };

  void * arg[1] = { & tmp_0 };
  int    num[1] = {   num_0 };

  return TPI_Run_many( 1 , sub , arg , num );
}

template<class W0, class W1>
inline
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 )
{
  typedef WorkerMethod<W0> WM_0 ;
  typedef WorkerMethod<W1> WM_1 ;

  WM_0 tmp1( worker_0 , method_0 );
  WM_1 tmp1( worker_1 , method_1 );

  TPI_parallel_subprogram sub[2] = {
    reinterpret_cast<TPI_parallel_subprogram>(& WM_0::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_1::run) };

  void * arg[2] = { & tmp_0 , & tmp_1 };
  int    num[2] = {   num_0 ,   num_1 };

  return TPI_Run_many( 2 , sub , arg , num );
}

template<class W0, class W1, class W3>
inline
int Run( W0 & worker_0 , void (W0::*method_0)(ThreadPool) , int num_0 ,
         W1 & worker_1 , void (W1::*method_1)(ThreadPool) , int num_1 ,
         W2 & worker_2 , void (W2::*method_2)(ThreadPool) , int num_2 )
{
  typedef WorkerMethod<W0> WM_0 ;
  typedef WorkerMethod<W1> WM_1 ;
  typedef WorkerMethod<W2> WM_2 ;

  WM_0 tmp1( worker_0 , method_0 );
  WM_1 tmp1( worker_1 , method_1 );
  WM_2 tmp1( worker_2 , method_2 );

  TPI_parallel_subprogram sub[3] = {
    reinterpret_cast<TPI_parallel_subprogram>(& WM_0::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_1::run) ,
    reinterpret_cast<TPI_parallel_subprogram>(& WM_2::run) };

  void * arg[3] = { & tmp_0 , & tmp_1 , & tmp_2 };
  int    num[3] = {   num_0 ,   num_1 ,   num_2 };

  return TPI_Run_many( 3 , sub , arg , num );
}

//----------------------------------------------------------------------

}

#endif

