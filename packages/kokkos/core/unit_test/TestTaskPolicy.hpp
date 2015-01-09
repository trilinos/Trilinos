/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/


#ifndef KOKKOS_UNITTEST_TASKPOLICY_HPP
#define KOKKOS_UNITTEST_TASKPOLICY_HPP

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <Kokkos_TaskPolicy.hpp>

namespace TestTaskPolicy {

//----------------------------------------------------------------------------

template< class ExecSpace >
struct FibChild {

  typedef long value_type ;

  Kokkos::TaskPolicy<ExecSpace> policy ;
  const value_type n ;
  int has_nested ;

  FibChild( const Kokkos::TaskPolicy<ExecSpace> & arg_policy
          , const value_type arg_n )
    : policy(arg_policy,2) /* default dependence capacity = 2 */
    , n( arg_n ), has_nested(0) {}

  inline
  void apply( value_type & result )
    {
      if ( n < 2 ) {
        result = n ;
      }
      else {
        if ( has_nested ) {
          const Kokkos::Future<long,ExecSpace> fib_1 = policy.get_dependence(this,0);
          const Kokkos::Future<long,ExecSpace> fib_2 = policy.get_dependence(this,1);

          result = fib_1.get() + fib_2.get();
        }
        else {
          // Spawn new children and respawn myself to sum their results:
          has_nested = 1 ;

          Kokkos::respawn( policy
                         , this
                         , Kokkos::spawn( policy , FibChild(policy,n-1) )
                         , Kokkos::spawn( policy , FibChild(policy,n-2) )
                         );
        }
      }
    }
};

template< class ExecSpace >
struct FibChild2 {

  typedef long value_type ;

  Kokkos::TaskPolicy<ExecSpace> policy ;
  const value_type n ;
  int has_nested ;

  FibChild2( const Kokkos::TaskPolicy<ExecSpace> & arg_policy
           , const value_type arg_n )
    : policy(arg_policy,2) /* default dependence capacity = 2 */
    , n( arg_n ), has_nested(0) {}

  inline
  void apply( value_type & result )
    {
      if ( ! has_nested ) {
        if ( n < 2 ) {
          result = n ;
        }
        else if ( n < 4 ) {
          // Spawn new children and respawn myself to sum their results:
          // result = Fib(n-1) + Fib(n-2)
          has_nested = 2 ;
          // Kokkos::respawn implements the following steps:
          policy.clear_dependence( this );
          policy.add_dependence( this , Kokkos::spawn( policy , FibChild2(policy,n-1) ) );
          policy.add_dependence( this , Kokkos::spawn( policy , FibChild2(policy,n-2) ) );
          policy.respawn( this );
        }
        else {
          // Spawn new children and respawn myself to sum their results:
          // result = Fib(n-1) + Fib(n-2)
          // result = ( Fib(n-2) + Fib(n-3) ) + ( Fib(n-3) + Fib(n-4) )
          // result = ( ( Fib(n-3) + Fib(n-4) ) + Fib(n-3) ) + ( Fib(n-3) + Fib(n-4) )
          // result = 3 * Fib(n-3) + 2 * Fib(n-4)
          has_nested = 4 ;
          // Kokkos::respawn implements the following steps:
          policy.clear_dependence( this );
          policy.add_dependence( this , Kokkos::spawn( policy , FibChild2(policy,n-3) ) );
          policy.add_dependence( this , Kokkos::spawn( policy , FibChild2(policy,n-4) ) );
          policy.respawn( this );
        }
     }
     else {
        const Kokkos::Future<long,ExecSpace> fib_a = policy.get_dependence(this,0);
        const Kokkos::Future<long,ExecSpace> fib_b = policy.get_dependence(this,1);

        result = ( has_nested == 2 ) ? fib_a.get() + fib_b.get()
                                     : 3 * fib_a.get() + 2 * fib_b.get() ;
      }
    }
};

namespace {

long eval_fib( long n )
{
  return n < 2 ? n : eval_fib(n-2) + eval_fib(n-1);
}

}

template< class ExecSpace >
void test_fib( long n )
{
  Kokkos::TaskPolicy<ExecSpace> policy(2);

  Kokkos::Future<long,ExecSpace> f = Kokkos::spawn( policy , FibChild<ExecSpace>(policy,n) );

  Kokkos::wait( f );

  if ( f.get() != eval_fib(n) ) {
    std::cout << "Fib(" << n << ") = " << f.get();
    std::cout << " != " << eval_fib(n);
    std::cout << std::endl ;
  }
}

template< class ExecSpace >
void test_fib2( long n )
{
  Kokkos::TaskPolicy<ExecSpace> policy(2); // default dependence capacity

  Kokkos::Future<long,ExecSpace> f = Kokkos::spawn( policy , FibChild2<ExecSpace>(policy,n) );

  Kokkos::wait( f );

  if ( f.get() != eval_fib(n) ) {
    std::cout << "Fib2(" << n << ") = " << f.get();
    std::cout << " != " << eval_fib(n);
    std::cout << std::endl ;
  }
}

//----------------------------------------------------------------------------

template< class ExecSpace >
struct Norm2 {

  typedef double value_type ;

  const double * const m_x ;

  Norm2( const double * x ) : m_x(x) {}

  inline
  void init( double & val ) const { val = 0 ; }

  inline
  void operator()( int i , double & val ) const { val += m_x[i] * m_x[i] ; }

  void apply( double & dst ) const { dst = std::sqrt( dst ); }
};

template< class ExecSpace >
void test_norm2( const int n )
{
  Kokkos::TaskPolicy< ExecSpace > policy ;

  double * const x = new double[n];

  for ( int i = 0 ; i < n ; ++i ) x[i] = 1 ;

  Kokkos::RangePolicy<ExecSpace> r(0,n);

  Kokkos::Future<double,ExecSpace> f = Kokkos::spawn_reduce( policy , r , Norm2<ExecSpace>(x) );

  Kokkos::wait( f );

#if defined(PRINT)
  std::cout << "Norm2: " << f.get() << std::endl ;
#endif

  delete[] x ;
}

//----------------------------------------------------------------------------

template< class Space >
struct TaskDep {

  typedef int value_type ;
  typedef Kokkos::TaskPolicy< Space > policy_type ;

  const policy_type policy ;
  const int         input ;

  TaskDep( const policy_type & arg_p , const int arg_i )
    : policy( arg_p ), input( arg_i ) {}

  void apply( int & val )
  {
    val = input ;
    const int num = policy.get_dependence( this );

    for ( int i = 0 ; i < num ; ++i ) {
      Kokkos::Future<int,Space> f = policy.get_dependence( this , i );
      val += f.get();
    }
  }
};


template< class Space >
void test_task_dep( const int n )
{
  enum { NTEST = 64 };

  Kokkos::TaskPolicy< Space > policy ;

  Kokkos::Future<int,Space> f[ NTEST ];

  for ( int i = 0 ; i < NTEST ; ++i ) {
    // Create task in the "constructing" state with capacity for 'n+1' dependences
    f[i] = policy.create( TaskDep<Space>(policy,0) , n + 1 );

    if ( f[i].get_task_state() != Kokkos::TASK_STATE_CONSTRUCTING ) {
      Kokkos::Impl::throw_runtime_exception("get_task_state() != Kokkos::TASK_STATE_CONSTRUCTING");
    }

    // Only use 'n' dependences

    for ( int j = 0 ; j < n ; ++j ) {
      Kokkos::Future<int,Space> nested = policy.create( TaskDep<Space>(policy,j+1) );
      policy.spawn( nested );
      // Add dependence to a "constructing" task
      policy.add_dependence( f[i] , nested );
    }

    // Spawn task from the "constructing" to the "waiting" state
    policy.spawn( f[i] );
  }

  const int answer = n % 1 ? n * ( ( n + 1 ) / 2 ) : ( n / 2 ) * ( n + 1 );

  int error = 0 ;
  for ( int i = 0 ; i < NTEST ; ++i ) {
    Kokkos::wait( f[i] );
    if ( f[i].get_task_state() != Kokkos::TASK_STATE_COMPLETE ) {
      Kokkos::Impl::throw_runtime_exception("get_task_state() != Kokkos::TASK_STATE_COMPLETE");
    }
    if ( answer != f[i].get() && 0 == error ) {
      std::cout << "test_task_dep(" << n << ") ERROR at[" << i << "]"
                << " answer(" << answer << ") != result(" << f[i].get() << ")" << std::endl ;
    }
  }
}

//----------------------------------------------------------------------------

template< class ExecSpace >
struct FibArray {

  typedef long value_type ;

  Kokkos::TaskPolicy<ExecSpace> policy ;
  const value_type n ;

  FibArray( const Kokkos::TaskPolicy<ExecSpace> & arg_policy
          , const value_type arg_n )
    : policy(arg_policy) , n( arg_n ) {}

  inline
  void apply( value_type & result )
    {
      if ( n < 2 ) {
        result = n ;
      }
      else {
        const Kokkos::Future<long,ExecSpace> fib_1 = policy.get_dependence(this,0);
        const Kokkos::Future<long,ExecSpace> fib_2 = policy.get_dependence(this,1);

        result = fib_1.get() + fib_2.get();
      }
    }
};

template< class Space >
void test_future_array( const int n )
{
  Kokkos::TaskPolicy< Space > policy ;

  Kokkos::FutureArray< long , Space > fa( n+1 < 2 ? 2 : n+1 );

  fa[0] = Kokkos::spawn( policy , FibArray<Space>( policy , 0 ) );
  fa[1] = Kokkos::spawn( policy , FibArray<Space>( policy , 1 ) );
  
  for ( int i = 2 ; i <= n ; ++i ) {
    fa[i] = Kokkos::spawn( policy
                         , FibArray<Space>( policy , i )
                         , fa[i-1] , fa[i-2] );
  }

  Kokkos::wait( fa[n] );

  for ( int i = 0 ; i <= n ; ++i ) {
    if ( fa[i].get() != eval_fib(i) ) {
      std::cout << "test_future_array Fib(" << i << ") = " << fa[i].get();
      std::cout << " != " << eval_fib(i);
      std::cout << std::endl ;
    }
  }
}

//----------------------------------------------------------------------------

} // namespace TestTaskPolicy

#endif /* #ifndef KOKKOS_UNITTEST_TASKPOLICY_HPP */


