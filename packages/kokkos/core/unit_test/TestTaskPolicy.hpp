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

  Kokkos::TaskPolicy< ExecSpace > policy ;

  const value_type n ;
  int has_nested ;

  FibChild( const Kokkos::TaskPolicy<ExecSpace> & arg_policy
          , const value_type arg_n )
    : policy(arg_policy), n( arg_n ), has_nested(0) {}

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
          Kokkos::Future<long,ExecSpace> nested[2] ;
          // Spawn new children and respawn myself to sum their results:
          has_nested = 1 ;
          nested[0] = policy.spawn( FibChild(policy,n-1) );
          nested[1] = policy.spawn( FibChild(policy,n-2) );
          policy.respawn( this , nested , 2 );
        }
      }
    }
};

template< class ExecSpace >
struct FibChild2 {

  typedef long value_type ;

  Kokkos::TaskPolicy< ExecSpace > policy ;

  const value_type n ;
  int has_nested ;

  FibChild2( const Kokkos::TaskPolicy<ExecSpace> & arg_policy
           , const value_type arg_n )
    : policy(arg_policy), n( arg_n ), has_nested(0) {}

  inline
  void apply( value_type & result )
    {
      if ( ! has_nested ) {
        Kokkos::Future<long,ExecSpace> nest[2] ;
        if ( n < 2 ) {
          result = n ;
        }
        else if ( n < 4 ) {
          // Spawn new children and respawn myself to sum their results:
          // result = Fib(n-1) + Fib(n-2)
          has_nested = 2 ;
          nest[0] = policy.spawn( FibChild2(policy,n-1) );
          nest[1] = policy.spawn( FibChild2(policy,n-2) );
          policy.respawn( this , nest , 2 );
        }
        else {
          // Spawn new children and respawn myself to sum their results:
          // result = Fib(n-1) + Fib(n-2)
          // result = ( Fib(n-2) + Fib(n-3) ) + ( Fib(n-3) + Fib(n-4) )
          // result = ( ( Fib(n-3) + Fib(n-4) ) + Fib(n-3) ) + ( Fib(n-3) + Fib(n-4) )
          // result = 3 * Fib(n-3) + 2 * Fib(n-4)
          has_nested = 4 ;
          nest[0] = policy.spawn( FibChild2(policy,n-3) );
          nest[1] = policy.spawn( FibChild2(policy,n-4) );
          policy.respawn( this , nest , 2 );
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
  Kokkos::TaskPolicy< ExecSpace > policy ;

  Kokkos::Future<long,ExecSpace> f = Kokkos::spawn( policy , FibChild<ExecSpace>(policy,n) );

  Kokkos::wait( policy , f );

  if ( f.get() != eval_fib(n) ) {
    std::cout << "Fib(" << n << ") = " << f.get();
    std::cout << " != " << eval_fib(n);
    std::cout << std::endl ;
  }
}

template< class ExecSpace >
void test_fib2( long n )
{
  Kokkos::TaskPolicy< ExecSpace > policy ;

  Kokkos::Future<long,ExecSpace> f = Kokkos::spawn( policy , FibChild2<ExecSpace>(policy,n) );

  Kokkos::wait( policy , f );

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

  Kokkos::Future<double,ExecSpace> f = Kokkos::spawn( policy.reduce(r) , Norm2<ExecSpace>(x) );

  Kokkos::wait( policy , f );

#if defined(PRINT)
  std::cout << "Norm2: " << f.get() << std::endl ;
#endif

  delete[] x ;
}

//----------------------------------------------------------------------------

} // namespace TestTaskPolicy

#endif /* #ifndef KOKKOS_UNITTEST_TASKPOLICY_HPP */


