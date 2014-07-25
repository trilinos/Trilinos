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


#ifndef KOKKOS_TASK_UNIT_TEST_HPP
#define KOKKOS_TASK_UNIT_TEST_HPP

#include <cmath>
#include <iostream>
#include <Kokkos_Task.hpp>

namespace Test {

//----------------------------------------------------------------------------

template< class ExecSpace >
struct FibChild : public Kokkos::task_serial<ExecSpace,long> {

#if 0
  Kokkos::TaskExecPolicy< ExecSpace > policy ;
#endif

  const long n ;
  int has_nested ;

#if 0
  FibChild( const Kokkos::TaskExecPolicy<ExecSpace> & arg_policy , const long arg_n )
    : policy(arg_policy), n( arg_n ), has_nested(0) {}
#else
  FibChild( const long arg_n ) : n( arg_n ), has_nested(0) {}
#endif

  inline
  void apply( long & result )
    {
      if ( n < 2 ) {
        result = n ;
      }
      else {
        if ( has_nested ) {
          const Kokkos::Future<long,ExecSpace> fib_1 = Kokkos::task_dependence(this,0);
          const Kokkos::Future<long,ExecSpace> fib_2 = Kokkos::task_dependence(this,1);

          result = fib_1.get() + fib_2.get();
        }
        else {
          Kokkos::Future<long,ExecSpace> nested[2] ;
          // Spawn new children and respawn myself to sum their results:
          has_nested = 1 ;
#if 1
          nested[0] = Kokkos::spawn( FibChild(n-1) );
          nested[1] = Kokkos::spawn( FibChild(n-2) );
          Kokkos::respawn( this , nested , nested + 2 );
#else
          nested[0] = Kokkos::spawn( policy , FibChild(n-1) );
          nested[1] = Kokkos::spawn( policy , FibChild(n-2) );
          Kokkos::respawn( policy.depend(nested,2) , *this );
#endif
        }
      }
    }
};

template< class ExecSpace >
struct FibChild2 : public Kokkos::task_serial<ExecSpace,long> {

  const long n ;
  int nested ;

  FibChild2( const long arg_n ) : n( arg_n ), nested(0) {}

  inline
  void apply( long & result )
    {
      if ( ! nested ) {
        Kokkos::Future<long,ExecSpace> nest[2] ;
        if ( n < 2 ) {
          result = n ;
        }
        else if ( n < 4 ) {
          // Spawn new children and respawn myself to sum their results:
          // result = Fib(n-1) + Fib(n-2)
          nested = 2 ;
          nest[0] = Kokkos::spawn( FibChild2(n-1) );
          nest[1] = Kokkos::spawn( FibChild2(n-2) );
          Kokkos::respawn( this , nest , nest + 2 );
        }
        else {
          // Spawn new children and respawn myself to sum their results:
          // result = Fib(n-1) + Fib(n-2)
          // result = ( Fib(n-2) + Fib(n-3) ) + ( Fib(n-3) + Fib(n-4) )
          // result = ( ( Fib(n-3) + Fib(n-4) ) + Fib(n-3) ) + ( Fib(n-3) + Fib(n-4) )
          // result = 3 * Fib(n-3) + 2 * Fib(n-4)
          nested = 4 ;
          nest[0] = Kokkos::spawn( FibChild2(n-3) );
          nest[1] = Kokkos::spawn( FibChild2(n-4) );
          Kokkos::respawn( this , nest , nest + 2 );
        }
     }
     else {
        const Kokkos::Future<long,ExecSpace> fib_a = Kokkos::task_dependence(this,0);
        const Kokkos::Future<long,ExecSpace> fib_b = Kokkos::task_dependence(this,1);

        result = ( nested == 2 ) ? fib_a.get() + fib_b.get()
                                 : 3 * fib_a.get() + 2 * fib_b.get() ;
      }
    }
};

long eval_fib( long n )
{
  return n < 2 ? n : eval_fib(n-2) + eval_fib(n-1);
}

template< class ExecSpace >
void test_fib( long n )
{
#if 0
  Kokkos::TaskExecPolicy< ExecSpace > policy ;

  Kokkos::Future<long> f = Kokkos::spawn( policy , FibChild(policy,n) );
#endif

  Kokkos::Future<long,ExecSpace> f = Kokkos::spawn( FibChild<ExecSpace>(n) );

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
  Kokkos::Future<long,ExecSpace> f = Kokkos::spawn( FibChild2<ExecSpace>(n) );

  Kokkos::wait( f );

  if ( f.get() != eval_fib(n) ) {
    std::cout << "Fib2(" << n << ") = " << f.get();
    std::cout << " != " << eval_fib(n);
    std::cout << std::endl ;
  }
}

//----------------------------------------------------------------------------

template< class ExecSpace >
struct Norm2 : public Kokkos::task_reduce<ExecSpace,double,size_t> {

  const double * const m_x ;

  Norm2( const int n , const double * x ) : Kokkos::task_reduce<ExecSpace,double,size_t>(n), m_x(x) {}

  void operator()( int i , double & val ) const { val += m_x[i] * m_x[i] ; }

  void final( double & dst ) const { dst = std::sqrt( dst ); }
};

template< class ExecSpace >
void test_norm2( const int n )
{
  double * const x = new double[n];

  for ( int i = 0 ; i < n ; ++i ) x[i] = 1 ;

  Kokkos::Future<double,ExecSpace> f = Kokkos::spawn( Norm2<ExecSpace>(n,x) );

#if 0

  Kokkos::Future<double,ExecSpace> f =
    Kokkos::spawn< Kokkos::task_reduce<double,ExecSpace,size_t> >
      ( n
      , [=]( int i , double & val ) { val += x[i] * x[i] ; } /* data parallel operator */
      , []( double & val ) { val = std::sqrt( val ); }       /* task serial operator */
      );

  Kokkos::parallel_reduce<ExecSpace,double,int>( n , Norm2(x) , result );
    // Compose task-functor from input functor type ...

#endif

  Kokkos::wait( f );

#if defined(PRINT)
  std::cout << "Norm2: " << f.get() << std::endl ;
#endif

  delete[] x ;
}

//----------------------------------------------------------------------------

} // namespace Test

#endif /* #ifndef KOKKOS_TASK_UNIT_TEST_HPP */

