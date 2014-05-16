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

template< class Device >
struct FibChild : public Kokkos::task_serial<Device,long> {

#if 0
  Kokkos::TaskPool< FibChild > fib_task_pool ;
#endif

  const long n ;
  int has_nested ;

  FibChild( const long arg_n ) : n( arg_n ), has_nested(0) {}

  inline
  void apply( long & result )
    {
      if ( n < 2 ) {
        result = n ;
      }
      else {
        if ( has_nested ) {
#if 1
          const Kokkos::Future<Device,long> fib_1 = Kokkos::task_dependence(this,0);
          const Kokkos::Future<Device,long> fib_2 = Kokkos::task_dependence(this,1);
#else
          const Kokkos::Future<Device,long> fib_1 = fib_task_pool.dependence(this,0);
          const Kokkos::Future<Device,long> fib_2 = fib_task_pool.dependence(this,1);
#endif

          result = fib_1.get() + fib_2.get();
        }
        else {
          Kokkos::Future<Device,long> nested[2] ;
          // Spawn new children and respawn myself to sum their results:
          has_nested = 1 ;
#if 1
          nested[0] = Kokkos::spawn( FibChild(n-1) );
          nested[1] = Kokkos::spawn( FibChild(n-2) );
          Kokkos::respawn( this , nested , nested + 2 );
#else
          nested[0] = fib_task_pool.spawn( n-1 );
          nested[1] = fib_task_pool.spawn( n-2 );
          fib_task_pool.respawn( this , nested , nested + 2 );
#endif
        }
      }
    }
};

template< class Device >
struct FibChild2 : public Kokkos::task_serial<Device,long> {

  const long n ;
  int nested ;

  FibChild2( const long arg_n ) : n( arg_n ), nested(0) {}

  inline
  void apply( long & result )
    {
      if ( ! nested ) {
        Kokkos::Future<Device,long> nest[2] ;
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
        const Kokkos::Future<Device,long> fib_a = Kokkos::task_dependence(this,0);
        const Kokkos::Future<Device,long> fib_b = Kokkos::task_dependence(this,1);

        result = ( nested == 2 ) ? fib_a.get() + fib_b.get()
                                 : 3 * fib_a.get() + 2 * fib_b.get() ;
      }
    }
};

long eval_fib( long n )
{
  return n < 2 ? n : eval_fib(n-2) + eval_fib(n-1);
}

template< class Device >
void test_fib( long n )
{
  Kokkos::Future<Device,long> f = Kokkos::spawn( FibChild<Device>(n) );

  Kokkos::wait( f );

  if ( f.get() != eval_fib(n) ) {
    std::cout << "Fib(" << n << ") = " << f.get();
    std::cout << " != " << eval_fib(n);
    std::cout << std::endl ;
  }
}

template< class Device >
void test_fib2( long n )
{
  Kokkos::Future<Device,long> f = Kokkos::spawn( FibChild2<Device>(n) );

  Kokkos::wait( f );

  if ( f.get() != eval_fib(n) ) {
    std::cout << "Fib2(" << n << ") = " << f.get();
    std::cout << " != " << eval_fib(n);
    std::cout << std::endl ;
  }
}

//----------------------------------------------------------------------------

template< class Device >
struct Norm2 : public Kokkos::task_reduce<Device,double,size_t> {

  const double * const m_x ;

  Norm2( const int n , const double * x ) : Kokkos::task_reduce<Device,double,size_t>(n), m_x(x) {}

  void operator()( int i , double & val ) const { val += m_x[i] * m_x[i] ; }

  void final( double & dst ) const { dst = std::sqrt( dst ); }
};

template< class Device >
void test_norm2( const int n )
{
  double * const x = new double[n];

  for ( int i = 0 ; i < n ; ++i ) x[i] = 1 ;

  Kokkos::Future<Device,double> f = Kokkos::spawn( Norm2<Device>(n,x) );

#if 0

  Kokkos::Future<Device,double> f =
    Kokkos::spawn< Kokkos::task_reduce<Device,double,size_t> >
      ( n
      , [=]( int i , double & val ) { val += x[i] * x[i] ; } /* data parallel operator */
      , []( double & val ) { val = std::sqrt( val ); }       /* task serial operator */
      );

  Kokkos::parallel_reduce<Device,double,int>( n , Norm2(x) , result );
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

