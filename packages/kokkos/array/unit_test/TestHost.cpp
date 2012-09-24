/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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

#include <gtest/gtest.h>

#include <KokkosArray_View.hpp>

#include <KokkosArray_CrsArray.hpp>

#include <KokkosArray_Host.hpp>

//----------------------------------------------------------------------------

#include <TestViewImpl.hpp>

#include <KokkosArray_Host_macros.hpp>
#include <TestMemoryTracking.hpp>
#include <TestViewAPI.hpp>

#include <TestCrsArray.hpp>
#include <TestReduce.hpp>
#include <TestMultiReduce.hpp>

#include <KokkosArray_Clear_macros.hpp>

namespace Test {

class host : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    const size_t node_count = KokkosArray::Host::detect_node_count();
/*
    const size_t node_core_count = KokkosArray::Host::detect_node_core_count();
    std::cout << "KokkosArray::Host node_count(" << node_count
              << ") X node_core_count(" << node_core_count
              << ")" << std::endl ;
*/
    KokkosArray::Host::initialize( node_count , 4 );
  }

  static void TearDownTestCase()
  {
    KokkosArray::Host::finalize();
  }
};

TEST_F( host, memory_tracking) {
  TestMemoryTracking();
}

TEST_F( host, view_impl) {
  test_view_impl< KokkosArray::Host >();
}

TEST_F( host, view_api) {
  TestViewAPI< double , KokkosArray::Host >();
}


TEST_F( host, crsarray) {
  TestCrsArray< KokkosArray::Host >();
}

TEST_F( host, long_reduce) {
  TestReduce< long ,   KokkosArray::Host >( 1000000 );
}

TEST_F( host, double_reduce) {
  TestReduce< double ,   KokkosArray::Host >( 1000000 );
}

TEST_F( host, long_multi_reduce) {
  TestReduceMulti< long , KokkosArray::Host >( 1000000 , 7 );
}

TEST_F( host , view_remap )
{
  enum { N0 = 3 , N1 = 2 , N2 = 8 , N3 = 9 };

  typedef KokkosArray::View< double*[N1][N2][N3] ,
                             KokkosArray::LayoutRight ,
                             KokkosArray::Host > output_type ;

  typedef KokkosArray::View< int**[N2][N3] ,
                             KokkosArray::LayoutLeft ,
                             KokkosArray::Host > input_type ;

  typedef KokkosArray::View< int*[N0][N2][N3] ,
                             KokkosArray::LayoutLeft ,
                             KokkosArray::Host > diff_type ;

  output_type output( "output" , N0 );
  input_type  input ( "input" , N0 , N1 );
  diff_type   diff  ( "diff" , N0 );

  int value = 0 ;
  for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
  for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
  for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
  for ( size_t i0 = 0 ; i0 < N0 ; ++i0 ) {
    input(i0,i1,i2,i3) = ++value ;
  }}}}

  // KokkosArray::deep_copy( diff , input ); // throw with incompatible shape
  KokkosArray::deep_copy( output , input );

  value = 0 ;
  for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
  for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
  for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
  for ( size_t i0 = 0 ; i0 < N0 ; ++i0 ) {
    ++value ;
    ASSERT_EQ( value , ((int) output(i0,i1,i2,i3) ) );
  }}}}
}

//----------------------------------------------------------------------------

struct HostFunctor {

  typedef int value_type ;

  volatile int & flag ;

  static void join( volatile int & update , const volatile int & input )
    { update += input ; }

  void operator()( const value_type & value ) const
    { flag += value + 1 ; }

  HostFunctor( int & f ) : flag(f) {}

  void operator()( KokkosArray::Impl::HostThread & thread ) const
    {
      int value = 0 ;
      thread.barrier();
      thread.barrier();
      thread.reduce< HostFunctor >( value , *this );
      thread.reduce< HostFunctor >( value , *this );
      thread.barrier();
    }
};

TEST_F( host , host_thread )
{
  const int N = 1000 ;
  int flag = 0 ;

  for ( int i = 0 ; i < 1000 ; ++i ) {
    KokkosArray::Impl::HostParallelLaunch< HostFunctor >( HostFunctor( flag ) );
  }

  ASSERT_EQ( flag , N * 2 );
}

//----------------------------------------------------------------------------

} // namespace test

