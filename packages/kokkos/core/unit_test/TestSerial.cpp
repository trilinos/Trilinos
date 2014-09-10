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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_hwloc.hpp>

#include <Kokkos_View.hpp>
#include <impl/Kokkos_ViewTileLeft.hpp>
#include <impl/Kokkos_Serial_TaskPolicy.hpp>

#include <Kokkos_CrsArray.hpp>

//----------------------------------------------------------------------------

#include <TestViewImpl.hpp>

#include <TestMemoryTracking.hpp>
#include <TestViewAPI.hpp>
#include <TestAtomic.hpp>
#include <TestTile.hpp>
#include <TestRange.hpp>
#include <TestTeam.hpp>
#include <TestCrsArray.hpp>
#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestAggregate.hpp>
#include <TestCompilerMacros.hpp>
#include <TestTaskPolicy.hpp>
#include <TestCXX11.hpp>
#include <TestTeamVector.hpp>

namespace Test {

class serial : public ::testing::Test {
protected:
  static void SetUpTestCase() {}
  static void TearDownTestCase() {}
};

TEST_F( serial, memory_tracking) {
  TestMemoryTracking();
}

TEST_F( serial, view_impl) {
  test_view_impl< Kokkos::Serial >();
}

TEST_F( serial, view_api) {
  TestViewAPI< double , Kokkos::Serial >();
}

TEST_F( serial , range_tag )
{
  TestRange< Kokkos::Serial >::test_for(1000);
  TestRange< Kokkos::Serial >::test_reduce(1000);
  TestRange< Kokkos::Serial >::test_scan(1000);
}

TEST_F( serial, crsarray) {
  TestCrsArray< Kokkos::Serial >();
}

TEST_F( serial, long_reduce) {
  TestReduce< long ,   Kokkos::Serial >( 1000000 );
}

TEST_F( serial, double_reduce) {
  TestReduce< double ,   Kokkos::Serial >( 1000000 );
}

TEST_F( serial, long_reduce_dynamic ) {
  TestReduceDynamic< long ,   Kokkos::Serial >( 1000000 );
}

TEST_F( serial, double_reduce_dynamic ) {
  TestReduceDynamic< double ,   Kokkos::Serial >( 1000000 );
}

TEST_F( serial, long_reduce_dynamic_view ) {
  TestReduceDynamicView< long ,   Kokkos::Serial >( 1000000 );
}

TEST_F( serial , scan )
{
  TestScan< Kokkos::Serial >::test_range( 1 , 1000 );
  TestScan< Kokkos::Serial >( 10 );
  TestScan< Kokkos::Serial >( 10000 );
}

TEST_F( serial , team_long_reduce) {
  TestReduceTeam< long ,   Kokkos::Serial >( 100000 );
}

TEST_F( serial , team_double_reduce) {
  TestReduceTeam< double ,   Kokkos::Serial >( 100000 );
}

TEST_F( serial , team_shared_request) {
  TestSharedTeam< Kokkos::Serial >();
}

TEST_F( serial  , team_scan )
{
  TestScanTeam< Kokkos::Serial >( 10 );
  TestScanTeam< Kokkos::Serial >( 10000 );
}


TEST_F( serial , view_remap )
{
  enum { N0 = 3 , N1 = 2 , N2 = 8 , N3 = 9 };

  typedef Kokkos::View< double*[N1][N2][N3] ,
                             Kokkos::LayoutRight ,
                             Kokkos::Serial > output_type ;

  typedef Kokkos::View< int**[N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Serial > input_type ;

  typedef Kokkos::View< int*[N0][N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Serial > diff_type ;

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

  // Kokkos::deep_copy( diff , input ); // throw with incompatible shape
  Kokkos::deep_copy( output , input );

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

TEST_F( serial , view_aggregate )
{
  TestViewAggregate< Kokkos::Serial >();
}

//----------------------------------------------------------------------------

TEST_F( serial , atomics )
{
  const int loop_count = 1e6 ;

  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Serial>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Serial>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Serial>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Serial>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Serial>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Serial>(100,3) ) );
}

//----------------------------------------------------------------------------

TEST_F( serial, tile_1x1)
{
  static const size_t dim = 9;
  typedef Kokkos::LayoutTileLeft<1,1> tile_layout;
  typedef ReduceTileErrors< Kokkos::Serial, tile_layout > functor_type;

  functor_type::array_type array("",dim,dim);
  ptrdiff_t errors = 0 ;
  Kokkos::parallel_reduce(dim, functor_type(array) , errors );
  EXPECT_EQ( errors, 0u);
}

TEST_F( serial, tile_2x2)
{
  static const size_t dim = 9;
  typedef Kokkos::LayoutTileLeft<2,2> tile_layout;
  typedef ReduceTileErrors< Kokkos::Serial, tile_layout > functor_type;

  functor_type::array_type array("",dim,dim);
  ptrdiff_t errors = 0 ;
  Kokkos::parallel_reduce(dim, functor_type(array) , errors );
  EXPECT_EQ( errors, ptrdiff_t(0) );
}

TEST_F( serial, tile_4x4)
{
  static const size_t dim = 9;
  typedef Kokkos::LayoutTileLeft<4,4> tile_layout;
  typedef ReduceTileErrors< Kokkos::Serial, tile_layout > functor_type;

  functor_type::array_type array("",dim,dim);
  ptrdiff_t errors = 0 ;
  Kokkos::parallel_reduce(dim, functor_type(array) , errors );
  EXPECT_EQ( errors, ptrdiff_t(0) );
}

TEST_F( serial, tile_8x8)
{
  static const size_t dim = 9;
  typedef Kokkos::LayoutTileLeft<8,8> tile_layout;
  typedef ReduceTileErrors< Kokkos::Serial, tile_layout > functor_type;

  functor_type::array_type array("",dim,dim);
  ptrdiff_t errors = 0 ;
  Kokkos::parallel_reduce(dim, functor_type(array) , errors );
  EXPECT_EQ( errors, ptrdiff_t(0) );
}

TEST_F( serial, tile_16x16)
{
  static const size_t dim = 9;
  typedef Kokkos::LayoutTileLeft<16,16> tile_layout;
  typedef ReduceTileErrors< Kokkos::Serial, tile_layout > functor_type;

  functor_type::array_type array("",dim,dim);
  ptrdiff_t errors = 0 ;
  Kokkos::parallel_reduce(dim, functor_type(array) , errors );
  EXPECT_EQ( errors, ptrdiff_t(0) );
}

//----------------------------------------------------------------------------

TEST_F( serial , compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::Serial >() ) );
}

//----------------------------------------------------------------------------

TEST_F( serial , task_policy )
{
  TestTaskPolicy::test_norm2< Kokkos::Serial >( 1000 );
  for ( long i = 0 ; i < 30 ; ++i ) TestTaskPolicy::test_fib< Kokkos::Serial >(i);
  for ( long i = 0 ; i < 40 ; ++i ) TestTaskPolicy::test_fib2< Kokkos::Serial >(i);
}

//----------------------------------------------------------------------------
#if defined( KOKKOS_HAVE_CXX11 ) && defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_SERIAL )
TEST_F( serial , cxx11 )
{
  if ( Kokkos::Impl::is_same< Kokkos::DefaultExecutionSpace , Kokkos::Serial >::value ) {
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >(1) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >(2) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >(3) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Serial >(4) ) );
  }
}
#endif

#if defined (KOKKOS_HAVE_CXX11)
TEST_F( serial , team_vector )
{
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(1) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(2) ) );
  ASSERT_TRUE( ( TestTeamVector::Test< Kokkos::Serial >(3) ) );
}
#endif
} // namespace test

