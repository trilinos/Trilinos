/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <gtest/gtest.h>

#include <impl/KokkosArray_MemoryView.hpp>
#include <KokkosArray_Value.hpp>
#include <KokkosArray_MultiVector.hpp>
#include <KokkosArray_Array.hpp>
#include <KokkosArray_MDArray.hpp>
#include <KokkosArray_PrefixSum.hpp>
#include <KokkosArray_CrsArray.hpp>

#include <KokkosArray_View.hpp>

#include <KokkosArray_Host.hpp>

//----------------------------------------------------------------------------

#include <KokkosArray_Host_macros.hpp>
#include <impl/KokkosArray_IndexMapLeft_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

//----------------------------------------------------------------------------

#include <TestViewImpl.hpp>

#include <KokkosArray_Host_macros.hpp>
#include <TestMemoryTracking.hpp>
#include <TestViewAPI.hpp>

#include <TestMemoryManagement.hpp>
#include <TestValue.hpp>
#include <TestMultiVector.hpp>
#include <TestArray.hpp>
#include <TestMDArray.hpp>
#include <TestMDArrayIndexMap.hpp>
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
    KokkosArray::Host::initialize( node_count , 2 );
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


TEST_F( host, memory_management_double) {
  TestMemoryManagement< double, KokkosArray::Host >();
}

TEST_F( host, memory_management_int) {
  TestMemoryManagement< int, KokkosArray::Host >();
}

TEST_F( host, value_view_double) {
  TestValue< double, KokkosArray::Host >();
}

TEST_F( host, value_view_int) {
  TestValue< int, KokkosArray::Host >();
}

TEST_F( host, multivector_double) {
  TestMultiVector< double, KokkosArray::Host >();
}

TEST_F( host, multivector_int) {
  TestMultiVector< int, KokkosArray::Host >();
}

TEST_F( host, crsarray) {
  TestCrsArray< KokkosArray::Host >();
}

TEST_F( host, mdarray_view_double) {
  TestMDArray< double, KokkosArray::Host >();
}

TEST_F( host, mdarray_view_int) {
  TestMDArray< int, KokkosArray::Host >();
}

TEST_F( host, mdarray_index_map) {
  TestMDArrayIndexMap< KokkosArray::Host >();
}

TEST_F( host, array_view_double) {
  TestArray< double, KokkosArray::Host >();
}

TEST_F( host, array_view_int) {
  TestArray< int, KokkosArray::Host >();
}

// TEST_F( host, array_index_map) {
//   TestArrayIndexMap< KokkosArray::Host >();
// }

TEST_F( host, long_reduce) {
  TestReduce< long ,   KokkosArray::Host >( 1000000 );
}

TEST_F( host, double_reduce) {
  TestReduce< double ,   KokkosArray::Host >( 1000000 );
}

TEST_F( host, long_multi_reduce) {
  TestReduceMulti< long , KokkosArray::Host >( 1000000 , 7 );
}


} // namespace test
