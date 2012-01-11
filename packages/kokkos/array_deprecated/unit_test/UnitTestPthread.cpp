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

#include <Kokkos_DevicePthread.hpp>
#include <Kokkos_DevicePthread_ValueView.hpp>
#include <Kokkos_DevicePthread_MultiVectorView.hpp>
#include <Kokkos_DevicePthread_MDArrayView.hpp>
#include <Kokkos_DevicePthread_ParallelFor.hpp>
#include <Kokkos_DevicePthread_ParallelReduce.hpp>

//----------------------------------------------------------------------------

#include <Kokkos_DevicePthread_macros.hpp>
#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
#include <UnitTestMDArrayIndexMap.hpp>
#include <UnitTestReduce.hpp>
#include <UnitTestMultiReduce.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Test {

class pthread : public ::testing::Test {
  protected:
    static void SetUpTestCase() {
      Kokkos::DevicePthread::initialize( 4 );
    }
    static void TearDownTestCase() {
      Kokkos::DevicePthread::finalize();
    }
};



TEST_F( pthread, memory_management_double) {
  UnitTestDeviceMemoryManagement< double, Kokkos::DevicePthread >();
}

TEST_F( pthread, memory_management_int) {
  UnitTestDeviceMemoryManagement< int, Kokkos::DevicePthread >();
}

TEST_F( pthread, value_view_double) {
  UnitTestValueView< double, Kokkos::DevicePthread >();
}

TEST_F( pthread, value_view_int) {
  UnitTestValueView< int, Kokkos::DevicePthread >();
}

TEST_F( pthread, multi_vector_view_double) {
  UnitTestMultiVectorView< double, Kokkos::DevicePthread >();
}

TEST_F( pthread, multi_vector_view_int) {
  UnitTestMultiVectorView< int, Kokkos::DevicePthread >();
}

TEST_F( pthread, mdarray_view_double) {
  UnitTestMDArrayView< double, Kokkos::DevicePthread >();
}

TEST_F( pthread, mdarray_view_int) {
  UnitTestMDArrayView< int, Kokkos::DevicePthread >();
}

TEST_F( pthread, mdarray_deep_copy_double) {
  UnitTestMDArrayDeepCopy< double, Kokkos::DevicePthread >();
}

TEST_F( pthread, mdarray_deep_copy_int) {
  UnitTestMDArrayDeepCopy< int, Kokkos::DevicePthread >();
}

TEST_F( pthread, mdarray_index_map) {
  UnitTestMDArrayIndexMap< Kokkos::DevicePthread >();
}

TEST_F( pthread, long_reduce) {
  UnitTestReduce< long ,   Kokkos::DevicePthread >( 1000000 );
}

TEST_F( pthread, double_reduce) {
  UnitTestReduce< double ,   Kokkos::DevicePthread >( 1000000 );
}

TEST_F( pthread, long_multi_reduce) {
  UnitTestReduceMulti< long , Kokkos::DevicePthread >( 1000000 , 7 );
}

} // namespace test

