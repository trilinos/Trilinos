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

#include <iostream>

#include <KokkosArray_View.hpp>

#include <KokkosArray_CrsArray.hpp>

#include <KokkosArray_Host.hpp>
#include <KokkosArray_Cuda.hpp>

//----------------------------------------------------------------------------

#include <TestViewImpl.hpp>

#include <KokkosArray_Cuda_macros.hpp>

#include <TestViewAPI.hpp>
#include <TestCrsArray.hpp>

#include <TestReduce.hpp>
#include <TestMultiReduce.hpp>

#include <KokkosArray_Clear_macros.hpp>

namespace Test {

void test_device_cuda_init() {
  KokkosArray::Cuda::initialize();
}

void test_device_cuda_view_impl()
{
  test_view_impl< KokkosArray::Cuda >();
}

void test_device_cuda_view_api()
{
  TestViewAPI< double , KokkosArray::Cuda >();
}

void test_device_cuda_crsarray() {
  TestCrsArray< KokkosArray::Cuda >();
}

void test_device_cuda_reduce() {
  TestReduce< long ,   KokkosArray::Cuda >( 10000000 );
  TestReduce< long ,   KokkosArray::Cuda >( 1000000 );
  TestReduce< double , KokkosArray::Cuda >( 1000000 );
}

void test_device_cuda_multi_reduce() {
  TestReduceMulti< long , KokkosArray::Cuda >( 1000000 , 7 );
}

}

