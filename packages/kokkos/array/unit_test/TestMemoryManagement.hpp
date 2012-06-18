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

#ifndef KOKKOS_MACRO_DEVICE
#error "KOKKOS_MACRO_DEVICE undefined"
#endif

#include <stdexcept>
#include <sstream>
#include <iostream>
#include <impl/KokkosArray_Preprocessing_macros.hpp>

/*--------------------------------------------------------------------------*/

namespace {

template< typename T, class > class TestMemoryManagement ;

template<typename T>
class TestMemoryManagement< T, KOKKOS_MACRO_DEVICE >
{
public:
  typedef KOKKOS_MACRO_DEVICE device ;

  typedef KokkosArray::Impl::MemoryManager< device::memory_space > Manager ;
  typedef KokkosArray::Impl::MemoryView< T , device::memory_space > dView ;

  TestMemoryManagement() { run_test(); }

  void run_test()
  {
    std::ostringstream s_init ;
    std::ostringstream s_clear ;

    Manager::print_memory_view( s_init );

    dView dx , dy , dz ;

    const bool error_constructing_views =
       ( dx || dx.test_support_view_count() ||
         dy || dy.test_support_view_count() ||
         dz || dz.test_support_view_count() );

    ASSERT_FALSE(error_constructing_views);

    dx.allocate( 10 , "dx" );

    const bool error_initializing_views =
       ( ! dx || 1 != dx.test_support_view_count() );

    ASSERT_FALSE(error_initializing_views);

    dy = dx ;

    const bool error_copying_views =
       ( ! dx || 2 != dx.test_support_view_count() ||
         ! dy || 2 != dy.test_support_view_count() ||
         dx != dy );

    ASSERT_FALSE(error_copying_views);

    dx = dView();

    bool error_clearing_views =
       ( dx || dx.test_support_view_count() ||
         ! dy || 1 != dy.test_support_view_count() );

    ASSERT_FALSE(error_clearing_views);

    dy = dView();

    Manager::print_memory_view( s_clear );

    error_clearing_views =
       ( dx || dx.test_support_view_count() ||
         dy || dy.test_support_view_count() ||
         dz || dz.test_support_view_count() ||
         s_init.str() != s_clear.str() );

    ASSERT_FALSE(error_clearing_views);

  }
};

}

/*--------------------------------------------------------------------------*/

