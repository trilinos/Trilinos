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

// Force OMP atomics

#define KOKKOS_ATOMICS_USE_OMP31
#include <Kokkos_Atomic.hpp>

#include <Kokkos_OpenMP.hpp>
#include <Kokkos_hwloc.hpp>

#include <Kokkos_UnorderedMap.hpp>

//----------------------------------------------------------------------------
#include <TestUnorderedMapInsert.hpp>

namespace Test {

class openmp : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    const std::pair<unsigned,unsigned> core_top =  Kokkos::hwloc::get_core_topology();
    const unsigned core_size = Kokkos::hwloc::get_core_capacity();

    const unsigned gang_count        = core_top.first ;
    const unsigned gang_worker_count = core_top.second ;

    Kokkos::OpenMP::initialize( gang_count , gang_worker_count*core_size );
  }

  static void TearDownTestCase()
  {
    Kokkos::OpenMP::finalize();

    omp_set_num_threads(0);

    ASSERT_EQ( 1 , omp_get_max_threads() );
  }
};

TEST_F( openmp, unordered_map_insert_10000000_9000000_9000000) {
  //9 million inserts of the same key
  uint32_t num_nodes      = 10000000;
  uint32_t num_inserts    =  9000000;
  uint32_t num_duplicates =  9000000;
  test_unordered_map_insert< Kokkos::OpenMP >(num_nodes,num_inserts,num_duplicates);
}

TEST_F( openmp, unordered_map_insert_10000000_9000000_3000000) {
  //9 million inserts of the same key
  uint32_t num_nodes      = 10000000;
  uint32_t num_inserts    =  9000000;
  uint32_t num_duplicates =  3000000;
  test_unordered_map_insert< Kokkos::OpenMP >(num_nodes,num_inserts,num_duplicates);
}

TEST_F( openmp, unordered_map_insert_10000000_9000000_30) {
  //9 million inserts of the same key
  uint32_t num_nodes      = 10000000;
  uint32_t num_inserts    =  9000000;
  uint32_t num_duplicates =  30;
  test_unordered_map_insert< Kokkos::OpenMP >(num_nodes,num_inserts,num_duplicates);
}

TEST_F( openmp, unordered_map_insert_10000000_9000000_3) {
  //9 million inserts of the same key
  uint32_t num_nodes      = 10000000;
  uint32_t num_inserts    =  9000000;
  uint32_t num_duplicates =  3;
  test_unordered_map_insert< Kokkos::OpenMP >(num_nodes,num_inserts,num_duplicates);
}

TEST_F( openmp, unordered_map_insert_10000000_9000000_1) {
  //9 million inserts of the same key
  uint32_t num_nodes      = 10000000;
  uint32_t num_inserts    =  9000000;
  uint32_t num_duplicates =  1;
  test_unordered_map_insert< Kokkos::OpenMP >(num_nodes,num_inserts,num_duplicates);
}

TEST_F( openmp, unordered_map_insert_10000000_2000000_1) {
  //9 million inserts of the same key
  uint32_t num_nodes      = 10000000;
  uint32_t num_inserts    = 20000000;
  uint32_t num_duplicates = 1;
  test_unordered_map_insert< Kokkos::OpenMP >(num_nodes,num_inserts,num_duplicates);
}

} // namespace test

