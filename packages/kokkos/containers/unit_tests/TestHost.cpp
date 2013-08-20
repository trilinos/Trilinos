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

#include <Kokkos_Host.hpp>
#include <Kokkos_hwloc.hpp>

#include <Kokkos_UnorderedMap.hpp>

#include <iomanip>


//----------------------------------------------------------------------------
#include <TestUnorderedMap.hpp>

namespace Test {

class host : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    std::cout << std::setprecision(5) << std::scientific;
    Kokkos::Host::initialize();
  }

  static void TearDownTestCase()
  {
    Kokkos::Host::finalize();
  }
};

#define HOST_INSERT_TEST( name, num_nodes, num_inserts, num_duplicates, repeat )                                \
  TEST_F( host, unordered_map_insert_##name##_##num_nodes##_##num_inserts##_##num_duplicates##_##repeat##x) {   \
    for (int i=0; i<repeat; ++i)                                                                                \
      test_insert_##name<Kokkos::Host>(num_nodes,num_inserts,num_duplicates);                                   \
  }

#define HOST_FAILED_INSERT_TEST( num_nodes, repeat )                            \
  TEST_F( host, unordered_map_failed_insert_##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      test_failed_insert<Kokkos::Host>(num_nodes);                             \
  }

#define HOST_ASSIGNEMENT_TEST( num_nodes, repeat )                             \
  TEST_F( host, unordered_map_assignment_operators_##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      test_assignement_operators<Kokkos::Host>(num_nodes);                     \
  }

#define HOST_DEEP_COPY( num_nodes, repeat )                             \
  TEST_F( host, unordered_map_deep_copy##num_nodes##_##repeat##x) {       \
    for (int i=0; i<repeat; ++i)                                               \
      test_deep_copy<Kokkos::Host>(num_nodes);                     \
  }

HOST_INSERT_TEST(close,               100000, 90000, 100, 1)
HOST_INSERT_TEST(far,                 100000, 90000, 100, 1)
HOST_INSERT_TEST(mark_pending_delete, 100000, 90000, 100, 1)
HOST_FAILED_INSERT_TEST( 10000, 1 )
HOST_ASSIGNEMENT_TEST( 10000, 1 )
HOST_DEEP_COPY( 10000, 1 )



#undef HOST_INSERT_TEST
#undef HOST_FAILED_INSERT_TEST
#undef HOST_ASSIGNEMENT_TEST
#undef HOST_DEEP_COPY

} // namespace Test

