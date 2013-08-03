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
#include <TestUnorderedMapInsert.hpp>

namespace Test {

class host : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    std::cout << std::setprecision(5) << std::scientific;
#if 1
    const std::pair<unsigned,unsigned> core_top =  Kokkos::hwloc::get_core_topology();
    const unsigned core_size = Kokkos::hwloc::get_core_capacity();

    const unsigned gang_count        = core_top.first ;
    const unsigned gang_worker_count = core_top.second * core_size ;

    Kokkos::Host::initialize( gang_count , gang_worker_count );
#else
    Kokkos::Host::initialize(1,1);
#endif

  }

  static void TearDownTestCase()
  {
    Kokkos::Host::finalize();
  }
};

#define HOST_INSERT_TEST( name, num_nodes, num_inserts, num_duplicates, repeat )                     \
  TEST_F( host, unordered_map_insert_##name##_##num_nodes##_##num_inserts##_##num_duplicates##_##repeat) {   \
    map_test_times test_times;                                                                       \
    for (int i=0; i<repeat; ++i)                                                                     \
      test_insert_##name<Kokkos::Host>(num_nodes,num_inserts,num_duplicates,test_times);             \
    std::cout << "Test Times" << std::endl;                                                          \
    std::cout << "      Construct: " << test_times.construct / repeat << std::endl;                  \
    std::cout << "  Santity Check: " << test_times.santity_check / repeat << std::endl;              \
    std::cout << "         Insert: " << test_times.insert / repeat << std::endl;                     \
    std::cout << "           Find: " << test_times.find / repeat << std::endl;                       \
  }

HOST_INSERT_TEST(close,     10000,     6000, 100, 10000)
HOST_INSERT_TEST(close,    100000,    90000, 100, 5000)
HOST_INSERT_TEST(close,   1000000,   900000, 100, 500)
HOST_INSERT_TEST(close,  10000000,  9000000, 100, 50)
HOST_INSERT_TEST(close, 100000000, 90000000, 100, 5)
HOST_INSERT_TEST(close, 100000000, 90000000,  10, 5)

HOST_INSERT_TEST(far,     10000,     6000, 100, 10000)
HOST_INSERT_TEST(far,    100000,    90000, 100, 5000)
HOST_INSERT_TEST(far,   1000000,   900000, 100, 500)
HOST_INSERT_TEST(far,  10000000,  9000000, 100, 50)
HOST_INSERT_TEST(far, 100000000, 90000000, 100, 5)
HOST_INSERT_TEST(far, 100000000, 90000000,  10, 5)

#undef HOST_INSERT_TEST

} // namespace Test

