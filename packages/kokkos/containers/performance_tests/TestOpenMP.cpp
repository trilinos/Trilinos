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

#include <iomanip>

#include <TestGlobal2LocalIds.hpp>

namespace Performance {

class openmp : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    std::cout << std::setprecision(5) << std::scientific;

    std::pair<unsigned, unsigned> team_league(1,4);
    if (Kokkos::hwloc::available()) {
      const std::pair<unsigned,unsigned> core_top =  Kokkos::hwloc::get_core_topology();
      team_league.first  = core_top.first ;
      team_league.second = core_top.second;
    }

    std::cout << "Threads: " << team_league.first << "x" << team_league.second << std::endl;

    Kokkos::OpenMP::initialize( team_league.first , team_league.second );
  }

  static void TearDownTestCase()
  {
    Kokkos::OpenMP::finalize();

    omp_set_num_threads(0);

    ASSERT_EQ( 1 , omp_get_max_threads() );
  }
};

TEST_F( openmp, global_2_local)
{
  std::cout << "OpenMP" << std::endl;
  std::cout << "size, create, generate, fill, find" << std::endl;
  for (unsigned i=Performance::begin_id_size; i<=Performance::end_id_size; i *= Performance::id_step)
    test_global_to_local_ids<Kokkos::OpenMP>(i);
}

} // namespace test

