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

#include <impl/KokkosArray_Timer.hpp>

#include <KokkosArray_Host.hpp>
#include <KokkosArray_hwloc.hpp>

#include <PerfTestHexGrad.hpp>
#include <PerfTestBlasKernels.hpp>
#include <PerfTestGramSchmidt.hpp>
#include <PerfTestDriver.hpp>

//------------------------------------------------------------------------

namespace Test {

class host : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    const std::pair<unsigned,unsigned> core_topo = KokkosArray::hwloc::get_core_topology();
    const unsigned                     core_size = KokkosArray::hwloc::get_core_capacity();
    const unsigned thread_count = std::min( 8u , core_topo.first *
                                                 core_topo.second *
                                                 core_size );

    KokkosArray::Host::initialize( core_topo.first , thread_count / core_topo.first );
  }

  static void TearDownTestCase()
  {
    KokkosArray::Host::finalize();
  }
};

TEST_F( host, hexgrad ) {
  EXPECT_NO_THROW(run_test_hexgrad< KokkosArray::Host>( 10, 20 ));
}

TEST_F( host, gramschmidt ) {
  EXPECT_NO_THROW(run_test_gramschmidt< KokkosArray::Host>( 10, 20 ));
}

} // namespace Test


