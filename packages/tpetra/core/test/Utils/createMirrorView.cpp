/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Details_createMirrorView.hpp"

namespace { // (anonymous)

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( TpetraUtils, create_mirror_view_from_raw_host_array )
  {
    using Tpetra::Details::create_mirror_view_from_raw_host_array;

    auto comm = Tpetra::TestingUtilities::getDefaultComm ();
    // Create a Map just to ensure that Kokkos gets initialized and finalized correctly.
    const Tpetra::Map<> map (comm->getSize (), 1, 0, comm);
    Tpetra::Map<>::device_type outputDevice;

    const int size = 42;
    {
      std::vector<int> x (size);
      for (int k = 0; k < size; ++k) {
        x[k] = k+1;
      }
      // Test the case where the input pointer is nonconst.
      auto x_view = create_mirror_view_from_raw_host_array (outputDevice,
                                                            x.data (),
                                                            x.size ());
      auto x_host = Kokkos::create_mirror_view (x_view);
      Kokkos::deep_copy (x_host, x_view);
      for (int k = 0; k < size; ++k) {
        TEST_EQUALITY(x_host(k), x[k]);
      }
    }
    {
      std::vector<int> x (size);
      for (int k = 0; k < size; ++k) {
        x[k] = k+1;
      }
      // Test the case where the input pointer is const.
      const int* x_raw = x.data ();
      auto x_view = create_mirror_view_from_raw_host_array (outputDevice,
                                                            x_raw,
                                                            x.size ());
      auto x_host = Kokkos::create_mirror_view (x_view);
      Kokkos::deep_copy (x_host, x_view);
      for (int k = 0; k < size; ++k) {
        TEST_EQUALITY(x_host(k), x[k]);
      }
    }
  }

} // namespace (anonymous)


