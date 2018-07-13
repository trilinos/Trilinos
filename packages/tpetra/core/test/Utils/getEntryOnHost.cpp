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
#include "Tpetra_Details_getEntryOnHost.hpp"

namespace { // (anonymous)

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST( TpetraUtils, get_entry_on_host )
  {
    using Tpetra::Details::getEntryOnHost;
    typedef Tpetra::Map<> map_type;
    typedef map_type::device_type device_type;

    auto comm = Tpetra::TestingUtilities::getDefaultComm ();
    // Create a Map just to ensure that Kokkos gets initialized and
    // finalized correctly.
    const map_type map (comm->getSize (), 1, 0, comm);

    typedef Kokkos::View<int*, device_type> dev_view_type;
    typedef dev_view_type::HostMirror::execution_space host_exec_space;
    // Don't just use HostMirror's memory space, because that could be
    // the same memory space (in the case of CudaUVMSpace).
    typedef Kokkos::Device<host_exec_space, Kokkos::HostSpace> host_device_type;
    // Same array layout means we can deep_copy between host_view_type
    // and dev_view_type.
    typedef Kokkos::View<int*, dev_view_type::array_layout, host_device_type> host_view_type;

    const int size = 42;
    host_view_type x_h ("x_h", size);
    for (int i = 0; i < size; ++i) {
      x_h(i) = i+1; // no entries are zero
    }
    dev_view_type x_d ("x_d", size);
    Kokkos::deep_copy (x_d, x_h);

    // Make sure that x_h and x_d really are distinct.  Otherwise,
    // getEntryOnHost might not be doing what we expect.
    for (int i = 0; i < size; ++i) {
      x_h(i) = -(i+1);
    }

    for (int i = 0; i < size; ++i) {
      const int curEnt = getEntryOnHost (x_d, i);
      TEST_EQUALITY( curEnt, i+1 );
    }
  }

} // namespace (anonymous)


