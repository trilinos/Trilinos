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

  using map_type = Tpetra::Map<>;
  using device_type = map_type::device_type;
  using device_view_type = Kokkos::View<int*, device_type>;

  using Tpetra::Details::getEntryOnHost;
  using Tpetra::Details::getEntriesOnHost;

  void ensureKokkosIsInitializedCorrectly() {
    auto comm = Tpetra::TestingUtilities::getDefaultComm();
    const map_type map(comm->getSize(), 1, 0, comm);
  }

  TEUCHOS_UNIT_TEST( TpetraUtils, get_entry_on_host )
  {
    ensureKokkosIsInitializedCorrectly();

    const int SIZE = 42;

    device_view_type deviceView("deviceView", SIZE);
    Kokkos::parallel_for("fill deviceView", SIZE, KOKKOS_LAMBDA(const int& i) {
          deviceView(i) = i+1;
        });

    for (int i=0; i<SIZE; ++i) {
      TEST_EQUALITY( getEntryOnHost(deviceView, i), i+1 );
    }
  }

  TEUCHOS_UNIT_TEST( TpetraUtils, get_two_entries_on_host )
  {
    ensureKokkosIsInitializedCorrectly();

    const int SIZE = 42;

    device_view_type deviceView("deviceView", SIZE);
    Kokkos::parallel_for("fill deviceView", SIZE, KOKKOS_LAMBDA(const int& i) {
          deviceView(i) = i+1;
        });

    for (int i=0; i<SIZE-1; ++i) {
      auto entries = getEntriesOnHost(deviceView, i, 2);
      TEST_EQUALITY( entries(0), i+1 );
      TEST_EQUALITY( entries(1), i+2 );
    }
  }

  TEUCHOS_UNIT_TEST( TpetraUtils, get_all_entries_on_host )
  {
    ensureKokkosIsInitializedCorrectly();

    const int SIZE = 42;

    device_view_type deviceView("deviceView", SIZE);
    Kokkos::parallel_for("fill deviceView", SIZE, KOKKOS_LAMBDA(const int& i) {
          deviceView(i) = i+1;
        });

    auto entries = getEntriesOnHost(deviceView, 0, SIZE);
    for (int i=0; i<SIZE-1; ++i) {
      TEST_EQUALITY( entries(i), i+1 );
    }
  }

} // namespace (anonymous)


