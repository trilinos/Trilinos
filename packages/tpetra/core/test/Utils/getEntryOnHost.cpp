// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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


