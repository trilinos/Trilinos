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


