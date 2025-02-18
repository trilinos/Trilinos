// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_WrappedDualView.hpp"
#include "Tpetra_Access.hpp"
#include "Tpetra_Core.hpp"
#include "Kokkos_DualView.hpp"
#include "Teuchos_UnitTestHarness.hpp"

namespace {

TEUCHOS_UNIT_TEST(WrappedDualView, InitFromNull) {

  // Test that wrapped dual views return dual views with
  // equal host and device counts, regardless of how they are constructed.
  // Note that constructing a DualView with a null Kokkos device view
  // and its host mirror returns use_count = 2 on host and 0 on device.

  std::cout << "\n" << std::endl;

  using device_t =  Kokkos::Device<Kokkos::DefaultExecutionSpace,
                                   Kokkos::DefaultExecutionSpace::memory_space>;

  using view_t = Kokkos::View<double*, device_t>;
  using view_mirror_t = typename view_t::HostMirror;

  using dualview_t = Kokkos::DualView<double*, device_t>;
  using wrapped_t = Tpetra::Details::WrappedDualView<dualview_t>;

  {
    // Here's why we need special logic in WrappedDualView constructor
    // when the input device view is NULL.

    view_t viewNull;
    view_mirror_t viewNull_mirror =
       create_mirror_view_and_copy(typename view_t::host_mirror_space(),
                                   viewNull);

    dualview_t dvNull(viewNull, viewNull_mirror);
    size_t use_h = dvNull.h_view.use_count();
    size_t use_d = dvNull.d_view.use_count();

    std::cout << "Null DualView:     "
              << "host.use_count = " << use_h << ";  "
              << "device.use_count = " << use_d 
              << std::endl;
    // For UVM or serial builds, use_h == use_d == 0.
    // But for non-UVM CUDA builds, use_h == 2 and use_d == 0.
    // This difference is bad for WrappedDualView.
    // Thus, WrappedDualView's constructor needs to check for a 
    // null device view before creating the HostMirror.
  }

  {
    // Happy case:  device view is non-null
    wrapped_t wrapped;
    {
      view_t v("viewFour", 4);
      wrapped = wrapped_t(v);
    }
    size_t use_h = wrapped.getHostView(Tpetra::Access::ReadOnly).use_count();
    size_t use_d = wrapped.getDeviceView(Tpetra::Access::ReadOnly).use_count();
    std::cout << "Wrapped "
              << wrapped.getDeviceView(Tpetra::Access::ReadOnly).label()
              << ":  host use_count = " << use_h
              << ";  device use_count = " << use_d << std::endl;
    TEST_EQUALITY(use_h, use_d);
  }

  {
    // Happy case:  device view is non-null, even though length zero
    wrapped_t wrapped;
    {
      view_t v("viewZero", 0);
      wrapped = wrapped_t(v);
    }
    size_t use_h = wrapped.getHostView(Tpetra::Access::ReadOnly).use_count();
    size_t use_d = wrapped.getDeviceView(Tpetra::Access::ReadOnly).use_count();
    std::cout << "Wrapped "
              << wrapped.getDeviceView(Tpetra::Access::ReadOnly).label()
              << ":  host use_count = " << use_h
              << ";  device use_count = " << use_d << std::endl;
    TEST_EQUALITY(use_h, use_d);
  }

  {
    // Unhappy case:  null device view won't work with HostMirror in
    //                DualView constructor in WrappedDualView constructor
    wrapped_t wrapped;
    {
      view_t v;
      wrapped = wrapped_t(v);
    }
    size_t use_h = wrapped.getHostView(Tpetra::Access::ReadOnly).use_count();
    size_t use_d = wrapped.getDeviceView(Tpetra::Access::ReadOnly).use_count();
    std::cout << "Wrapped nullview"
              << wrapped.getDeviceView(Tpetra::Access::ReadOnly).label()
              << ":  host use_count = " << use_h
              << ";  device use_count = " << use_d << std::endl;
    TEST_EQUALITY(use_h, use_d);
  }
}

} // namespace

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard scopeGuard(&argc, &argv);
  const int errCode = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  return errCode;
}
