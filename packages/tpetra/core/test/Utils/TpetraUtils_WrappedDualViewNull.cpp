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
