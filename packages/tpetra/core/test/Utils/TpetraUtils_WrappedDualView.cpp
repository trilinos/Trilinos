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

#include <Tpetra_Details_WrappedDualView.hpp>
#include <Tpetra_Access.hpp>
#include <Tpetra_Core.hpp>

#include <Kokkos_DualView.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

namespace {

using DualViewType = Kokkos::DualView<int*, Kokkos::Cuda>;
using WrappedDualViewType = Tpetra::Details::WrappedDualView<DualViewType>;

class WrappedDualViewFixture {
public:
  WrappedDualViewFixture()
    : comm(Teuchos::DefaultComm<int>::getComm()),
      commSize(comm->getSize()),
      viewSize(commSize*4),
      dualView("dualView", viewSize)
  {
    for (int i=0; i<viewSize; i++) {
      dualView.view_host()(i) = 0;
    }
    dualView.modify_host();
    dualView.sync_device();
  }

  DualViewType getDualView() {
    return dualView;
  }

  void fillDualViewOnHost() {
    auto hostView = dualView.view_host();
    fillViewOnHost(hostView);
    dualView.modify_host();
  }

  void fillDualViewOnDevice() {
    auto deviceView = dualView.view_device();
    fillViewOnDevice(deviceView);
    dualView.modify_device();
  }

  bool valuesInitializedToZero() {
    auto hostView = dualView.view_host();
    auto deviceView = dualView.view_device();
    return valuesCorrectOnHost(hostView, 0) && valuesCorrectOnDevice(deviceView, 0);
  }

  template <typename ViewType>
  void fillViewOnHost(ViewType view) {
    for (int i=0; i<viewSize; i++) {
      view(i) = i;
    }
  }

  template <typename ViewType>
  void multiplyOnHost(ViewType view, int multiplier) {
    for (int i=0; i<viewSize; i++) {
      view(i) = multiplier*view(i);
    }
  }

  template <typename ViewType>
  bool valuesCorrectOnHost(ViewType view, int multiplier = 1) {
    bool result = true;
    for (int i=0; i<viewSize; i++) {
      result &= (view(i) == multiplier*i);
    }
    return result;
  }

  template <typename ViewType>
  void fillViewOnDevice(ViewType view) {
    Kokkos::parallel_for("fill on device", viewSize, KOKKOS_LAMBDA(const int& i) {
          view(i) = i;
        });
  }

  template <typename ViewType>
  void multiplyOnDevice(ViewType view, int multiplier) {
    Kokkos::parallel_for("multiply on device", viewSize, KOKKOS_LAMBDA(const int& i) {
          view(i) = multiplier*view(i);
        });
  }

  template <typename ViewType>
  bool valuesCorrectOnDevice(ViewType view, int multiplier = 1) {
    int result = 0;
    Kokkos::parallel_reduce("check on device", viewSize,
        KOKKOS_LAMBDA(const int& i, int& localResult) {
          localResult = (view(i) == multiplier*i) ? 0 : 1;
        }, result);
    return (result == 0);
  }

private:
  Teuchos::RCP<const Teuchos::Comm<int>> comm;
  int commSize;
  int viewSize;
  DualViewType dualView;
};

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostReadOnly) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHost();

  const WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostReadWrite) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHost();

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadWrite);
  fixture.multiplyOnHost(hostView, 2);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostView, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostWriteOnly) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::WriteOnly);
  fixture.fillViewOnHost(hostView);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceReadOnly) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnDevice();

  const WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceReadWrite) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnDevice();

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadWrite);
  fixture.multiplyOnDevice(deviceView, 2);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceView, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceWriteOnly) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::WriteOnly);
  fixture.fillViewOnDevice(deviceView);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, hostThrowsIfDeviceViewAlive_ReadOnly) {
  WrappedDualViewFixture fixture;

  const WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadOnly);
  TEST_THROW(wrappedView.getHostView(Tpetra::Access::ReadOnly), std::runtime_error);
}

TEUCHOS_UNIT_TEST(WrappedDualView, hostThrowsIfDeviceViewAlive_ReadWrite) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadWrite);
  TEST_THROW(wrappedView.getHostView(Tpetra::Access::ReadWrite), std::runtime_error);
}

TEUCHOS_UNIT_TEST(WrappedDualView, hostThrowsIfDeviceViewAlive_WriteOnly) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::WriteOnly);
  TEST_THROW(wrappedView.getHostView(Tpetra::Access::WriteOnly), std::runtime_error);
}

TEUCHOS_UNIT_TEST(WrappedDualView, deviceThrowsIfHostViewAlive_ReadOnly) {
  WrappedDualViewFixture fixture;

  const WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadOnly);
  TEST_THROW(wrappedView.getDeviceView(Tpetra::Access::ReadOnly), std::runtime_error);
}

TEUCHOS_UNIT_TEST(WrappedDualView, deviceThrowsIfHostViewAlive_ReadWrite) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadWrite);
  TEST_THROW(wrappedView.getDeviceView(Tpetra::Access::ReadWrite), std::runtime_error);
}

TEUCHOS_UNIT_TEST(WrappedDualView, deviceThrowsIfHostViewAlive_WriteOnly) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::WriteOnly);
  TEST_THROW(wrappedView.getDeviceView(Tpetra::Access::WriteOnly), std::runtime_error);
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostReadOnly_syncToHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());
  fixture.fillDualViewOnDevice();

  const WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostReadWrite_syncToHost_modifyOnHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());
  fixture.fillDualViewOnDevice();

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::ReadWrite);
    TEST_ASSERT(fixture.valuesCorrectOnHost(hostView));
    fixture.multiplyOnHost(hostView, 2);
  }

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceView, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostWriteOnly_clearSyncState_modifyOnHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::WriteOnly);
    fixture.fillViewOnDevice(deviceView);
    fixture.multiplyOnDevice(deviceView, 2);
  }

  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::WriteOnly);
    fixture.fillViewOnHost(hostView);
  }

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceReadOnly_syncToDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());
  fixture.fillDualViewOnHost();

  const WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceReadWrite_syncToDevice_modifyOnDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());
  fixture.fillDualViewOnHost();

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadWrite);
    TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceView));
    fixture.multiplyOnDevice(deviceView, 2);
  }

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostView, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceWriteOnly_clearSyncState_modifyOnDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::WriteOnly);
    fixture.fillViewOnHost(hostView);
    fixture.multiplyOnHost(hostView, 2);
  }

  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::WriteOnly);
    fixture.fillViewOnDevice(deviceView);
  }

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostView));
}

}

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard scopeGuard(&argc, &argv);
  const int errCode = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  return errCode;
}
