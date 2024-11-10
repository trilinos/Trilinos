// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Kokkos_StaticCrsGraph.hpp"

#include <Tpetra_Details_WrappedDualView.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>

#include <Kokkos_DualView.hpp>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StackedTimer.hpp>

namespace {

using DeviceType = Tpetra::Map<>::device_type;

using DualViewType = Kokkos::DualView<int*, DeviceType>;
using WrappedDualViewType = Tpetra::Details::WrappedDualView<DualViewType>;

using HostViewType = typename DualViewType::t_host;
using DeviceViewType = typename DualViewType::t_dev;
using ConstDeviceViewType = typename DualViewType::t_dev::const_type;

using ConstDualViewType = Kokkos::DualView<const int*, DeviceType>;
using WrappedConstDualViewType = Tpetra::Details::WrappedDualView<ConstDualViewType>;

class WrappedDualViewFixture {
public:
  static constexpr bool deviceMemoryIsHostAccessible = Kokkos::SpaceAccessibility<Kokkos::Serial, typename DeviceType::memory_space>::accessible;

  WrappedDualViewFixture()
    : viewSize(16),
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

  ConstDualViewType getConstDualView() {
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

  void fillDualViewOnHostDevice() {
    fillDualViewOnHost();
    dualView.sync_device();
  }

  bool valuesInitializedToZero() {
    auto hostView = dualView.view_host();
    auto deviceView = dualView.view_device();
    return valuesCorrectOnHost(hostView, 0) && valuesCorrectOnDevice(deviceView, 0);
  }

  template <typename ViewType>
  void fillViewOnHost(ViewType view) {
    fillViewOnHost(view, 0, viewSize);
  }

  template <typename ViewType>
  void fillViewOnHost(ViewType view, int startIndex, int length) {
    for (int i=0; i<length; i++) {
      int value = i + startIndex;
      view(i) = value;
    }
  }

  template <typename ViewType>
  void multiplyOnHost(ViewType view, int multiplier) {
    for (unsigned i=0; i<view.size(); i++) {
      view(i) = multiplier*view(i);
    }
  }

  template <typename ViewType>
  bool valuesCorrectOnHost(ViewType view, int multiplier = 1) {
    return valuesCorrectOnHost(view, 0, viewSize, multiplier);
  }

  template <typename ViewType>
  bool valuesCorrectOnHost(ViewType view, int startIndex, int length, int multiplier = 1) {
    bool result = (static_cast<int>(view.size()) == length);
    for (int i=0; i<length && result; i++) {
      int value = multiplier*(i + startIndex);
      result &= (view(i) == value);
    }
    return result;
  }

  template <typename ViewType>
  void fillViewOnDevice(ViewType view) {
    fillViewOnDevice(view, 0, viewSize);
  }

  template <typename ViewType>
  void fillViewOnDevice(ViewType view, int startIndex, int length) {
    Kokkos::parallel_for("fill on device", length, KOKKOS_LAMBDA(const int& i) {
          int value = i + startIndex;
          view(i) = value;
        });
  }

  template <typename ViewType>
  void multiplyOnDevice(ViewType view, int multiplier) {
    Kokkos::parallel_for("multiply on device", view.size(), KOKKOS_LAMBDA(const int& i) {
          view(i) = multiplier*view(i);
        });
  }

  template <typename ViewType>
  bool valuesCorrectOnDevice(ViewType view, int multiplier = 1) {
    return valuesCorrectOnDevice(view, 0, viewSize, multiplier);
  }

  template <typename ViewType>
  bool valuesCorrectOnDevice(ViewType view, int startIndex, int length, int multiplier = 1) {
    int result = 0;
    if (static_cast<int>(view.size()) != length) {
      result++;
    }
    else {
      Kokkos::parallel_reduce("check on device", length,
          KOKKOS_LAMBDA(const int& i, int& localResult) {
            int value = multiplier*(i + startIndex);
            localResult = (view(i) == value) ? 0 : 1;
          }, result);
    }
    return (result == 0);
  }

  template <typename ViewType>
  bool extentCorrect(ViewType view) {
    return (view.extent(0) == dualView.extent(0));
  }

  int getViewSize() {
    return viewSize;
  }

private:
  int viewSize;
  DualViewType dualView;
};


}
