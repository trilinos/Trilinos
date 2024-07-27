// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Kokkos_StaticCrsGraph.hpp"

#include "TpetraUtils_WrappedDualViewFixture.hpp"

#include <Tpetra_Details_WrappedDualView.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Access.hpp>
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

TEUCHOS_UNIT_TEST(WrappedDualView, defaultConstructorAvailable) {
  const WrappedDualViewType wrappedView;
}

TEUCHOS_UNIT_TEST(WrappedDualView, deviceViewConstructor) {
  WrappedDualViewFixture fixture;
  WrappedDualViewType wrappedView;

  {
    DeviceViewType deviceView("device view", fixture.getViewSize());
    fixture.fillViewOnDevice(deviceView);

    wrappedView = WrappedDualViewType(deviceView);
  }

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, extent) {
  WrappedDualViewFixture fixture;
  const WrappedDualViewType wrappedView(fixture.getDualView());

  TEST_ASSERT(fixture.extentCorrect(wrappedView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostReadOnly_constData) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();

  const WrappedConstDualViewType wrappedView(fixture.getConstDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceReadOnly_constData) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();

  const WrappedConstDualViewType wrappedView(fixture.getConstDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostSubviewReadOnly_constData) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();

  const WrappedConstDualViewType wrappedView(fixture.getConstDualView());

  int startIndex = 4;
  int length = 6;
  auto hostSubview = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubview, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceSubviewReadOnly_constData) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();

  const WrappedConstDualViewType wrappedView(fixture.getConstDualView());

  int startIndex = 4;
  int length = 6;
  auto deviceSubview = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubview, startIndex, length));
}

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

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostOverwriteAll) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
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

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceOverwriteAll) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
  fixture.fillViewOnDevice(deviceView);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, hostThrowsIfDeviceViewAlive_ReadOnly) {
  WrappedDualViewFixture fixture;

  const WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadOnly);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getHostView(Tpetra::Access::ReadOnly));
  }
  else {
    TEST_THROW(wrappedView.getHostView(Tpetra::Access::ReadOnly), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, hostThrowsIfDeviceViewAlive_ReadWrite) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::ReadWrite);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getHostView(Tpetra::Access::ReadWrite));
  }
  else {
    TEST_THROW(wrappedView.getHostView(Tpetra::Access::ReadWrite), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, hostThrowsIfDeviceViewAlive_OverwriteAll) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getHostView(Tpetra::Access::OverwriteAll));
  }
  else {
    TEST_THROW(wrappedView.getHostView(Tpetra::Access::OverwriteAll), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, deviceThrowsIfHostViewAlive_ReadOnly) {
  WrappedDualViewFixture fixture;

  const WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadOnly);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getDeviceView(Tpetra::Access::ReadOnly));
  }
  else {
    TEST_THROW(wrappedView.getDeviceView(Tpetra::Access::ReadOnly), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, deviceThrowsIfHostViewAlive_ReadWrite) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadWrite);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getDeviceView(Tpetra::Access::ReadWrite));
  }
  else {
    TEST_THROW(wrappedView.getDeviceView(Tpetra::Access::ReadWrite), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, deviceThrowsIfHostViewAlive_OverwriteAll) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getDeviceView(Tpetra::Access::OverwriteAll));
  }
  else {
    TEST_THROW(wrappedView.getDeviceView(Tpetra::Access::OverwriteAll), std::runtime_error);
  }
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

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostOverwriteAll_clearSyncState_modifyOnHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceView);
    fixture.multiplyOnDevice(deviceView, 2);
  }

  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
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

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceOverwriteAll_clearSyncState_modifyOnDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostView);
    fixture.multiplyOnHost(hostView, 2);
  }

  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceView);
  }

  auto hostView = wrappedView.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostView));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostSubviewReadOnly) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHost();

  const WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 3;
  int length = 4;
  auto hostSubview = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubview, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostSubviewReadWrite) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHost();

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 2;
  int length = 5;
  auto hostSubview = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadWrite);
  fixture.multiplyOnHost(hostSubview, 2);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubview, startIndex, length, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostSubviewOverwriteAll) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 5;
  int length = 5;
  auto hostSubview = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::OverwriteAll);
  fixture.fillViewOnHost(hostSubview, startIndex, length);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubview, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceSubviewReadOnly) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnDevice();

  const WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 3;
  int length = 4;
  auto deviceSubview = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubview, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceSubviewReadWrite) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnDevice();

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 2;
  int length = 5;
  auto deviceSubview = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadWrite);
  fixture.multiplyOnDevice(deviceSubview, 2);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubview, startIndex, length, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceSubviewOverwriteAll) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 5;
  int length = 5;
  auto deviceSubview = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::OverwriteAll);
  fixture.fillViewOnDevice(deviceSubview, startIndex, length);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubview, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostSubviewReadOnly_syncToHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());
  fixture.fillDualViewOnDevice();

  const WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 4;
  int length = 3;
  auto hostSubview = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubview, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostSubviewReadWrite_syncToHost_modifyOnHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());
  fixture.fillDualViewOnDevice();

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 5;
  int length = 2;
  {
    auto hostSubview = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadWrite);
    TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubview, startIndex, length));
    fixture.multiplyOnHost(hostSubview, 2);
  }

  auto deviceSubview = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubview, startIndex, length, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostSubviewOverwriteAll_syncToHost_modifyOnHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceView);
    fixture.multiplyOnDevice(deviceView, 2);
  }

  int startIndex = 0;
  int length = 4;
  {
    auto hostSubview = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostSubview, startIndex, length);
  }

  int startIndexUnchanged = length;
  int lengthUnchanged = fixture.getViewSize() - length;

  {
    auto hostSubviewFromDevice = wrappedView.getHostSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewFromDevice, startIndexUnchanged, lengthUnchanged, 2));
  }

  auto deviceSubviewUnchanged = wrappedView.getDeviceSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewUnchanged, startIndexUnchanged, lengthUnchanged, 2));

  auto deviceSubviewChanged = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewChanged, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceSubviewReadOnly_syncToDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());
  fixture.fillDualViewOnHost();

  const WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 4;
  int length = 3;
  auto deviceSubview = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubview, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceSubviewReadWrite_syncToDevice_modifyOnDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());
  fixture.fillDualViewOnHost();

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 5;
  int length = 2;
  {
    auto deviceSubview = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadWrite);
    TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubview, startIndex, length));
    fixture.multiplyOnDevice(deviceSubview, 2);
  }

  auto hostSubview = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadWrite);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubview, startIndex, length, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceSubviewOverwriteAll_syncToDevice_modifyOnDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostView);
    fixture.multiplyOnHost(hostView, 2);
  }

  int startIndex = 0;
  int length = 4;
  {
    auto deviceSubview = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceSubview, startIndex, length);
  }

  int startIndexUnchanged = length;
  int lengthUnchanged = fixture.getViewSize() - length;

  {
    auto deviceSubviewFromHost = wrappedView.getDeviceSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewFromHost, startIndexUnchanged, lengthUnchanged, 2));
  }

  auto hostSubviewUnchanged = wrappedView.getHostSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewUnchanged, startIndexUnchanged, lengthUnchanged, 2));

  auto hostSubviewChanged = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewChanged, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, hostSubviewThrowsIfNonOverlappingDeviceViewAlive_ReadOnly) {
  WrappedDualViewFixture fixture;

  const WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndexA = 0;
  int lengthA = 4;

  int startIndexB = 8;
  int lengthB = 4;

  auto deviceView = wrappedView.getDeviceSubview(startIndexA, lengthA, Tpetra::Access::ReadOnly);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getHostSubview(startIndexB, lengthB, Tpetra::Access::ReadOnly));
  }
  else {
    TEST_THROW(wrappedView.getHostSubview(startIndexB, lengthB, Tpetra::Access::ReadOnly), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, hostSubviewThrowsIfNonOverlappingDeviceViewAlive_ReadWrite) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndexA = 0;
  int lengthA = 4;

  int startIndexB = 8;
  int lengthB = 4;

  auto deviceView = wrappedView.getDeviceSubview(startIndexA, lengthA, Tpetra::Access::ReadWrite);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getHostSubview(startIndexB, lengthB, Tpetra::Access::ReadWrite));
  }
  else {
    TEST_THROW(wrappedView.getHostSubview(startIndexB, lengthB, Tpetra::Access::ReadWrite), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, hostSubviewThrowsIfNonOverlappingDeviceViewAlive_OverwriteAll) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndexA = 0;
  int lengthA = 4;

  int startIndexB = 8;
  int lengthB = 4;

  auto deviceView = wrappedView.getDeviceSubview(startIndexA, lengthA, Tpetra::Access::OverwriteAll);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getHostSubview(startIndexB, lengthB, Tpetra::Access::OverwriteAll));
  }
  else {
    TEST_THROW(wrappedView.getHostSubview(startIndexB, lengthB, Tpetra::Access::OverwriteAll), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, deviceSubviewThrowsIfNonOverlappingHostViewAlive_ReadOnly) {
  WrappedDualViewFixture fixture;

  const WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndexA = 0;
  int lengthA = 4;

  int startIndexB = 8;
  int lengthB = 4;

  auto hostView = wrappedView.getHostSubview(startIndexA, lengthA, Tpetra::Access::ReadOnly);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getDeviceSubview(startIndexB, lengthB, Tpetra::Access::ReadOnly));
  }
  else {
    TEST_THROW(wrappedView.getDeviceSubview(startIndexB, lengthB, Tpetra::Access::ReadOnly), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, deviceSubviewThrowsIfNonOverlappingHostViewAlive_ReadWrite) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndexA = 0;
  int lengthA = 4;

  int startIndexB = 8;
  int lengthB = 4;

  auto hostView = wrappedView.getHostSubview(startIndexA, lengthA, Tpetra::Access::ReadWrite);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getDeviceSubview(startIndexB, lengthB, Tpetra::Access::ReadWrite));
  }
  else {
    TEST_THROW(wrappedView.getDeviceSubview(startIndexB, lengthB, Tpetra::Access::ReadWrite), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, deviceSubviewThrowsIfNonOverlappingHostViewAlive_OverwriteAll) {
  WrappedDualViewFixture fixture;

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndexA = 0;
  int lengthA = 4;

  int startIndexB = 8;
  int lengthB = 4;

  auto hostView = wrappedView.getHostSubview(startIndexA, lengthA, Tpetra::Access::OverwriteAll);
  if (fixture.deviceMemoryIsHostAccessible) {
    TEST_NOTHROW(wrappedView.getDeviceSubview(startIndexB, lengthB, Tpetra::Access::OverwriteAll));
  }
  else {
    TEST_THROW(wrappedView.getDeviceSubview(startIndexB, lengthB, Tpetra::Access::OverwriteAll), std::runtime_error);
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, aliasedSubviewConstructor) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();
  const WrappedDualViewType wrappedView(fixture.getDualView());

  int startIndex = 4;
  int length = 6;
  const WrappedDualViewType wrappedSubview(wrappedView, startIndex, length);

  TEST_EQUALITY(wrappedSubview.extent(0), static_cast<size_t>(length));
  {
    auto hostSubview = wrappedSubview.getHostView(Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubview, startIndex, length));
  }
  {
    auto deviceSubview = wrappedSubview.getDeviceView(Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubview, startIndex, length));
  }
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostViewWrappedSubviewOverwriteAll_syncToHost_modifyOnHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceView);
    fixture.multiplyOnDevice(deviceView, 2);
  }

  int startIndex = 0;
  int length = 4;
  WrappedDualViewType wrappedSubview(wrappedView, startIndex, length);
  {
    auto hostSubview = wrappedSubview.getHostView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostSubview, startIndex, length);
  }

  int startIndexUnchanged = length;
  int lengthUnchanged = fixture.getViewSize() - length;

  {
    auto hostSubviewFromDevice = wrappedView.getHostSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewFromDevice, startIndexUnchanged, lengthUnchanged, 2));
  }

  auto deviceSubviewUnchanged = wrappedView.getDeviceSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewUnchanged, startIndexUnchanged, lengthUnchanged, 2));

  auto deviceSubviewChanged = wrappedSubview.getDeviceView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewChanged, startIndex, length));

  auto deviceSubviewChangedOriginal = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewChangedOriginal, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceViewWrappedSubviewOverwriteAll_syncToDevice_modifyOnDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostView);
    fixture.multiplyOnHost(hostView, 2);
  }

  int startIndex = 0;
  int length = 4;
  WrappedDualViewType wrappedSubview(wrappedView, startIndex, length);
  {
    auto deviceSubview = wrappedSubview.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceSubview, startIndex, length);
  }

  int startIndexUnchanged = length;
  int lengthUnchanged = fixture.getViewSize() - length;

  {
    auto deviceSubviewFromHost = wrappedView.getDeviceSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewFromHost, startIndexUnchanged, lengthUnchanged, 2));
  }

  auto hostSubviewUnchanged = wrappedView.getHostSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewUnchanged, startIndexUnchanged, lengthUnchanged, 2));

  auto hostSubviewChanged = wrappedSubview.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewChanged, startIndex, length));

  auto hostSubviewChangedOriginal = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewChangedOriginal, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostViewIntermediateWrappedSubviewOverwriteAll_syncToHost_modifyOnHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceView);
    fixture.multiplyOnDevice(deviceView, 2);
  }

  WrappedDualViewType intermediateWrappedSubview(wrappedView, 0, 8);

  int startIndex = 0;
  int length = 4;
  WrappedDualViewType wrappedSubview(intermediateWrappedSubview, startIndex, length);
  {
    auto hostSubview = wrappedSubview.getHostView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostSubview, startIndex, length);
  }

  int startIndexUnchanged = length;
  int lengthUnchanged = fixture.getViewSize() - length;

  {
    auto hostSubviewFromDevice = wrappedView.getHostSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewFromDevice, startIndexUnchanged, lengthUnchanged, 2));
  }

  auto deviceSubviewUnchanged = wrappedView.getDeviceSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewUnchanged, startIndexUnchanged, lengthUnchanged, 2));

  auto deviceSubviewChanged = wrappedSubview.getDeviceView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewChanged, startIndex, length));

  auto deviceSubviewChangedOriginal = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewChangedOriginal, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceViewIntermediateWrappedSubviewOverwriteAll_syncToDevice_modifyOnDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostView);
    fixture.multiplyOnHost(hostView, 2);
  }

  WrappedDualViewType intermediateWrappedSubview(wrappedView, 0, 8);

  int startIndex = 0;
  int length = 4;
  WrappedDualViewType wrappedSubview(intermediateWrappedSubview, startIndex, length);
  {
    auto deviceSubview = wrappedSubview.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceSubview, startIndex, length);
  }

  int startIndexUnchanged = length;
  int lengthUnchanged = fixture.getViewSize() - length;

  {
    auto deviceSubviewFromHost = wrappedView.getDeviceSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewFromHost, startIndexUnchanged, lengthUnchanged, 2));
  }

  auto hostSubviewUnchanged = wrappedView.getHostSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewUnchanged, startIndexUnchanged, lengthUnchanged, 2));

  auto hostSubviewChanged = wrappedSubview.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewChanged, startIndex, length));

  auto hostSubviewChangedOriginal = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewChangedOriginal, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostSubviewIntermediateWrappedSubviewOverwriteAll_syncToHost_modifyOnHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceView);
    fixture.multiplyOnDevice(deviceView, 2);
  }

  WrappedDualViewType intermediateWrappedSubview(wrappedView, 0, 8);
  int startIndex = 0;
  int length = 4;
  {
    auto hostSubview = intermediateWrappedSubview.getHostSubview(startIndex, length, Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostSubview, startIndex, length);
  }

  int startIndexUnchanged = length;
  int lengthUnchanged = fixture.getViewSize() - length;

  {
    auto hostSubviewFromDevice = wrappedView.getHostSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewFromDevice, startIndexUnchanged, lengthUnchanged, 2));
  }

  auto deviceSubviewUnchanged = wrappedView.getDeviceSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewUnchanged, startIndexUnchanged, lengthUnchanged, 2));

  auto deviceSubviewChanged = intermediateWrappedSubview.getDeviceSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewChanged, startIndex, length));

  auto deviceSubviewChangedOriginal = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewChangedOriginal, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceSubviewIntermediateWrappedSubviewOverwriteAll_syncToDevice_modifyOnDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostView);
    fixture.multiplyOnHost(hostView, 2);
  }

  WrappedDualViewType intermediateWrappedSubview(wrappedView, 0, 8);
  int startIndex = 0;
  int length = 4;
  {
    auto deviceSubview = intermediateWrappedSubview.getDeviceSubview(startIndex, length, Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceSubview, startIndex, length);
  }

  int startIndexUnchanged = length;
  int lengthUnchanged = fixture.getViewSize() - length;

  {
    auto deviceSubviewFromHost = wrappedView.getDeviceSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewFromHost, startIndexUnchanged, lengthUnchanged, 2));
  }

  auto hostSubviewUnchanged = wrappedView.getHostSubview(startIndexUnchanged, lengthUnchanged, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewUnchanged, startIndexUnchanged, lengthUnchanged, 2));

  auto hostSubviewChanged = intermediateWrappedSubview.getHostSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewChanged, startIndex, length));

  auto hostSubviewChangedOriginal = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewChangedOriginal, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceViewWrappedSubviewInMiddleOverwriteAll_syncToDevice_modifyOnDevice) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceView);
    fixture.multiplyOnDevice(deviceView, 2);
  }

  int startIndex = 3;
  int length = 3;
  WrappedDualViewType wrappedSubview(wrappedView, startIndex, length);
  {
    auto hostSubview = wrappedSubview.getHostView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostSubview, startIndex, length);
  }

  int startIndexFront = 0;
  int lengthFront = startIndex;

  int startIndexBack = startIndex + length;
  int lengthBack = fixture.getViewSize() - startIndexBack;

  {
    auto hostSubviewFromDeviceFront = wrappedView.getHostSubview(startIndexFront, lengthFront, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewFromDeviceFront, startIndexFront, lengthFront, 2));

    auto hostSubviewFromDeviceBack = wrappedView.getHostSubview(startIndexBack, lengthBack, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewFromDeviceBack, startIndexBack, lengthBack, 2));
  }

  auto deviceSubviewFront = wrappedView.getDeviceSubview(startIndexFront, lengthFront, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewFront, startIndexFront, lengthFront, 2));

  auto deviceSubviewBack = wrappedView.getDeviceSubview(startIndexBack, lengthBack, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewBack, startIndexBack, lengthBack, 2));

  auto deviceSubviewChanged = wrappedSubview.getDeviceView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewChanged, startIndex, length));

  auto deviceSubviewChangedOriginal = wrappedView.getDeviceSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewChangedOriginal, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostViewWrappedSubviewInMiddleOverwriteAll_syncToHost_modifyOnHost) {
  WrappedDualViewFixture fixture;
  TEST_ASSERT(fixture.valuesInitializedToZero());

  WrappedDualViewType wrappedView(fixture.getDualView());

  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnHost(hostView);
    fixture.multiplyOnHost(hostView, 2);
  }

  int startIndex = 3;
  int length = 3;
  WrappedDualViewType wrappedSubview(wrappedView, startIndex, length);
  {
    auto deviceSubview = wrappedSubview.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.fillViewOnDevice(deviceSubview, startIndex, length);
  }

  int startIndexFront = 0;
  int lengthFront = startIndex;

  int startIndexBack = startIndex + length;
  int lengthBack = fixture.getViewSize() - startIndexBack;

  {
    auto deviceSubviewFromDeviceFront = wrappedView.getDeviceSubview(startIndexFront, lengthFront, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewFromDeviceFront, startIndexFront, lengthFront, 2));

    auto deviceSubviewFromDeviceBack = wrappedView.getDeviceSubview(startIndexBack, lengthBack, Tpetra::Access::ReadOnly);
    TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewFromDeviceBack, startIndexBack, lengthBack, 2));
  }

  auto hostSubviewFront = wrappedView.getHostSubview(startIndexFront, lengthFront, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewFront, startIndexFront, lengthFront, 2));

  auto hostSubviewBack = wrappedView.getHostSubview(startIndexBack, lengthBack, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewBack, startIndexBack, lengthBack, 2));

  auto hostSubviewChanged = wrappedSubview.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewChanged, startIndex, length));

  auto hostSubviewChangedOriginal = wrappedView.getHostSubview(startIndex, length, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewChangedOriginal, startIndex, length));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostTwoSubviews_ReadOnly) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();

  WrappedDualViewType wrappedView(fixture.getDualView());
  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.multiplyOnDevice(deviceView, 2);
  }

  int startFirstHalf = 0;
  int lengthHalf = fixture.getViewSize()/2;
  int startSecondHalf = lengthHalf;

  const WrappedDualViewType wrappedSubview(wrappedView, startFirstHalf, lengthHalf);
  auto hostSubviewFirstHalf = wrappedSubview.getHostView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewFirstHalf, startFirstHalf, lengthHalf, 2));

  auto hostSubviewSecondHalf = wrappedView.getHostSubview(startSecondHalf, lengthHalf, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewSecondHalf, startSecondHalf, lengthHalf, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceTwoSubviews_ReadOnly) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();

  WrappedDualViewType wrappedView(fixture.getDualView());
  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
    fixture.multiplyOnHost(hostView, 2);
  }

  int startFirstHalf = 0;
  int lengthHalf = fixture.getViewSize()/2;
  int startSecondHalf = lengthHalf;

  const WrappedDualViewType wrappedSubview(wrappedView, startFirstHalf, lengthHalf);
  auto deviceSubviewFirstHalf = wrappedSubview.getDeviceView(Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewFirstHalf, startFirstHalf, lengthHalf, 2));

  auto deviceSubviewSecondHalf = wrappedView.getDeviceSubview(startSecondHalf, lengthHalf, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewSecondHalf, startSecondHalf, lengthHalf, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostSubviewOfSubviewAndSubview_ReadOnly) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();

  int startFirstHalf = 0;
  int lengthHalf = fixture.getViewSize()/2;
  int startSecondHalf = lengthHalf;
  int lengthQuarter = lengthHalf/2;

  WrappedDualViewType wrappedView(fixture.getDualView());
  WrappedDualViewType wrappedSubview(wrappedView, startFirstHalf, lengthHalf);
  {
    auto deviceView = wrappedView.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.multiplyOnDevice(deviceView, 2);
  }

  auto hostSubviewFirstHalf = wrappedSubview.getHostSubview(startFirstHalf, lengthQuarter, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewFirstHalf, startFirstHalf, lengthQuarter, 2));

  auto hostSubviewSecondHalf = wrappedView.getHostSubview(startSecondHalf, lengthHalf, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewSecondHalf, startSecondHalf, lengthHalf, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceSubviewOfSubviewAndSubview_ReadOnly) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();

  int startFirstHalf = 0;
  int lengthHalf = fixture.getViewSize()/2;
  int startSecondHalf = lengthHalf;
  int lengthQuarter = lengthHalf/2;

  WrappedDualViewType wrappedView(fixture.getDualView());
  WrappedDualViewType wrappedSubview(wrappedView, startFirstHalf, lengthHalf);
  {
    auto hostView = wrappedView.getHostView(Tpetra::Access::OverwriteAll);
    fixture.multiplyOnHost(hostView, 2);
  }

  auto deviceSubviewFirstHalf = wrappedSubview.getDeviceSubview(startFirstHalf, lengthQuarter, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewFirstHalf, startFirstHalf, lengthQuarter, 2));

  auto deviceSubviewSecondHalf = wrappedView.getDeviceSubview(startSecondHalf, lengthHalf, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewSecondHalf, startSecondHalf, lengthHalf, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessHostTwoSubviewsOfSubview_ReadOnly) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startSubview = 4;
  int lengthSubview = 8;
  WrappedDualViewType wrappedSubview(wrappedView, startSubview, lengthSubview);
  {
    auto deviceView = wrappedSubview.getDeviceView(Tpetra::Access::OverwriteAll);
    fixture.multiplyOnDevice(deviceView, 2);
    TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceView, startSubview, lengthSubview, 2));
  }

  int startFirstHalf = 0;
  int lengthHalf = lengthSubview/2;
  int startSecondHalf = lengthHalf;

  auto hostSubviewFirstHalf = wrappedSubview.getHostSubview(startFirstHalf, lengthHalf, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewFirstHalf, startFirstHalf+startSubview, lengthHalf, 2));

  auto hostSubviewSecondHalf = wrappedSubview.getHostSubview(startSecondHalf, lengthHalf, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnHost(hostSubviewSecondHalf, startSecondHalf+startSubview, lengthHalf, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, accessDeviceTwoSubviewsOfSubview_ReadOnly) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();

  WrappedDualViewType wrappedView(fixture.getDualView());

  int startSubview = 4;
  int lengthSubview = 8;
  WrappedDualViewType wrappedSubview(wrappedView, startSubview, lengthSubview);
  {
    auto hostView = wrappedSubview.getHostView(Tpetra::Access::OverwriteAll);
    fixture.multiplyOnHost(hostView, 2);
    TEST_ASSERT(fixture.valuesCorrectOnHost(hostView, startSubview, lengthSubview, 2));
  }

  int startFirstHalf = 0;
  int lengthHalf = lengthSubview/2;
  int startSecondHalf = lengthHalf;

  auto deviceSubviewFirstHalf = wrappedSubview.getDeviceSubview(startFirstHalf, lengthHalf, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewFirstHalf, startFirstHalf+startSubview, lengthHalf, 2));

  auto deviceSubviewSecondHalf = wrappedSubview.getDeviceSubview(startSecondHalf, lengthHalf, Tpetra::Access::ReadOnly);
  TEST_ASSERT(fixture.valuesCorrectOnDevice(deviceSubviewSecondHalf, startSecondHalf+startSubview, lengthHalf, 2));
}

TEUCHOS_UNIT_TEST(WrappedDualView, attemptConstructUnmanaged) {
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();
  WrappedDualViewType wrappedView(fixture.getDualView());
  auto owningView = wrappedView.getDeviceView(Tpetra::Access::ReadWrite);
  static_assert(WrappedDualViewType::t_dev::rank == 1,
      "This test requires WrappedDualViewType to be rank 1. If this breaks, use a custom type here.");

  //Although this view doesn't have Unmanaged memory traits in its type,
  //it behaves as if it did (does not do reference counting), and has use_count() == 0
  typename WrappedDualViewType::t_dev unmanagedView(owningView.data(), owningView.extent(0));
  //This should throw - WrappedDualView must be able to take ownership
  //of the device memory from user, but the user's view does not own it
  try
  {
    WrappedDualViewType cannotConstructThis(unmanagedView);
    TEST_ASSERT(false);
  }
  catch(std::exception&)
  {}
}

}

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard scopeGuard(&argc, &argv);
  const int errCode = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  return errCode;
}
