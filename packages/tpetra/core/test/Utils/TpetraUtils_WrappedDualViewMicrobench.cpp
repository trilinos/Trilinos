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

TEUCHOS_UNIT_TEST(WrappedDualView, hostViewMicrobench) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;
  using Teuchos::StackedTimer;
  using std::cout;
  using std::endl;
  /*
  Write a loop over many multiple calls to getLocalViewHost(Access::...)
  Make sure the sync does not happen; want to measure only the overhead
     Compare to just grabbing the host pointer from a DualView

       for each access tag {
         auto tmp = getLocalViewHost(Access::ReadOnly);
         start timer
         for many many iterations
           auto blah = getLocalViewHost(Access::...);
         stop timer
       }
  */
  WrappedDualViewFixture fixture;
  fixture.fillDualViewOnHostDevice();
  RCP<StackedTimer> timer = rcp(new StackedTimer("hostView"));
  TimeMonitor::setStackedTimer(timer);
  size_t iterations = 4000000;
  std::string dvTimer = "hostView: raw dualView";
  std::string roTimer = "hostView: ReadOnly";
  std::string oaTimer = "hostView: OverwriteAll";
  std::string rwTimer = "hostView: ReadWrite";
  
  //get communicator
  RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  auto dualV = fixture.getDualView();
  WrappedDualViewType wdv(fixture.getDualView());

  //raw dualView
  {
    auto warmUp = dualV.view_host();
    TimeMonitor t(*TimeMonitor::getNewTimer(dvTimer));
    for (size_t i = 0; i < iterations; i++) {
      auto tmp = dualV.view_host();
    }
  }

  //ReadOnly
  {
    auto warmUp = wdv.getHostView(Tpetra::Access::ReadOnly);
    TimeMonitor t(*TimeMonitor::getNewTimer(roTimer));
    for (size_t i = 0; i < iterations; i++) {
      auto tmp = wdv.getHostView(Tpetra::Access::ReadOnly);
    }
  }

  //OverwriteAll
  {
    auto warmUp = wdv.getHostView(Tpetra::Access::OverwriteAll);
    TimeMonitor t(*TimeMonitor::getNewTimer(oaTimer));
    for (size_t i = 0; i < iterations; i++) {
      auto tmp = wdv.getHostView(Tpetra::Access::OverwriteAll);
    }
  }

  //ReadWrite
  {
    auto warmUp = wdv.getHostView(Tpetra::Access::ReadWrite);
    TimeMonitor t(*TimeMonitor::getNewTimer(rwTimer));
    for (size_t i = 0; i < iterations; i++) {
      auto tmp = wdv.getHostView(Tpetra::Access::ReadWrite);
    }
  }
  StackedTimer::OutputOptions options;
  options.print_warnings = false;
  timer->report(std::cout, comm, options);
}

}

int main(int argc, char* argv[]) {
  Tpetra::ScopeGuard scopeGuard(&argc, &argv);
  const int errCode = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  return errCode;
}
