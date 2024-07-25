// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_StatTimeMonitor.hpp"

namespace TSQR {

  StatTimeMonitor::StatTimeMonitor (Teuchos::Time& timer, TimeStats& stats)
    : timer_ (timer), stats_ (stats)
  {
    if (timer_.isRunning())
      timer_.stop(); // Just for sanity
    if (timer_.numCalls() == 0)
      stats_.init();
    timer_.start (true);
  }

  StatTimeMonitor::~StatTimeMonitor ()
  {
    const double curTime = timer_.stop();
    stats_.update (curTime);
  }

} // namespace TSQR
