// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_TrivialTimer.hpp"
#include "Tsqr_verifyTimerConcept.hpp"

namespace TSQR {

TrivialTimer::TrivialTimer(const std::string& theName, bool doStart)
  : name_(theName)
  , counter_(0)
  , isRunning_(false) {
  if (doStart)
    start();
}

void TrivialTimer::verifyConcept() {
  TSQR::Test::verifyTimerConcept<TrivialTimer>();
}

void TrivialTimer::start(bool reset) {
  isRunning_ = true;
}

double
TrivialTimer::stop() {
  isRunning_ = false;
  ++counter_;
  return static_cast<double>(counter_);
}

}  // namespace TSQR
