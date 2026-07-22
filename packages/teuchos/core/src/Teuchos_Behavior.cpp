// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_BEHAVIOR_HPP
#define TEUCHOS_BEHAVIOR_HPP

#include "Teuchos_Behavior.hpp"
#include "Teuchos_EnvVariables.hpp"
#include "Teuchos_TestForException.hpp"

namespace Teuchos {

constexpr const std::string_view FENCE_TIMERS = "TEUCHOS_FENCE_TIMERS";

constexpr bool fenceTimersDefault() {
#ifdef HAVE_TEUCHOS_TIMER_KOKKOS_FENCE
  return true;
#else
  return false;
#endif // HAVE_TEUCHOS_TIMER_KOKKOS_FENCE
}


bool Behavior::fenceTimers() {
  constexpr bool defaultValue = Teuchos::fenceTimersDefault();

  static bool value_ = defaultValue;
  static bool initialized_ = false;
  return Teuchos::idempotentlyGetEnvironmentVariable(value_, initialized_, Teuchos::FENCE_TIMERS,
                                                     defaultValue);
}

} // namespace Teuchos

#endif
