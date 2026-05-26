// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_EnvVariables.hpp"
#include "MueLu_Behavior.hpp"
#include <array>
#include <map>

// environ should be available on posix platforms
#if not(defined(WIN) && (_MSC_VER >= 1900))
// needs to be in the global namespace
extern char **environ;
#endif

namespace MueLu {

namespace BehaviorDetails {

constexpr const std::string_view DEBUG = "MUELU_DEBUG";

}  // namespace BehaviorDetails

namespace {  // (anonymous)

constexpr bool debugDefault() {
#ifdef HAVE_MUELU_DEBUG
  return true;
#else
  return false;
#endif  // HAVE_MUELU_DEBUG
}

}  // namespace

bool Behavior::debug() {
  constexpr bool defaultValue = debugDefault();

  static bool value_       = defaultValue;
  static bool initialized_ = false;
  return Teuchos::idempotentlyGetEnvironmentVariable(
      value_, initialized_, BehaviorDetails::DEBUG, defaultValue);
}

}  // namespace MueLu
