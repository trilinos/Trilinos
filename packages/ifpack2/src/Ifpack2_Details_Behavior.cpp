// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_EnvVariables.hpp"
#include "Ifpack2_Details_Behavior.hpp"

namespace Ifpack2 {
namespace Details {

namespace {

constexpr const char DEBUG[] = "IFPACK2_DEBUG";

constexpr bool debugDefault() {
#ifdef HAVE_IFPACK2_DEBUG
  return true;
#else
  return false;
#endif  // HAVE_IFPACK2_DEBUG
}

}  // namespace

bool Behavior::debug() {
  constexpr bool defaultValue = debugDefault();

  static bool value_       = defaultValue;
  static bool initialized_ = false;
  return Teuchos::idempotentlyGetEnvironmentVariable(
      value_, initialized_, DEBUG, defaultValue);
}

}  // namespace Details
}  // namespace Ifpack2