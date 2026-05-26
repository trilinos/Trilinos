// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_BEHAVIOR_HPP
#define MUELU_BEHAVIOR_HPP

#include <stddef.h>
#include <string>

namespace MueLu {
class Behavior {
 public:
  /// \brief Whether MueLu is in debug mode.
  ///
  /// "Debug mode" means that MueLu does extra error checks that may
  /// require more MPI communication or local computation.  It may
  /// also produce more detailed error messages, and more copious
  /// debug output.
  static bool debug();
};

}  // namespace MueLu

#endif
