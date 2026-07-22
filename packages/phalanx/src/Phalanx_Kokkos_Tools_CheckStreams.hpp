// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PHALANX_KOKKOS_TOOLS_VERIFYSTREAM_HPP
#define PHALANX_KOKKOS_TOOLS_VERIFYSTREAM_HPP

namespace PHX {
  /// Function that sets kokkos-tools callbacks to make sure that the
  /// default stream is not being using in any kokkos calls. Checks
  /// are only active for cuda and hip backends as the id for the
  /// serial backend is always the same.
  void set_enforce_no_default_stream_use();

  /// Function that unsets kokkos-tools callbacks.
  void unset_enforce_no_default_stream_use();
}

#endif
