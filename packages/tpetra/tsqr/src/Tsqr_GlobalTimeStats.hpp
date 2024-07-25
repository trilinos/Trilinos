// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_GlobalTimeStats_hpp
#define __TSQR_GlobalTimeStats_hpp

#include "Tsqr_TimeStats.hpp"

namespace TSQR {

  // Forward declaration
  template< class Scalar >
  class MessengerBase;

  /// Produce global time statistics out of all the local ones.
  ///
  /// \param comm [in] Encapsulation of the interprocess communicator
  /// \param localStats [in] Local (to this process) time statistics
  ///
  /// \return Global (over all processes) time statistics
  TimeStats
  globalTimeStats (MessengerBase<double>& comm,
                   const TimeStats& localStats);

} // namespace TSQR

#endif // __TSQR_GlobalTimeStats_hpp
