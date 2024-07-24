// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_GlobalTimeStats.hpp"
#include "Tsqr_MessengerBase.hpp"
#include <stdexcept>

namespace TSQR {

  TimeStats
  globalTimeStats (MessengerBase<double>& comm,
                   const TimeStats& localStats)
  {
    TimeStats globalStats;

    // Casting int to double is exact.
    const double localCount = static_cast<double> (localStats.count ());

    // The counts on all MPI processes should be equal if we are going
    // to try to compare them.  We test for this by computing the min
    // and max counts (which should be equal) over all processes.
    const double minCount = comm.globalMin (localCount);
    const double maxCount = comm.globalMax (localCount);
    if (minCount != maxCount) {
      throw std::logic_error ("Global stats don\'t make "
                              "sense, because counts differ");
    }
    // minCount == maxCount, so we can use either one of them.  The
    // cast back from double is exact, because the double originally
    // came from an int.
    const int newCount = static_cast<int> (minCount);

    // Casting int to double is exact.
    const double P = static_cast<double> (comm.size ());
    const double newMin = comm.globalMin (localStats.min ());
    const double newMax = comm.globalMax (localStats.max ());
    const double newMean = comm.globalSum (localStats.mean () / P);

    // Note that this is not the sum of the totals of all the
    // processes, but rather the "global total."  I've chosen to
    // define that as the max of the totals.
    const double newTotal = comm.globalMax (localStats.total ());

    return TimeStats (newCount, newMin, newMax, newMean, newTotal);
  }
} // namespace TSQR
