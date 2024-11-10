// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Tsqr_TimeStats.hpp"
#include <limits>

namespace TSQR {

  TimeStats::TimeStats () { init (); }

  void
  TimeStats::init () {
    min_ = std::numeric_limits< double >::infinity();
    max_ = -std::numeric_limits< double >::infinity();
    mean_ = double (0);
    total_ = double (0);
    count_ = int (0);
  }

  void
  TimeStats::update (const double curTime) {
    total_ += curTime;
    count_++;

    if (curTime < min_)
      min_ = curTime;
    if (curTime > max_)
      max_ = curTime;

    // Mean(1:n) = ((n-1) / n) * Mean(1:n-1) + x(n) / n.
    //
    // Casting int to double is exact.
    const double scale = double(count_ - 1) / double(count_);
    mean_ = scale * mean_ + curTime / double(count_);
  }

  void
  TimeStats::print (std::ostream& out,
                    const bool humanReadable,
                    const std::string& label,
                    const std::string& labelLabel,
                    const bool printHeaders) const
  {
    using std::endl;

    if (humanReadable) {
      const char prefix[] = "-- ";
      out << label << ":" << endl;
      if (count() == 0) {
        out << prefix << "No values collected" << endl;
      }
      else if (count() == 1) {
        out << prefix << "One value collected: " << min() << endl;
      }
      else {
        out << prefix << "Count: " << count() << endl
            << prefix << "Min:   " << min() << endl
            << prefix << "Mean:  " << mean() << endl
            << prefix << "Max:   " << max() << endl
            << prefix << "Total: " << total() << endl;
      }
    }
    else {
      // "%" identifies this line as a "comment" line to filter out.
      // First print field labels on one line, then print field
      // values on the next line.
      if (printHeaders) {
        out << "%" << labelLabel
            << "," << "count"
            << "," << "min"
            << "," << "mean"
            << "," << "max"
            << "," << "total"
            << endl;
      }
      out << label
          << "," << count()
          << "," << min()
          << "," << mean()
          << "," << max()
          << "," << total()
          << endl;
    }
  }

} // namespace TSQR
