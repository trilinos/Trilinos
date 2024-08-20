// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_TimeStats_hpp
#define __TSQR_TimeStats_hpp

#include <iostream>
#include <string>

namespace TSQR {

  /// \brief Collect running timing statistics.
  ///
  /// TimeStats collects running statistics on a particular timing,
  /// which is collected count() times.  When you get a new timing,
  /// call update().
  class TimeStats {
  public:
    //! Default constructor: resets statistics to initial values.
    TimeStats ();

    //! Reset the statistics.
    void init ();

    //! Add a new data point and update the running statistics.
    void update (const double curTime);

    /// \brief Print to the given output stream \c out.
    ///
    /// \param out [in/out] Output stream to which to print
    /// \param humanReadable [in] Whether to print in a format easy
    ///   for humans to read, or easy for automatic parsing
    /// \param label [in] If not humanReadable, then print this string
    ///   as a row identifier at the beginning of the row
    /// \param labelLabel [in] If not humanReadable, then use this
    ///   as the column header for the "label" (first) column
    /// \param printHeaders [in] If not humanReadable, then print
    ///   column headers, preceded by a "%" so that the parser will
    ///   ignore the line
    void
    print (std::ostream& out,
           const bool humanReadable,
           const std::string& label,
           const std::string& labelLabel,
           const bool printHeaders) const;

    //! Min value seen thus far (+Inf if no data have been collected).
    double min () const { return min_; }

    //! Max value seen thus far (-Inf if no data have been collected).
    double max () const { return max_; }

    //! Arithmetic mean thus far (0 if no data has been collected)
    double mean () const { return mean_; }

    //! Total thus far (0 if no data has been collected)
    double total () const { return total_; }

    //! Count of data points collected thus far
    int count () const { return count_; }

    /// \brief Construct a TimeStats object from its constituent data.
    ///
    /// This is useful for computing global statistics over many MPI
    /// processes, or for otherwise combining different TimeStats
    /// objects.
    ///
    /// \note This design is suboptimal, because it makes it hard for
    ///   new statistics to be added to the class.
    TimeStats (const int newCount,
               const double newMin,
               const double newMax,
               const double newMean,
               const double newTotal) :
      min_ (newMin),
      max_ (newMax),
      mean_ (newMean),
      total_ (newTotal),
      count_ (newCount)
    {}

  private:
    double min_, max_, mean_, total_;
    int count_;
  };



} // namespace TSQR

#endif // __TSQR_TimeStats_hpp
