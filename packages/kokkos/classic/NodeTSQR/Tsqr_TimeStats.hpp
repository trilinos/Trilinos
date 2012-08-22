/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#ifndef __TSQR_TimeStats_hpp
#define __TSQR_TimeStats_hpp

#include <Teuchos_RCP.hpp>
#include <ostream>
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \brief Collect running statistics
  ///
  /// TimeStats collects running statistics on a particular timing,
  /// which is collected count() times.  When you get a new timing,
  /// call update().
  class TimeStats {
  public:
    TimeStats();
    ~TimeStats() {}

    /// Reset the statistics
    ///
    void init ();

    /// Add a new data point and update the running statistics
    ///
    void update (const double curTime);

    /// Print to out
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

    /// Minimum value seen thus far (+Inf if no data has been
    /// collected)
    double min() const { return min_; }

    /// Maximum value seen thus far (-Inf if no data has been
    /// collected)
    double max() const { return max_; }

    /// Arithmetic mean thus far (0 if no data has been collected)
    ///
    double mean() const { return mean_; }

    /// Total thus far (0 if no data has been collected)
    ///
    double total() const { return total_; }

    /// Count of data points collected thus far
    ///
    int count() const { return count_; }

    /// Construct a TimeStats object from its constituent data.  This
    /// is useful for computing global statistics over many MPI
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
