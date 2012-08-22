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

#include <Tsqr_TimeStats.hpp>
#include <limits>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  TimeStats::TimeStats() { init(); }

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

    if (humanReadable)
      {
	const char prefix[] = "-- ";
	out << label << ":" << endl;
	if (count() == 0)
	  out << prefix << "No values collected" << endl;
	else if (count() == 1)
	  out << prefix << "One value collected: " << min() << endl;
	else
	  {
	    out << prefix << "Count: " << count() << endl
		<< prefix << "Min:   " << min() << endl
		<< prefix << "Mean:  " << mean() << endl
		<< prefix << "Max:   " << max() << endl
		<< prefix << "Total: " << total() << endl;
	  }
      }
    else
      {
	// "%" identifies this line as a "comment" line to filter out.
	// First print field labels on one line, then print field
	// values on the next line.
	if (printHeaders)
	  {
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
