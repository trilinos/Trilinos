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

#include <Tsqr_GlobalTimeStats.hpp>
#include <Tsqr_MessengerBase.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  TimeStats
  globalTimeStats (const Teuchos::RCP< MessengerBase< double > >& comm,
		   const TimeStats& localStats)
  {
    TimeStats globalStats;

    // Casting int to double is exact.
    const double localCount = static_cast<double> (localStats.count());

    // The counts on all MPI processes should be equal if we are going
    // to try to compare them.  We test for this by computing the min
    // and max counts (which should be equal) over all processes.
    const double minCount = comm->globalMin (localCount);
    const double maxCount = comm->globalMax (localCount);
    if (minCount != maxCount)
      throw std::logic_error ("Global stats don\'t make "
			      "sense, because counts differ");
    // minCount == maxCount, so we can use either one of them.  The
    // cast back from double is exact, because the double originally
    // came from an int.
    const int newCount = static_cast<int> (minCount);
      
    // Casting int to double is exact.
    const double P = static_cast<double> (comm->size());
    const double newMin = comm->globalMin (localStats.min());
    const double newMax = comm->globalMax (localStats.max());
    const double newMean = comm->globalSum (localStats.mean() / P);

    // Note that this is not the sum of the totals of all the
    // processes, but rather the "global total."  I've chosen to
    // define that as the max of the totals.
    const double newTotal = comm->globalMax (localStats.total());

    return TimeStats (newCount, newMin, newMax, newMean, newTotal);
  }


} // namespace TSQR
