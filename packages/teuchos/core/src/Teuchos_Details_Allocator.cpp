// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
// @HEADER

#include <Teuchos_Details_Allocator.hpp>
#include <sstream>

namespace Teuchos {
namespace Details {

void
AllocationLogger::
logAllocation (std::ostream& out,
               const size_type numEntries,
               const size_type numBytes,
               const char typeName[],
               const bool verbose)
{
  using std::endl;
  curAllocInBytes_ += numBytes;
  if (curAllocInBytes_ > maxAllocInBytes_) {
    maxAllocInBytes_ = curAllocInBytes_;
  }

  if (verbose) {
    // Identify this as a Teuchos allocation.
    out << "Teuchos,alloc," << numEntries << "," << typeName << "," << numBytes << endl;
  }
}

void
AllocationLogger::
logDeallocation (std::ostream& out,
                 const size_type numEntries,
                 const size_type numBytes,
                 const char typeName[],
                 const bool verbose)
{
  using std::endl;
  curAllocInBytes_ -= numBytes;

  if (verbose) {
    // First field identifies this as a Teuchos allocation.  Use the
    // same number of characters for "deall"(ocation) as
    // "alloc"(ation) above, so that the columns line up nicely for
    // human reading.  Print deallocations as negative allocations.
    // This makes it easy for a sed or awk script to compute totals.
    out << "Teuchos,deall,-" << numEntries << "," << typeName
        << ",-" << numBytes << endl;
  }
}

AllocationLogger::size_type
AllocationLogger::curAllocInBytes () { return curAllocInBytes_; }

AllocationLogger::size_type
AllocationLogger::maxAllocInBytes () { return maxAllocInBytes_; }

void
AllocationLogger::resetAllocationCounts ()
{
  curAllocInBytes_ = 0;
  maxAllocInBytes_ = 0;
}

AllocationLogger::size_type AllocationLogger::curAllocInBytes_ = 0;
AllocationLogger::size_type AllocationLogger::maxAllocInBytes_ = 0;

} // namespace Details
} // namespace Teuchos
