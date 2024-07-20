// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
