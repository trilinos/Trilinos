// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_DETAILS_COPYANDPERMUTEBLOCKMULTIVECTOR_DEF_HPP
#define TPETRA_DETAILS_COPYANDPERMUTEBLOCKMULTIVECTOR_DEF_HPP

#include "Tpetra_BlockMultiVector.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include <iostream>
#include <sstream>

/// \file Tpetra_Details_packBlockMultiVector_def.hpp
/// \brief Definition of function that packs the entries of a
///   Tpetra::BlockMultiVector for communication.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {
namespace Details {

template<class SC, class LO, class GO, class NT>
std::pair<int, std::unique_ptr<std::string>>
copyAndPermuteBlockMultiVector
  (BlockMultiVector<SC, LO, GO, NT>& tgt,
   const BlockMultiVector<SC, LO, GO, NT>& src,
   const size_t numSameIDs,
   const Kokkos::DualView<
     const LO*,
     typename BlockMultiVector<SC, LO, GO, NT>::buffer_device_type
   >& permuteToLIDs,
   const Kokkos::DualView<
     const LO*,
     typename BlockMultiVector<SC, LO, GO, NT>::buffer_device_type
   >& permuteFromLIDs)
{
  using std::endl;
  const char funcName[] = "Tpetra::Details::copyAndPermuteBlockMultiVector";
  std::unique_ptr<std::string> errMsg;
  const bool debug = ::Tpetra::Details::Behavior::debug ("BlockMultiVector");
  const int myRank = src.getMap ()->getComm ()->getRank ();

  std::unique_ptr<std::string> prefix;
  if (debug) {
    std::ostringstream os;
    os << "Proc " << myRank << ": " << funcName << ": ";
    prefix = std::unique_ptr<std::string> (new std::string (os.str ()));
    os << "Start" << endl;
    std::cerr << os.str ();
  }
  if (permuteToLIDs.need_sync_host () || permuteFromLIDs.need_sync_host ()) {
    std::ostringstream os;
    os << "Proc " << myRank << ": " << funcName << ": permuteToLIDs and/or "
      "permuteFromLIDs need sync to host.  This should never happen.";
    errMsg = std::unique_ptr<std::string> (new std::string (os.str ()));
    return { -1, std::move (errMsg) };
  }

  const size_t numPermuteLIDs = static_cast<size_t> (permuteToLIDs.extent (0));
  if (static_cast<size_t> (permuteFromLIDs.extent (0)) != numPermuteLIDs) {
    std::ostringstream os;
    os << "Proc " << myRank << ": " << funcName << ": "
      "permuteToLIDs.extent(0) = " << numPermuteLIDs << " != "
      "permuteFromLIDs.extent(0) = " << permuteFromLIDs.extent (0) << ".";
    errMsg = std::unique_ptr<std::string> (new std::string (os.str ()));
    return { -2, std::move (errMsg) };
  }

  const size_t numVecs = static_cast<size_t> (src.getNumVectors ());
  if (debug) {
    std::ostringstream os;
    os << *prefix << "Sames: numSameIDs=" << numSameIDs
       << ", numVecs=" << numVecs << endl;
    std::cerr << os.str ();
  }
  // FIXME (mfh 23 Apr 2014) This implementation is sequential and
  // assumes UVM.
  for (size_t j = 0; j < numVecs; ++j) {
    for (size_t lclRow = 0; lclRow < numSameIDs; ++lclRow) {
      deep_copy (tgt.getLocalBlock (lclRow, j),
                 src.getLocalBlock (lclRow, j));
    }
  }

  if (debug) {
    std::ostringstream os;
    os << *prefix << "Permutes: numPermuteLIDs=" << numPermuteLIDs << endl;
    std::cerr << os.str ();
  }
  auto permuteToLIDs_h = permuteToLIDs.view_host ();
  auto permuteFromLIDs_h = permuteFromLIDs.view_host ();

  // FIXME (mfh 20 June 2012) For an Export with duplicate GIDs on the
  // same process, this merges their values by replacement of the last
  // encountered GID, not by the specified merge rule (such as ADD).
  for (size_t j = 0; j < numVecs; ++j) {
    for (size_t k = numSameIDs; k < numPermuteLIDs; ++k) {
      deep_copy (tgt.getLocalBlock (permuteToLIDs_h[k], j),
                 src.getLocalBlock (permuteFromLIDs_h[k], j));
    }
  }

  if (debug) {
    std::ostringstream os;
    os << *prefix << "Done" << endl;
    std::cerr << os.str ();
  }
  return { 0, std::move (errMsg) };
}

} // namespace Details
} // namespace Tpetra

#define TPETRA_DETAILS_COPYANDPERMUTEBLOCKMULTIVECTOR_INSTANT( SC, LO, GO, NT )  \
namespace Details { \
  template std::pair<int, std::unique_ptr<std::string>> \
  copyAndPermuteBlockMultiVector<SC, LO, GO, NT> ( \
    BlockMultiVector<SC, LO, GO, NT>&, \
    const BlockMultiVector<SC, LO, GO, NT>&, \
    const size_t, \
    const Kokkos::DualView< \
      const LO*, \
      BlockMultiVector<SC, LO, GO, NT>::buffer_device_type \
    >&, \
    const Kokkos::DualView< \
      const LO*, \
      BlockMultiVector<SC, LO, GO, NT>::buffer_device_type \
    >& ); \
}

#endif // TPETRA_DETAILS_COPYANDPERMUTEBLOCKMULTIVECTOR_DEF_HPP
