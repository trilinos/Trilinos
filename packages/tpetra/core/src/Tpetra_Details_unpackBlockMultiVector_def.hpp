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

#ifndef TPETRA_DETAILS_UNPACKBLOCKMULTIVECTOR_DEF_HPP
#define TPETRA_DETAILS_UNPACKBLOCKMULTIVECTOR_DEF_HPP

#include "Tpetra_BlockMultiVector.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Experimental_BlockView.hpp"
#include "Kokkos_ArithTraits.hpp"
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
unpackBlockMultiVector (BlockMultiVector<SC, LO, GO, NT>& tgt,
                        const Kokkos::DualView<
                          const LO*,
                          typename BlockMultiVector<SC, LO, GO, NT>::buffer_device_type
                        >& importLIDs,
                        Kokkos::DualView<
                          typename BlockMultiVector<SC, LO, GO, NT>::impl_scalar_type*,
                          typename BlockMultiVector<SC, LO, GO, NT>::buffer_device_type
                        > imports,
                        Kokkos::DualView<
                          size_t*,
                          typename BlockMultiVector<SC, LO, GO, NT>::buffer_device_type
                        > numPacketsPerLID,
                        const size_t constantNumPackets,
                        Distributor& distor,
                        const CombineMode combineMode)
{
  using std::endl;
  using BMV = BlockMultiVector<SC, LO, GO, NT>;
  using IST = typename BMV::impl_scalar_type;
  using KAT = Kokkos::ArithTraits<IST>;
  using packet_type = IST;
  using buffer_device_type = typename BMV::buffer_device_type;
  using imports_dv_type = Kokkos::DualView<packet_type*, buffer_device_type>;
  using host_device_type = typename imports_dv_type::t_host::device_type;
  using src_little_vec_type = Kokkos::View<const IST*, host_device_type,
    Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
  const char funcName[] = "Tpetra::Details::unpackBlockMultiVector";
  std::unique_ptr<std::string> errMsg;
  const bool debug = ::Tpetra::Details::Behavior::debug ("BlockMultiVector");

  std::unique_ptr<std::string> prefix;
  if (debug) {
    const int myRank = tgt.getMap ()->getComm ()->getRank ();
    std::ostringstream os;
    os << "Proc " << myRank << ": " << funcName << ": ";
    prefix = std::unique_ptr<std::string> (new std::string (os.str ()));
    os << "Start" << endl;
    std::cerr << os.str ();
  }

  if (combineMode != ADD && combineMode != REPLACE &&
      combineMode != INSERT && combineMode != ABSMAX &&
      combineMode != ZERO) {
    std::ostringstream os;
    os << funcName << ": Invalid CombineMode: " << combineMode << ".  "
      "Valid CombineMode values: ADD, REPLACE, INSERT, ABSMAX, and ZERO.";
    errMsg = std::unique_ptr<std::string> (new std::string (os.str ()));
    return { -1, std::move (errMsg) };
  }
  if (combineMode == ZERO) {
    return { 0, std::move (errMsg) };
  }

  if (importLIDs.need_sync_host ()) {
    const int myRank = tgt.getMap ()->getComm ()->getRank ();
    std::ostringstream os;
    os << "Proc " << myRank << ": " << funcName << ": importLIDs "
      "needs sync to host.  This should never happen.";
    errMsg = std::unique_ptr<std::string> (new std::string (os.str ()));
    return { -2, std::move (errMsg) };
  }

  // Number of things to pack per LID is the block size.
  // Input LIDs correspond to the mesh points, not the DOFs.
  const size_t numMeshLIDs = static_cast<size_t> (importLIDs.extent (0));
  const size_t blockSize = static_cast<size_t> (tgt.getBlockSize ());
  const size_t numVecs = static_cast<size_t> (tgt.getNumVectors ());

  const size_t requiredImportsSize = numMeshLIDs * blockSize * numVecs;
  if (static_cast<size_t> (imports.extent (0)) < requiredImportsSize) {
    const int myRank = tgt.getMap ()->getComm ()->getRank ();
    std::ostringstream os;
    os << "Proc " << myRank << ": " << funcName << ": imports.extent(0) = "
       << imports.extent (0) << " < requiredImportsSize = "
       << requiredImportsSize << ".";
    errMsg = std::unique_ptr<std::string> (new std::string (os.str ()));
    return { -3, std::move (errMsg) };
  }

  auto importLIDs_h = importLIDs.view_host ();
  if (imports.need_sync_host ()) {
    imports.sync_host ();
  }
  auto imports_h = imports.view_host ();

  size_t curImportPos = 0;
  for (size_t meshLidIndex = 0; meshLidIndex < numMeshLIDs; ++meshLidIndex) {
    for (size_t j = 0; j < numVecs; ++j, curImportPos += blockSize) {
      const LO meshLid = importLIDs_h[meshLidIndex];
      const IST* const curImportPtr = &imports_h[curImportPos];

      src_little_vec_type X_src (curImportPtr, blockSize);
      auto X_dst = tgt.getLocalBlock (meshLid, j);

      if (combineMode == INSERT || combineMode == REPLACE) {
        deep_copy (X_dst, X_src);
      }
      else if (combineMode == ADD) {
        using ::Tpetra::Experimental::AXPY;
        AXPY (static_cast<IST> (KAT::one ()), X_src, X_dst);
      }
      else if (combineMode == ABSMAX) {
        ::Tpetra::Experimental::Impl::absMax (X_dst, X_src);
      }
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

#define TPETRA_DETAILS_UNPACKBLOCKMULTIVECTOR_INSTANT( SC, LO, GO, NT )  \
namespace Details { \
  template std::pair<int, std::unique_ptr<std::string>> \
  unpackBlockMultiVector<SC, LO, GO, NT> ( \
    BlockMultiVector<SC, LO, GO, NT>&, \
    const Kokkos::DualView< \
      const LO*, \
      BlockMultiVector<SC, LO, GO, NT>::buffer_device_type \
    >&, \
    Kokkos::DualView< \
      BlockMultiVector<SC, LO, GO, NT>::impl_scalar_type*, \
      BlockMultiVector<SC, LO, GO, NT>::buffer_device_type \
    >, \
    Kokkos::DualView< \
      size_t*, \
      BlockMultiVector<SC, LO, GO, NT>::buffer_device_type \
    >, \
    const size_t, \
    Distributor&, \
    const CombineMode); \
}

#endif // TPETRA_DETAILS_UNPACKBLOCKMULTIVECTOR_DEF_HPP
