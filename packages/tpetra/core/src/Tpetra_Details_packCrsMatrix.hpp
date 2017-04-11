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

#ifndef TPETRA_DETAILS_PACKCRSMATRIX_HPP
#define TPETRA_DETAILS_PACKCRSMATRIX_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Exceptions.hpp"

/// \file Tpetra_Details_packCrsMatrix.hpp
/// \brief Functions for packing the KokkosSparse::CrsMatrix

namespace Tpetra {
//
// Users must never rely on anything in the Details namespace.
//
namespace Details {


template<class SC, class LO, class GO, class NT, const bool classic>
void
allocatePackSpace(Teuchos::Array<char>& exports,
                  size_t& totalNumEntries,
                  const Teuchos::ArrayView<const LO>& exportLIDs,
                  const typename Tpetra::CrsMatrix<SC,LO,GO,NT,classic>::local_matrix_type& lclMatrix)
{
  //const char tfecfFuncName[] = "allocatePackSpace: ";

  typedef typename Tpetra::CrsMatrix<SC,LO,GO,NT,classic>::impl_scalar_type IST;

  // The number of export LIDs must fit in LocalOrdinal, assuming
  // that the LIDs are distinct and valid on the calling process.
  const LO numExportLIDs = static_cast<LO>(exportLIDs.size());

  // Count the total number of matrix entries to send.
  totalNumEntries = 0;
  for (LO i = 0; i < numExportLIDs; ++i) {
    const LO lclRow = exportLIDs[i];
    size_t curNumEntries;
    curNumEntries = static_cast<size_t>(lclMatrix.graph.row_map[lclRow+1] -
                                        lclMatrix.graph.row_map[lclRow]);

    // FIXME (mfh 25 Jan 2015) We should actually report invalid row
    // indices as an error.  Just consider them nonowned for now.
    if(curNumEntries == Teuchos::OrdinalTraits<size_t>::invalid()) {
      curNumEntries = 0;
    }
    totalNumEntries += curNumEntries;
  }

  // FIXME (mfh 24 Feb 2013, 24 Mar 2017) This code is only correct
  // if sizeof(IST) is a meaningful representation of the amount of
  // data in a Scalar instance.  (LO and GO are always built-in
  // integer types.)
  //
  // Allocate the exports array.  It does NOT need padding for
  // alignment, since we use memcpy to write to / read from send /
  // receive buffers.
  const size_t allocSize = static_cast<size_t>(numExportLIDs) * sizeof(LO)
                          + totalNumEntries * (sizeof(IST) + sizeof(GO));
  if (static_cast<size_t>(exports.size()) < allocSize) {
    exports.resize(allocSize);
  }
}

template<class SC, class LO, class GO, class NT, const bool classic>
bool
packCrsMatrixRow(char* const numEntOut,
        char* const valOut,
        char* const indOut,
        const size_t numEnt,
        const LO lclRow,
        const typename Tpetra::Map<LO,GO,NT>& colMap,
        const typename Tpetra::CrsMatrix<SC,LO,GO,NT,classic>::local_matrix_type& lclMatrix)
{
  typedef typename Tpetra::CrsMatrix<SC,LO,GO,NT,classic>::impl_scalar_type IST;
  typedef typename Tpetra::CrsMatrix<SC,LO,GO,NT,classic>::local_matrix_type::size_type offset_type;
  typedef Kokkos::pair<offset_type, offset_type> pair_type;

  const LO numEntLO = static_cast<LO>(numEnt);
  memcpy(numEntOut, &numEntLO, sizeof(LO));

  if (numEnt == 0) {
    return true; // nothing more to pack
  }

#ifdef HAVE_TPETRA_DEBUG
  if (lclRow >= lclMatrix.numRows() ||
      (static_cast<size_t>(lclRow + 1) >=
        static_cast<size_t>(lclMatrix.graph.row_map.dimension_0()))) {
#else // NOT HAVE_TPETRA_DEBUG
  if (lclRow >= lclMatrix.numRows()) {
#endif // HAVE_TPETRA_DEBUG
    // It's bad if this is not a valid local row index.  One thing
    // we can do is just pack the flag invalid value for the column
    // indices.  That makes sure that the receiving process knows
    // something went wrong.
    const GO flag = Tpetra::Details::OrdinalTraits<GO>::invalid();
    for (size_t k = 0; k < numEnt; ++k) {
      memcpy(indOut + k * sizeof(GO), &flag, sizeof(GO));
    }
    // The values don't actually matter, but we might as well pack
    // something here.
    const IST zero = Kokkos::ArithTraits<IST>::zero();
    for (size_t k = 0; k < numEnt; ++k) {
      memcpy(valOut + k * sizeof(GO), &zero, sizeof(GO));
    }
    return false;
  }

  // FIXME (mfh 24 Mar 2017) Everything here assumes UVM.  If we
  // want to fix that, we need to write a pack kernel for the whole
  // matrix, that runs on device.

  // Since the matrix is locally indexed on the calling process, we
  // have to use its column Map (which it _must_ have in this case)
  // to convert to global indices.
  const offset_type rowBeg = lclMatrix.graph.row_map[lclRow];
  const offset_type rowEnd = lclMatrix.graph.row_map[lclRow + 1];

  auto indIn = Kokkos::subview(lclMatrix.graph.entries,
                                pair_type(rowBeg, rowEnd));
  auto valIn = Kokkos::subview(lclMatrix.values,
                                pair_type(rowBeg, rowEnd));

  // Copy column indices one at a time, so that we don't need
  // temporary storage.
  for (size_t k = 0; k < numEnt; ++k) {
    const GO gblIndIn = colMap.getGlobalElement(indIn[k]);
    memcpy(indOut + k * sizeof(GO), &gblIndIn, sizeof(GO));
  }
  memcpy(valOut, valIn.ptr_on_device(), numEnt * sizeof(IST));
  return true;
}

template<class SC, class LO, class GO, class NT, const bool classic>
void
packCrsMatrix(const Teuchos::ArrayView<const LO>& exportLIDs,
     Teuchos::Array<char>& exports,
     const Teuchos::ArrayView<size_t>& numPacketsPerLID,
     size_t& constantNumPackets,
     Distributor& distor,
     const typename Tpetra::Map<LO,GO,NT>& colMap,
     const typename Tpetra::CrsMatrix<SC,LO,GO,NT,classic>::local_matrix_type& lclMatrix)
{
  typedef typename Tpetra::CrsMatrix<SC,LO,GO,NT,classic>::impl_scalar_type IST;

  const size_t numExportLIDs = static_cast<size_t>(exportLIDs.size());
  TEUCHOS_TEST_FOR_EXCEPTION
    (numExportLIDs != static_cast<size_t>(numPacketsPerLID.size()),
      std::invalid_argument, "exportLIDs.size() = " << numExportLIDs
      << " != numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");

  // Setting this to zero tells the caller to expect a possibly
  // different ("nonconstant") number of packets per local index
  // (i.e., a possibly different number of entries per row).
  constantNumPackets = 0;

  // The pack buffer 'exports' enters this method possibly
  // unallocated.  Do the first two parts of "Count, allocate, fill,
  // compute."
  size_t totalNumEntries = 0;
  allocatePackSpace<SC,LO,GO,NT,classic>(exports, totalNumEntries, exportLIDs, lclMatrix);
  const size_t bufSize = static_cast<size_t>(exports.size());

  // Compute the number of "packets" (in this case, bytes) per
  // export LID (in this case, local index of the row to send), and
  // actually pack the data.
  //
  // FIXME (mfh 24 Feb 2013, 25 Jan 2015) This code is only correct
  // if sizeof(Scalar) is a meaningful representation of the amount
  // of data in a Scalar instance.  (LO and GO are always built-in
  // integer types.)

  // Variables for error reporting in the loop.
  size_t firstBadIndex = 0; // only valid if outOfBounds == true.
  size_t firstBadOffset = 0;   // only valid if outOfBounds == true.
  size_t firstBadNumBytes = 0; // only valid if outOfBounds == true.
  bool outOfBounds = false;
  bool packErr = false;

  char* const exportsRawPtr = exports.getRawPtr();
  size_t offset = 0; // current index into 'exports' array.

  // Since the graph is static, we can go straight to lclMatrix
  // for the data.  This makes packing faster, and also makes
  // thread-parallel packing possible: see #800.
  for (size_t i = 0; i < numExportLIDs; ++i) {
    const LO lclRow = exportLIDs[i];

    // FIXME (mfh 24 Mar 2017) Everything here assumes UVM.  If
    // we want to fix that, we need to write a pack kernel for
    // the whole matrix, that runs on device.
    size_t numEnt = static_cast<size_t>(lclMatrix.graph.row_map[lclRow+1] -
                                        lclMatrix.graph.row_map[lclRow]);

    // Only pack this row's data if it has a nonzero number of
    // entries.  We can do this because receiving processes get the
    // number of packets, and will know that zero packets means zero
    // entries.
    if (numEnt == 0) {
      numPacketsPerLID[i] = 0;
    }
    else {
      char* const numEntBeg = exportsRawPtr + offset;
      char* const numEntEnd = numEntBeg + sizeof(LO);
      char* const valBeg = numEntEnd;
      char* const valEnd = valBeg + numEnt * sizeof(SC);
      char* const indBeg = valEnd;
      const size_t numBytes = sizeof(LO) +
        numEnt * (sizeof(IST) + sizeof(GO));
      if (offset > bufSize || offset + numBytes > bufSize) {
        firstBadIndex = i;
        firstBadOffset = offset;
        firstBadNumBytes = numBytes;
        outOfBounds = true;
        break;
      }

      packErr = ! packCrsMatrixRow<SC,LO,GO,NT,classic>(numEntBeg, valBeg, indBeg,
                                                        numEnt, lclRow, colMap, lclMatrix);
      if (packErr) {
        firstBadIndex = i;
        firstBadOffset = offset;
        firstBadNumBytes = numBytes;
        break;
      }
      // numPacketsPerLID[i] is the number of "packets" in the
      // current local row i.  Packet=char (really "byte") so use
      // the number of bytes of the packed data for that row.
      numPacketsPerLID[i] = numBytes;
      offset += numBytes;
    }
  }

}

} // namespace Details

} // namespace Tpetra

#endif // TPETRA_DETAILS_PACKCRSMATRIX_HPP
