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

#ifndef TPETRA_DETAILS_UNPACKCRSMATRIX_HPP
#define TPETRA_DETAILS_UNPACKCRSMATRIX_HPP

#include "TpetraCore_config.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Kokkos_Core.hpp"
#include <memory>
#include <string>

/// \file Tpetra_Details_unpackCrsMatrix.hpp
/// \brief Functions for unpacking the entries of a Tpetra::CrsMatrix
///   for communication, in the case where it is valid to go to the
///   KokkosSparse::CrsMatrix (local sparse matrix data structure)
///   directly.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// Forward declaration of Distributor
class Distributor;
#endif // DOXYGEN_SHOULD_SKIP_THIS

//
// Users must never rely on anything in the Details namespace.
//
namespace Details {

/// \brief combine values in a single CRS matrix row
///
/// \warning The allowed \c combineMode are:
///   ADD, REPLACE, and ABSMAX. INSERT is not allowed.
template<class LocalMatrixType, class LocalMapType>
typename LocalMatrixType::ordinal_type
combineCrsMatrixValues(LocalMatrixType& lclMatrix,
                       LocalMapType& lclColMap,
                       const typename LocalMatrixType::ordinal_type lclRow,
                       const typename LocalMatrixType::ordinal_type numEnt,
                       const typename LocalMatrixType::value_type vals[],
                       const typename LocalMapType::global_ordinal_type cols[],
                       const Tpetra::CombineMode combineMode,
                       const bool atomic)
{
  typedef typename LocalMatrixType::ordinal_type LO;

  // INSERT doesn't make sense for a static graph, since you
  // aren't allowed to change the structure of the graph.
  // However, all the other combine modes work.

  LO numValid = 0; // number of valid input column indices
  if (combineMode == ADD) {
    for (int k=0; k<numEnt; k++) {
      LO lclColInd = lclColMap.getLocalElement(cols[k]);
      numValid += lclMatrix.sumIntoValues(lclRow, &lclColInd, 1, &vals[k], false, atomic);
    }
  }
  else if (combineMode == REPLACE) {
    for (int k=0; k<numEnt; k++) {
      LO lclColInd = lclColMap.getLocalElement(cols[k]);
      numValid += lclMatrix.replaceValues(lclRow, &lclColInd, 1, &vals[k], false, atomic);
    }
  }
  else if (combineMode == ABSMAX) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument,
        "ABSMAX combine mode is not yet implemented for a matrix that has a "
        "static graph (i.e., was constructed with the CrsMatrix constructor "
        "that takes a const CrsGraph pointer).");
  }
  else if (combineMode == INSERT) {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument,
        "INSERT combine mode is not allowed if the matrix has a static graph "
        "(i.e., was constructed with the CrsMatrix constructor that takes a "
        "const CrsGraph pointer).");
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error, "Invalid combine mode; should never get "
        "here!  Please report this bug to the Tpetra developers.");
  }

  return numValid;

}


/// \brief Reduction result for UnpackCrsMatrixAndCombineFunctor below.
///
/// The reduction result finds the offset and number of bytes associated with
/// the first out of bounds error or unpacking error occurs.
template<class LO>
struct UnpackCrsMatrixError {

  LO bad_index; // only valid if outOfBounds == true.
  size_t num_ent;
  size_t offset;   // only valid if outOfBounds == true.
  size_t expected_num_bytes;
  size_t num_bytes; // only valid if outOfBounds == true.
  size_t buf_size;
  bool wrong_num_bytes_error;
  bool out_of_bounds_error;
  bool unpacking_error;

  UnpackCrsMatrixError():
    bad_index(Tpetra::Details::OrdinalTraits<LO>::invalid()),
    num_ent(0),
    offset(0),
    expected_num_bytes(0),
    num_bytes(0),
    buf_size(0),
    wrong_num_bytes_error(false),
    out_of_bounds_error(false),
    unpacking_error(false) {}

  bool success() const
  {
    // Any possible error would result in bad_index being changed from
    // `invalid` to the index associated with the error.
    return bad_index == Tpetra::Details::OrdinalTraits<LO>::invalid();
  }

  std::string summary() const
  {
    std::ostringstream os;
    if (wrong_num_bytes_error || out_of_bounds_error || unpacking_error) {
      if (wrong_num_bytes_error) {
        os << "At index i = " << bad_index
           << ", expectedNumBytes > numBytes.";
      }
      else if (out_of_bounds_error) {
        os << "First invalid offset into 'imports' "
           << "unpack buffer at index i = " << bad_index << ".";
      }
      else if (unpacking_error) {
        os << "First error in unpackRow() at index i = "
           << bad_index << ".";
      }
      os << "  importLIDs[i]: " << bad_index
         << ", bufSize: " << buf_size
         << ", offset: " << offset
         << ", numBytes: " << num_bytes
         << ", expectedNumBytes: " << expected_num_bytes
         << ", numEnt: " << num_ent;
    }
    return os.str();
  }
};


/// \brief Unpacks and combines a single row of the CrsMatrix.
///
/// Data (bytes) describing the row of the CRS matrix are "unpacked"
/// from a single (concatenated) char* in to the row of the matrix
///
/// \tparam LocalMatrixType the specialization of the KokkosSparse::CrsMatrix
///   local matrix
/// \tparam LocalMapType the type of the local column map
template<class CountType, class OffsetType, class LocalMatrixType, class LocalMapType>
struct UnpackCrsMatrixAndCombineFunctor {

  typedef CountType count_type;
  typedef OffsetType offset_type;
  typedef LocalMatrixType local_matrix_type;
  typedef LocalMapType local_map_type;
  typedef typename local_matrix_type::value_type IST;
  typedef typename local_matrix_type::ordinal_type LO;
  typedef typename local_map_type::global_ordinal_type GO;

  typedef Kokkos::View<const count_type*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> num_packets_per_lid_view_type;
  typedef Kokkos::View<offset_type*, Kokkos::HostSpace> offsets_view_type;
  typedef Kokkos::View<const char*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> imports_view_type;
  typedef Kokkos::View<const LO*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged> import_lids_view_type;

  typedef UnpackCrsMatrixError<LO> value_type;

  static_assert (std::is_same<LO, typename LocalMatrixType::ordinal_type>::value,
                 "LocalMapType::local_ordinal_type and "
                 "LocalMatrixType::ordinal_type must be the same.");

  num_packets_per_lid_view_type num_packets_per_lid_;
  offsets_view_type offsets_;
  imports_view_type imports_;
  import_lids_view_type import_lids_;
  LocalMatrixType local_matrix_;
  LocalMapType local_col_map_;
  Tpetra::CombineMode combine_mode_;
  bool atomic_;

  UnpackCrsMatrixAndCombineFunctor(num_packets_per_lid_view_type num_packets_per_lid,
      offsets_view_type offsets, imports_view_type imports,
      import_lids_view_type import_lids,
      local_matrix_type local_matrix, local_map_type local_col_map,
      Tpetra::CombineMode combine_mode, bool atomic) :
    num_packets_per_lid_(num_packets_per_lid), offsets_(offsets),
    imports_(imports), import_lids_(import_lids), local_matrix_(local_matrix),
    local_col_map_(local_col_map), combine_mode_(combine_mode), atomic_(atomic)
  {}

  KOKKOS_INLINE_FUNCTION void init(value_type& dst) const
  {
    dst.bad_index = Tpetra::Details::OrdinalTraits<LO>::invalid();
    dst.num_ent = 0;
    dst.offset = 0;
    dst.expected_num_bytes = 0;
    dst.num_bytes = 0;
    dst.buf_size = 0;
    dst.wrong_num_bytes_error = false;
    dst.out_of_bounds_error = false;
    dst.unpacking_error = false;
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst, const volatile value_type& src) const
  {
    // The dst object should reflect the first bad index and all other
    // associated error codes and data.  Thus, we need only check if the src
    // object shows an error and if its associated `bad_index` is less than the
    // dst.bad_index (if dst shows errors).
    LO invalid = Tpetra::Details::OrdinalTraits<LO>::invalid();
    if (src.bad_index != invalid) {
      // An error in the src, check if whether:
      //   1. The dst shows errors
      //   2. If dst does show errors, if src bad_index is less than its bad
      //      index
      if (dst.bad_index == invalid || src.bad_index < dst.bad_index) {
        dst.bad_index = src.bad_index;
        dst.num_ent = src.num_ent;
        dst.offset = src.offset;
        dst.expected_num_bytes = src.expected_num_bytes;
        dst.num_bytes = src.num_bytes;
        dst.buf_size = src.buf_size;
        dst.wrong_num_bytes_error = src.wrong_num_bytes_error;
        dst.out_of_bounds_error = src.out_of_bounds_error;
        dst.unpacking_error = src.unpacking_error;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO i, value_type& dst) const
  {
    const LO local_row = import_lids_[i];
    const size_t num_bytes = num_packets_per_lid_(i);

    if (num_bytes > 0) { // there is actually something in the row

      const size_t buf_size = imports_.size();
      const char* const num_ent_beg = imports_.ptr_on_device() + offsets_(i);
      const char* const num_ent_end = num_ent_beg + sizeof(LO);

      // Now we know how many entries to expect in the received data
      // for this row.
      LO num_ent = 0;
      memcpy(&num_ent, num_ent_beg, sizeof(LO));

      const char* const val_beg = num_ent_end;
      const char* const val_end = val_beg + static_cast<size_t>(num_ent) * sizeof(IST);
      const char* const ind_beg = val_end;
      const size_t expected_num_bytes = sizeof(LO) +
        static_cast<size_t>(num_ent) * (sizeof(IST) + sizeof(GO));

      if (expected_num_bytes > num_bytes) {
        dst.wrong_num_bytes_error = true;
      }
      else if (offsets_(i) > buf_size || offsets_(i) + num_bytes > buf_size) {
        dst.out_of_bounds_error = true;
      }
      else {

        // Unpack this row
        IST* vals = new IST[num_ent];
        memcpy(vals, val_beg, num_ent * sizeof(IST));

        GO* cols = new GO[num_ent];
        memcpy(cols, ind_beg, num_ent * sizeof(GO));

        // FIXME (mfh 23 Mar 2017) It would make sense to use the return
        // value here as more than just a "did it succeed" Boolean test.

        // FIXME (mfh 23 Mar 2017) CrsMatrix_NonlocalSumInto_Ignore test
        // expects this method to ignore incoming entries that do not
        // exist on the process that owns those rows.  We would like to
        // distinguish between "errors" resulting from ignored entries,
        // vs. actual errors.

        // Combine values
        auto num_modified =
          combineCrsMatrixValues(local_matrix_, local_col_map_,
              local_row, num_ent, vals, cols, combine_mode_, atomic_);

        delete [] vals;
        delete [] cols;
      }

      if (dst.wrong_num_bytes_error || dst.out_of_bounds_error || dst.unpacking_error) {
        dst.bad_index = i;
        dst.num_ent = num_ent;
        dst.offset = offsets_(i);
        dst.expected_num_bytes = expected_num_bytes;
        dst.num_bytes = num_bytes;
        dst.buf_size = buf_size;
      }
    }
  }

};

/// \brief Unpack the imported column indices and values, and combine into matrix.
///
/// \warning The allowed \c combineMode are:
///   ADD, REPLACE, and ABSMAX. INSERT is not allowed.
template<class LocalMatrixType, class LocalMapType>
bool
unpackCrsMatrixAndCombine(
    LocalMatrixType& lclMatrix,
    const LocalMapType& lclColMap,
    std::unique_ptr<std::string>& errStr,
    const Teuchos::ArrayView<const typename LocalMatrixType::ordinal_type>& importLIDs,
    const Teuchos::ArrayView<const char>& imports,
    const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
    size_t constantNumPackets,
    const int myRank,
    Distributor & /* distor */,
    CombineMode combineMode,
    const bool atomic)
{
  using Kokkos::View;
  using Kokkos::HostSpace;
  using Kokkos::MemoryUnmanaged;
  using ::Tpetra::Details::computeOffsetsFromCounts;
  typedef LocalMatrixType local_matrix_type;
  typedef LocalMapType local_map_type;
  typedef typename LocalMapType::local_ordinal_type LO;

  static_assert (std::is_same<LO, typename LocalMatrixType::ordinal_type>::value,
                 "LocalMapType::local_ordinal_type and "
                 "LocalMatrixType::ordinal_type must be the same.");

  // Get the number of packets on each row.
  typedef size_t count_type;
  typedef View<const LO*, HostSpace, MemoryUnmanaged> LIVT;  // Local indices view type
  typedef View<const count_type*, HostSpace, MemoryUnmanaged> CVT; // Counts view type
  typedef View<const char*, HostSpace, MemoryUnmanaged> IPVT; // Imports view type

  typedef typename CVT::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, LO> range_type;

  const size_t numImportLIDs = static_cast<size_t>(importLIDs.size());
  if (numImportLIDs != static_cast<size_t>(numPacketsPerLID.size())) {
    std::ostringstream os;
    os << "importLIDs.size() (" << numImportLIDs << ") != "
       << "numPacketsPerLID.size() (" << numPacketsPerLID.size() << ").";
    if (errStr.get() == NULL) {
      errStr = std::unique_ptr<std::string>(new std::string());
    }
    *errStr = os.str();
    return false;
  }

  CVT num_packets_per_lid(numPacketsPerLID.getRawPtr(), numPacketsPerLID.size());
  LIVT import_lids(importLIDs.getRawPtr(), importLIDs.size());
  IPVT imports_v(imports.getRawPtr(), imports.size());

  // Get the offsets
  typedef size_t offset_type;
  typedef View<offset_type*, HostSpace> OVT; // Offsets view type
  OVT offsets("unpackCrsMatrixAndCombine: offsets", numImportLIDs+1);
  computeOffsetsFromCounts(offsets, num_packets_per_lid);

  // Now do the actual unpack!
  typedef UnpackCrsMatrixAndCombineFunctor<count_type, offset_type,
          local_matrix_type, local_map_type> unpack_functor_type;

  typename unpack_functor_type::value_type result;
  unpack_functor_type unpack_functor(num_packets_per_lid, offsets,
      imports_v, import_lids, lclMatrix, lclColMap, combineMode, atomic);
  Kokkos::parallel_reduce(range_type(0, numImportLIDs), unpack_functor, result);

  if (!result.success()) {
    std::ostringstream os;
    os << "Proc " << myRank << ": packCrsMatrix failed.  "
       << result.summary()
       << std::endl;
    if (errStr.get () != NULL) {
      errStr = std::unique_ptr<std::string> (new std::string ());
      *errStr = os.str ();
    }
  }

  return result.success();

}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_UNPACKCRSMATRIX_HPP
