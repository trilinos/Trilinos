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
  bool invalid_combine_mode_error;

  KOKKOS_INLINE_FUNCTION UnpackCrsMatrixError():
    bad_index(Tpetra::Details::OrdinalTraits<LO>::invalid()),
    num_ent(0),
    offset(0),
    expected_num_bytes(0),
    num_bytes(0),
    buf_size(0),
    wrong_num_bytes_error(false),
    out_of_bounds_error(false),
    invalid_combine_mode_error(false) {}

  bool success() const
  {
    // Any possible error would result in bad_index being changed from
    // `invalid` to the index associated with the error.
    return bad_index == Tpetra::Details::OrdinalTraits<LO>::invalid();
  }

  std::string summary() const
  {
    std::ostringstream os;
    if (wrong_num_bytes_error ||
        out_of_bounds_error ||
        invalid_combine_mode_error) {
      if (wrong_num_bytes_error) {
        os << "At index i = " << bad_index
           << ", expectedNumBytes > numBytes.";
      }
      else if (out_of_bounds_error) {
        os << "First invalid offset into 'imports' "
           << "unpack buffer at index i = " << bad_index << ".";
      }
      else if (invalid_combine_mode_error) {
        os << "Invalid combine mode; code should never get "
              "here!  Please report this bug to the Tpetra developers.";
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
/// \tparam NumPacketsPerLIDType the specialization of the Kokkos::View of counts
/// \tparam OffsetsType the specialization of the Kokkos::View of offsets
/// \tparam ImportsType the specialization of the Kokkos::View of imports
/// \tparam ImportLIDsType the specialization of the Kokkos::View of import
//local IDs
/// \tparam LocalMatrixType the specialization of the KokkosSparse::CrsMatrix
///   local matrix
/// \tparam LocalMapType the type of the local column map
template<class NumPacketsPerLIDType, class OffsetsType,
  class ImportsType, class ImportLIDsType,
  class LocalMatrixType, class LocalMapType>
struct UnpackCrsMatrixAndCombineFunctor {

  typedef NumPacketsPerLIDType num_packets_per_lid_type;
  typedef OffsetsType offsets_type;
  typedef ImportsType imports_type;
  typedef ImportLIDsType import_lids_type;
  typedef LocalMatrixType local_matrix_type;
  typedef LocalMapType local_map_type;
  typedef typename local_matrix_type::value_type IST;
  typedef typename local_matrix_type::ordinal_type LO;
  typedef typename local_map_type::global_ordinal_type GO;
  typedef UnpackCrsMatrixError<LO> value_type;

  static_assert (Kokkos::Impl::is_view<NumPacketsPerLIDType>::value,
                 "NumPacketsPerLIDType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<OffsetsType>::value,
                 "OffsetsType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<ImportsType>::value,
                 "ImportsType must be a Kokkos::View.");
  static_assert (Kokkos::Impl::is_view<ImportLIDsType>::value,
                 "ImportLIDsType must be a Kokkos::View.");
  static_assert (std::is_same<LO, typename LocalMatrixType::ordinal_type>::value,
                 "LocalMapType::local_ordinal_type and "
                 "LocalMatrixType::ordinal_type must be the same.");

  NumPacketsPerLIDType num_packets_per_lid_;
  OffsetsType offsets_;
  ImportsType imports_;
  ImportLIDsType import_lids_;
  LocalMatrixType local_matrix_;
  LocalMapType local_col_map_;
  Tpetra::CombineMode combine_mode_;
  bool atomic_;

  UnpackCrsMatrixAndCombineFunctor(
      const num_packets_per_lid_type& num_packets_per_lid,
      const offsets_type& offsets,
      const imports_type& imports,
      const import_lids_type& import_lids,
      local_matrix_type& local_matrix,
      const local_map_type& local_col_map,
      const Tpetra::CombineMode combine_mode,
      const bool atomic) :
    num_packets_per_lid_(num_packets_per_lid),
    offsets_(offsets),
    imports_(imports),
    import_lids_(import_lids),
    local_matrix_(local_matrix),
    local_col_map_(local_col_map),
    combine_mode_(combine_mode),
    atomic_(atomic)
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
    dst.invalid_combine_mode_error = false;
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
        dst.invalid_combine_mode_error = src.invalid_combine_mode_error;
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

        // FIXME (mfh 23 Mar 2017) CrsMatrix_NonlocalSumInto_Ignore test
        // expects this method to ignore incoming entries that do not
        // exist on the process that owns those rows.  We would like to
        // distinguish between "errors" resulting from ignored entries,
        // vs. actual errors.

        // Combine a single column/value at a time to avoid temporary storage
        LO num_modified = 0;
        if (combine_mode_ == ADD) {
          for (size_t k = 0; k < static_cast<size_t>(num_ent); ++k) {
            IST val = 0;
            GO col = 0;
            memcpy(&val, val_beg + k * sizeof(IST), sizeof(IST));
            memcpy(&col, ind_beg + k * sizeof(GO), sizeof(GO));
            LO local_col_idx = local_col_map_.getLocalElement(col);
            num_modified += local_matrix_.sumIntoValues(local_row,
                &local_col_idx, 1, &val, false, atomic_);
          }
        }
        else if (combine_mode_ == REPLACE) {
          for (size_t k = 0; k < static_cast<size_t>(num_ent); ++k) {
            IST val = 0;
            GO col = 0;
            memcpy(&val, val_beg + k * sizeof(IST), sizeof(IST));
            memcpy(&col, ind_beg + k * sizeof(GO), sizeof(GO));
            LO local_col_idx = local_col_map_.getLocalElement(col);
            num_modified += local_matrix_.replaceValues(local_row,
                &local_col_idx, 1, &val, false, atomic_);
          }
        }
        else {
          dst.invalid_combine_mode_error = true;
        }
      }

      if (dst.wrong_num_bytes_error ||
          dst.out_of_bounds_error ||
          dst.invalid_combine_mode_error) {
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
  using ::Tpetra::Details::computeOffsetsFromCounts;
  using Tpetra::Details::create_mirror_view_from_raw_host_array;
  typedef LocalMatrixType local_matrix_type;
  typedef LocalMapType local_map_type;
  typedef typename LocalMapType::local_ordinal_type LO;

  static_assert (std::is_same<LO, typename LocalMatrixType::ordinal_type>::value,
                 "LocalMapType::local_ordinal_type and "
                 "LocalMatrixType::ordinal_type must be the same.");

  typedef typename LocalMatrixType::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, LO> range_type;

  const size_t numImportLIDs = static_cast<size_t>(importLIDs.size());

  {
    // Check for correct input
    int qa_errors = 0;
    std::ostringstream os;
    if (combineMode == ABSMAX) {
      qa_errors += 1;
      os << "ABSMAX combine mode is not yet implemented for a matrix that has a "
         << "static graph (i.e., was constructed with the CrsMatrix constructor "
         << "that takes a const CrsGraph pointer).\n";
    }
    else if (combineMode == INSERT) {
      qa_errors += 1;
      os << "INSERT combine mode is not allowed if the matrix has a static graph "
         << "(i.e., was constructed with the CrsMatrix constructor that takes a "
         << "const CrsGraph pointer).\n";
    }
    else if (!(combineMode == ADD || combineMode == REPLACE)) {
      qa_errors += 1;
      // Unknown combine mode!
      os << "Invalid combine mode; should never get "
         << "here!  Please report this bug to the Tpetra developers.\n";
    }

    // Check that sizes of input objects are consistent.
    if (numImportLIDs != static_cast<size_t>(numPacketsPerLID.size())) {
      qa_errors += 1;
      os << "importLIDs.size() (" << numImportLIDs << ") != "
         << "numPacketsPerLID.size() (" << numPacketsPerLID.size() << ").";
    }

    if (qa_errors) {
      if (errStr.get() == NULL) {
        errStr = std::unique_ptr<std::string>(new std::string());
      }
      *errStr = os.str();
      return false;
    }
  } // end QA error checking

  // numPacketsPerLID, importLIDs, and imports are input, so we have to copy
  // them to device.  Since unpacking is done directly in to the local matrix
  // (lclMatrix), not copying needs to be performed after unpacking.
  typename LocalMatrixType::device_type outputDevice;
  auto num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array(outputDevice,
        numPacketsPerLID.getRawPtr(), numPacketsPerLID.size(),
        true, "num_packets_per_lid_d");

  auto import_lids_d =
    create_mirror_view_from_raw_host_array(outputDevice,
        importLIDs.getRawPtr(), importLIDs.size(),
        true, "import_lids_d");

  auto imports_d =
    create_mirror_view_from_raw_host_array(outputDevice,
        imports.getRawPtr(), imports.size(),
        true, "imports_d");

  // Get the offsets
  View<size_t*, device_type>
    offsets_d("offsets_d", numImportLIDs+1);
  computeOffsetsFromCounts(offsets_d, num_packets_per_lid_d);

  // Now do the actual unpack!
  typedef UnpackCrsMatrixAndCombineFunctor<
    decltype(num_packets_per_lid_d),
    decltype(offsets_d),
    decltype(imports_d),
    decltype(import_lids_d),
    local_matrix_type,
    local_map_type> unpack_functor_type;
  unpack_functor_type unpack_functor(num_packets_per_lid_d, offsets_d,
      imports_d, import_lids_d, lclMatrix, lclColMap, combineMode, atomic);

  typename unpack_functor_type::value_type result;
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
