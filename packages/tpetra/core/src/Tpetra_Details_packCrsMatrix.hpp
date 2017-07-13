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

#include "TpetraCore_config.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_createMirrorView.hpp"
#include "Kokkos_Core.hpp"
#include <memory>
#include <string>

/// \file Tpetra_Details_packCrsMatrix.hpp
/// \brief Functions for packing the entries of a Tpetra::CrsMatrix
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

template<class OutputOffsetsViewType,
         class CountsViewType,
         class InputOffsetsViewType,
         class InputLocalRowIndicesViewType,
         const bool debug =
#ifdef HAVE_TPETRA_DEBUG
         true
#else
         false
#endif // HAVE_TPETRA_DEBUG
         >
class NumPacketsAndOffsetsFunctor {
public:
  typedef typename OutputOffsetsViewType::non_const_value_type output_offset_type;
  typedef typename CountsViewType::non_const_value_type count_type;
  typedef typename InputOffsetsViewType::non_const_value_type input_offset_type;
  typedef typename InputLocalRowIndicesViewType::non_const_value_type local_row_index_type;
  // output Views drive where execution happens.
  typedef typename OutputOffsetsViewType::device_type device_type;
  static_assert (std::is_same<typename CountsViewType::device_type::execution_space,
                   typename device_type::execution_space>::value,
                 "OutputOffsetsViewType and CountsViewType must have the same execution space.");
  static_assert (Kokkos::Impl::is_view<OutputOffsetsViewType>::value,
                 "OutputOffsetsViewType must be a Kokkos::View.");
  static_assert (std::is_same<typename OutputOffsetsViewType::value_type, output_offset_type>::value,
                 "OutputOffsetsViewType must be a nonconst Kokkos::View.");
  static_assert (std::is_integral<output_offset_type>::value,
                 "The type of each entry of OutputOffsetsViewType must be a built-in integer type.");
  static_assert (Kokkos::Impl::is_view<CountsViewType>::value,
                 "CountsViewType must be a Kokkos::View.");
  static_assert (std::is_same<typename CountsViewType::value_type, output_offset_type>::value,
                 "CountsViewType must be a nonconst Kokkos::View.");
  static_assert (std::is_integral<count_type>::value,
                 "The type of each entry of CountsViewType must be a built-in integer type.");
  static_assert (Kokkos::Impl::is_view<InputOffsetsViewType>::value,
                 "InputOffsetsViewType must be a Kokkos::View.");
  static_assert (std::is_integral<input_offset_type>::value,
                 "The type of each entry of InputOffsetsViewType must be a built-in integer type.");
  static_assert (Kokkos::Impl::is_view<InputLocalRowIndicesViewType>::value,
                 "InputLocalRowIndicesViewType must be a Kokkos::View.");
  static_assert (std::is_integral<local_row_index_type>::value,
                 "The type of each entry of InputLocalRowIndicesViewType must be a built-in integer type.");

  NumPacketsAndOffsetsFunctor (const OutputOffsetsViewType& outputOffsets,
                               const CountsViewType& counts,
                               const InputOffsetsViewType& rowOffsets,
                               const InputLocalRowIndicesViewType& lclRowInds,
                               const count_type sizeOfLclCount,
                               const count_type sizeOfValue,
                               const count_type sizeOfGblColInd) :
    outputOffsets_ (outputOffsets),
    counts_ (counts),
    rowOffsets_ (rowOffsets),
    lclRowInds_ (lclRowInds),
    sizeOfLclCount_ (sizeOfLclCount),
    sizeOfValue_ (sizeOfValue),
    sizeOfGblColInd_ (sizeOfGblColInd),
    error_ ("error") // don't forget this, or you'll get segfaults!
  {
    if (debug) {
      const size_t numRowsToPack = static_cast<size_t> (lclRowInds_.dimension_0 ());

      if (numRowsToPack != static_cast<size_t> (counts_.dimension_0 ())) {
        std::ostringstream os;
        os << "lclRowInds.dimension_0() = " << numRowsToPack
           << " != counts.dimension_0() = " << counts_.dimension_0 ()
           << ".";
        throw std::invalid_argument (os.str ());
      }
      if (static_cast<size_t> (numRowsToPack + 1) != static_cast<size_t> (outputOffsets_.dimension_0 ())) {
        std::ostringstream os;
        os << "lclRowInds.dimension_0() + 1 = " << (numRowsToPack + 1)
           << " != outputOffsets.dimension_0() = " << outputOffsets_.dimension_0 ()
           << ".";
        throw std::invalid_argument (os.str ());
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const local_row_index_type& curInd,
              output_offset_type& update,
              const bool final) const
  {
    if (debug) {
      if (curInd < static_cast<local_row_index_type> (0)) {
        error_ () = 1;
        return;
      }
    }

    if (final) {
      if (debug) {
        if (curInd >= static_cast<local_row_index_type> (outputOffsets_.dimension_0 ())) {
          error_ () = 2;
          return;
        }
      }
      outputOffsets_(curInd) = update;
    }

    if (curInd < static_cast<local_row_index_type> (counts_.dimension_0 ())) {
      const auto lclRow = lclRowInds_(curInd);
      if (static_cast<size_t> (lclRow + 1) >= static_cast<size_t> (rowOffsets_.dimension_0 ()) ||
          static_cast<local_row_index_type> (lclRow) < static_cast<local_row_index_type> (0)) {
        error_ () = 3;
        return;
      }
      // count_type could differ from the type of each row offset.
      // For example, row offsets might each be 64 bits, but if their
      // difference always fits in 32 bits, we may then safely use a
      // 32-bit count_type.
      const count_type count =
        static_cast<count_type> (rowOffsets_(lclRow+1) - rowOffsets_(lclRow));

      // We pack first the number of entries in the row, then that many
      // global column indices, then that many values.  However, if the
      // number of entries in the row is zero, we pack nothing.
      const count_type numBytes = (count == 0) ?
        static_cast<count_type> (0) :
        sizeOfLclCount_ + count * (sizeOfGblColInd_ + sizeOfValue_);

      if (final) {
        counts_(curInd) = numBytes;
      }
      update += numBytes;
    }
  }

  // mfh 31 May 2017: Don't need init or join.  If you have join, MUST
  // have join both with and without volatile!  Otherwise intrawarp
  // joins are really slow on GPUs.

  //! Host function for getting the error.
  int getError () const {
    auto error_h = Kokkos::create_mirror_view (error_);
    Kokkos::deep_copy (error_h, error_);
    return error_h ();
  }

private:
  OutputOffsetsViewType outputOffsets_;
  CountsViewType counts_;
  typename InputOffsetsViewType::const_type rowOffsets_;
  typename InputLocalRowIndicesViewType::const_type lclRowInds_;
  count_type sizeOfLclCount_;
  count_type sizeOfValue_;
  count_type sizeOfGblColInd_;
  Kokkos::View<int, device_type> error_;
};

template<class OutputOffsetsViewType,
         class CountsViewType,
         class InputOffsetsViewType,
         class InputLocalRowIndicesViewType>
std::pair<typename CountsViewType::non_const_value_type, bool>
computeNumPacketsAndOffsets (std::unique_ptr<std::ostringstream>& errStr,
                             const OutputOffsetsViewType& outputOffsets,
                             const CountsViewType& counts,
                             const InputOffsetsViewType& rowOffsets,
                             const InputLocalRowIndicesViewType& lclRowInds,
                             const typename CountsViewType::non_const_value_type sizeOfLclCount,
                             const typename CountsViewType::non_const_value_type sizeOfValue,
                             const typename CountsViewType::non_const_value_type sizeOfGblColInd)
{
  typedef NumPacketsAndOffsetsFunctor<OutputOffsetsViewType,
    CountsViewType, typename InputOffsetsViewType::const_type,
    typename InputLocalRowIndicesViewType::const_type> functor_type;
  typedef typename CountsViewType::non_const_value_type count_type;
  typedef typename OutputOffsetsViewType::size_type size_type;
  typedef typename OutputOffsetsViewType::execution_space execution_space;
  typedef typename functor_type::local_row_index_type LO;
  typedef Kokkos::RangePolicy<execution_space, LO> range_type;
  const char prefix[] = "computeNumPacketsAndOffsets: ";

  const count_type numRowsToPack = lclRowInds.dimension_0 ();
  if (numRowsToPack == 0) {
    return {static_cast<count_type> (0), true}; // nothing to pack, but no error
  }
  else {
    if (rowOffsets.dimension_0 () <= static_cast<size_type> (1)) {
      if (errStr.get () == NULL) {
        errStr = std::unique_ptr<std::ostringstream> (new std::ostringstream ());
      }
      std::ostringstream& os = *errStr;
      os << prefix
         << "There is at least one row to pack, but the matrix has no rows.  "
        "lclRowInds.dimension_0() = " << numRowsToPack << ", but "
        "rowOffsets.dimension_0() = " << rowOffsets.dimension_0 () << " <= 1."
         << std::endl;
      return {static_cast<count_type> (0), false};
    }
    if (outputOffsets.dimension_0 () != static_cast<size_type> (numRowsToPack + 1)) {
      if (errStr.get () == NULL) {
        errStr = std::unique_ptr<std::ostringstream> (new std::ostringstream ());
      }
      std::ostringstream& os = *errStr;
      os << prefix
         << "Output dimension does not match number of rows to pack.  "
         << "outputOffsets.dimension_0() = " << outputOffsets.dimension_0 ()
         << " != lclRowInds.dimension_0() + 1 = "
         << static_cast<size_type> (numRowsToPack + 1) << "." << std::endl;
      return {static_cast<count_type> (0), false};
    }

    functor_type f (outputOffsets, counts, rowOffsets,
                    lclRowInds, sizeOfLclCount,
                    sizeOfValue, sizeOfGblColInd);
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(static_cast<size_t> (outputOffsets.dimension_0 ()) !=
                            static_cast<size_t> (numRowsToPack + 1));
    TEUCHOS_TEST_FOR_EXCEPT(static_cast<size_t> (counts.dimension_0 ()) !=
                            static_cast<size_t> (numRowsToPack));
    TEUCHOS_TEST_FOR_EXCEPT(static_cast<size_t> (outputOffsets.dimension_0 ()) !=
                            static_cast<size_t> (numRowsToPack + 1));
#endif // HAVE_TPETRA_DEBUG

    Kokkos::parallel_scan (range_type (0, numRowsToPack+1), f);

    // At least in debug mode, this functor checks for errors.
    const int errCode = f.getError ();
    if (errCode != 0) {
      if (errStr.get () == NULL) {
        errStr = std::unique_ptr<std::ostringstream> (new std::ostringstream ());
      }
      std::ostringstream& os = *errStr;
      os << prefix
         << "NumPacketsAndOffsetsFunctor reported error code " << errCode
         << " != 0." << std::endl;
      return {0, false};
    }

#if 0
    size_t total = 0;
    for (LO k = 0; k < numRowsToPack; ++k) {
      total += counts[k];
    }
    if (outputOffsets(numRowsToPack) != total) {
      if (errStr.get () == NULL) {
        errStr = std::unique_ptr<std::ostringstream> (new std::ostringstream ());
      }
      std::ostringstream& os = *errStr;
      os << prefix
         << "outputOffsets(numRowsToPack=" << numRowsToPack << ") "
         << outputOffsets(numRowsToPack) << " != sum of counts = "
         << total << "." << std::endl;
      if (numRowsToPack != 0) {
        // Only print the array if it's not too long.
        if (numRowsToPack < static_cast<LO> (10)) {
          os << "outputOffsets: [";
          for (LO i = 0; i <= numRowsToPack; ++i) {
            os << outputOffsets(i);
            if (static_cast<LO> (i + 1) <= numRowsToPack) {
              os << ",";
            }
          }
          os << "]" << std::endl;
          os << "counts: [";
          for (LO i = 0; i < numRowsToPack; ++i) {
            os << counts(i);
            if (static_cast<LO> (i + 1) < numRowsToPack) {
              os << ",";
            }
          }
          os << "]" << std::endl;
        }
        else {
          os << "outputOffsets(" << (numRowsToPack-1) << ") = "
             << outputOffsets(numRowsToPack-1) << "." << std::endl;
        }
      }
      return {outputOffsets(numRowsToPack), false};
    }
#endif // HAVE_TPETRA_DEBUG

    // Get last entry of outputOffsets, which is the sum of the entries
    // of counts.  Don't assume UVM.
    auto outputOffsets_last = Kokkos::subview (outputOffsets, numRowsToPack);
    auto outputOffsets_last_h = Kokkos::create_mirror_view (outputOffsets_last);
    return {static_cast<count_type> (outputOffsets_last_h ()), true};
  }
}

/// \brief Reduction result for PackCrsMatrixFunctor below.
///
/// The reduction result finds the offset and number of bytes associated with
/// the first out of bounds error or packing error occurs.
template<class LO>
struct PackCrsMatrixError {

  LO bad_index; // only valid if outOfBounds == true.
  size_t bad_offset;   // only valid if outOfBounds == true.
  size_t bad_num_bytes; // only valid if outOfBounds == true.
  bool out_of_bounds_error;
  bool packing_error;

  KOKKOS_INLINE_FUNCTION PackCrsMatrixError () :
    bad_index (Tpetra::Details::OrdinalTraits<LO>::invalid ()),
    bad_offset (0),
    bad_num_bytes (0),
    out_of_bounds_error (false),
    packing_error (false)
  {}

  bool success() const
  {
    // Any possible error would result in bad_index being changed from
    // `invalid` to the index associated with the error.
    return bad_index == Tpetra::Details::OrdinalTraits<LO>::invalid();
  }

  std::string summary() const
  {
    std::ostringstream os;
    os << "First bad index: " << bad_index
       << ", first bad offset: " << bad_offset
       << ", first bad number of bytes: " << bad_num_bytes
       << ", out of bounds error?: " << (out_of_bounds_error ? "true" : "false");
    return os.str();
  }
};

template<class NumPacketsPerLidViewType,
         class OffsetsViewType,
         class ExportsViewType,
         class ExportLidsViewType,
         class LocalMatrixType,
         class LocalMapType>
struct PackCrsMatrixFunctor {
  typedef NumPacketsPerLidViewType num_packets_per_lid_view_type;
  typedef OffsetsViewType offsets_view_type;
  typedef ExportsViewType exports_view_type;
  typedef ExportLidsViewType export_lids_view_type;

  typedef typename num_packets_per_lid_view_type::non_const_value_type count_type;
  typedef typename offsets_view_type::non_const_value_type offset_type;
  typedef typename LocalMatrixType::value_type IST;
  typedef typename LocalMatrixType::ordinal_type LO;
  typedef typename LocalMapType::global_ordinal_type GO;
  typedef PackCrsMatrixError<LO> value_type;

  static_assert (std::is_same<LO, typename LocalMatrixType::ordinal_type>::value,
                 "LocalMapType::local_ordinal_type and "
                 "LocalMatrixType::ordinal_type must be the same.");

  num_packets_per_lid_view_type num_packets_per_lid_;
  offsets_view_type offsets_;
  exports_view_type exports_;
  export_lids_view_type export_lids_;
  LocalMatrixType local_matrix_;
  LocalMapType local_col_map_;

  PackCrsMatrixFunctor (const num_packets_per_lid_view_type& num_packets_per_lid,
                        const offsets_view_type& offsets,
                        const exports_view_type& exports,
                        const export_lids_view_type& export_lids,
                        const LocalMatrixType& local_matrix,
                        const LocalMapType& local_col_map) :
    num_packets_per_lid_ (num_packets_per_lid),
    offsets_ (offsets),
    exports_ (exports),
    export_lids_ (export_lids),
    local_matrix_ (local_matrix),
    local_col_map_ (local_col_map)
  {}

  KOKKOS_INLINE_FUNCTION void init(value_type& dst) const
  {
    dst.bad_index = Tpetra::Details::OrdinalTraits<LO>::invalid();
    dst.bad_offset = 0;
    dst.bad_num_bytes = 0;
    dst.out_of_bounds_error = false;
    dst.packing_error = false;
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst, const volatile value_type& src) const
  {
    // The dst object should reflect the first (least) bad index and
    // all other associated error codes and data.  Thus, we need only
    // check if the src object shows an error and if its associated
    // `bad_index` is less than the dst.bad_index (if dst shows
    // errors).
    LO invalid = Tpetra::Details::OrdinalTraits<LO>::invalid();
    if (src.bad_index != invalid) {
      // An error in the src; check if
      //   1. The dst shows errors
      //   2. If dst does show errors, if src bad_index is less than
      //      its bad index
      if (dst.bad_index == invalid || src.bad_index < dst.bad_index) {
        dst.bad_index = src.bad_index;
        dst.bad_offset = src.bad_offset;
        dst.bad_num_bytes = src.bad_num_bytes;
        dst.out_of_bounds_error = src.out_of_bounds_error;
        dst.packing_error = src.packing_error;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const LO i, value_type& dst) const
  {
    const LO export_lid = export_lids_[i];
    const size_t buf_size = exports_.size();
    const size_t num_ent =
      static_cast<size_t> (local_matrix_.graph.row_map[export_lid+1]
                           - local_matrix_.graph.row_map[export_lid]);

    // Only pack this row's data if it has a nonzero number of
    // entries.  We can do this because receiving processes get the
    // number of packets, and will know that zero packets means zero
    // entries.
    if (num_ent != 0) {
      char* const num_ent_beg = exports_.ptr_on_device() + offsets_(i);
      char* const num_ent_end = num_ent_beg + sizeof (LO);
      char* const val_beg = num_ent_end;
      char* const val_end = val_beg + num_ent * sizeof (IST);
      char* const ind_beg = val_end;
      const size_t num_bytes = num_packets_per_lid_(i);

      if ((offsets_(i) > buf_size || offsets_(i) + num_bytes > buf_size)) {
        dst.out_of_bounds_error = true;
      }
      else {
        dst.packing_error = ! packCrsMatrixRow(local_matrix_, local_col_map_,
            num_ent_beg, val_beg, ind_beg, num_ent, export_lid);
      }
      if (dst.out_of_bounds_error || dst.packing_error) {
        dst.bad_index = i;
        dst.bad_offset = offsets_(i);
        dst.bad_num_bytes = num_bytes;
      }
    }
  }
};

/// \brief Packs a single row of the CrsMatrix.
///
/// Data (bytes) describing the row of the CRS matrix are "packed"
/// (concatenated) in to a single char* as
///
///   LO number of entries  |
///   GO column indices      > -- number of entries | column indices | values --
///   SC values             |
///
/// \tparam LocalMatrixType the specialization of the KokkosSparse::CrsMatrix
///   local matrix
/// \tparam LocalMapType the type of the local column map
template<class LocalMatrixType, class LocalMapType>
KOKKOS_FUNCTION bool
packCrsMatrixRow (const LocalMatrixType& lclMatrix,
                  const LocalMapType& lclColMap,
                  char* const numEntOut,
                  char* const valOut,
                  char* const indOut,
                  const size_t numEnt, // number of entries in row
                  const typename LocalMatrixType::ordinal_type lclRow)
{
  using Kokkos::subview;
  typedef LocalMatrixType local_matrix_type;
  typedef LocalMapType local_map_type;
  typedef typename local_matrix_type::value_type IST;
  typedef typename local_matrix_type::ordinal_type LO;
  typedef typename local_matrix_type::size_type offset_type;
  typedef typename local_map_type::global_ordinal_type GO;
  typedef Kokkos::pair<offset_type, offset_type> pair_type;

  const LO numEntLO = static_cast<LO> (numEnt);
  // As of CUDA 6, it's totally fine to use memcpy in a CUDA device
  // function.  It does what one would expect.
  memcpy (numEntOut, &numEntLO, sizeof (LO));

  if (numEnt == 0) {
    return true; // nothing more to pack
  }

#ifdef HAVE_TPETRA_DEBUG
  if (lclRow >= lclMatrix.numRows () ||
      (static_cast<size_t> (lclRow + 1) >=
       static_cast<size_t> (lclMatrix.graph.row_map.dimension_0 ()))) {
#else // NOT HAVE_TPETRA_DEBUG
  if (lclRow >= lclMatrix.numRows ()) {
#endif // HAVE_TPETRA_DEBUG
    // It's bad if this is not a valid local row index.  One thing
    // we can do is just pack the flag invalid value for the column
    // indices.  That makes sure that the receiving process knows
    // something went wrong.
    const GO flagInd = Tpetra::Details::OrdinalTraits<GO>::invalid ();
    for (size_t k = 0; k < numEnt; ++k) {
      // As of CUDA 6, it's totally fine to use memcpy in a CUDA
      // device function.  It does what one would expect.
      memcpy (indOut + k * sizeof (GO), &flagInd, sizeof (GO));
    }
    // The values don't actually matter, but we might as well pack
    // something here.
    const IST zero = Kokkos::ArithTraits<IST>::zero ();
    for (size_t k = 0; k < numEnt; ++k) {
      // As of CUDA 6, it's totally fine to use memcpy in a CUDA
      // device function.  It does what one would expect.
      memcpy (valOut + k * sizeof (IST), &zero, sizeof (IST));
    }
    return false;
  }

  // Since the matrix is locally indexed on the calling process, we
  // have to use its column Map (which it _must_ have in this case)
  // to convert to global indices.
  const offset_type rowBeg = lclMatrix.graph.row_map[lclRow];
  const offset_type rowEnd = lclMatrix.graph.row_map[lclRow + 1];

  auto indIn = subview (lclMatrix.graph.entries, pair_type (rowBeg, rowEnd));
  auto valIn = subview (lclMatrix.values, pair_type (rowBeg, rowEnd));

  // Copy column indices one at a time, so that we don't need
  // temporary storage.
  for (size_t k = 0; k < numEnt; ++k) {
    const GO gblIndIn = lclColMap.getGlobalElement (indIn[k]);
    // As of CUDA 6, it's totally fine to use memcpy in a CUDA
    // device function.  It does what one would expect.
    memcpy (indOut + k * sizeof (GO), &gblIndIn, sizeof (GO));
  }
  // As of CUDA 6, it's totally fine to use memcpy in a CUDA device
  // function.  It does what one would expect.
  memcpy (valOut, valIn.ptr_on_device (), numEnt * sizeof (IST));
  return true;
}




/// \brief Pack specified entries of the given local sparse matrix for
///   communication.
///
/// \warning This is an implementation detail of Tpetra::CrsMatrix.
///
/// \param errStr [in/out] If an error occurs on any participating
///   process, allocate this if it is null, then fill the string with
///   local error reporting.  This is purely local to the process that
///   reports the error.  The caller is responsible for synchronizing
///   across processes.
///
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local matrix.
///
/// \param exportLIDs [in] Local indices of the rows to pack.
///
/// \param lclMatrix [in] The local sparse matrix to pack.
///
/// \return true if no errors occurred on the calling process, else
///   false.  This is purely local to the process that discovered the
///   error.  The caller is responsible for synchronizing across
///   processes.
template<class LocalMatrixType, class LocalMapType>
bool
packCrsMatrix (const LocalMatrixType& lclMatrix,
               const LocalMapType& lclColMap,
               std::unique_ptr<std::string>& errStr,
               Teuchos::Array<char>& exports,
               const Teuchos::ArrayView<size_t>& numPacketsPerLID,
               size_t& constantNumPackets,
               const Teuchos::ArrayView<const typename LocalMatrixType::ordinal_type>& exportLIDs,
               const int myRank,
               Distributor& /* dist */)
{
  using Kokkos::View;
  using Kokkos::HostSpace;
  using Kokkos::MemoryUnmanaged;
  using ::Tpetra::Details::computeOffsetsFromCounts;
  typedef typename LocalMapType::local_ordinal_type LO;
  typedef typename LocalMapType::global_ordinal_type GO;
  typedef typename LocalMatrixType::value_type IST;
  typedef typename LocalMatrixType::device_type device_type;
  typedef typename device_type::execution_space execution_space;
  typedef typename Kokkos::RangePolicy<execution_space, LO> range_type;
  const char prefix[] = "Tpetra::Details::packCrsMatrix: ";

  static_assert (std::is_same<LO, typename LocalMatrixType::ordinal_type>::value,
                 "LocalMapType::local_ordinal_type and "
                 "LocalMatrixType::ordinal_type must be the same.");
  // Setting this to zero tells the caller to expect a possibly
  // different ("nonconstant") number of packets per local index
  // (i.e., a possibly different number of entries per row).
  constantNumPackets = 0;

  const size_t numExportLIDs = static_cast<size_t> (exportLIDs.size ());
  if (numExportLIDs != static_cast<size_t> (numPacketsPerLID.size ())) {
    std::ostringstream os;
    os << prefix << "exportLIDs.size() = " << numExportLIDs
       << " != numPacketsPerLID.size() = " << numPacketsPerLID.size ()
       << "." << std::endl;
    if (errStr.get () == NULL) {
      errStr = std::unique_ptr<std::string> (new std::string (os.str ()));
    }
    else {
      *errStr = *errStr + os.str ();
    }
    return false;
  }
  if (numExportLIDs != 0 && numPacketsPerLID.getRawPtr () == NULL) {
    std::ostringstream os;
    os << prefix << "numExportLIDs = " << numExportLIDs << " != 0, but "
       << "numPacketsPerLID.getRawPtr() = " << numPacketsPerLID.getRawPtr ()
       << " != NULL." << std::endl;
    if (errStr.get () == NULL) {
      errStr = std::unique_ptr<std::string> (new std::string (os.str ()));
    }
    else {
      *errStr = *errStr + os.str ();
    }
    return false;
  }

  if (numExportLIDs == 0) {
    exports.resize (0);
    return true; // nothing to pack
  }

  typename LocalMatrixType::device_type outputDevice;
  using Tpetra::Details::create_mirror_view_from_raw_host_array;
  // This is an output array, so we don't have to copy to device here.
  // However, we'll have to remember to copy back to host when done.
  auto numPktPerLid_d =
    create_mirror_view_from_raw_host_array (outputDevice,
                                            numPacketsPerLID.getRawPtr (),
                                            numPacketsPerLID.size (),
                                            false,
                                            "numPktPerLid");
  // This is an input array, so we have to copy to device here.
  // However, we never need to copy it back to host.
  auto packLids_d =
    create_mirror_view_from_raw_host_array (outputDevice,
                                            exportLIDs.getRawPtr (),
                                            exportLIDs.size (),
                                            true,
                                            "packLids");
  // Array of offsets into the pack buffer.
  Kokkos::View<size_t*, device_type> packOffsets_d ("packOffsets",
                                                    numExportLIDs + 1);
  // Compute number of packets per LID (row to send), as well as
  // corresponding offsets (the prefix sum of the packet counts).
  std::unique_ptr<std::ostringstream> errStrm;
  std::pair<size_t, bool> countResult =
    computeNumPacketsAndOffsets (errStrm,
                                 packOffsets_d,
                                 numPktPerLid_d,
                                 lclMatrix.graph.row_map,
                                 packLids_d,
                                 sizeof (LO),
                                 sizeof (IST),
                                 sizeof (GO));

  if (! countResult.second) {
    if (errStr.get () == NULL) {
      errStr = std::unique_ptr<std::string> (new std::string (errStrm->str ()));
    }
    else {
      *errStr = *errStr + errStrm->str ();
    }
    return false;
  }

  // The counts were an output of computeNumPacketsAndOffsets, so we
  // have to copy them back to host.
  typename decltype (numPktPerLid_d)::HostMirror numPktPerLid_h (numPacketsPerLID.getRawPtr (),
                                                                 numPacketsPerLID.size ());
  Kokkos::deep_copy (numPktPerLid_h, numPktPerLid_d);

  // Resize the output pack buffer if needed.
  if (countResult.first > static_cast<size_t> (exports.size ())) {
    exports.resize (countResult.first);
  }
  // Make device version of output pack buffer.  This is an output
  // array, so we don't have to copy to device here.
  auto packBuf_d =
    create_mirror_view_from_raw_host_array (outputDevice,
                                            exports.getRawPtr (),
                                            countResult.first,
                                            false,
                                            "packBuf");
  // Now do the actual pack!
  typedef PackCrsMatrixFunctor<
    decltype (numPktPerLid_d),
    decltype (packOffsets_d),
    decltype (packBuf_d),
    decltype (packLids_d),
    LocalMatrixType,
    LocalMapType> pack_functor_type;
  pack_functor_type packer (numPktPerLid_d, packOffsets_d, packBuf_d,
                            packLids_d, lclMatrix, lclColMap);
  typename pack_functor_type::value_type result;
  Kokkos::parallel_reduce (range_type (0, numExportLIDs), packer, result);

  if (! result.success ()) {
    std::ostringstream os;
    os << "Proc " << myRank << ": packCrsMatrix failed.  "
       << result.summary()
       << std::endl;
    if (errStr.get () != NULL) {
      errStr = std::unique_ptr<std::string> (new std::string ());
      *errStr = os.str ();
    }
    return false;
  }

  // Copy pack result back to host, if needed.
  typename decltype (packBuf_d)::HostMirror packBuf_h (exports.getRawPtr (),
                                                       countResult.first);
  Kokkos::deep_copy (packBuf_h, packBuf_d);
  return true; // if we got this far, we succeeded
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PACKCRSMATRIX_HPP
