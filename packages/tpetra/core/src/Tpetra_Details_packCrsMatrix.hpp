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
#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
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
///
/// Data (bytes) describing the row of the CRS matrix are "packed"
/// (concatenated) in to a (view of) char* object in the following order:
///
///   1. number of entries (LocalOrdinal)
///   2. global column indices (GlobalOrdinal)
///   3. proces IDs (optional, int)
///   4. row values (Scalar)
///
/// The functions in this file are companions to
/// Tpetra_Details_unpackCrsMatrix.hpp, i.e., Tpetra_Details_unpackCrsMatrix.hpp
/// implements the reverse of the packing order described above to ensure proper
/// unpacking.

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// Forward declaration of Distributor
class Distributor;
#endif // DOXYGEN_SHOULD_SKIP_THIS

//
// Users must never rely on anything in the Details namespace.
//
namespace Details {

/// \brief Compute the number of packets and offsets for the pack procedure
///
/// \tparam OutputOffsetsViewType the type of the output offsets view
/// \tparam CountsViewType the type of the counts view
/// \tparam InputOffsetsViewType the type of the input offsets view
/// \tparam InputLocalRowIndicesViewType the type of the local row indices view
/// \tparam InputLocalRowPidsViewType the type of the local process IDs view
template<class OutputOffsetsViewType,
         class CountsViewType,
         class InputOffsetsViewType,
         class InputLocalRowIndicesViewType,
         class InputLocalRowPidsViewType,
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
  typedef typename InputLocalRowPidsViewType::non_const_value_type local_row_pid_type;
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
                               const InputLocalRowPidsViewType& lclRowPids,
                               const count_type sizeOfLclCount,
                               const count_type sizeOfGblColInd,
                               const count_type sizeOfPid,
                               const count_type sizeOfValue) :
    outputOffsets_ (outputOffsets),
    counts_ (counts),
    rowOffsets_ (rowOffsets),
    lclRowInds_ (lclRowInds),
    lclRowPids_ (lclRowPids),
    sizeOfLclCount_ (sizeOfLclCount),
    sizeOfGblColInd_ (sizeOfGblColInd),
    sizeOfPid_ (sizeOfPid),
    sizeOfValue_ (sizeOfValue),
    error_ ("error") // don't forget this, or you'll get segfaults!
  {
    if (debug) {
      const size_t numRowsToPack = static_cast<size_t> (lclRowInds_.dimension_0 ());

      if (numRowsToPack != static_cast<size_t> (counts_.dimension_0 ())) {
        std::ostringstream os;
        os << "lclRowInds.dimension_0() = " << numRowsToPack
           << " != counts.dimension_0() = " << counts_.dimension_0 ()
           << ".";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str ());
      }
      if (static_cast<size_t> (numRowsToPack + 1) != static_cast<size_t> (outputOffsets_.dimension_0 ())) {
        std::ostringstream os;
        os << "lclRowInds.dimension_0() + 1 = " << (numRowsToPack + 1)
           << " != outputOffsets.dimension_0() = " << outputOffsets_.dimension_0 ()
           << ".";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str ());
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
      // global column indices, then that many pids (if any), then that many
      // values.  However, if the number of entries in the row is zero, we pack
      // nothing.
      const count_type numBytes = (count == 0) ?
        static_cast<count_type> (0) :
        sizeOfLclCount_ + count * (sizeOfGblColInd_ +
                                   (lclRowPids_.size() > 0 ? sizeOfPid_ : 0) +
                                   sizeOfValue_);

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
  typename InputLocalRowPidsViewType::const_type lclRowPids_;
  count_type sizeOfLclCount_;
  count_type sizeOfGblColInd_;
  count_type sizeOfPid_;
  count_type sizeOfValue_;
  Kokkos::View<int, device_type> error_;
};

/// \brief Compute the number of packets and offsets for the pack procedure
///
/// \tparam OutputOffsetsViewType the type of the output offsets view
/// \tparam CountsViewType the type of the counts view
/// \tparam InputOffsetsViewType the type of the input offsets view
/// \tparam InputLocalRowIndicesViewType the type of the local row indices view
/// \tparam InputLocalRowPidsViewType the type of the local process IDs view
///
/// This is the high level interface to the NumPacketsAndOffsetsFunctor functor
template<class OutputOffsetsViewType,
         class CountsViewType,
         class InputOffsetsViewType,
         class InputLocalRowIndicesViewType,
         class InputLocalRowPidsViewType>
typename CountsViewType::non_const_value_type
computeNumPacketsAndOffsets (const OutputOffsetsViewType& outputOffsets,
                             const CountsViewType& counts,
                             const InputOffsetsViewType& rowOffsets,
                             const InputLocalRowIndicesViewType& lclRowInds,
                             const InputLocalRowPidsViewType& lclRowPids,
                             const typename CountsViewType::non_const_value_type sizeOfLclCount,
                             const typename CountsViewType::non_const_value_type sizeOfGblColInd,
                             const typename CountsViewType::non_const_value_type sizeOfPid,
                             const typename CountsViewType::non_const_value_type sizeOfValue)
{
  typedef NumPacketsAndOffsetsFunctor<OutputOffsetsViewType,
    CountsViewType, typename InputOffsetsViewType::const_type,
    typename InputLocalRowIndicesViewType::const_type,
    typename InputLocalRowPidsViewType::const_type> functor_type;
  typedef typename CountsViewType::non_const_value_type count_type;
  typedef typename OutputOffsetsViewType::size_type size_type;
  typedef typename OutputOffsetsViewType::execution_space execution_space;
  typedef typename functor_type::local_row_index_type LO;
  typedef Kokkos::RangePolicy<execution_space, LO> range_type;
  const char prefix[] = "computeNumPacketsAndOffsets: ";


  count_type count = 0;
  const count_type numRowsToPack = lclRowInds.dimension_0 ();

  if (numRowsToPack == 0) {
    return count;
  }
  else {

    TEUCHOS_TEST_FOR_EXCEPTION(
        rowOffsets.dimension_0() <= static_cast<size_type>(1),
        std::invalid_argument,
        prefix << "There is at least one row to pack, but the matrix has no rows.  "
        << "lclRowInds.dimension_0() = " << numRowsToPack << ", but "
        << "rowOffsets.dimension_0() = " << rowOffsets.dimension_0 () << " <= 1.");

    TEUCHOS_TEST_FOR_EXCEPTION(
        outputOffsets.dimension_0() != static_cast<size_type> (numRowsToPack + 1),
        std::invalid_argument,
        prefix << "Output dimension does not match number of rows to pack.  "
        << "outputOffsets.dimension_0() = " << outputOffsets.dimension_0 ()
        << " != lclRowInds.dimension_0() + 1 = "
        << static_cast<size_type> (numRowsToPack + 1) << ".");

    functor_type f (outputOffsets, counts, rowOffsets,
                    lclRowInds, lclRowPids, sizeOfLclCount,
                    sizeOfGblColInd, sizeOfPid, sizeOfValue);
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
    TEUCHOS_TEST_FOR_EXCEPTION(errCode != 0, std::runtime_error,
        prefix << "NumPacketsAndOffsetsFunctor reported error code " << errCode
        << " != 0.");

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
      count = outputOffsets(numRowsToPack);
      return {false, errStr};
    }
#endif // HAVE_TPETRA_DEBUG

    // Get last entry of outputOffsets, which is the sum of the entries
    // of counts.  Don't assume UVM.
    auto outputOffsets_last = Kokkos::subview (outputOffsets, numRowsToPack);
    auto outputOffsets_last_h = Kokkos::create_mirror_view (outputOffsets_last);

    return static_cast<count_type>(outputOffsets_last_h());
  }
}

/// \brief Packs a single row of the CrsMatrix.
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam ColumnMap the type of the local column map
///
/// Data (bytes) describing the row of the CRS matrix are "packed"
/// (concatenated) in to a single (view of) char* in the following order:
///
///   1. LO number of entries
///   2. GO column indices
///   3. int proces IDs
///   4. SC values
///
template<class ST, class ColumnMap>
KOKKOS_FUNCTION
Kokkos::pair<int, size_t>
packCrsMatrixRow(const ColumnMap& col_map,
                 const typename PackTraits<typename ColumnMap::local_ordinal_type, typename ColumnMap::device_type>::output_buffer_type& exports,
                 const typename PackTraits<typename ColumnMap::local_ordinal_type, typename ColumnMap::device_type>::input_array_type& lids_in,
                 const typename PackTraits<int, typename ColumnMap::device_type>::input_array_type& pids_in,
                 const typename PackTraits<ST, typename ColumnMap::device_type>::input_array_type& vals_in,
                 const size_t offset,
                 const size_t num_ent,
                 const size_t num_bytes_per_value)
{
    using Kokkos::subview;
    typedef typename ColumnMap::local_ordinal_type LO;
    typedef typename ColumnMap::global_ordinal_type GO;
    typedef typename ColumnMap::device_type DT;
    typedef Kokkos::pair<int, size_t> return_type;

    // NOTE (mfh 02 Feb 2015) This assumes that output_buffer_type is
    // the same, no matter what type we're packing.  It's a reasonable
    // assumption, given that we go through the trouble of PackTraits
    // just so that we can pack data of different types in the same
    // buffer.
    typedef typename PackTraits<LO, DT>::output_buffer_type output_buffer_type;
    typedef typename output_buffer_type::size_type size_type;
    typedef typename Kokkos::pair<size_type, size_type> pair_type;

    if (num_ent == 0) {
      // Empty rows always take zero bytes, to ensure sparsity.
      return return_type(0, 0);
    }

    const LO num_ent_LO = static_cast<LO>(num_ent); // packValueCount wants this
    const size_t num_ent_beg = offset;
    const size_t num_ent_len = PackTraits<LO, DT>::packValueCount(num_ent_LO);

    const GO gid = 0; // packValueCount wants this
    const size_t gids_beg = num_ent_beg + num_ent_len;
    const size_t gids_len = num_ent * PackTraits<GO, DT>::packValueCount(gid);

    const bool pack_pids = pids_in.size() > 0;
    const int pid = 0; // packValueCount wants this
    const size_t pids_beg = gids_beg + gids_len;
    const size_t pids_len = (pack_pids) ? num_ent * PackTraits<int, DT>::packValueCount(pid) : 0;

    const size_t vals_beg = gids_beg + gids_len + pids_len;
    const size_t vals_len = num_ent * num_bytes_per_value;

    output_buffer_type num_ent_out =
      subview(exports, pair_type(num_ent_beg, num_ent_beg + num_ent_len));
    output_buffer_type gids_out =
      subview(exports, pair_type(gids_beg, gids_beg + gids_len));
    output_buffer_type pids_out;
    if (pack_pids) {
      pids_out = subview(exports, pair_type(pids_beg, pids_beg + pids_len));
    }
    output_buffer_type vals_out =
      subview(exports, pair_type(vals_beg, vals_beg + vals_len));

    size_t num_bytes_out = 0;
    int error_code = 0;
    num_bytes_out += PackTraits<LO, DT>::packValue(num_ent_out, num_ent_LO);

    {
      // Copy column indices one at a time, so that we don't need
      // temporary storage.
      for (size_t k=0; k<num_ent; k++) {
        GO gid = col_map.getGlobalElement(lids_in[k]);
        num_bytes_out += PackTraits<GO,DT>::packValue(gids_out, k, gid);
      }

      // Copy PIDs one at a time, so that we don't need temporary storage.
      if (pack_pids) {
        for (size_t k=0; k<num_ent; k++) {
          int pid = pids_in[lids_in[k]];
          num_bytes_out += PackTraits<int,DT>::packValue(pids_out, k, pid);
        }
      }

      // Copy the values
      auto p = PackTraits<ST, DT>::packArray(vals_out, vals_in, num_ent);
      error_code += p.first;
      num_bytes_out += p.second;
    }

    if (error_code != 0) {
      return return_type(10, num_bytes_out);
    }

    const size_t expected_num_bytes = num_ent_len + gids_len + pids_len + vals_len;
    if (num_bytes_out != expected_num_bytes) {
      return return_type(11, num_bytes_out);
    }
    return return_type(0, num_bytes_out);
}

template<class LocalMatrix, class LocalMap>
struct PackCrsMatrixFunctor {
  typedef LocalMatrix local_matrix_type;
  typedef LocalMap local_map_type;
  typedef typename local_matrix_type::value_type ST;
  typedef typename local_map_type::local_ordinal_type LO;
  typedef typename local_map_type::global_ordinal_type GO;
  typedef typename local_matrix_type::device_type DT;

  typedef typename PackTraits<size_t,DT>::input_array_type num_packets_per_lid_view_type;
  typedef typename PackTraits<size_t,DT>::input_array_type offsets_view_type;
  typedef typename PackTraits<size_t,DT>::output_buffer_type exports_view_type;
  typedef typename PackTraits<LO,DT>::input_array_type export_lids_view_type;
  typedef typename PackTraits<int,DT>::input_array_type source_pids_view_type;

  typedef typename num_packets_per_lid_view_type::non_const_value_type count_type;
  typedef typename offsets_view_type::non_const_value_type offset_type;

  typedef Kokkos::pair<int, LO> value_type;

  static_assert (std::is_same<LO, typename local_matrix_type::ordinal_type>::value,
                 "local_map_type::local_ordinal_type and "
                 "local_matrix_type::ordinal_type must be the same.");

  LO InvalidIndex = OrdinalTraits<LO>::invalid();

  local_matrix_type local_matrix;
  local_map_type local_col_map;
  exports_view_type exports;
  num_packets_per_lid_view_type num_packets_per_lid;
  export_lids_view_type export_lids;
  source_pids_view_type source_pids;
  offsets_view_type offsets;
  size_t num_bytes_per_value;

  PackCrsMatrixFunctor (const local_matrix_type& local_matrix_in,
                        const local_map_type& local_col_map_in,
                        const exports_view_type& exports_in,
                        const num_packets_per_lid_view_type& num_packets_per_lid_in,
                        const export_lids_view_type& export_lids_in,
                        const source_pids_view_type& source_pids_in,
                        const offsets_view_type& offsets_in, 
                        const size_t num_bytes_per_value_in) :
    local_matrix (local_matrix_in),
    local_col_map (local_col_map_in),
    exports (exports_in),
    num_packets_per_lid (num_packets_per_lid_in),
    export_lids (export_lids_in),
    source_pids (source_pids_in),
    offsets (offsets_in),
    num_bytes_per_value (num_bytes_per_value_in)
  {}

  KOKKOS_INLINE_FUNCTION void init(value_type& dst) const
  {
    dst = Kokkos::make_pair(0, InvalidIndex);
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst, const volatile value_type& src) const
  {
    // `dst` should reflect the first (least) bad index and
    // all other associated error codes and data.  Thus, we need only
    // check if the `src` object shows an error and if its associated
    // bad index is less than `dst`'s bad index.
    if (src.second != InvalidIndex) {
      // An error in the src; check if
      //   1. `dst` shows errors
      //   2. If `dst` does show errors, if src's bad index is less than
      //      *this' bad index
      if (dst.second == InvalidIndex || src.second < dst.second) {
        dst = src;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const LO i, value_type& dst) const
  {
    using Kokkos::View;
    typedef typename PackTraits<LO,DT>::output_buffer_type output_buffer_type;
    typedef typename output_buffer_type::size_type size_type;
    typedef typename Kokkos::pair<size_type, size_type> pair_type;

    const size_t offset = offsets(i);
    const LO export_lid = export_lids[i];
    const size_t buf_size = exports.size();
    const size_t num_bytes = num_packets_per_lid(i);
    const size_t num_ent =
      static_cast<size_t> (local_matrix.graph.row_map[export_lid+1]
                         - local_matrix.graph.row_map[export_lid]);

    // Only pack this row's data if it has a nonzero number of
    // entries.  We can do this because receiving processes get the
    // number of packets, and will know that zero packets means zero
    // entries.
    if (num_ent == 0) {
      return;
    }

#ifdef HAVE_TPETRA_DEBUG
    if (export_lid >= local_matrix.numRows() ||
        (static_cast<size_t>(export_lid + 1) >=
         static_cast<size_t>(local_matrix.graph.row_map.dimension_0())))
#else // NOT HAVE_TPETRA_DEBUG
    if (export_lid >= local_matrix.numRows ())
#endif // HAVE_TPETRA_DEBUG
    {
      // Invalid row
      dst = Kokkos::make_pair(1, i);
      return;
    }

    if ((offset > buf_size || offset + num_bytes > buf_size)) {
      // Out of bounds
      dst = Kokkos::make_pair(2, i);
      return;
    }

    // We can now pack this row

    // Since the matrix is locally indexed on the calling process, we
    // have to use its column Map (which it _must_ have in this case)
    // to convert to global indices.
    const size_type row_beg = local_matrix.graph.row_map[export_lid];
    const size_type row_end = local_matrix.graph.row_map[export_lid + 1];

    auto vals_in = subview(local_matrix.values, pair_type(row_beg, row_end));
    auto lids_in = subview(local_matrix.graph.entries, pair_type(row_beg, row_end));

    auto p = packCrsMatrixRow<ST,local_map_type>(local_col_map, exports,
        lids_in, source_pids, vals_in, offset, num_ent, num_bytes_per_value);

    int error_code_this_row = p.first;
    size_t num_bytes_packed_this_row = p.second;
    if (error_code_this_row != 0) {
      // Bad pack
      dst = Kokkos::make_pair(error_code_this_row, i);
    } else if (num_bytes_packed_this_row != num_bytes) {
      dst = Kokkos::make_pair(3, i);
    }
  }
};

/// \brief Perform the pack operation for the matrix
///
/// \tparam LocalMatrix the specialization of the KokkosSparse::CrsMatrix
///   local matrix
/// \tparam LocalMap the type of the local column map
///
/// This is a higher level interface to the PackCrsMatrixFunctor
template<class LocalMatrix, class LocalMap>
void
do_pack(const LocalMatrix& local_matrix,
        const LocalMap& local_map,
        const typename PackTraits<size_t,typename LocalMatrix::device_type>::output_buffer_type& exports,
        const typename PackTraits<size_t,typename LocalMatrix::device_type>::input_array_type& num_packets_per_lid,
        const typename PackTraits<typename LocalMap::local_ordinal_type,typename LocalMatrix::device_type>::input_array_type& export_lids,
        const typename PackTraits<int,typename LocalMatrix::device_type>::input_array_type& source_pids,
        const typename PackTraits<size_t,typename LocalMatrix::device_type>::input_array_type& offsets,
        const size_t num_bytes_per_value)
{
  typedef typename LocalMap::local_ordinal_type LO;
  typedef typename LocalMatrix::device_type DT;
  typedef Kokkos::RangePolicy<typename DT::execution_space, LO> range_type;

  const char prefix[] = "Tpetra::Details::do_pack: ";

  // Now do the actual pack!
  typedef PackCrsMatrixFunctor<LocalMatrix, LocalMap> pack_functor_type;
  pack_functor_type f(local_matrix, local_map, exports, num_packets_per_lid,
      export_lids, source_pids, offsets, num_bytes_per_value);

  typename pack_functor_type::value_type x;
  Kokkos::parallel_reduce(range_type(0, num_packets_per_lid.dimension_0()), f, x);
  auto x_h = x.to_std_pair();
  TEUCHOS_TEST_FOR_EXCEPTION(x_h.first != 0, std::runtime_error,
      prefix << "PackCrsMatrixFunctor reported error code "
             << x_h.first << " for the first bad row " << x.second);

  return;
}

/// \brief Pack specified entries of the given local sparse matrix for
///   communication.
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam NT The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \warning This is an implementation detail of Tpetra::CrsMatrix.
///
/// \param sourceMatrix [in] the CrsMatrix source
///
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param num_packets_per_lid [out] Entry k gives the number of bytes
///   packed for row export_lids[k] of the local matrix.
///
/// \param export_lids [in] Local indices of the rows to pack.
///
/// \param export_pids [in] Process IDs for each row
///
/// \param constant_num_packets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
template<typename ST, typename LO, typename GO, typename NT>
void
packCrsMatrixImpl(const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
                  Kokkos::DualView<char*, typename NT::device_type>& exports,
                  const Kokkos::View<size_t*, typename NT::device_type>& num_packets_per_lid,
                  const Kokkos::View<const LO*, typename NT::device_type>& export_lids,
                  const Kokkos::View<const int*, typename NT::device_type>& export_pids,
                  size_t& constant_num_packets,
                  Distributor& /* dist */)
{

  using Kokkos::View;

  typedef typename NT::device_type DT;
  typedef typename DT::execution_space execution_space;
  typedef typename Kokkos::DualView<char*, DT> exports_view_type;
  typedef typename exports_view_type::t_dev::memory_space dev_memory_space;

  const char prefix[] = "Tpetra::Details::packCrsMatrixImpl: ";

  auto local_matrix = sourceMatrix.getLocalMatrix();
  auto local_col_map = sourceMatrix.getColMap()->getLocalMap();

  // Setting this to zero tells the caller to expect a possibly
  // different ("nonconstant") number of packets per local index
  // (i.e., a possibly different number of entries per row).
  constant_num_packets = 0;

  const size_t num_export_lids = static_cast<size_t>(export_lids.dimension_0());
  if (num_export_lids != static_cast<size_t>(num_packets_per_lid.dimension_0())) {
    std::ostringstream os;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        prefix << "num_export_lids.dimension_0() = " << num_export_lids
               << " != num_packets_per_lid.dimension_0() = " << num_packets_per_lid.dimension_0()
               << ".");
  }

  if (num_export_lids != 0 && num_packets_per_lid.ptr_on_device() == NULL) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        prefix << "num_export_lids = " << num_export_lids << " != 0, but "
               << "num_packets_per_lid.ptr_on_device() = " << num_packets_per_lid.ptr_on_device()
               << " != NULL.");
  }

  LO lid = 0;
  const size_t num_bytes_per_lid = PackTraits<LO,DT>::packValueCount(lid);

  GO gid = 0;
  const size_t num_bytes_per_gid = PackTraits<GO,DT>::packValueCount(gid);

  int pid = 0;
  const size_t num_bytes_per_pid = PackTraits<int,DT>::packValueCount(pid);

  size_t num_bytes_per_value = 0;
  if (PackTraits<ST,DT>::compileTimeSize) {
    ST val; // assume that ST is default constructible
    num_bytes_per_value = PackTraits<ST,DT>::packValueCount(val);
  }
  else {
    // Since the packed data come from the source matrix, we can use the source
    // matrix to get the number of bytes per Scalar value stored in the matrix.
    // This assumes that all Scalar values in the source matrix require the same
    // number of bytes.  If the source matrix has no entries on the calling
    // process, then we hope that some process does have some idea how big
    // a Scalar value is.  Of course, if no processes have any entries, then no
    // values should be packed (though this does assume that in our packing
    // scheme, rows with zero entries take zero bytes).
    size_t num_bytes_per_value_l = 0;
    if (local_matrix.values.dimension_0() > 0) {
      const ST& val = local_matrix.values(0);
      num_bytes_per_value_l = PackTraits<ST,DT>::packValueCount(val);
    }
    Teuchos::reduceAll<int, size_t>(*(sourceMatrix.getComm()),
                                    Teuchos::REDUCE_MAX,
                                    num_bytes_per_value_l,
                                    Teuchos::outArg(num_bytes_per_value));
  }

  if (num_export_lids == 0) {
    // FIXME (26 Apr 2016) Fences around (UVM) allocations only
    // temporarily needed for #227 debugging.  Should be able to
    // remove them after that's fixed.
    execution_space::fence ();
    exports = Kokkos::DualView<char*, DT>("exports", 0);
    execution_space::fence ();
    return;
  }

  // Array of offsets into the pack buffer.
  Kokkos::View<size_t*, DT> offsets("offsets", num_export_lids + 1);

  // Compute number of packets per LID (row to send), as well as
  // corresponding offsets (the prefix sum of the packet counts).
  size_t count =
    computeNumPacketsAndOffsets (offsets, num_packets_per_lid,
                                 local_matrix.graph.row_map, export_lids, export_pids,
                                 num_bytes_per_lid, num_bytes_per_gid,
                                 num_bytes_per_pid, num_bytes_per_value);

  // Resize the output pack buffer if needed.
  if (count > static_cast<size_t>(exports.dimension_0())) {
    // FIXME (26 Apr 2016) Fences around (UVM) allocations only
    // temporarily needed for #227 debugging.  Should be able to
    // remove them after that's fixed.
    execution_space::fence ();
    exports = Kokkos::DualView<char*, DT>("exports", count);
    execution_space::fence ();
  }

  exports.template modify<dev_memory_space>();

  do_pack(local_matrix, local_col_map,
          exports.template view<dev_memory_space>(), num_packets_per_lid,
          export_lids, export_pids, offsets, num_bytes_per_value);

  // If we got this far, we succeeded. Copy pack result back to host, if needed.
  exports.template sync<Kokkos::HostSpace>();

  return;
}

/// \brief Pack specified entries of the given local sparse matrix for
///   communication.
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam NT The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \param sourceMatrix [in] the CrsMatrix source
///
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local matrix.
///
/// \param exportLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
///
/// \param distor [in] The distributor (not used)
///
/// This is the public interface to the pack machinery
/// converts passed Teuchos::ArrayView objects to Kokkos::View objects (and
/// copies back in to the Teuchos::ArrayView objects, if needed).  When
/// CrsMatrix migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<typename ST, typename LO, typename GO, typename NT>
void
packCrsMatrix (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
               Teuchos::Array<char>& exports,
               const Teuchos::ArrayView<size_t>& numPacketsPerLID,
               const Teuchos::ArrayView<const LO>& exportLIDs,
               size_t& constantNumPackets,
               Distributor& distor)
{
  typedef typename CrsMatrix<ST,LO,GO,NT>::local_matrix_type local_matrix_type;
  typedef typename local_matrix_type::device_type device_type;
  typename local_matrix_type::device_type outputDevice;

  // Convert all Teuchos::Array to Kokkos::View

  // This is an output array, so we don't have to copy to device here.
  // However, we'll have to remember to copy back to host when done.
  auto num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array(outputDevice, numPacketsPerLID.getRawPtr(),
                                           numPacketsPerLID.size(), false,
                                           "num_packets_per_lid");

  // This is an input array, so we have to copy to device here.
  // However, we never need to copy it back to host.
  auto export_lids_d =
    create_mirror_view_from_raw_host_array(outputDevice, exportLIDs.getRawPtr(),
                                           exportLIDs.size(), true,
                                           "export_lids");

  // Create an empty array of PIDs
  Kokkos::View<int*, device_type> export_pids_d("export_pids", 0);

  // Dual view of Exports
  typedef typename Kokkos::DualView<char*, device_type> exports_view_type;
  typedef typename exports_view_type::t_dev::memory_space dev_memory_space;
  Kokkos::DualView<char*, device_type> exports_dv("exports", 0);

  packCrsMatrixImpl<ST,LO,GO,NT>(sourceMatrix, exports_dv,
                                 num_packets_per_lid_d,
                                 export_lids_d, export_pids_d,
                                 constantNumPackets, distor);

  // The counts are an output of packCrsMatrixImpl, so we
  // have to copy them back to host.
  typename decltype(num_packets_per_lid_d)::HostMirror num_packets_per_lid_h(
      numPacketsPerLID.getRawPtr(), numPacketsPerLID.size());
  Kokkos::deep_copy(num_packets_per_lid_h, num_packets_per_lid_d);

  // The exports are an output of packCrsMatrixImpl, so we
  // have to copy them back to host.
  auto exports_d = exports_dv.template view<dev_memory_space>();
  exports.resize(exports_d.dimension_0());
  typename decltype(exports_d)::HostMirror exports_h(
      exports.getRawPtr(), exports.size());
  Kokkos::deep_copy(exports_h, exports_d);

}

/// \brief Pack specified entries of the given local sparse matrix for
///   communication.
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam NT The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \param sourceMatrix [in] the CrsMatrix source
///
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local matrix.
///
/// \param exportLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
///
/// \param distor [in] The distributor (not used)
///
/// This is the public interface to the pack machinery
/// converts passed Teuchos::ArrayView objects to Kokkos::View objects (and
/// copies back in to the Teuchos::ArrayView objects, if needed).  When
/// CrsMatrix migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<typename ST, typename LO, typename GO, typename NT>
void
packCrsMatrixWithOwningPIDs (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
                             Kokkos::DualView<char*, typename NT::device_type>& exports_dv,
                             const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                             const Teuchos::ArrayView<const LO>& exportLIDs,
                             const Teuchos::ArrayView<const int>& sourcePIDs,
                             size_t& constantNumPackets,
                             Distributor &distor)
{
  typedef typename CrsMatrix<ST,LO,GO,NT>::local_matrix_type local_matrix_type;
  typename local_matrix_type::device_type outputDevice;

  // Convert all Teuchos::Array to Kokkos::View

  // This is an output array, so we don't have to copy to device here.
  // However, we'll have to remember to copy back to host when done.
  auto num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array(outputDevice, numPacketsPerLID.getRawPtr(),
                                           numPacketsPerLID.size(), false,
                                           "num_packets_per_lid");

  // This is an input array, so we have to copy to device here.
  // However, we never need to copy it back to host.
  auto export_lids_d =
    create_mirror_view_from_raw_host_array(outputDevice, exportLIDs.getRawPtr(),
                                           exportLIDs.size(), true,
                                           "export_lids");

  // This is an input array, so we have to copy to device here.
  // However, we never need to copy it back to host.
  auto export_pids_d =
    create_mirror_view_from_raw_host_array(outputDevice, sourcePIDs.getRawPtr(),
                                           sourcePIDs.size(), true,
                                           "export_pids");

  packCrsMatrixImpl<ST,LO,GO,NT>(sourceMatrix, exports_dv,
                                 num_packets_per_lid_d,
                                 export_lids_d, export_pids_d,
                                 constantNumPackets, distor);

  // The counts are an output of packCrsMatrixImpl, so we
  // have to copy them back to host.
  typename decltype(num_packets_per_lid_d)::HostMirror num_packets_per_lid_h(
      numPacketsPerLID.getRawPtr(), numPacketsPerLID.size());
  Kokkos::deep_copy(num_packets_per_lid_h, num_packets_per_lid_d);

}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PACKCRSMATRIX_HPP
