// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_PACKCRSMATRIX_DEF_HPP
#define TPETRA_DETAILS_PACKCRSMATRIX_DEF_HPP

#include "TpetraCore_config.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_castAwayConstDualView.hpp"
#include "Tpetra_Details_createMirrorView.hpp"
#include "Tpetra_Details_getEntryOnHost.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_Details_Profiling.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

/// \file Tpetra_Details_packCrsMatrix_def.hpp
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
/// Tpetra_Details_unpackCrsMatrixAndCombine.hpp, i.e.,
/// Tpetra_Details_unpackCrsMatrixAndCombine.hpp implements the
/// reverse of the packing order described above to ensure proper
/// unpacking.

namespace Tpetra {

//
// Users must never rely on anything in the Details namespace.
//
namespace Details {

namespace PackCrsMatrixImpl {
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
  static_assert (Kokkos::is_view<OutputOffsetsViewType>::value,
                 "OutputOffsetsViewType must be a Kokkos::View.");
  static_assert (std::is_same<typename OutputOffsetsViewType::value_type, output_offset_type>::value,
                 "OutputOffsetsViewType must be a nonconst Kokkos::View.");
  static_assert (std::is_integral<output_offset_type>::value,
                 "The type of each entry of OutputOffsetsViewType must be a built-in integer type.");
  static_assert (Kokkos::is_view<CountsViewType>::value,
                 "CountsViewType must be a Kokkos::View.");
  static_assert (std::is_same<typename CountsViewType::value_type, output_offset_type>::value,
                 "CountsViewType must be a nonconst Kokkos::View.");
  static_assert (std::is_integral<count_type>::value,
                 "The type of each entry of CountsViewType must be a built-in integer type.");
  static_assert (Kokkos::is_view<InputOffsetsViewType>::value,
                 "InputOffsetsViewType must be a Kokkos::View.");
  static_assert (std::is_integral<input_offset_type>::value,
                 "The type of each entry of InputOffsetsViewType must be a built-in integer type.");
  static_assert (Kokkos::is_view<InputLocalRowIndicesViewType>::value,
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
      const size_t numRowsToPack = static_cast<size_t> (lclRowInds_.extent (0));

      if (numRowsToPack != static_cast<size_t> (counts_.extent (0))) {
        std::ostringstream os;
        os << "lclRowInds.extent(0) = " << numRowsToPack
           << " != counts.extent(0) = " << counts_.extent (0)
           << ".";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, os.str ());
      }
      if (static_cast<size_t> (numRowsToPack + 1) !=
          static_cast<size_t> (outputOffsets_.extent (0))) {
        std::ostringstream os;
        os << "lclRowInds.extent(0) + 1 = " << (numRowsToPack + 1)
           << " != outputOffsets.extent(0) = " << outputOffsets_.extent (0)
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
        if (curInd >= static_cast<local_row_index_type> (outputOffsets_.extent (0))) {
          error_ () = 2;
          return;
        }
      }
      outputOffsets_(curInd) = update;
    }

    if (curInd < static_cast<local_row_index_type> (counts_.extent (0))) {
      const auto lclRow = lclRowInds_(curInd);
      if (static_cast<size_t> (lclRow + 1) >= static_cast<size_t> (rowOffsets_.extent (0)) ||
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

      // We pack first the number of entries in the row, then that
      // many global column indices, then that many pids (if any),
      // then that many values.  However, if the number of entries in
      // the row is zero, we pack nothing.
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
    // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
    // Note: In the UVM case, this would otherwise be a no-op
    // and thus not fence, so the value might not be correct on return
    // In the non-UVM case, create_mirror_view will block for the allocation
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
  const count_type numRowsToPack = lclRowInds.extent (0);

  if (numRowsToPack == 0) {
    return count;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION
      (rowOffsets.extent (0) <= static_cast<size_type> (1),
       std::invalid_argument, prefix << "There is at least one row to pack, "
       "but the matrix has no rows.  lclRowInds.extent(0) = " <<
       numRowsToPack << ", but rowOffsets.extent(0) = " <<
       rowOffsets.extent (0) << " <= 1.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (outputOffsets.extent (0) !=
       static_cast<size_type> (numRowsToPack + 1), std::invalid_argument,
       prefix << "Output dimension does not match number of rows to pack.  "
       << "outputOffsets.extent(0) = " << outputOffsets.extent (0)
       << " != lclRowInds.extent(0) + 1 = "
       << static_cast<size_type> (numRowsToPack + 1) << ".");
    TEUCHOS_TEST_FOR_EXCEPTION
      (counts.extent (0) != numRowsToPack, std::invalid_argument,
       prefix << "counts.extent(0) = " << counts.extent (0)
       << " != numRowsToPack = " << numRowsToPack << ".");

    functor_type f (outputOffsets, counts, rowOffsets,
                    lclRowInds, lclRowPids, sizeOfLclCount,
                    sizeOfGblColInd, sizeOfPid, sizeOfValue);
    Kokkos::parallel_scan (range_type (0, numRowsToPack + 1), f);

    // At least in debug mode, this functor checks for errors.
    const int errCode = f.getError ();
    TEUCHOS_TEST_FOR_EXCEPTION
      (errCode != 0, std::runtime_error, prefix << "parallel_scan error code "
       << errCode << " != 0.");

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
    using Tpetra::Details::getEntryOnHost;
    return static_cast<count_type> (getEntryOnHost (outputOffsets,
                                                    numRowsToPack));
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
template<class ST, class ColumnMap, class BufferDeviceType>
KOKKOS_FUNCTION
Kokkos::pair<int, size_t>
packCrsMatrixRow (const ColumnMap& col_map,
                  const Kokkos::View<char*, BufferDeviceType>& exports,
                  const typename PackTraits<typename ColumnMap::local_ordinal_type>::input_array_type& lids_in,
                  const typename PackTraits<int>::input_array_type& pids_in,
                  const typename PackTraits<ST>::input_array_type& vals_in,
                  const size_t offset,
                  const size_t num_ent,
                  const size_t num_bytes_per_value,
                  const bool pack_pids)
{
  using Kokkos::subview;
  using LO = typename ColumnMap::local_ordinal_type;
  using GO = typename ColumnMap::global_ordinal_type;
  using return_type = Kokkos::pair<int, size_t>;

  if (num_ent == 0) {
    // Empty rows always take zero bytes, to ensure sparsity.
    return return_type (0, 0);
  }

  const LO num_ent_LO = static_cast<LO> (num_ent); // packValueCount wants this
  const size_t num_ent_beg = offset;
  const size_t num_ent_len = PackTraits<LO>::packValueCount (num_ent_LO);

  const size_t gids_beg = num_ent_beg + num_ent_len;
  const size_t gids_len = num_ent * PackTraits<GO>::packValueCount (GO (0));

  const size_t pids_beg = gids_beg + gids_len;
  const size_t pids_len = pack_pids ?
    num_ent * PackTraits<int>::packValueCount (int (0)) :
    static_cast<size_t> (0);

  const size_t vals_beg = gids_beg + gids_len + pids_len;
  const size_t vals_len = num_ent * num_bytes_per_value;

  char* const num_ent_out = exports.data () + num_ent_beg;
  char* const gids_out = exports.data () + gids_beg;
  char* const pids_out = pack_pids ? exports.data () + pids_beg : NULL;
  char* const vals_out = exports.data () + vals_beg;

  size_t num_bytes_out = 0;
  int error_code = 0;
  num_bytes_out += PackTraits<LO>::packValue (num_ent_out, num_ent_LO);

  {
    // Copy column indices one at a time, so that we don't need
    // temporary storage.
    for (size_t k = 0; k < num_ent; ++k) {
      const LO lid = lids_in[k];
      const GO gid = col_map.getGlobalElement (lid);
      num_bytes_out += PackTraits<GO>::packValue (gids_out, k, gid);
    }
    // Copy PIDs one at a time, so that we don't need temporary storage.
    if (pack_pids) {
      for (size_t k = 0; k < num_ent; ++k) {
        const LO lid = lids_in[k];
        const int pid = pids_in[lid];
        num_bytes_out += PackTraits<int>::packValue (pids_out, k, pid);
      }
    }
    const auto p =
      PackTraits<ST>::packArray (vals_out, vals_in.data (), num_ent);
    error_code += p.first;
    num_bytes_out += p.second;
  }

  if (error_code != 0) {
    return return_type (10, num_bytes_out);
  }

  const size_t expected_num_bytes =
    num_ent_len + gids_len + pids_len + vals_len;
  if (num_bytes_out != expected_num_bytes) {
    return return_type (11, num_bytes_out);
  }
  return return_type (0, num_bytes_out);
}

template<class LocalMatrix, class LocalMap, class BufferDeviceType>
struct PackCrsMatrixFunctor {
  typedef LocalMatrix local_matrix_device_type;
  typedef LocalMap local_map_type;
  typedef typename local_matrix_device_type::value_type ST;
  typedef typename local_map_type::local_ordinal_type LO;
  typedef typename local_map_type::global_ordinal_type GO;
  typedef typename local_matrix_device_type::device_type DT;

  typedef Kokkos::View<const size_t*, BufferDeviceType>
    num_packets_per_lid_view_type;
  typedef Kokkos::View<const size_t*, BufferDeviceType> offsets_view_type;
  typedef Kokkos::View<char*, BufferDeviceType> exports_view_type;
  using export_lids_view_type = typename PackTraits<LO>::input_array_type;
  using source_pids_view_type = typename PackTraits<int>::input_array_type;

  typedef typename num_packets_per_lid_view_type::non_const_value_type
    count_type;
  typedef typename offsets_view_type::non_const_value_type
    offset_type;
  typedef Kokkos::pair<int, LO> value_type;

  static_assert (std::is_same<LO, typename local_matrix_device_type::ordinal_type>::value,
                 "local_map_type::local_ordinal_type and "
                 "local_matrix_device_type::ordinal_type must be the same.");

  local_matrix_device_type local_matrix;
  local_map_type local_col_map;
  exports_view_type exports;
  num_packets_per_lid_view_type num_packets_per_lid;
  export_lids_view_type export_lids;
  source_pids_view_type source_pids;
  offsets_view_type offsets;
  size_t num_bytes_per_value;
  bool pack_pids;

  PackCrsMatrixFunctor (const local_matrix_device_type& local_matrix_in,
                        const local_map_type& local_col_map_in,
                        const exports_view_type& exports_in,
                        const num_packets_per_lid_view_type& num_packets_per_lid_in,
                        const export_lids_view_type& export_lids_in,
                        const source_pids_view_type& source_pids_in,
                        const offsets_view_type& offsets_in,
                        const size_t num_bytes_per_value_in,
                        const bool pack_pids_in) :
    local_matrix (local_matrix_in),
    local_col_map (local_col_map_in),
    exports (exports_in),
    num_packets_per_lid (num_packets_per_lid_in),
    export_lids (export_lids_in),
    source_pids (source_pids_in),
    offsets (offsets_in),
    num_bytes_per_value (num_bytes_per_value_in),
    pack_pids (pack_pids_in)
  {
    const LO numRows = local_matrix_in.numRows ();
    const LO rowMapDim =
      static_cast<LO> (local_matrix.graph.row_map.extent (0));
    TEUCHOS_TEST_FOR_EXCEPTION
      (numRows != 0 && rowMapDim != numRows + static_cast<LO> (1),
       std::logic_error, "local_matrix.graph.row_map.extent(0) = "
       << rowMapDim << " != numRows (= " << numRows << " ) + 1.");
  }

  KOKKOS_INLINE_FUNCTION void init (value_type& dst) const
  {
    using ::Tpetra::Details::OrdinalTraits;
    dst = Kokkos::make_pair (0, OrdinalTraits<LO>::invalid ());
  }

  KOKKOS_INLINE_FUNCTION void
  join (value_type& dst, const value_type& src) const
  {
    // `dst` should reflect the first (least) bad index and all other
    // associated error codes and data, so prefer keeping it.
    if (src.first != 0 && dst.first == 0) {
      dst = src;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const LO i, value_type& dst) const
  {
    const size_t offset = offsets[i];
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

    if (export_lid >= local_matrix.numRows ()) {
      if (dst.first != 0) { // keep only the first error
        dst = Kokkos::make_pair (1, i); // invalid row
      }
      return;
    }
    else if ((offset > buf_size || offset + num_bytes > buf_size)) {
      if (dst.first != 0) { // keep only the first error
        dst = Kokkos::make_pair (2, i); // out of bounds
      }
      return;
    }

    // We can now pack this row

    // Since the matrix is locally indexed on the calling process, we
    // have to use its column Map (which it _must_ have in this case)
    // to convert to global indices.
    const auto row_beg = local_matrix.graph.row_map[export_lid];
    const auto row_end = local_matrix.graph.row_map[export_lid + 1];
    auto vals_in = subview (local_matrix.values,
                            Kokkos::make_pair (row_beg, row_end));
    auto lids_in = subview (local_matrix.graph.entries,
                            Kokkos::make_pair (row_beg, row_end));
    typedef local_map_type LMT;
    typedef BufferDeviceType BDT;
    auto p = packCrsMatrixRow<ST, LMT, BDT> (local_col_map, exports, lids_in,
                                             source_pids, vals_in, offset,
                                             num_ent, num_bytes_per_value,
                                             pack_pids);
    int error_code_this_row = p.first;
    size_t num_bytes_packed_this_row = p.second;
    if (error_code_this_row != 0) {
      if (dst.first != 0) { // keep only the first error
        dst = Kokkos::make_pair (error_code_this_row, i); // bad pack
      }
    }
    else if (num_bytes_packed_this_row != num_bytes) {
      if (dst.first != 0) { // keep only the first error
        dst = Kokkos::make_pair (3, i);
      }
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
template<class LocalMatrix, class LocalMap, class BufferDeviceType>
void
do_pack (const LocalMatrix& local_matrix,
         const LocalMap& local_map,
         const Kokkos::View<char*, BufferDeviceType>& exports,
         const typename PackTraits<size_t>::input_array_type& num_packets_per_lid,
         const typename PackTraits<typename LocalMap::local_ordinal_type>::input_array_type& export_lids,
         const typename PackTraits<int>::input_array_type& source_pids,
         const Kokkos::View<const size_t*, BufferDeviceType>& offsets,
         const size_t num_bytes_per_value,
         const bool pack_pids)
{
  using LO = typename LocalMap::local_ordinal_type;
  using DT = typename LocalMatrix::device_type;
  using range_type = Kokkos::RangePolicy<typename DT::execution_space, LO>;
  const char prefix[] = "Tpetra::Details::do_pack: ";

  if (export_lids.extent (0) != 0) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (static_cast<size_t> (offsets.extent (0)) !=
       static_cast<size_t> (export_lids.extent (0) + 1),
       std::invalid_argument, prefix << "offsets.extent(0) = "
       << offsets.extent (0) << " != export_lids.extent(0) (= "
       << export_lids.extent (0) << ") + 1.");
    TEUCHOS_TEST_FOR_EXCEPTION
      (export_lids.extent (0) != num_packets_per_lid.extent (0),
       std::invalid_argument, prefix << "export_lids.extent(0) = " <<
       export_lids.extent (0) << " != num_packets_per_lid.extent(0) = "
       << num_packets_per_lid.extent (0) << ".");
    // If exports has nonzero length at this point, then the matrix
    // has at least one entry to pack.  Thus, if packing process
    // ranks, we had better have at least one process rank to pack.
    TEUCHOS_TEST_FOR_EXCEPTION
      (pack_pids && exports.extent (0) != 0 &&
       source_pids.extent (0) == 0, std::invalid_argument, prefix <<
       "pack_pids is true, and exports.extent(0) = " <<
       exports.extent (0)  << " != 0, meaning that we need to pack at "
       "least one matrix entry, but source_pids.extent(0) = 0.");
  }

  using pack_functor_type =
    PackCrsMatrixFunctor<LocalMatrix, LocalMap, BufferDeviceType>;
  pack_functor_type f (local_matrix, local_map, exports,
                       num_packets_per_lid, export_lids,
                       source_pids, offsets, num_bytes_per_value,
                       pack_pids);

  typename pack_functor_type::value_type result;
  range_type range (0, num_packets_per_lid.extent (0));
  Kokkos::parallel_reduce (range, f, result);

  if (result.first != 0) {
    // We can't deep_copy from AnonymousSpace Views, so we can't print
    // out any information from them in case of error.
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::runtime_error, prefix << "PackCrsMatrixFunctor "
       "reported error code " << result.first << " for the first "
       "bad row " << result.second << ".");
  }
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
/// \param export_pids [in] Process ranks for the column indices in each packed row.
///
/// \param constant_num_packets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
template<typename ST, typename LO, typename GO, typename NT, typename BufferDeviceType>
void
packCrsMatrix (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
               Kokkos::DualView<char*, BufferDeviceType>& exports,
               const Kokkos::View<size_t*, BufferDeviceType>& num_packets_per_lid,
               const Kokkos::View<const LO*, BufferDeviceType>& export_lids,
               const Kokkos::View<const int*, typename NT::device_type>& export_pids,
               size_t& constant_num_packets,
               const bool pack_pids)
{
  ::Tpetra::Details::ProfilingRegion region_pack_crs_matrix(
    "Tpetra::Details::PackCrsMatrixImpl::packCrsMatrix",
    "Import/Export"
  );
  using Kokkos::View;
  typedef BufferDeviceType DT;
  typedef Kokkos::DualView<char*, BufferDeviceType> exports_view_type;
  const char prefix[] = "Tpetra::Details::PackCrsMatrixImpl::packCrsMatrix: ";
  constexpr bool debug = false;

  auto local_matrix = sourceMatrix.getLocalMatrixDevice ();
  auto local_col_map = sourceMatrix.getColMap ()->getLocalMap ();

  // Setting this to zero tells the caller to expect a possibly
  // different ("nonconstant") number of packets per local index
  // (i.e., a possibly different number of entries per row).
  constant_num_packets = 0;

  const size_t num_export_lids =
    static_cast<size_t> (export_lids.extent (0));
  TEUCHOS_TEST_FOR_EXCEPTION
    (num_export_lids !=
     static_cast<size_t> (num_packets_per_lid.extent (0)),
     std::invalid_argument, prefix << "num_export_lids.extent(0) = "
     << num_export_lids << " != num_packets_per_lid.extent(0) = "
     << num_packets_per_lid.extent (0) << ".");
  if (num_export_lids != 0) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (num_packets_per_lid.data () == NULL, std::invalid_argument,
       prefix << "num_export_lids = "<< num_export_lids << " != 0, but "
       "num_packets_per_lid.data() = "
       << num_packets_per_lid.data () << " == NULL.");
  }

  const size_t num_bytes_per_lid = PackTraits<LO>::packValueCount (LO (0));
  const size_t num_bytes_per_gid = PackTraits<GO>::packValueCount (GO (0));
  const size_t num_bytes_per_pid = PackTraits<int>::packValueCount (int (0));

  size_t num_bytes_per_value = 0;
  if (PackTraits<ST>::compileTimeSize) {
    // Assume ST is default constructible; packValueCount wants an instance.
    num_bytes_per_value = PackTraits<ST>::packValueCount (ST ());
  }
  else {
    // Since the packed data come from the source matrix, we can use
    // the source matrix to get the number of bytes per Scalar value
    // stored in the matrix.  This assumes that all Scalar values in
    // the source matrix require the same number of bytes.  If the
    // source matrix has no entries on the calling process, then we
    // hope that some process does have some idea how big a Scalar
    // value is.  Of course, if no processes have any entries, then no
    // values should be packed (though this does assume that in our
    // packing scheme, rows with zero entries take zero bytes).
    size_t num_bytes_per_value_l = 0;
    if (local_matrix.values.extent(0) > 0) {
      const ST& val = local_matrix.values(0);
      num_bytes_per_value_l = PackTraits<ST>::packValueCount (val);
    }
    using Teuchos::reduceAll;
    reduceAll<int, size_t> (* (sourceMatrix.getComm ()),
                            Teuchos::REDUCE_MAX,
                            num_bytes_per_value_l,
                            Teuchos::outArg (num_bytes_per_value));
  }

  if (num_export_lids == 0) {
    exports = exports_view_type ("exports", 0);
    return;
  }

  // Array of offsets into the pack buffer.
  Kokkos::View<size_t*, DT> offsets ("offsets", num_export_lids + 1);

  // Compute number of packets per LID (row to send), as well as
  // corresponding offsets (the prefix sum of the packet counts).
  const size_t count =
    computeNumPacketsAndOffsets (offsets, num_packets_per_lid,
                                 local_matrix.graph.row_map, export_lids,
                                 export_pids,
                                 num_bytes_per_lid, num_bytes_per_gid,
                                 num_bytes_per_pid, num_bytes_per_value);

  // Resize the output pack buffer if needed.
  if (count > static_cast<size_t> (exports.extent (0))) {
    exports = exports_view_type ("exports", count);
    if (debug) {
      std::ostringstream os;
      os << "*** exports resized to " << count << std::endl;
      std::cerr << os.str ();
    }
  }
  if (debug) {
    std::ostringstream os;
    os << "*** count: " << count << ", exports.extent(0): "
       << exports.extent (0) << std::endl;
    std::cerr << os.str ();
  }

  // If exports has nonzero length at this point, then the matrix has
  // at least one entry to pack.  Thus, if packing process ranks, we
  // had better have at least one process rank to pack.
  TEUCHOS_TEST_FOR_EXCEPTION
    (pack_pids && exports.extent (0) != 0 &&
     export_pids.extent (0) == 0, std::invalid_argument, prefix <<
     "pack_pids is true, and exports.extent(0) = " <<
     exports.extent (0)  << " != 0, meaning that we need to pack at least "
     "one matrix entry, but export_pids.extent(0) = 0.");

  typedef typename std::decay<decltype (local_matrix)>::type
    local_matrix_device_type;
  typedef typename std::decay<decltype (local_col_map)>::type
    local_map_type;

  exports.modify_device ();
  auto exports_d = exports.view_device ();
  do_pack<local_matrix_device_type, local_map_type, DT>
    (local_matrix, local_col_map, exports_d, num_packets_per_lid,
     export_lids, export_pids, offsets, num_bytes_per_value,
     pack_pids);
  // If we got this far, we succeeded.
}

} // namespace PackCrsMatrixImpl

template<typename ST, typename LO, typename GO, typename NT>
void
packCrsMatrix (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
               Teuchos::Array<char>& exports,
               const Teuchos::ArrayView<size_t>& numPacketsPerLID,
               const Teuchos::ArrayView<const LO>& exportLIDs,
               size_t& constantNumPackets)
{

  using device_type = typename CrsMatrix<ST,LO,GO,NT>::local_matrix_device_type::device_type;
  using buffer_device_type = typename DistObject<char, LO, GO, NT>::buffer_device_type;
  using host_exec_space = typename Kokkos::View<size_t*, device_type>::HostMirror::execution_space;
  using device_exec_space = typename device_type::execution_space;
  using host_dev_type = Kokkos::Device<host_exec_space, Kokkos::HostSpace>;

  // Convert all Teuchos::Array to Kokkos::View

  // This is an output array, so we don't have to copy to device here.
  // However, we'll have to remember to copy back to host when done.
  Kokkos::View<size_t*, buffer_device_type> num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array (buffer_device_type (),
                                            numPacketsPerLID.getRawPtr (),
                                            numPacketsPerLID.size (), false,
                                            "num_packets_per_lid");
  // FIXME (mfh 05 Feb 2019) We should just pass the exportLIDs
  // DualView through here, instead of recreating a device View from a
  // host ArrayView that itself came from a DualView.
  //
  // This is an input array, so we have to copy to device here.
  // However, we never need to copy it back to host.
  Kokkos::View<const LO*, buffer_device_type> export_lids_d =
    create_mirror_view_from_raw_host_array (buffer_device_type (),
                                            exportLIDs.getRawPtr (),
                                            exportLIDs.size (), true,
                                            "export_lids");

  Kokkos::View<int*, device_type> export_pids_d; // output arg
  Kokkos::DualView<char*, buffer_device_type> exports_dv; // output arg
  constexpr bool pack_pids = false;
  PackCrsMatrixImpl::packCrsMatrix<ST, LO, GO, NT, buffer_device_type> (
      sourceMatrix, exports_dv, num_packets_per_lid_d, export_lids_d,
      export_pids_d, constantNumPackets, pack_pids);

  // The counts are an output of PackCrsMatrixImpl::packCrsMatrix, so we have to
  // copy them back to host.
  Kokkos::View<size_t*, host_dev_type> num_packets_per_lid_h
    (numPacketsPerLID.getRawPtr (),
     numPacketsPerLID.size ());
  // DEEP_COPY REVIEW - DEVICE-TO-HOST
  Kokkos::deep_copy (device_exec_space(), num_packets_per_lid_h, num_packets_per_lid_d);

  // FIXME (mfh 23 Aug 2017) If we're forced to use a DualView for
  // exports_dv above, then we have two host copies for exports_h.

  // The exports are an output of PackCrsMatrixImpl::packCrsMatrix, so we have
  // to copy them back to host.
  if (static_cast<size_t> (exports.size ()) !=
      static_cast<size_t> (exports_dv.extent (0))) {
    exports.resize (exports_dv.extent (0));
  }
  Kokkos::View<char*, host_dev_type> exports_h (exports.getRawPtr (),
                                                exports.size ());
  // DEEP_COPY REVIEW - DEVICE-TO-HOST
  Kokkos::deep_copy (device_exec_space(), exports_h, exports_dv.d_view);
}

template<typename ST, typename LO, typename GO, typename NT>
void
packCrsMatrixNew(
  const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
  Kokkos::DualView<char*, typename DistObject<char, LO, GO, NT>::buffer_device_type>& exports,
  const Kokkos::DualView<size_t*, typename DistObject<char, LO, GO, NT>::buffer_device_type>& numPacketsPerLID,
  const Kokkos::DualView<const LO*, typename DistObject<char, LO, GO, NT>::buffer_device_type>& exportLIDs,
  size_t& constantNumPackets)
{
  using device_type = typename CrsMatrix<ST, LO, GO, NT>::device_type;
  using buffer_device_type = typename DistObject<char, LO, GO, NT>::buffer_device_type;

  // Create an empty array of PIDs, since the interface needs it.
  Kokkos::View<int*, device_type> exportPIDs_d ("exportPIDs", 0);
  constexpr bool pack_pids = false;

  // Write-only device access
  auto numPacketsPerLID_nc = numPacketsPerLID; // const DV& -> DV
  numPacketsPerLID_nc.clear_sync_state ();
  numPacketsPerLID_nc.modify_device ();
  auto numPacketsPerLID_d = numPacketsPerLID.view_device ();

  // Read-only device access
  TEUCHOS_ASSERT( ! exportLIDs.need_sync_device () );
  auto exportLIDs_d = exportLIDs.view_device ();

  ::Tpetra::Details::ProfilingRegion region_pack_crs_matrix_new(
    "Tpetra::Details::packCrsMatrixNew",
    "Import/Export"
  );
  PackCrsMatrixImpl::packCrsMatrix<ST,LO,GO,NT,buffer_device_type> (
      sourceMatrix, exports, numPacketsPerLID_d, exportLIDs_d,
      exportPIDs_d, constantNumPackets, pack_pids);
}

template<typename ST, typename LO, typename GO, typename NT>
void
packCrsMatrixWithOwningPIDs (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
                             Kokkos::DualView<char*, typename DistObject<char, LO, GO, NT>::buffer_device_type>& exports_dv,
                             const Teuchos::ArrayView<size_t>& numPacketsPerLID,
                             const Teuchos::ArrayView<const LO>& exportLIDs,
                             const Teuchos::ArrayView<const int>& sourcePIDs,
                             size_t& constantNumPackets)
{
  typedef typename CrsMatrix<ST,LO,GO,NT>::local_matrix_device_type local_matrix_device_type;
  typedef typename DistObject<char, LO, GO, NT>::buffer_device_type buffer_device_type;
  typedef typename Kokkos::DualView<char*, buffer_device_type>::t_host::execution_space host_exec_space;
  typedef Kokkos::Device<host_exec_space, Kokkos::HostSpace> host_dev_type;

  typename local_matrix_device_type::device_type outputDevice;
  typedef typename NT::execution_space execution_space;


  const bool verbose = ::Tpetra::Details::Behavior::verbose ();
  std::unique_ptr<std::string> prefix;
  if (verbose) {
    const int myRank = [&] () {
      auto map = sourceMatrix.getMap ();
      if (map.get () == nullptr) {
        return -1;
      }
      auto comm = map->getComm ();
      if (comm.get () == nullptr) {
        return -2;
      }
      return comm->getRank ();
    } ();
    std::ostringstream os;
    os << "Proc " << myRank << ": packCrsMatrixWithOwningPIDs: ";
    prefix = std::unique_ptr<std::string> (new std::string (os.str ()));

    std::ostringstream os2;
    os2 << *prefix << "start" << std::endl;
    std::cerr << os2.str ();
  }

  // Convert all Teuchos::Array to Kokkos::View

  // This is an output array, so we don't have to copy to device here.
  // However, we'll have to remember to copy back to host when done.
  auto num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array (buffer_device_type (),
                                            numPacketsPerLID.getRawPtr (),
                                            numPacketsPerLID.size (), false,
                                            "num_packets_per_lid");

  // This is an input array, so we have to copy to device here.
  // However, we never need to copy it back to host.
  auto export_lids_d =
    create_mirror_view_from_raw_host_array (buffer_device_type (),
                                            exportLIDs.getRawPtr (),
                                            exportLIDs.size (), true,
                                            "export_lids");
  // This is an input array, so we have to copy to device here.
  // However, we never need to copy it back to host.
  auto export_pids_d =
    create_mirror_view_from_raw_host_array (outputDevice,
                                            sourcePIDs.getRawPtr (),
                                            sourcePIDs.size (), true,
                                            "export_pids");
  constexpr bool pack_pids = true;
  try {
    PackCrsMatrixImpl::packCrsMatrix
      (sourceMatrix, exports_dv, num_packets_per_lid_d, export_lids_d,
       export_pids_d, constantNumPackets, pack_pids);
  }
  catch (std::exception& e) {
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "PackCrsMatrixImpl::packCrsMatrix threw: "
         << e.what () << std::endl;
      std::cerr << os.str ();
    }
    throw;
  }
  catch (...) {
    if (verbose) {
      std::ostringstream os;
      os << *prefix << "PackCrsMatrixImpl::packCrsMatrix threw an exception "
        "not a subclass of std::exception" << std::endl;
      std::cerr << os.str ();
    }
    throw;
  }

  if (numPacketsPerLID.size () != 0) {
    try {
      // The counts are an output of PackCrsMatrixImpl::packCrsMatrix,
      // so we have to copy them back to host.
      Kokkos::View<size_t*, host_dev_type> num_packets_per_lid_h
        (numPacketsPerLID.getRawPtr (), numPacketsPerLID.size ());
      // DEEP_COPY REVIEW - DEVICE-TO-HOST
      Kokkos::deep_copy (execution_space(), num_packets_per_lid_h, num_packets_per_lid_d);
    }
    catch (std::exception& e) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Kokkos::deep_copy threw: " << e.what () << std::endl;
        std::cerr << os.str ();
      }
      throw;
    }
    catch (...) {
      if (verbose) {
        std::ostringstream os;
        os << *prefix << "Kokkos::deep_copy threw an exception not a subclass "
          "of std::exception" << std::endl;
        std::cerr << os.str ();
      }
      throw;
    }
  }

  if (verbose) {
    std::ostringstream os;
    os << *prefix << "done" << std::endl;
    std::cerr << os.str ();
  }
}

} // namespace Details
} // namespace Tpetra

#define TPETRA_DETAILS_PACKCRSMATRIX_INSTANT( ST, LO, GO, NT ) \
  template void \
  Details::packCrsMatrix<ST, LO, GO, NT> (const CrsMatrix<ST, LO, GO, NT>&, \
    Teuchos::Array<char>&, \
    const Teuchos::ArrayView<size_t>&, \
    const Teuchos::ArrayView<const LO>&, \
    size_t&); \
  template void \
  Details::packCrsMatrixNew<ST, LO, GO, NT> (const CrsMatrix<ST, LO, GO, NT>&, \
    Kokkos::DualView<char*, DistObject<char, LO, GO, NT>::buffer_device_type>&, \
    const Kokkos::DualView<size_t*, DistObject<char, LO, GO, NT>::buffer_device_type>&, \
    const Kokkos::DualView<const LO*, DistObject<char, LO, GO, NT>::buffer_device_type>&, \
    size_t&); \
  template void \
  Details::packCrsMatrixWithOwningPIDs<ST, LO, GO, NT> (const CrsMatrix<ST, LO, GO, NT>&, \
    Kokkos::DualView<char*, DistObject<char, LO, GO, NT>::buffer_device_type>&, \
    const Teuchos::ArrayView<size_t>&, \
    const Teuchos::ArrayView<const LO>&, \
    const Teuchos::ArrayView<const int>&, \
    size_t&);

#endif // TPETRA_DETAILS_PACKCRSMATRIX_DEF_HPP
