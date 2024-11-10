// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_PACKCRSGRAPH_DEF_HPP
#define TPETRA_DETAILS_PACKCRSGRAPH_DEF_HPP

#include "TpetraCore_config.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Tpetra_Details_Behavior.hpp"
#include "Tpetra_Details_castAwayConstDualView.hpp"
#include "Tpetra_Details_createMirrorView.hpp"
#include "Tpetra_Details_getEntryOnHost.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_CrsGraph_decl.hpp"
#include <memory>
#include <string>

/// \file Tpetra_Details_packCrsGraph_def.hpp
/// \brief Functions for packing the entries of a Tpetra::CrsGraph
///   for communication, in the case where it is valid to go to the
///   KokkosSparse::CrsGraph (local sparse graph data structure)
///   directly.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.
///
/// Data (bytes) describing the row of the CRS graph are "packed"
/// (concatenated) in to a (view of) GO* object in the following order:
///
///   1. number of entries (LocalOrdinal)
///   2. global column indices (GlobalOrdinal)
///   3. proces IDs (optional, int)
///
/// The functions in this file are companions to
/// Tpetra_Details_unpackCrsGraphAndCombine.hpp, i.e.,
/// Tpetra_Details_unpackCrsGraphAndCombine.hpp implements the
/// reverse of the packing order described above to ensure proper
/// unpacking.

namespace Tpetra {

//
// Users must never rely on anything in the Details namespace.
//
namespace Details {

namespace PackCrsGraphImpl {
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
class NumPacketsAndOffsetsFunctor{
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

  NumPacketsAndOffsetsFunctor(const OutputOffsetsViewType& outputOffsets,
                               const CountsViewType& counts,
                               const InputOffsetsViewType& rowOffsets,
                               const InputLocalRowIndicesViewType& lclRowInds,
                               const InputLocalRowPidsViewType& lclRowPids) :
    outputOffsets_ (outputOffsets),
    counts_ (counts),
    rowOffsets_ (rowOffsets),
    lclRowInds_ (lclRowInds),
    lclRowPids_ (lclRowPids),
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

      // We pack first the global column indices and then pids (if any),
      // However, if the number of entries in the row is zero, we pack nothing.
      const count_type numEntToPack = (count == 0)
                                    ? static_cast<count_type>(0)
                                    : count * (1 + (lclRowPids_.size() > 0 ? 1 : 0));

      if (final) {
        counts_(curInd) = numEntToPack;
      }
      update += numEntToPack;
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
computeNumPacketsAndOffsets(const OutputOffsetsViewType& outputOffsets,
                            const CountsViewType& counts,
                            const InputOffsetsViewType& rowOffsets,
                            const InputLocalRowIndicesViewType& lclRowInds,
                            const InputLocalRowPidsViewType& lclRowPids)
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
       "but the graph has no rows.  lclRowInds.extent(0) = " <<
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

    functor_type f (outputOffsets, counts, rowOffsets, lclRowInds, lclRowPids);
    Kokkos::parallel_scan ("Tpetra::Details::computeNumPacketsAndOffsets::scan", range_type (0, numRowsToPack + 1), f);

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

/// \brief Packs a single row of the CrsGraph.
///
/// \tparam ColumnMap the type of the local column map
///
/// Data (bytes) describing the row of the CRS graph are "packed"
/// (concatenated) in to a single (view of) GO* in the following order:
///
///   1. GO column indices
///   2. int proces IDs
///
template<class Packet,
         class LocalMapType,
         class BufferDeviceType,
         class InputLidsType,
         class InputPidsType>
KOKKOS_FUNCTION
size_t
packRow(const LocalMapType& col_map,
        const Kokkos::View<Packet*, BufferDeviceType>& exports,
        const InputLidsType& lids_in,
        const InputPidsType& pids_in,
        const size_t offset,
        const size_t num_ent,
        const bool pack_pids)
{
  using LO = typename LocalMapType::local_ordinal_type;
  using GO = typename LocalMapType::global_ordinal_type;

  if (num_ent == 0) {
    // Empty rows always take zero bytes, to ensure sparsity.
    return static_cast<size_t>(0);
  }

  size_t num_ent_packed = num_ent;
  if (pack_pids) {
    num_ent_packed += num_ent;
  }

  // Copy column indices one at a time, so that we don't need
  // temporary storage.
  for (size_t k = 0; k < num_ent; ++k) {
    const LO lid = lids_in[k];
    const GO gid = col_map.getGlobalElement (lid);
    exports(offset+k) = gid;
  }
  // Copy PIDs one at a time, so that we don't need temporary storage.
  if (pack_pids) {
    for (size_t k = 0; k < num_ent; ++k) {
      const LO lid = lids_in[k];
      const int pid = pids_in[lid];
      exports(offset+num_ent+k) = static_cast<GO>(pid);
    }
  }

  return num_ent_packed;
}

template<class Packet,
         class LocalGraph,
         class LocalMap,
         class BufferDeviceType>
struct PackCrsGraphFunctor {
  using local_graph_type = LocalGraph;
  using local_map_type = LocalMap;
  using LO = typename local_map_type::local_ordinal_type;
  using GO = typename local_map_type::global_ordinal_type;

  using num_packets_per_lid_view_type =
    Kokkos::View<const size_t*, BufferDeviceType>;
  using offsets_view_type = Kokkos::View<const size_t*, BufferDeviceType>;
  using exports_view_type = Kokkos::View<Packet*, BufferDeviceType>;
  using export_lids_view_type =
    typename PackTraits<LO>::input_array_type;
  using source_pids_view_type =
    typename PackTraits<int>::input_array_type;

  using count_type =
    typename num_packets_per_lid_view_type::non_const_value_type;
  using offset_type = typename offsets_view_type::non_const_value_type;
  using value_type = Kokkos::pair<int, LO>;

  static_assert (std::is_same<LO, typename local_graph_type::data_type>::value,
                 "local_map_type::local_ordinal_type and "
                 "local_graph_type::data_type must be the same.");

  local_graph_type local_graph;
  local_map_type local_col_map;
  exports_view_type exports;
  num_packets_per_lid_view_type num_packets_per_lid;
  export_lids_view_type export_lids;
  source_pids_view_type source_pids;
  offsets_view_type offsets;
  bool pack_pids;

  PackCrsGraphFunctor(const local_graph_type& local_graph_in,
                      const local_map_type& local_col_map_in,
                      const exports_view_type& exports_in,
                      const num_packets_per_lid_view_type& num_packets_per_lid_in,
                      const export_lids_view_type& export_lids_in,
                      const source_pids_view_type& source_pids_in,
                      const offsets_view_type& offsets_in,
                      const bool pack_pids_in) :
    local_graph (local_graph_in),
    local_col_map (local_col_map_in),
    exports (exports_in),
    num_packets_per_lid (num_packets_per_lid_in),
    export_lids (export_lids_in),
    source_pids (source_pids_in),
    offsets (offsets_in),
    pack_pids (pack_pids_in)
  {
    const LO numRows = local_graph_in.numRows ();
    const LO rowMapDim =
      static_cast<LO> (local_graph.row_map.extent (0));
    TEUCHOS_TEST_FOR_EXCEPTION
      (numRows != 0 && rowMapDim != numRows + static_cast<LO> (1),
       std::logic_error, "local_graph.row_map.extent(0) = "
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
    const size_t num_packets_this_lid = num_packets_per_lid(i);
    const size_t num_ent =
      static_cast<size_t> (local_graph.row_map[export_lid+1]
                         - local_graph.row_map[export_lid]);

    // Only pack this row's data if it has a nonzero number of
    // entries.  We can do this because receiving processes get the
    // number of packets, and will know that zero packets means zero
    // entries.
    if (num_ent == 0) {
      return;
    }

    if (export_lid >= static_cast<LO>(local_graph.numRows())) {
      if (dst.first != 0) { // keep only the first error
        dst = Kokkos::make_pair (1, i); // invalid row
      }
      return;
    }
    else if ((offset > buf_size || offset + num_packets_this_lid > buf_size)) {
      if (dst.first != 0) { // keep only the first error
        dst = Kokkos::make_pair (2, i); // out of bounds
      }
      return;
    }

    // We can now pack this row

    // Since the graph is locally indexed on the calling process, we
    // have to use its column Map (which it _must_ have in this case)
    // to convert to global indices.
    const auto row_beg = local_graph.row_map[export_lid];
    const auto row_end = local_graph.row_map[export_lid + 1];
    auto lids_in = Kokkos::subview (local_graph.entries,
                                    Kokkos::make_pair (row_beg, row_end));
    size_t num_ent_packed_this_row =
      packRow (local_col_map, exports, lids_in,
               source_pids, offset, num_ent, pack_pids);
    if (num_ent_packed_this_row != num_packets_this_lid) {
      if (dst.first != 0) { // keep only the first error
        dst = Kokkos::make_pair (3, i);
      }
    }
  }
};

/// \brief Perform the pack operation for the graph
///
/// \tparam LocalGraph the specialization of the KokkosSparse::CrsGraph
///   local graph
/// \tparam LocalMap the type of the local column map
///
/// This is a higher level interface to the PackCrsGraphFunctor
template<class Packet,
         class LocalGraph,
         class LocalMap,
         class BufferDeviceType>
void
do_pack(const LocalGraph& local_graph,
        const LocalMap& local_map,
        const Kokkos::View<Packet*, BufferDeviceType>& exports,
        const typename PackTraits<
            size_t
        >::input_array_type& num_packets_per_lid,
        const typename PackTraits<
          typename LocalMap::local_ordinal_type
        >::input_array_type& export_lids,
        const typename PackTraits<
          int
        >::input_array_type& source_pids,
        const Kokkos::View<const size_t*, BufferDeviceType>& offsets,
        const bool pack_pids)
{
  using LO = typename LocalMap::local_ordinal_type;
  using execution_space = typename LocalGraph::device_type::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, LO>;
  const char prefix[] = "Tpetra::Details::PackCrsGraphImpl::do_pack: ";

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
    // If exports has nonzero length at this point, then the graph
    // has at least one entry to pack.  Thus, if packing process
    // ranks, we had better have at least one process rank to pack.
    TEUCHOS_TEST_FOR_EXCEPTION
      (pack_pids && exports.extent (0) != 0 &&
       source_pids.extent (0) == 0, std::invalid_argument, prefix <<
       "pack_pids is true, and exports.extent(0) = " <<
       exports.extent (0)  << " != 0, meaning that we need to pack at "
       "least one graph entry, but source_pids.extent(0) = 0.");
  }

  using pack_functor_type =
    PackCrsGraphFunctor<Packet, LocalGraph, LocalMap,
                        BufferDeviceType>;
  pack_functor_type f (local_graph, local_map, exports,
                       num_packets_per_lid, export_lids,
                       source_pids, offsets, pack_pids);

  typename pack_functor_type::value_type result;
  range_type range (0, num_packets_per_lid.extent (0));
  Kokkos::parallel_reduce ("Tpetra::Details::computeNumPacketsAndOffsets::reduce",range, f, result);

  if (result.first != 0) {
    // We can't deep_copy from AnonymousSpace Views, so we can't
    // print out any information from them in case of error.
    std::ostringstream os;
    if (result.first == 1) { // invalid local row index
      os << "invalid local row index";
    }
    else if (result.first == 2) { // invalid offset
      os << "invalid offset";
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::runtime_error, prefix << "PackCrsGraphFunctor "
       "reported error code " << result.first << " (" << os.str ()
       << ") for the first bad row " << result.second << ".");
  }
}

/// \brief Pack specified entries of the given local sparse graph for
///   communication.
///
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam NT The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \warning This is an implementation detail of Tpetra::CrsGraph.
///
/// \param sourceGraph [in] the CrsGraph source
///
/// \param exports [in/out] Output pack buffer; resized if needed.
///
/// \param num_packets_per_lid [out] Entry k gives the number of bytes
///   packed for row export_lids[k] of the local graph.
///
/// \param export_lids [in] Local indices of the rows to pack.
///
/// \param export_pids [in] Process ranks for the column indices in each packed row.
///
/// \param constant_num_packets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
template<typename LO, typename GO, typename NT>
void
packCrsGraph
(const CrsGraph<LO,GO,NT>& sourceGraph,
 Kokkos::DualView<
   typename CrsGraph<LO,GO,NT>::packet_type*,
   typename CrsGraph<LO,GO,NT>::buffer_device_type
 >& exports,
 const Kokkos::View<
   size_t*,
   typename CrsGraph<LO,GO,NT>::buffer_device_type
 >& num_packets_per_lid,
 const Kokkos::View<
   const LO*,
   typename CrsGraph<LO, GO, NT>::buffer_device_type
 >& export_lids,
 const Kokkos::View<
   const int*,
   typename CrsGraph<LO, GO, NT>::buffer_device_type
 >& export_pids,
 size_t& constant_num_packets,
 const bool pack_pids)
{
  using Kokkos::View;
  using crs_graph_type = CrsGraph<LO, GO, NT>;
  using packet_type = typename crs_graph_type::packet_type;
  using buffer_device_type = typename crs_graph_type::buffer_device_type;
  using exports_view_type = Kokkos::DualView<packet_type*, buffer_device_type>;
  using local_graph_device_type = typename crs_graph_type::local_graph_device_type;
  using local_map_type = typename Tpetra::Map<LO, GO, NT>::local_map_type;
  const char prefix[] = "Tpetra::Details::packCrsGraph: ";
  constexpr bool debug = false;

  local_graph_device_type local_graph = sourceGraph.getLocalGraphDevice ();
  local_map_type local_col_map = sourceGraph.getColMap ()->getLocalMap ();

  // Setting this to zero tells the caller to expect a possibly
  // different ("nonconstant") number of packets per local index
  // (i.e., a possibly different number of entries per row).
  constant_num_packets = 0;

  const size_t num_export_lids (export_lids.extent (0));
  TEUCHOS_TEST_FOR_EXCEPTION
    (num_export_lids != size_t (num_packets_per_lid.extent (0)),
     std::invalid_argument, prefix << "num_export_lids.extent(0) = "
     << num_export_lids << " != num_packets_per_lid.extent(0) = "
     << num_packets_per_lid.extent (0) << ".");
  if (num_export_lids != 0) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (num_packets_per_lid.data () == nullptr, std::invalid_argument,
       prefix << "num_export_lids = "<< num_export_lids << " != 0, but "
       "num_packets_per_lid.data() = "
       << num_packets_per_lid.data () << " == NULL.");
  }

  if (num_export_lids == 0) {
    exports = exports_view_type ("exports", 0);
    return;
  }

  // Array of offsets into the pack buffer.
  View<size_t*, buffer_device_type> offsets ("offsets", num_export_lids + 1);

  // Compute number of packets per LID (row to send), as well as
  // corresponding offsets (the prefix sum of the packet counts).
  const size_t count =
    computeNumPacketsAndOffsets(offsets, num_packets_per_lid,
                                local_graph.row_map, export_lids, export_pids);

  // Resize the output pack buffer if needed.
  if (count > size_t (exports.extent (0))) {
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

  // If exports has nonzero length at this point, then the graph has
  // at least one entry to pack.  Thus, if packing process ranks, we
  // had better have at least one process rank to pack.
  TEUCHOS_TEST_FOR_EXCEPTION
    (pack_pids && exports.extent (0) != 0 &&
     export_pids.extent (0) == 0, std::invalid_argument, prefix <<
     "pack_pids is true, and exports.extent(0) = " <<
     exports.extent (0)  << " != 0, meaning that we need to pack at least "
     "one graph entry, but export_pids.extent(0) = 0.");

  exports.modify_device ();
  auto exports_d = exports.view_device ();
  do_pack<packet_type, local_graph_device_type, local_map_type, buffer_device_type>
    (local_graph, local_col_map, exports_d, num_packets_per_lid,
     export_lids, export_pids, offsets, pack_pids);
  // If we got this far, we succeeded.
}

} // namespace PackCrsGraphImpl

template<typename LO, typename GO, typename NT>
void
packCrsGraph (const CrsGraph<LO, GO, NT>& sourceGraph,
              Teuchos::Array<typename CrsGraph<LO,GO,NT>::packet_type>& exports,
              const Teuchos::ArrayView<size_t>& numPacketsPerLID,
              const Teuchos::ArrayView<const LO>& exportLIDs,
              size_t& constantNumPackets)
{
  using Kokkos::HostSpace;
  using Kokkos::MemoryUnmanaged;
  using Kokkos::View;
  using crs_graph_type = CrsGraph<LO, GO, NT>;
  using packet_type = typename crs_graph_type::packet_type;
  using BDT = typename crs_graph_type::buffer_device_type;

  // Convert all Teuchos::Array to Kokkos::View

  // This is an output array, so we don't have to copy to device here.
  // However, we'll have to remember to copy back to host when done.
  BDT outputDevice;
  View<size_t*, BDT> num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array (outputDevice,
                                            numPacketsPerLID.getRawPtr (),
                                            numPacketsPerLID.size (), false,
                                            "num_packets_per_lid");
  // This is an input array, so we have to copy to device here.
  // However, we never need to copy it back to host.
  View<const LO*, BDT> export_lids_d =
    create_mirror_view_from_raw_host_array (outputDevice,
                                            exportLIDs.getRawPtr (),
                                            exportLIDs.size (), true,
                                            "export_lids");
  View<const int*, BDT> export_pids_d;
  Kokkos::DualView<packet_type*, BDT> exports_dv;
  constexpr bool pack_pids = false;

  static_assert
    (std::is_same<
       typename decltype (num_packets_per_lid_d)::non_const_value_type,
       size_t>::value,
     "num_packets_per_lid_d's non_const_value_type should be size_t.");
  static_assert
    (std::is_same<
       typename decltype (num_packets_per_lid_d)::device_type,
       BDT>::value,
     "num_packets_per_lid_d's BDT should be size_t.");
  static_assert
    (std::is_same<
       typename decltype (export_lids_d)::device_type,
       BDT>::value,
     "export_lids_d's device_type should be BDT.");
  static_assert
    (std::is_same<
       typename decltype (export_pids_d)::non_const_value_type,
       int>::value,
     "export_pids_d's non_const_value_type should be int.");
  static_assert
    (std::is_same<
       typename decltype (export_pids_d)::device_type,
       BDT>::value,
     "export_pids_d's device_type should be BDT.");

  PackCrsGraphImpl::packCrsGraph
    (sourceGraph, exports_dv, num_packets_per_lid_d, export_lids_d,
     export_pids_d, constantNumPackets, pack_pids);

  // The counts are an output of packCrsGraph, so we have to copy
  // them back to host.
  View<size_t*, HostSpace, MemoryUnmanaged>
    num_packets_per_lid_h (numPacketsPerLID.getRawPtr (),
                           numPacketsPerLID.size ());

  // DEEP_COPY REVIEW - DEVICE-TO-HOST
  using execution_space = typename BDT::execution_space;
  Kokkos::deep_copy (execution_space(), num_packets_per_lid_h, num_packets_per_lid_d);

  // FIXME (mfh 23 Aug 2017) If we're forced to use a DualView for
  // exports_dv above, then we have two host copies for exports_h.

  // The exports are an output of packCrsGraph, so we have to
  // copy them back to host.
  if (static_cast<size_t> (exports.size ()) !=
      static_cast<size_t> (exports_dv.extent (0))) {
    exports.resize (exports_dv.extent (0));
  }
  View<packet_type*, HostSpace, MemoryUnmanaged>
    exports_h (exports.getRawPtr (), exports.size ());
  // DEEP_COPY REVIEW - DEVICE-TO-HOST
  Kokkos::deep_copy (execution_space(), exports_h, exports_dv.d_view);
  execution_space().fence();
}

/// \brief Pack specified entries of the given local sparse graph for
///   communication ("new" DistObject interface version).
template<typename LO, typename GO, typename NT>
void
packCrsGraphNew (const CrsGraph<LO,GO,NT>& sourceGraph,
                 const Kokkos::DualView<
                   const LO*,
                   typename CrsGraph<LO,GO,NT>::buffer_device_type
                 >& export_lids,
                 const Kokkos::DualView<
                   const int*,
                   typename CrsGraph<LO,GO,NT>::buffer_device_type
                 >& export_pids,
                 Kokkos::DualView<
                   typename CrsGraph<LO,GO,NT>::packet_type*,
                   typename CrsGraph<LO,GO,NT>::buffer_device_type>& exports,
                 Kokkos::DualView<
                   size_t*,
                   typename CrsGraph<LO,GO,NT>::buffer_device_type
                 > num_packets_per_lid,
                 size_t& constant_num_packets,
                 const bool pack_pids)
{
  using Kokkos::View;
  using crs_graph_type = CrsGraph<LO,GO,NT>;
  using BDT = typename crs_graph_type::buffer_device_type;
  using PT = typename crs_graph_type::packet_type;
  using exports_dual_view_type = Kokkos::DualView<PT*, BDT>;
  using LGT = typename crs_graph_type::local_graph_device_type;
  using LMT = typename crs_graph_type::map_type::local_map_type;
  const char prefix[] = "Tpetra::Details::packCrsGraphNew: ";

  const LGT local_graph = sourceGraph.getLocalGraphDevice ();
  const LMT local_col_map = sourceGraph.getColMap ()->getLocalMap ();

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
  TEUCHOS_TEST_FOR_EXCEPTION
    (num_export_lids != 0 &&
     num_packets_per_lid.view_device ().data () == nullptr,
     std::invalid_argument, prefix << "num_export_lids = "<< num_export_lids
     << " != 0, but num_packets_per_lid.view_device().data() = nullptr.");

  if (num_export_lids == 0) {
    exports = exports_dual_view_type ();
    return;
  }

  // Array of offsets into the pack buffer.
  using offsets_type = Kokkos::View<size_t*, BDT>;
  offsets_type offsets ("offsets", num_export_lids + 1);

  // Compute number of packets per LID (row to send), as well as
  // corresponding offsets (the prefix sum of the packet counts).
  num_packets_per_lid.clear_sync_state ();
  num_packets_per_lid.modify_device ();
  using PackCrsGraphImpl::computeNumPacketsAndOffsets;
  const size_t count =
    computeNumPacketsAndOffsets (offsets, num_packets_per_lid.view_device (),
                                 local_graph.row_map,
                                 export_lids.view_device (),
                                 export_pids.view_device ());

  // Resize the output pack buffer if needed.
  if (count > static_cast<size_t> (exports.extent (0))) {
    exports = exports_dual_view_type ("exports", count);
  }

  // If exports has nonzero length at this point, then the graph has
  // at least one entry to pack.  Thus, if packing process ranks, we
  // had better have at least one process rank to pack.
  TEUCHOS_TEST_FOR_EXCEPTION
    (pack_pids && exports.extent (0) != 0 &&
     export_pids.extent (0) == 0, std::invalid_argument, prefix <<
     "pack_pids is true, and exports.extent(0) = " <<
     exports.extent (0)  << " != 0, meaning that we need to pack at least "
     "one graph entry, but export_pids.extent(0) = 0.");

  exports.modify_device ();
  using PackCrsGraphImpl::do_pack;
  do_pack<PT, LGT, LMT, BDT> (local_graph, local_col_map,
                              exports.view_device (),
                              num_packets_per_lid.view_device (),
                              export_lids.view_device (),
                              export_pids.view_device (),
                              offsets, pack_pids);
}

template<typename LO, typename GO, typename NT>
void
packCrsGraphWithOwningPIDs
(const CrsGraph<LO, GO, NT>& sourceGraph,
 Kokkos::DualView<
   typename CrsGraph<LO, GO, NT>::packet_type*,
   typename CrsGraph<LO, GO, NT>::buffer_device_type
 >& exports_dv,
 const Teuchos::ArrayView<size_t>& numPacketsPerLID,
 const Teuchos::ArrayView<const LO>& exportLIDs,
 const Teuchos::ArrayView<const int>& sourcePIDs,
 size_t& constantNumPackets)
{
  using Kokkos::HostSpace;
  using Kokkos::MemoryUnmanaged;
  using Kokkos::View;
  using crs_graph_type = CrsGraph<LO, GO, NT>;
  using buffer_device_type = typename crs_graph_type::buffer_device_type;

  // Convert all Teuchos::Array to Kokkos::View

  // This is an output array, so we don't have to copy to device here.
  // However, we'll have to remember to copy back to host when done.
  View<size_t*, buffer_device_type> num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array (buffer_device_type (),
                                            numPacketsPerLID.getRawPtr (),
                                            numPacketsPerLID.size (), false,
                                            "num_packets_per_lid");

  // This is an input array, so we have to copy to device here.
  // However, we never need to copy it back to host.
  View<const LO*, buffer_device_type> export_lids_d =
    create_mirror_view_from_raw_host_array (buffer_device_type (),
                                            exportLIDs.getRawPtr (),
                                            exportLIDs.size (), true,
                                            "export_lids");
  // This is an input array, so we have to copy to device here.
  // However, we never need to copy it back to host.
  View<const int*, buffer_device_type> export_pids_d =
    create_mirror_view_from_raw_host_array (buffer_device_type (),
                                            sourcePIDs.getRawPtr (),
                                            sourcePIDs.size (), true,
                                            "export_pids");
  constexpr bool pack_pids = true;
  PackCrsGraphImpl::packCrsGraph
    (sourceGraph, exports_dv, num_packets_per_lid_d, export_lids_d,
     export_pids_d, constantNumPackets, pack_pids);

  // The counts are an output of packCrsGraph, so we
  // have to copy them back to host.
  View<size_t*, HostSpace, MemoryUnmanaged> num_packets_per_lid_h
    (numPacketsPerLID.getRawPtr (), numPacketsPerLID.size ());
  // DEEP_COPY REVIEW - DEVICE-TO-HOST
  using execution_space = typename buffer_device_type::execution_space;
  Kokkos::deep_copy (execution_space(),
    num_packets_per_lid_h, num_packets_per_lid_d);
  execution_space().fence();
}

} // namespace Details
} // namespace Tpetra

#define TPETRA_DETAILS_PACKCRSGRAPH_INSTANT( LO, GO, NT ) \
  template void \
  Details::packCrsGraph<LO, GO, NT> ( \
    const CrsGraph<LO, GO, NT>&, \
    Teuchos::Array<CrsGraph<LO,GO,NT>::packet_type>&, \
    const Teuchos::ArrayView<size_t>&, \
    const Teuchos::ArrayView<const LO>&, \
    size_t&); \
  template void \
  Details::packCrsGraphNew<LO, GO, NT> ( \
    const CrsGraph<LO, GO, NT>&, \
    const Kokkos::DualView< \
      const LO*, \
      CrsGraph<LO,GO,NT>::buffer_device_type>&, \
    const Kokkos::DualView< \
      const int*, \
      CrsGraph<LO,GO,NT>::buffer_device_type>&, \
    Kokkos::DualView< \
      CrsGraph<LO,GO,NT>::packet_type*, \
      CrsGraph<LO,GO,NT>::buffer_device_type>&, \
    Kokkos::DualView< \
      size_t*, \
      CrsGraph<LO,GO,NT>::buffer_device_type>, \
    size_t&, \
    const bool); \
  template void \
  Details::packCrsGraphWithOwningPIDs<LO, GO, NT> ( \
    const CrsGraph<LO, GO, NT>&, \
    Kokkos::DualView<CrsGraph<LO,GO,NT>::packet_type*, CrsGraph<LO,GO,NT>::buffer_device_type>&, \
    const Teuchos::ArrayView<size_t>&, \
    const Teuchos::ArrayView<const LO>&, \
    const Teuchos::ArrayView<const int>&, \
    size_t&);

#endif // TPETRA_DETAILS_PACKCRSGRAPH_DEF_HPP
