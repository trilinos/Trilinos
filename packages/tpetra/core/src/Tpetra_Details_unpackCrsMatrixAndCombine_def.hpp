// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_DETAILS_UNPACKCRSMATRIXANDCOMBINE_DEF_HPP
#define TPETRA_DETAILS_UNPACKCRSMATRIXANDCOMBINE_DEF_HPP

#include <memory>
#include <string>
#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Tpetra_Details_castAwayConstDualView.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_createMirrorView.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_Details_getEntryOnHost.hpp"
#include "Tpetra_Details_DefaultTypes.hpp"
#include "Tpetra_Details_ExecutionSpaces.hpp"

/// \file Tpetra_Details_unpackCrsMatrixAndCombine_def.hpp
/// \brief Definition of functions for unpacking the entries of a
///   Tpetra::CrsMatrix for communication, in the case where it is
///   valid to go to the KokkosSparse::CrsMatrix (local sparse matrix
///   data structure) directly.
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
/// Tpetra_Details_packCrsMatrix.hpp, i.e., Tpetra_Details_packCrsMatrix.hpp
/// implements the packing order described above to ensure proper unpacking.

namespace Tpetra {

//
// Users must never rely on anything in the Details namespace.
//
namespace Details {

namespace UnpackAndCombineCrsMatrixImpl {

/// \brief Unpack a single row of a CrsMatrix
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
template<class ST, class LO, class GO>
KOKKOS_FUNCTION int
unpackRow(const typename PackTraits<GO>::output_array_type& gids_out,
          const typename PackTraits<int>::output_array_type& pids_out,
          const typename PackTraits<ST>::output_array_type& vals_out,
          const char imports[],
          const size_t offset,
          const size_t /* num_bytes */,
          const size_t num_ent,
          const size_t bytes_per_value)
{
  if (num_ent == 0) {
    // Empty rows always take zero bytes, to ensure sparsity.
    return 0;
  }
  bool unpack_pids = pids_out.size() > 0;

  const size_t num_ent_beg = offset;
  const size_t num_ent_len = PackTraits<LO>::packValueCount (LO (0));

  const size_t gids_beg = num_ent_beg + num_ent_len;
  const size_t gids_len =
    num_ent * PackTraits<GO>::packValueCount (GO (0));

  const size_t pids_beg = gids_beg + gids_len;
  const size_t pids_len = unpack_pids ?
    size_t (num_ent * PackTraits<int>::packValueCount (int (0))) :
    size_t (0);

  const size_t vals_beg = gids_beg + gids_len + pids_len;
  const size_t vals_len = num_ent * bytes_per_value;

  const char* const num_ent_in = imports + num_ent_beg;
  const char* const gids_in = imports + gids_beg;
  const char* const pids_in = unpack_pids ? imports + pids_beg : nullptr;
  const char* const vals_in = imports + vals_beg;

  size_t num_bytes_out = 0;
  LO num_ent_out;
  num_bytes_out += PackTraits<LO>::unpackValue (num_ent_out, num_ent_in);
  if (static_cast<size_t> (num_ent_out) != num_ent) {
    return 20; // error code
  }

  {
    Kokkos::pair<int, size_t> p;
    p = PackTraits<GO>::unpackArray (gids_out.data (), gids_in, num_ent);
    if (p.first != 0) {
      return 21; // error code
    }
    num_bytes_out += p.second;

    if (unpack_pids) {
      p = PackTraits<int>::unpackArray (pids_out.data (), pids_in, num_ent);
      if (p.first != 0) {
        return 22; // error code
      }
      num_bytes_out += p.second;
    }

    p = PackTraits<ST>::unpackArray (vals_out.data (), vals_in, num_ent);
    if (p.first != 0) {
      return 23; // error code
    }
    num_bytes_out += p.second;
  }

  const size_t expected_num_bytes = num_ent_len + gids_len + pids_len + vals_len;
  if (num_bytes_out != expected_num_bytes) {
    return 24; // error code
  }
  return 0; // no errors
} //unpackRow

/// \brief Unpacks and combines a single row of the CrsMatrix.
///
/// \tparam LocalMatrix KokkosSparse::CrsMatrix specialization.
/// \tparam LocalMap Type of the "local" column map
/// \tparam BufferDeviceType Type of the "buffer device type."
///   See Trilinos GitHub Issue #1088 for details.
///
/// Data (bytes) describing the row of the CRS matrix are "unpacked"
/// from a single (concatenated) (view of) char* directly into the
/// row of the matrix.
template<class LocalMatrix, class LocalMap, class BufferDeviceType>
struct UnpackCrsMatrixAndCombineFunctor {
  typedef LocalMatrix local_matrix_type;
  typedef LocalMap local_map_type;

  typedef typename local_matrix_type::value_type ST;
  typedef typename local_map_type::local_ordinal_type LO;
  typedef typename local_map_type::global_ordinal_type GO;
  typedef typename local_map_type::device_type DT;
  typedef typename DT::execution_space XS;

  typedef Kokkos::View<const size_t*, BufferDeviceType>
    num_packets_per_lid_type;
  typedef Kokkos::View<const size_t*, DT> offsets_type;
  typedef Kokkos::View<const char*, BufferDeviceType> input_buffer_type;
  typedef Kokkos::View<const LO*, BufferDeviceType> import_lids_type;

  typedef Kokkos::View<int, DT> error_type;
  using member_type = typename Kokkos::TeamPolicy<XS>::member_type;

  static_assert (std::is_same<LO, typename local_matrix_type::ordinal_type>::value,
                 "LocalMap::local_ordinal_type and "
                 "LocalMatrix::ordinal_type must be the same.");

  local_matrix_type local_matrix;
  local_map_type local_col_map;
  input_buffer_type imports;
  num_packets_per_lid_type num_packets_per_lid;
  import_lids_type import_lids;
  Kokkos::View<const LO*[2], DT> batch_info;
  offsets_type offsets;
  Tpetra::CombineMode combine_mode;
  size_t batch_size;
  size_t bytes_per_value;
  bool atomic;
  error_type error_code;

  UnpackCrsMatrixAndCombineFunctor(
      const local_matrix_type& local_matrix_in,
      const local_map_type& local_col_map_in,
      const input_buffer_type& imports_in,
      const num_packets_per_lid_type& num_packets_per_lid_in,
      const import_lids_type& import_lids_in,
      const Kokkos::View<const LO*[2], DT>& batch_info_in,
      const offsets_type& offsets_in,
      const Tpetra::CombineMode combine_mode_in,
      const size_t batch_size_in,
      const size_t bytes_per_value_in,
      const bool atomic_in) :
    local_matrix (local_matrix_in),
    local_col_map (local_col_map_in),
    imports (imports_in),
    num_packets_per_lid (num_packets_per_lid_in),
    import_lids (import_lids_in),
    batch_info (batch_info_in),
    offsets (offsets_in),
    combine_mode (combine_mode_in),
    batch_size (batch_size_in),
    bytes_per_value (bytes_per_value_in),
    atomic (atomic_in),
    error_code("error")
  {}

  KOKKOS_INLINE_FUNCTION
  void operator()(member_type team_member) const
  {
    using Kokkos::View;
    using Kokkos::subview;
    using Kokkos::MemoryUnmanaged;

    const LO batch = team_member.league_rank();
    const LO lid_no = batch_info(batch, 0);
    const LO batch_no = batch_info(batch, 1);

    const size_t num_bytes = num_packets_per_lid(lid_no);

    // Only unpack data if there is a nonzero number of bytes.
    if (num_bytes == 0)
      return;

    // there is actually something in the row
    const LO import_lid = import_lids(lid_no);
    const size_t buf_size = imports.size();
    const size_t offset = offsets(lid_no);

    // Get the number of entries to expect in the received data for this row.
    LO num_ent_LO = 0;
    const char* const in_buf = imports.data() + offset;
    (void) PackTraits<LO>::unpackValue(num_ent_LO, in_buf);
    const size_t num_entries_in_row = static_cast<size_t>(num_ent_LO);

    // Count the number of bytes expected to unpack
    size_t expected_num_bytes = 0;
    {
      expected_num_bytes += PackTraits<LO>::packValueCount(LO(0));
      expected_num_bytes += num_entries_in_row * PackTraits<GO>::packValueCount(GO(0));
      expected_num_bytes += num_entries_in_row * PackTraits<ST>::packValueCount(ST());
    }

    if (expected_num_bytes > num_bytes)
    {
// FIXME_SYCL Enable again once a SYCL conforming printf implementation is available.
#ifndef KOKKOS_ENABLE_SYCL
      printf(
        "*** Error: UnpackCrsMatrixAndCombineFunctor: "
        "At row %d, the expected number of bytes (%d) != number of unpacked bytes (%d)\n",
        (int) lid_no, (int) expected_num_bytes, (int) num_bytes
      );
#endif
      Kokkos::atomic_compare_exchange_strong(error_code.data(), 0, 21);
      return;
    }

    if (offset > buf_size || offset + num_bytes > buf_size)
    {
// FIXME_SYCL Enable again once a SYCL conforming printf implementation is available.
#ifndef KOKKOS_ENABLE_SYCL
      printf(
        "*** Error: UnpackCrsMatrixAndCombineFunctor: "
        "At row %d, the offset (%d) > buffer size (%d)\n",
        (int) lid_no, (int) offset, (int) buf_size
      );
#endif
      Kokkos::atomic_compare_exchange_strong(error_code.data(), 0, 22);
      return;
    }

    // Determine the number of entries to unpack in this batch
    size_t num_entries_in_batch = 0;
    if (num_entries_in_row <= batch_size)
      num_entries_in_batch = num_entries_in_row;
    else if (num_entries_in_row >= (batch_no + 1) * batch_size)
      num_entries_in_batch = batch_size;
    else
      num_entries_in_batch = num_entries_in_row - batch_no * batch_size;

    const size_t bytes_per_lid = PackTraits<LO>::packValueCount(LO(0));
    const size_t num_ent_start = offset;
    const size_t num_ent_end = num_ent_start + bytes_per_lid;

    const size_t bytes_per_gid = PackTraits<GO>::packValueCount(GO(0));
    const size_t gids_start = num_ent_end;
    const size_t gids_end = gids_start + num_entries_in_row * bytes_per_gid;

    const size_t vals_start = gids_end;

    const size_t shift = batch_no * batch_size;
    const char* const num_ent_in = imports.data() + num_ent_start;
    const char* const gids_in = imports.data() + gids_start + shift * bytes_per_gid;
    const char* const vals_in = imports.data() + vals_start + shift * bytes_per_value;

    LO num_ent_out;
    (void)PackTraits<LO>::unpackValue(num_ent_out, num_ent_in);
    if (static_cast<size_t>(num_ent_out) != num_entries_in_row)
    {
// FIXME_SYCL Enable again once a SYCL conforming printf implementation is available.
#ifndef KOKKOS_ENABLE_SYCL
      printf(
        "*** Error: UnpackCrsMatrixAndCombineFunctor: "
        "At row %d, number of entries (%d) != number of entries unpacked (%d)\n",
        (int) lid_no, (int) num_entries_in_row, (int) num_ent_out
      );
#endif
      Kokkos::atomic_compare_exchange_strong(error_code.data(), 0, 23);
    }

    constexpr bool matrix_has_sorted_rows = true; // see #6282
    //Note BMK 6-22: this lambda must use capture-by-value [=] and not capture-by-ref [&].
    //By ref triggers compiler bug in CUDA 10.
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team_member, num_entries_in_batch),
      [=](const LO& j)
      {
        size_t distance = 0;

        GO gid_out;
        distance = j * bytes_per_gid;
        (void) PackTraits<GO>::unpackValue(gid_out, gids_in + distance);
        auto lid_out = local_col_map.getLocalElement(gid_out);

        // Column indices come in as global indices, in case the
        // source object's column Map differs from the target object's
        // (this's) column Map, and must be converted local index values

        // assume that ST is default constructible
        ST val_out;
        distance = j * bytes_per_value;
        (void) PackTraits<ST>::unpackValue(val_out, vals_in + distance);

        if (combine_mode == ADD) {
          // NOTE (mfh 20 Nov 2019) Must assume atomic is required, unless
          // different threads don't touch the same row (i.e., no
          // duplicates in incoming LIDs list).
          const bool use_atomic_updates = atomic;
          (void)local_matrix.sumIntoValues(
            import_lid,
            &lid_out,
            1,
            &val_out,
            matrix_has_sorted_rows,
            use_atomic_updates
          );
        } else if (combine_mode == REPLACE) {
          // NOTE (mfh 20 Nov 2019): It's never correct to use REPLACE
          // combine mode with multiple incoming rows that touch the same
          // target matrix entries, so we never need atomic updates.
          const bool use_atomic_updates = false;
          (void)local_matrix.replaceValues(
            import_lid,
            &lid_out,
            1,
            &val_out,
            matrix_has_sorted_rows,
            use_atomic_updates
          );
        } else {
          // should never get here
// FIXME_SYCL Enable again once a SYCL conforming printf implementation is available.
#ifndef KOKKOS_ENABLE_SYCL
          printf(
            "*** Error: UnpackCrsMatrixAndCombineFunctor: "
            "At row %d, an unknown error occurred during unpack\n", (int) lid_no
          );
#endif
          Kokkos::atomic_compare_exchange_strong(error_code.data(), 0, 31);
        }
      }
    );

    team_member.team_barrier();

  }

  //! Host function for getting the error.
  int error() const {
    auto error_code_h = Kokkos::create_mirror_view_and_copy(
      Kokkos::HostSpace(), error_code
    );
    return error_code_h();
  }

}; //UnpackCrsMatrixAndCombineFunctor 

struct MaxNumEntTag {};
struct TotNumEntTag {};

/// \brief Kokkos::parallel_reduce functor to determine the number of
///   entries (to unpack) in a KokkosSparse::CrsMatrix to pack.
///
/// \tparam LO The type of local indices.
/// \tparam DT The Kokkos device type.
/// \tparam BDT The "buffer device type."
///
/// The maximum number of entries is required for allocated device scratch space
template<class LO, class DT, class BDT>
class NumEntriesFunctor {
public:
  typedef Kokkos::View<const size_t*, BDT> num_packets_per_lid_type;
  typedef Kokkos::View<const size_t*, DT> offsets_type;
  typedef Kokkos::View<const char*, BDT> input_buffer_type;
  // This needs to be public, since it appears in the argument list of
  // public methods (see below).  Otherwise, build errors may happen.
  typedef size_t value_type;

private:
  num_packets_per_lid_type num_packets_per_lid;
  offsets_type offsets;
  input_buffer_type imports;

public:
  NumEntriesFunctor (const num_packets_per_lid_type num_packets_per_lid_in,
                     const offsets_type& offsets_in,
                     const input_buffer_type& imports_in) :
    num_packets_per_lid (num_packets_per_lid_in),
    offsets (offsets_in),
    imports (imports_in)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const MaxNumEntTag, const LO i, value_type& update) const {
    // Get how many entries to expect in the received data for this row.
    const size_t num_bytes = num_packets_per_lid(i);
    if (num_bytes > 0) {
      LO num_ent_LO = 0; // output argument of unpackValue
      const char* const in_buf = imports.data () + offsets(i);
      (void) PackTraits<LO>::unpackValue (num_ent_LO, in_buf);
      const size_t num_ent = static_cast<size_t> (num_ent_LO);

      update = (update < num_ent) ? num_ent : update;
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (const MaxNumEntTag,
        value_type& dst,
        const value_type& src) const
  {
    if (dst < src) dst = src;
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const TotNumEntTag, const LO i, value_type& tot_num_ent) const {
    // Get how many entries to expect in the received data for this row.
    const size_t num_bytes = num_packets_per_lid(i);
    if (num_bytes > 0) {
      LO num_ent_LO = 0; // output argument of unpackValue
      const char* const in_buf = imports.data () + offsets(i);
      (void) PackTraits<LO>::unpackValue (num_ent_LO, in_buf);
      tot_num_ent += static_cast<size_t> (num_ent_LO);
    }
  }
}; //NumEntriesFunctor 

/// \brief Maximum number of entries in any row of the packed matrix.
///
/// \tparam LO The type of local indices.
/// \tparam DT The Kokkos device type.
/// \tparam The "buffer device type."
///
/// This procedure is a higher-level interface to NumEntriesFunctor.
template<class LO, class DT, class BDT>
size_t
compute_maximum_num_entries (
  const Kokkos::View<const size_t*, BDT>& num_packets_per_lid,
  const Kokkos::View<const size_t*, DT>& offsets,
  const Kokkos::View<const char*, BDT>& imports)
{
  typedef typename DT::execution_space XS;
  typedef Kokkos::RangePolicy<XS, Kokkos::IndexType<LO>,
    MaxNumEntTag> range_policy;

  NumEntriesFunctor<LO, DT, BDT> functor (num_packets_per_lid, offsets,
                                          imports);
  const LO numRowsToUnpack =
    static_cast<LO> (num_packets_per_lid.extent (0));
  size_t max_num_ent = 0;
  Kokkos::parallel_reduce ("Max num entries in CRS",
                           range_policy (0, numRowsToUnpack),
                           functor, max_num_ent);
  return max_num_ent;
}

/// \brief Total number of entries in any row of the packed matrix.
///
/// \tparam LO The type of local indices.
/// \tparam DT The Kokkos device type.
/// \tparam BDT The "buffer device type."
///
/// This procedure is a high level interface to NumEntriesFunctor
template<class LO, class DT, class BDT>
size_t
compute_total_num_entries (
  const Kokkos::View<const size_t*, BDT>& num_packets_per_lid,
  const Kokkos::View<const size_t*, DT>& offsets,
  const Kokkos::View<const char*, BDT>& imports)
{
  typedef typename DT::execution_space XS;
  typedef Kokkos::RangePolicy<XS, Kokkos::IndexType<LO>, TotNumEntTag> range_policy;
  size_t tot_num_ent = 0;
  NumEntriesFunctor<LO, DT, BDT> functor (num_packets_per_lid, offsets,
                                          imports);
  const LO numRowsToUnpack =
    static_cast<LO> (num_packets_per_lid.extent (0));
  Kokkos::parallel_reduce ("Total num entries in CRS to unpack",
                           range_policy (0, numRowsToUnpack),
                           functor, tot_num_ent);
  return tot_num_ent;
}

template<class LO>
KOKKOS_INLINE_FUNCTION
size_t
unpackRowCount(const char imports[],
               const size_t offset,
               const size_t num_bytes)
{
  using PT = PackTraits<LO>;

  LO num_ent_LO = 0;
  if (num_bytes > 0) {
    const size_t p_num_bytes = PT::packValueCount(num_ent_LO);
    if (p_num_bytes > num_bytes) {
      return OrdinalTraits<size_t>::invalid();
    }
    const char* const in_buf = imports + offset;
    (void) PT::unpackValue(num_ent_LO, in_buf);
  }
  return static_cast<size_t>(num_ent_LO);
}

/// \brief Compute the index and batch number associated with each batch
///
/// batch_info(i, 0) is the local index of the ith batch
/// batch_info(i, 1) is the local batch number of the ith batch
template<class View1, class View2>
inline
bool
compute_batch_info(
  const View1& batches_per_lid,
  View2& batch_info
)
{
  using LO = typename View2::value_type;
  size_t batch = 0;
  for (size_t i=0; i<batches_per_lid.extent(0); i++)
  {
    for (size_t batch_no=0; batch_no<batches_per_lid(i); batch_no++)
    {
      batch_info(batch, 0) = static_cast<LO>(i);
      batch_info(batch, 1) = batch_no;
      batch++;
    }
  }
  return batch == batch_info.extent(0);
}

/// \brief Perform the unpack operation for the matrix
///
/// \tparam LocalMatrix the specialization of the KokkosSparse::CrsMatrix
///   local matrix
/// \tparam LocalMap the type of the local column map
///
/// This is a higher level interface to the UnpackCrsMatrixAndCombineFunctor
template<class LocalMatrix, class LocalMap, class BufferDeviceType>
void
unpackAndCombineIntoCrsMatrix(
    const LocalMatrix& local_matrix,
    const LocalMap& local_map,
    const Kokkos::View<const char*, BufferDeviceType>& imports,
    const Kokkos::View<const size_t*, BufferDeviceType>& num_packets_per_lid,
    const typename PackTraits<typename LocalMap::local_ordinal_type>::input_array_type import_lids,
    const Tpetra::CombineMode combine_mode)
{
  using ST = typename LocalMatrix::value_type;
  using LO = typename LocalMap::local_ordinal_type;
  using DT = typename LocalMap::device_type;
  using XS = typename DT::execution_space;
  const char prefix[] =
    "Tpetra::Details::UnpackAndCombineCrsMatrixImpl::"
    "unpackAndCombineIntoCrsMatrix: ";

  const size_t num_import_lids = static_cast<size_t>(import_lids.extent(0));
  if (num_import_lids == 0) {
    // Nothing to unpack
    return;
  }

  {
    // Check for correct input
    TEUCHOS_TEST_FOR_EXCEPTION(combine_mode == ABSMAX,
        std::invalid_argument,
        prefix << "ABSMAX combine mode is not yet implemented for a matrix that has a "
        "static graph (i.e., was constructed with the CrsMatrix constructor "
        "that takes a const CrsGraph pointer).");

    TEUCHOS_TEST_FOR_EXCEPTION(combine_mode == INSERT,
        std::invalid_argument,
        prefix << "INSERT combine mode is not allowed if the matrix has a static graph "
        "(i.e., was constructed with the CrsMatrix constructor that takes a "
        "const CrsGraph pointer).");

    // Unknown combine mode!
    TEUCHOS_TEST_FOR_EXCEPTION(!(combine_mode == ADD || combine_mode == REPLACE),
        std::invalid_argument,
        prefix << "Invalid combine mode; should never get "
        "here!  Please report this bug to the Tpetra developers.");

    // Check that sizes of input objects are consistent.
    bool bad_num_import_lids =
      num_import_lids != static_cast<size_t>(num_packets_per_lid.extent(0));
    TEUCHOS_TEST_FOR_EXCEPTION(bad_num_import_lids,
        std::invalid_argument,
        prefix << "importLIDs.size() (" << num_import_lids << ") != "
        "numPacketsPerLID.size() (" << num_packets_per_lid.extent(0) << ").");
  } // end QA error checking

  // Get the offsets
  Kokkos::View<size_t*, DT> offsets("offsets", num_import_lids+1);
  computeOffsetsFromCounts(offsets, num_packets_per_lid);

  // Determine the sizes of the unpack batches
  size_t max_num_ent = compute_maximum_num_entries<LO,DT>(num_packets_per_lid, offsets, imports);
  const size_t default_batch_size = Tpetra::Details::Behavior::hierarchicalUnpackBatchSize();
  const size_t batch_size = std::min(default_batch_size, max_num_ent);

  // To achieve some balance amongst threads, unpack each row in equal size batches
  size_t num_batches = 0;
  Kokkos::View<LO*[2], DT> batch_info("", num_batches);
  Kokkos::View<size_t*, DT> batches_per_lid("", num_import_lids);
  // Compute meta data that allows batch unpacking
  Kokkos::parallel_reduce(
    Kokkos::RangePolicy<XS, Kokkos::IndexType<size_t>>(0, num_import_lids),
    KOKKOS_LAMBDA(const size_t i, size_t& batches)
    {
      const size_t num_entries_in_row = unpackRowCount<LO>(
        imports.data(), offsets(i), num_packets_per_lid(i)
      );
      batches_per_lid(i) =
        (num_entries_in_row <= batch_size) ?
        1 :
        num_entries_in_row / batch_size + (num_entries_in_row % batch_size != 0);
      batches += batches_per_lid(i);
    },
    num_batches
  );
  Kokkos::resize(batch_info, num_batches);

  Kokkos::HostSpace host_space;
  auto batches_per_lid_h = Kokkos::create_mirror_view(host_space, batches_per_lid);
  // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
  Kokkos::deep_copy(XS(), batches_per_lid_h, batches_per_lid);

  auto batch_info_h = Kokkos::create_mirror_view(host_space, batch_info);

  (void) compute_batch_info(batches_per_lid_h, batch_info_h);
  // DEEP_COPY REVIEW - HOSTMIRROR-TO-DEVICE
  Kokkos::deep_copy(XS(), batch_info, batch_info_h);

  // FIXME (TJF SEP 2017)
  // The scalar type is not necessarily default constructible
  size_t bytes_per_value = PackTraits<ST>::packValueCount(ST());

  // Now do the actual unpack!
  const bool atomic = XS().concurrency() != 1;
  using functor = UnpackCrsMatrixAndCombineFunctor<LocalMatrix, LocalMap, BufferDeviceType>;
  functor f(
    local_matrix,
    local_map,
    imports,
    num_packets_per_lid,
    import_lids,
    batch_info,
    offsets,
    combine_mode,
    batch_size,
    bytes_per_value,
    atomic
  );

  using policy = Kokkos::TeamPolicy<XS, Kokkos::IndexType<LO>>;
  const size_t team_size = Tpetra::Details::Behavior::hierarchicalUnpackTeamSize();
  if (!Spaces::is_gpu_exec_space<XS>() || team_size == Teuchos::OrdinalTraits<size_t>::invalid())
  {
    Kokkos::parallel_for(policy(static_cast<LO>(num_batches), Kokkos::AUTO), f);
  }
  else
  {
    Kokkos::parallel_for(policy(static_cast<LO>(num_batches), static_cast<int>(team_size)), f);
  }

  auto error_code = f.error();
  TEUCHOS_TEST_FOR_EXCEPTION(
    error_code != 0,
    std::runtime_error,
    prefix << "UnpackCrsMatrixAndCombineFunctor reported error code " << error_code
  );
} //unpackAndCombineIntoCrsMatrix (Kokkos version)

template<class LocalMatrix, class BufferDeviceType>
size_t
unpackAndCombineWithOwningPIDsCount(
  const LocalMatrix& local_matrix,
  const typename PackTraits<typename LocalMatrix::ordinal_type>::input_array_type permute_from_lids,
  const Kokkos::View<const char*, BufferDeviceType, void, void>& imports,
  const Kokkos::View<const size_t*, BufferDeviceType, void, void>& num_packets_per_lid,
  const size_t num_same_ids)
{
  using Kokkos::parallel_reduce;
  typedef typename LocalMatrix::ordinal_type LO;
  typedef typename LocalMatrix::device_type device_type;
  typedef typename device_type::execution_space XS;
  typedef typename Kokkos::View<LO*, device_type>::size_type size_type;
  typedef Kokkos::RangePolicy<XS, Kokkos::IndexType<LO> > range_policy;
  typedef BufferDeviceType BDT;

  size_t count = 0;
  LO num_items;

  // Number of matrix entries to unpack (returned by this function).
  num_items = static_cast<LO>(num_same_ids);
  if (num_items) {
    size_t kcnt = 0;
    parallel_reduce(range_policy(0, num_items),
      KOKKOS_LAMBDA(const LO lid, size_t& update) {
        update += static_cast<size_t>(local_matrix.graph.row_map[lid+1]
                                     -local_matrix.graph.row_map[lid]);
      }, kcnt);
    count += kcnt;
  }

  // Count entries copied directly from the source matrix with permuting.
  num_items = static_cast<LO>(permute_from_lids.extent(0));
  if (num_items) {
    size_t kcnt = 0;
    parallel_reduce(range_policy(0, num_items),
      KOKKOS_LAMBDA(const LO i, size_t& update) {
        const LO lid = permute_from_lids(i);
        update += static_cast<size_t> (local_matrix.graph.row_map[lid+1]
                                     - local_matrix.graph.row_map[lid]);
      }, kcnt);
    count += kcnt;
  }

  {
    // Count entries received from other MPI processes.
    const size_type np = num_packets_per_lid.extent(0);
    Kokkos::View<size_t*, device_type> offsets("offsets", np+1);
    computeOffsetsFromCounts(offsets, num_packets_per_lid);
    count +=
      compute_total_num_entries<LO, device_type, BDT> (num_packets_per_lid,
                                                       offsets, imports);
  }

  return count;
} //unpackAndCombineWithOwningPIDsCount (Kokkos version)

/// \brief Setup row pointers for remotes
template<class LO, class DT, class BDT>
int
setupRowPointersForRemotes(
  const typename PackTraits<size_t>::output_array_type& tgt_rowptr,
  const typename PackTraits<LO>::input_array_type& import_lids,
  const Kokkos::View<const char*, BDT>& imports,
  const Kokkos::View<const size_t*, BDT>& num_packets_per_lid,
  const typename PackTraits<size_t>::input_array_type& offsets)
{
  using Kokkos::parallel_reduce;
  typedef typename DT::execution_space XS;
  typedef typename PackTraits<size_t>::input_array_type::size_type size_type;
  typedef Kokkos::RangePolicy<XS, Kokkos::IndexType<size_type> > range_policy;

  const size_t InvalidNum = OrdinalTraits<size_t>::invalid();
  const size_type N = num_packets_per_lid.extent(0);

  int errors = 0;
  parallel_reduce ("Setup row pointers for remotes",
    range_policy (0, N),
    KOKKOS_LAMBDA (const size_t i, int& k_error) {
      typedef typename std::remove_reference< decltype( tgt_rowptr(0) ) >::type atomic_incr_type;
      const size_t num_bytes = num_packets_per_lid(i);
      const size_t offset = offsets(i);
      const size_t num_ent = unpackRowCount<LO> (imports.data(), offset, num_bytes);
      if (num_ent == InvalidNum) {
        k_error += 1;
      }
      Kokkos::atomic_fetch_add (&tgt_rowptr (import_lids(i)), atomic_incr_type(num_ent));
    }, errors);
  return errors;
}

// Convert array of row lengths to a CRS pointer array
template<class DT>
void
makeCrsRowPtrFromLengths(
    const typename PackTraits<size_t>::output_array_type& tgt_rowptr,
    const Kokkos::View<size_t*,DT>& new_start_row)
{
  using Kokkos::parallel_scan;
  typedef typename DT::execution_space XS;
  typedef typename Kokkos::View<size_t*,DT>::size_type size_type;
  typedef Kokkos::RangePolicy<XS, Kokkos::IndexType<size_type> > range_policy;
  const size_type N = new_start_row.extent(0);
  parallel_scan(range_policy(0, N),
    KOKKOS_LAMBDA(const size_t& i, size_t& update, const bool& final) {
      auto cur_val = tgt_rowptr(i);
      if (final) {
        tgt_rowptr(i) = update;
        new_start_row(i) = tgt_rowptr(i);
      }
      update += cur_val;
    }
  );
}

template<class LocalMatrix, class LocalMap>
void
copyDataFromSameIDs(
    const typename PackTraits<typename LocalMap::global_ordinal_type>::output_array_type& tgt_colind,
    const typename PackTraits<int>::output_array_type& tgt_pids,
    const typename PackTraits<typename LocalMatrix::value_type>::output_array_type& tgt_vals,
    const Kokkos::View<size_t*, typename LocalMap::device_type>& new_start_row,
    const typename PackTraits<size_t>::output_array_type& tgt_rowptr,
    const typename PackTraits<int>::input_array_type& src_pids,
    const LocalMatrix& local_matrix,
    const LocalMap& local_col_map,
    const size_t num_same_ids,
    const int my_pid)
{
  using Kokkos::parallel_for;
  typedef typename LocalMap::device_type DT;
  typedef typename LocalMap::local_ordinal_type LO;
  typedef typename DT::execution_space XS;
  typedef Kokkos::RangePolicy<XS, Kokkos::IndexType<size_t> > range_policy;

  parallel_for(range_policy(0, num_same_ids),
    KOKKOS_LAMBDA(const size_t i) {
      typedef typename std::remove_reference< decltype( new_start_row(0) ) >::type atomic_incr_type;

      const LO src_lid    = static_cast<LO>(i);
      size_t src_row = local_matrix.graph.row_map(src_lid);

      const LO tgt_lid      = static_cast<LO>(i);
      const size_t tgt_row = tgt_rowptr(tgt_lid);

      const size_t nsr = local_matrix.graph.row_map(src_lid+1)
                       - local_matrix.graph.row_map(src_lid);
      Kokkos::atomic_fetch_add(&new_start_row(tgt_lid), atomic_incr_type(nsr));

      for (size_t j=local_matrix.graph.row_map(src_lid);
                  j<local_matrix.graph.row_map(src_lid+1); ++j) {
        LO src_col = local_matrix.graph.entries(j);
        tgt_vals(tgt_row + j - src_row)   = local_matrix.values(j);
        tgt_colind(tgt_row + j - src_row) = local_col_map.getGlobalElement(src_col);
        tgt_pids(tgt_row + j - src_row) = (src_pids(src_col) != my_pid) ? src_pids(src_col) : -1;
      }
    }
  );
}

template<class LocalMatrix, class LocalMap>
void
copyDataFromPermuteIDs(
    const typename PackTraits<typename LocalMap::global_ordinal_type>::output_array_type& tgt_colind,
    const typename PackTraits<int>::output_array_type& tgt_pids,
    const typename PackTraits<typename LocalMatrix::value_type>::output_array_type& tgt_vals,
    const Kokkos::View<size_t*,typename LocalMap::device_type>& new_start_row,
    const typename PackTraits<size_t>::output_array_type& tgt_rowptr,
    const typename PackTraits<int>::input_array_type& src_pids,
    const typename PackTraits<typename LocalMap::local_ordinal_type>::input_array_type& permute_to_lids,
    const typename PackTraits<typename LocalMap::local_ordinal_type>::input_array_type& permute_from_lids,
    const LocalMatrix& local_matrix,
    const LocalMap& local_col_map,
    const int my_pid)
{
  using Kokkos::parallel_for;
  typedef typename LocalMap::device_type DT;
  typedef typename LocalMap::local_ordinal_type LO;
  typedef typename DT::execution_space XS;
  typedef typename PackTraits<LO>::input_array_type::size_type size_type;
  typedef Kokkos::RangePolicy<XS, Kokkos::IndexType<size_type> > range_policy;

  const size_type num_permute_to_lids = permute_to_lids.extent(0);

  parallel_for(range_policy(0, num_permute_to_lids),
    KOKKOS_LAMBDA(const size_t i) {
      typedef typename std::remove_reference< decltype( new_start_row(0) ) >::type atomic_incr_type;

      const LO src_lid = permute_from_lids(i);
      const size_t src_row = local_matrix.graph.row_map(src_lid);

      const LO tgt_lid = permute_to_lids(i);
      const size_t tgt_row = tgt_rowptr(tgt_lid);

      size_t nsr = local_matrix.graph.row_map(src_lid+1)
                 - local_matrix.graph.row_map(src_lid);
      Kokkos::atomic_fetch_add(&new_start_row(tgt_lid), atomic_incr_type(nsr));

      for (size_t j=local_matrix.graph.row_map(src_lid);
                  j<local_matrix.graph.row_map(src_lid+1); ++j) {
        LO src_col = local_matrix.graph.entries(j);
        tgt_vals(tgt_row + j - src_row)   = local_matrix.values(j);
        tgt_colind(tgt_row + j - src_row) = local_col_map.getGlobalElement(src_col);
        tgt_pids(tgt_row + j - src_row) = (src_pids(src_col) != my_pid) ? src_pids(src_col) : -1;
      }
    }
  );
}

template<typename LocalMatrix, typename LocalMap, typename BufferDeviceType>
int
unpackAndCombineIntoCrsArrays2(
    const typename PackTraits<typename LocalMap::global_ordinal_type>::output_array_type& tgt_colind,
    const typename PackTraits<int>::output_array_type& tgt_pids,
    const typename PackTraits<typename LocalMatrix::value_type>::output_array_type& tgt_vals,
    const Kokkos::View<size_t*,typename LocalMap::device_type>& new_start_row,
    const typename PackTraits<size_t>::input_array_type& offsets,
    const typename PackTraits<typename LocalMap::local_ordinal_type>::input_array_type& import_lids,
    const Kokkos::View<const char*, BufferDeviceType, void, void>& imports,
    const Kokkos::View<const size_t*, BufferDeviceType, void, void>& num_packets_per_lid,
    const LocalMatrix& /* local_matrix */,
    const LocalMap /*& local_col_map*/,
    const int my_pid,
    const size_t bytes_per_value)
{
  using Kokkos::View;
  using Kokkos::subview;
  using Kokkos::MemoryUnmanaged;
  using Kokkos::parallel_reduce;
  using Kokkos::atomic_fetch_add;
  using Tpetra::Details::PackTraits;
  typedef typename LocalMap::device_type DT;
  typedef typename LocalMap::local_ordinal_type LO;
  typedef typename LocalMap::global_ordinal_type GO;
  typedef typename LocalMatrix::value_type ST;
  typedef typename DT::execution_space XS;
  typedef typename Kokkos::View<LO*, DT>::size_type size_type;
  typedef typename Kokkos::pair<size_type, size_type> slice;
  typedef Kokkos::RangePolicy<XS, Kokkos::IndexType<size_type> > range_policy;

  typedef View<int*,DT, MemoryUnmanaged> pids_out_type;
  typedef View<GO*, DT, MemoryUnmanaged> gids_out_type;
  typedef View<ST*, DT, MemoryUnmanaged> vals_out_type;

  const size_t InvalidNum = OrdinalTraits<size_t>::invalid();

  int errors = 0;
  const size_type num_import_lids = import_lids.size();

  // RemoteIDs: Loop structure following UnpackAndCombine
  parallel_reduce ("Unpack and combine into CRS",
    range_policy (0, num_import_lids),
    KOKKOS_LAMBDA (const size_t i, int& k_error) {
      typedef typename std::remove_reference< decltype( new_start_row(0) ) >::type atomic_incr_type;
      const size_t num_bytes = num_packets_per_lid(i);
      const size_t offset = offsets(i);
      if (num_bytes == 0) {
        // Empty buffer means that the row is empty.
        return;
      }
      size_t num_ent = unpackRowCount<LO>(imports.data(), offset, num_bytes);
      if (num_ent == InvalidNum) {
        k_error += 1;
        return;
      }
      const LO lcl_row = import_lids(i);
      const size_t start_row = atomic_fetch_add(&new_start_row(lcl_row), atomic_incr_type(num_ent));
      const size_t end_row = start_row + num_ent;

      gids_out_type gids_out = subview(tgt_colind, slice(start_row, end_row));
      vals_out_type vals_out = subview(tgt_vals, slice(start_row, end_row));
      pids_out_type pids_out = subview(tgt_pids, slice(start_row, end_row));

      k_error += unpackRow<ST,LO,GO>(gids_out, pids_out, vals_out,
                                     imports.data(), offset, num_bytes,
                                     num_ent, bytes_per_value);

      // Correct target PIDs.
      for (size_t j = 0; j < static_cast<size_t>(num_ent); ++j) {
        const int pid = pids_out(j);
        pids_out(j) = (pid != my_pid) ? pid : -1;
      }
    }, errors);

  return errors;
}

template<typename LocalMatrix, typename LocalMap, typename BufferDeviceType>
void
unpackAndCombineIntoCrsArrays(
    const LocalMatrix & local_matrix,
    const LocalMap & local_col_map,
    const typename PackTraits<typename LocalMap::local_ordinal_type>::input_array_type& import_lids,
    const Kokkos::View<const char*, BufferDeviceType, void, void>& imports,
    const Kokkos::View<const size_t*, BufferDeviceType, void, void>& num_packets_per_lid,
    const typename PackTraits<typename LocalMap::local_ordinal_type>::input_array_type& permute_to_lids,
    const typename PackTraits<typename LocalMap::local_ordinal_type>::input_array_type& permute_from_lids,
    const typename PackTraits<size_t>::output_array_type& tgt_rowptr,
    const typename PackTraits<typename LocalMap::global_ordinal_type>::output_array_type& tgt_colind,
    const typename PackTraits<typename LocalMatrix::value_type>::output_array_type& tgt_vals,
    const typename PackTraits<int>::input_array_type& src_pids,
    const typename PackTraits<int>::output_array_type& tgt_pids,
    const size_t num_same_ids,
    const size_t tgt_num_rows,
    const size_t tgt_num_nonzeros,
    const int my_tgt_pid,
    const size_t bytes_per_value)
{
  using Kokkos::View;
  using Kokkos::subview;
  using Kokkos::parallel_for;
  using Kokkos::MemoryUnmanaged;
  using Tpetra::Details::PackTraits;
  typedef typename LocalMap::device_type DT;
  typedef typename LocalMap::local_ordinal_type LO;
  typedef typename DT::execution_space XS;
  typedef typename Kokkos::View<LO*, DT>::size_type size_type;
  typedef Kokkos::RangePolicy<XS, Kokkos::IndexType<size_t> > range_policy;
  typedef BufferDeviceType BDT;

  const char prefix[] = "unpackAndCombineIntoCrsArrays: ";

  const size_t N = tgt_num_rows;

  // In the case of reduced communicators, the sourceMatrix won't have
  // the right "my_pid", so thus we have to supply it.
  const int my_pid = my_tgt_pid;

  // Zero the rowptr
  parallel_for(range_policy(0, N+1),
    KOKKOS_LAMBDA(const size_t i) {
      tgt_rowptr(i) = 0;
    }
  );

  // same IDs: Always first, always in the same place
  parallel_for(range_policy(0, num_same_ids),
    KOKKOS_LAMBDA(const size_t i) {
      const LO tgt_lid = static_cast<LO>(i);
      const LO src_lid = static_cast<LO>(i);
      tgt_rowptr(tgt_lid) = local_matrix.graph.row_map(src_lid+1)
                          - local_matrix.graph.row_map(src_lid);
    }
  );

  // Permute IDs: Still local, but reordered
  const size_type num_permute_to_lids = permute_to_lids.extent(0);
  parallel_for(range_policy(0, num_permute_to_lids),
    KOKKOS_LAMBDA(const size_t i) {
      const LO tgt_lid = permute_to_lids(i);
      const LO src_lid = permute_from_lids(i);
      tgt_rowptr(tgt_lid) = local_matrix.graph.row_map(src_lid+1)
                          - local_matrix.graph.row_map(src_lid);
    }
  );

  // Get the offsets from the number of packets per LID
  const size_type num_import_lids = import_lids.extent(0);
  View<size_t*, DT> offsets("offsets", num_import_lids+1);
  computeOffsetsFromCounts(offsets, num_packets_per_lid);

#ifdef HAVE_TPETRA_DEBUG
  {
    auto nth_offset_h = getEntryOnHost(offsets, num_import_lids);
    const bool condition =
      nth_offset_h != static_cast<size_t>(imports.extent (0));
    TEUCHOS_TEST_FOR_EXCEPTION
      (condition, std::logic_error, prefix
       << "The final offset in bytes " << nth_offset_h
       << " != imports.size() = " << imports.extent(0)
       << ".  Please report this bug to the Tpetra developers.");
  }
#endif // HAVE_TPETRA_DEBUG

  // Setup row pointers for remotes
  int k_error =
    setupRowPointersForRemotes<LO,DT,BDT>(tgt_rowptr,
      import_lids, imports, num_packets_per_lid, offsets);
  TEUCHOS_TEST_FOR_EXCEPTION(k_error != 0, std::logic_error, prefix
    << " Error transferring data to target row pointers.  "
       "Please report this bug to the Tpetra developers.");

  // If multiple processes contribute to the same row, we may need to
  // update row offsets.  This tracks that.
  View<size_t*, DT> new_start_row ("new_start_row", N+1);

  // Turn row length into a real CRS row pointer
  makeCrsRowPtrFromLengths(tgt_rowptr, new_start_row);

  // SameIDs: Copy the data over
  copyDataFromSameIDs(tgt_colind, tgt_pids, tgt_vals, new_start_row,
      tgt_rowptr, src_pids, local_matrix, local_col_map, num_same_ids, my_pid);

  copyDataFromPermuteIDs(tgt_colind, tgt_pids, tgt_vals, new_start_row,
      tgt_rowptr, src_pids, permute_to_lids, permute_from_lids,
      local_matrix, local_col_map, my_pid);

  if (imports.extent(0) <= 0) {
    return;
  }

  int unpack_err = unpackAndCombineIntoCrsArrays2(tgt_colind, tgt_pids,
      tgt_vals, new_start_row, offsets, import_lids, imports, num_packets_per_lid,
      local_matrix, local_col_map, my_pid, bytes_per_value);
  TEUCHOS_TEST_FOR_EXCEPTION(
    unpack_err != 0, std::logic_error, prefix << "unpack loop failed.  This "
    "should never happen.  Please report this bug to the Tpetra developers.");

  return;
}

} // namespace UnpackAndCombineCrsMatrixImpl

/// \brief Unpack the imported column indices and values, and combine into matrix.
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam Node The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \param sourceMatrix [in] the CrsMatrix source
///
/// \param imports [in] Input pack buffer
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local matrix.
///
/// \param importLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
///
/// \param distor [in] The distributor (not used)
///
/// \param combineMode [in] the mode to use for combining values
///
/// \param atomic [in] whether or not do atomic adds/replaces in to the matrix
///
/// \warning The allowed \c combineMode are:
///   ADD, REPLACE, and ABSMAX. INSERT is not allowed.
///
/// This is the public interface to the unpack and combine machinery and
/// converts passed Teuchos::ArrayView objects to Kokkos::View objects (and
/// copies back in to the Teuchos::ArrayView objects, if needed).  When
/// CrsMatrix migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<typename ST, typename LO, typename GO, typename Node>
void
unpackCrsMatrixAndCombine(
    const CrsMatrix<ST, LO, GO, Node>& sourceMatrix,
    const Teuchos::ArrayView<const char>& imports,
    const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
    const Teuchos::ArrayView<const LO>& importLIDs,
    size_t /* constantNumPackets */,
    CombineMode combineMode)
{
  using Kokkos::View;
  typedef typename Node::device_type device_type;
  typedef typename CrsMatrix<ST, LO, GO, Node>::local_matrix_device_type local_matrix_device_type;
  static_assert (std::is_same<device_type, typename local_matrix_device_type::device_type>::value,
                 "Node::device_type and LocalMatrix::device_type must be the same.");

  // Convert all Teuchos::Array to Kokkos::View.
  device_type outputDevice;

  // numPacketsPerLID, importLIDs, and imports are input, so we have to copy
  // them to device.  Since unpacking is done directly in to the local matrix
  // (lclMatrix), no copying needs to be performed after unpacking.
  auto num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array(outputDevice, numPacketsPerLID.getRawPtr(),
        numPacketsPerLID.size(), true, "num_packets_per_lid");

  auto import_lids_d =
    create_mirror_view_from_raw_host_array(outputDevice, importLIDs.getRawPtr(),
        importLIDs.size(), true, "import_lids");

  auto imports_d =
    create_mirror_view_from_raw_host_array(outputDevice, imports.getRawPtr(),
        imports.size(), true, "imports");

  auto local_matrix = sourceMatrix.getLocalMatrixDevice();
  auto local_col_map = sourceMatrix.getColMap()->getLocalMap();

//KDDKDD This loop doesn't appear to do anything; what is it?
//KDDKDD  for (int i=0; i<importLIDs.size(); i++)
//KDDKDD  {
//KDDKDD    auto lclRow = importLIDs[i];
//KDDKDD    Teuchos::ArrayView<const LO> A_indices;
//KDDKDD    Teuchos::ArrayView<const ST> A_values;
//KDDKDD    sourceMatrix.getLocalRowView(lclRow, A_indices, A_values);
//KDDKDD  }
  // Now do the actual unpack!
  UnpackAndCombineCrsMatrixImpl::unpackAndCombineIntoCrsMatrix(
      local_matrix, local_col_map, imports_d, num_packets_per_lid_d,
      import_lids_d, combineMode);

}

template<typename ST, typename LO, typename GO, typename NT>
void
unpackCrsMatrixAndCombineNew(
  const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
  Kokkos::DualView<char*,
    typename DistObject<char, LO, GO, NT>::buffer_device_type> imports,
  Kokkos::DualView<size_t*,
    typename DistObject<char, LO, GO, NT>::buffer_device_type> numPacketsPerLID,
  const Kokkos::DualView<const LO*,
    typename DistObject<char, LO, GO, NT>::buffer_device_type>& importLIDs,
  const size_t /* constantNumPackets */,
  const CombineMode combineMode)
{
  using Kokkos::View;
  using crs_matrix_type = CrsMatrix<ST, LO, GO, NT>;
  using dist_object_type = DistObject<char, LO, GO, NT>;
  using device_type = typename crs_matrix_type::device_type;
  using local_matrix_device_type = typename crs_matrix_type::local_matrix_device_type;
  using buffer_device_type = typename dist_object_type::buffer_device_type;

  static_assert
    (std::is_same<device_type, typename local_matrix_device_type::device_type>::value,
     "crs_matrix_type::device_type and local_matrix_device_type::device_type "
     "must be the same.");

  if (numPacketsPerLID.need_sync_device()) {
    numPacketsPerLID.sync_device ();
  }
  auto num_packets_per_lid_d = numPacketsPerLID.view_device ();

  TEUCHOS_ASSERT( ! importLIDs.need_sync_device () );
  auto import_lids_d = importLIDs.view_device ();

  if (imports.need_sync_device()) {
    imports.sync_device ();
  }
  auto imports_d = imports.view_device ();

  auto local_matrix = sourceMatrix.getLocalMatrixDevice ();
  auto local_col_map = sourceMatrix.getColMap ()->getLocalMap ();
  typedef decltype (local_col_map) local_map_type;

  UnpackAndCombineCrsMatrixImpl::unpackAndCombineIntoCrsMatrix<
      local_matrix_device_type,
      local_map_type,
      buffer_device_type
    > (local_matrix, local_col_map, imports_d, num_packets_per_lid_d,
       import_lids_d, combineMode);
}

/// \brief Special version of Tpetra::Details::unpackCrsMatrixAndCombine
///   that also unpacks owning process ranks.
///
/// Perform the count for unpacking the imported column indices pids,
/// and values, and combining them into matrix.  Return (a ceiling on)
/// the number of local stored entries ("nonzeros") in the matrix.  If
/// there are no shared rows in the sourceMatrix this count is exact.
///
/// Note: This routine also counts the copyAndPermute nonzeros in
/// addition to those that come in via import.
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam Node The Kokkos Node type.  See the documentation of Map
///   for requirements.
///
/// \param sourceMatrix [in] the CrsMatrix source
///
/// \param imports [in] Input pack buffer
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local matrix.
///
/// \param importLIDs [in] Local indices of the rows to pack.
///
/// \param constantNumPackets [out] Setting this to zero tells the caller
///   to expect a possibly /// different ("nonconstant") number of packets per local index
///   (i.e., a possibly different number of entries per row).
///
/// \param distor [in] The distributor (not used)
///
/// \param combineMode [in] the mode to use for combining values
///
/// \param numSameIds [in]
///
/// \param permuteToLIDs [in]
///
/// \param permuteFromLIDs [in]
///
/// \warning The allowed \c combineMode are:
///   ADD, REPLACE, and ABSMAX. INSERT is not allowed.
//
/// \warning This method is intended for expert developer use
///   only, and should never be called by user code.
///
/// Note: This is the public interface to the unpack and combine machinery and
/// converts passed Teuchos::ArrayView objects to Kokkos::View objects (and
/// copies back in to the Teuchos::ArrayView objects, if needed).  When
/// CrsMatrix migrates fully to adopting Kokkos::DualView objects for its storage
/// of data, this procedure could be bypassed.
template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
size_t
unpackAndCombineWithOwningPIDsCount (
    const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & sourceMatrix,
    const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
    const Teuchos::ArrayView<const char> &imports,
    const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
    size_t /* constantNumPackets */,
    CombineMode /* combineMode */,
    size_t numSameIDs,
    const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
    const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs)
{
  using Kokkos::MemoryUnmanaged;
  using Kokkos::View;
  typedef typename Node::device_type DT;
  const char prefix[] = "unpackAndCombineWithOwningPIDsCount: ";

  TEUCHOS_TEST_FOR_EXCEPTION
    (permuteToLIDs.size () != permuteFromLIDs.size (), std::invalid_argument,
     prefix << "permuteToLIDs.size() = " << permuteToLIDs.size () << " != "
     "permuteFromLIDs.size() = " << permuteFromLIDs.size() << ".");
  // FIXME (mfh 26 Jan 2015) If there are no entries on the calling
  // process, then the matrix is neither locally nor globally indexed.
  const bool locallyIndexed = sourceMatrix.isLocallyIndexed ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (! locallyIndexed, std::invalid_argument, prefix << "The input "
    "CrsMatrix 'sourceMatrix' must be locally indexed.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (importLIDs.size () != numPacketsPerLID.size (), std::invalid_argument,
     prefix << "importLIDs.size() = " << importLIDs.size () << " != "
     "numPacketsPerLID.size() = " << numPacketsPerLID.size () << ".");

  auto local_matrix = sourceMatrix.getLocalMatrixDevice ();

  using kokkos_device_type = Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>;

  Kokkos::View<LocalOrdinal const *, kokkos_device_type, void, void > permute_from_lids_d =
    create_mirror_view_from_raw_host_array (DT (),
                                            permuteFromLIDs.getRawPtr (),
                                            permuteFromLIDs.size (), true,
                                            "permute_from_lids");

  Kokkos::View<const char*, kokkos_device_type, void, void > imports_d =
    create_mirror_view_from_raw_host_array (DT (),
                                            imports.getRawPtr (),
                                            imports.size (), true,
                                            "imports");

  Kokkos::View<const size_t*, kokkos_device_type, void, void > num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array (DT (),
                                            numPacketsPerLID.getRawPtr (),
                                            numPacketsPerLID.size (), true,
                                            "num_packets_per_lid");

  return UnpackAndCombineCrsMatrixImpl::unpackAndCombineWithOwningPIDsCount(
      local_matrix, permute_from_lids_d, imports_d,
      num_packets_per_lid_d, numSameIDs);
} //unpackAndCombineWithOwningPIDsCount (Teuchos::Array version)

/// \brief unpackAndCombineIntoCrsArrays
///
/// \note You should call unpackAndCombineWithOwningPIDsCount first
///   and allocate all arrays accordingly, before calling this
///   function.
///
/// Note: The SourcePids vector (on input) should contain owning PIDs
/// for each column in the (source) ColMap, as from
/// Tpetra::Import_Util::getPids, with the "-1 for local" option being
/// used.
///
/// Note: The TargetPids vector (on output) will contain owning PIDs
/// for each entry in the matrix, with the "-1 for local" for locally
/// owned entries.

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
unpackAndCombineIntoCrsArrays (
    const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & sourceMatrix,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                        Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void > import_lids_d,
    const Kokkos::View<const char*, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void > imports_d,
    const Kokkos::View<const size_t*, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void > num_packets_per_lid_d,
    const size_t numSameIDs,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void > permute_to_lids_d,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void > permute_from_lids_d,
    size_t TargetNumRows,
    const int MyTargetPID,
    Kokkos::View<size_t*,typename Node::device_type> &crs_rowptr_d,
    Kokkos::View<GlobalOrdinal*,typename Node::device_type>     &crs_colind_d,
    Kokkos::View<typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type*,typename Node::device_type>& crs_vals_d,
    const Teuchos::ArrayView<const int>& SourcePids,
    Kokkos::View<int*,typename Node::device_type> &TargetPids)
{
  using execution_space = typename Node::execution_space;
  using Tpetra::Details::PackTraits;

  using Kokkos::View;
  using Kokkos::deep_copy;

  using Teuchos::ArrayView;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MAX;
  using Teuchos::reduceAll;

  typedef typename Node::device_type DT;

  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef typename matrix_type::impl_scalar_type ST;

  const char prefix[] = "Tpetra::Details::unpackAndCombineIntoCrsArrays_new: ";
#  ifdef HAVE_TPETRA_MMM_TIMINGS
   using Teuchos::TimeMonitor;
   Teuchos::RCP<TimeMonitor> tm;
#  endif

  using Kokkos::MemoryUnmanaged;

  TEUCHOS_TEST_FOR_EXCEPTION
    (permute_to_lids_d.size () != permute_from_lids_d.size (), std::invalid_argument,
     prefix << "permute_to_lids_d.size() = " << permute_to_lids_d.size () << " != "
     "permute_from_lids_d.size() = " << permute_from_lids_d.size() << ".");
  // FIXME (mfh 26 Jan 2015) If there are no entries on the calling
  // process, then the matrix is neither locally nor globally indexed.
  const bool locallyIndexed = sourceMatrix.isLocallyIndexed ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (! locallyIndexed, std::invalid_argument, prefix << "The input "
    "CrsMatrix 'sourceMatrix' must be locally indexed.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (((size_t)import_lids_d.size ()) != num_packets_per_lid_d.size (), std::invalid_argument,
     prefix << "import_lids_d.size() = " << import_lids_d.size () << " != "
     "num_packets_per_lid_d.size() = " << num_packets_per_lid_d.size () << ".");

  auto local_matrix = sourceMatrix.getLocalMatrixDevice ();

  // TargetNumNonzeros is number of nonzeros in local matrix.
# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("unpackAndCombineWithOwningPIDsCount"))));
# endif
  size_t TargetNumNonzeros =
     UnpackAndCombineCrsMatrixImpl::unpackAndCombineWithOwningPIDsCount(
      local_matrix, permute_from_lids_d, imports_d,
      num_packets_per_lid_d, numSameIDs);
# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::null;
# endif

# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("resize CRS pointers"))));
# endif
  Kokkos::resize(crs_rowptr_d,TargetNumRows+1);
  Kokkos::resize(crs_colind_d,TargetNumNonzeros);
  Kokkos::resize(crs_vals_d,TargetNumNonzeros);
# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::null;
# endif

  TEUCHOS_TEST_FOR_EXCEPTION(
    permute_to_lids_d.size () != permute_from_lids_d.size (), std::invalid_argument,
    prefix << "permuteToLIDs.size() = " << permute_to_lids_d.size ()
    << "!= permute_from_lids_d.size() = " << permute_from_lids_d.size () << ".");

  if (static_cast<size_t> (TargetPids.size ()) != TargetNumNonzeros) {
    Kokkos::resize(TargetPids,TargetNumNonzeros);
  }
  Kokkos::deep_copy(execution_space(), TargetPids, -1);

  // Grab pointers for sourceMatrix
  auto local_col_map = sourceMatrix.getColMap()->getLocalMap();

# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("create mirror views from inputs"))));
# endif
  // Convert input arrays to Kokkos::Views
  DT outputDevice;

  auto src_pids_d =
    create_mirror_view_from_raw_host_array(outputDevice, SourcePids.getRawPtr(),
        SourcePids.size(), true, "src_pids");

# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::null;
# endif

  size_t bytes_per_value = 0;
  if (PackTraits<ST>::compileTimeSize) {
    // assume that ST is default constructible
    bytes_per_value = PackTraits<ST>::packValueCount(ST());
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
    size_t bytes_per_value_l = 0;
    if (local_matrix.values.extent(0) > 0) {
      const ST& val = local_matrix.values(0);
      bytes_per_value_l = PackTraits<ST>::packValueCount(val);
    } else {
      const ST& val = crs_vals_d(0);
      bytes_per_value_l = PackTraits<ST>::packValueCount(val);
    }
    Teuchos::reduceAll<int, size_t>(*(sourceMatrix.getComm()),
                                    Teuchos::REDUCE_MAX,
                                    bytes_per_value_l,
                                    outArg(bytes_per_value));
  }

# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("unpackAndCombineIntoCrsArrays"))));
# endif
  UnpackAndCombineCrsMatrixImpl::unpackAndCombineIntoCrsArrays(
      local_matrix, local_col_map, import_lids_d, imports_d,
      num_packets_per_lid_d, permute_to_lids_d, permute_from_lids_d,
      crs_rowptr_d, crs_colind_d, crs_vals_d, src_pids_d, TargetPids,
      numSameIDs, TargetNumRows, TargetNumNonzeros, MyTargetPID,
      bytes_per_value);
# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::null;
# endif

  // Copy outputs back to host
# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("copy back to host"))));
# endif

  Kokkos::parallel_for("setLocalEntriesToPID", Kokkos::RangePolicy<typename DT::execution_space>(0,TargetPids.size()), KOKKOS_LAMBDA (const size_t i) {
    if (TargetPids(i) == -1) TargetPids(i) = MyTargetPID;
  });

} //unpackAndCombineIntoCrsArrays

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void
unpackAndCombineIntoCrsArrays (
    const CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & sourceMatrix,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                        Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void > import_lids_d,
    const Kokkos::View<const char*, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void > imports_d,
    const Kokkos::View<const size_t*, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void > num_packets_per_lid_d,
    const size_t numSameIDs,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void > permute_to_lids_d,
    const Kokkos::View<LocalOrdinal const *, 
          Kokkos::Device<typename Node::device_type::execution_space,
                         Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename Node::device_type>>,
          void, void > permute_from_lids_d,
    size_t TargetNumRows,
    const int MyTargetPID,
    Teuchos::ArrayRCP<size_t>& CRS_rowptr,
    Teuchos::ArrayRCP<GlobalOrdinal>& CRS_colind,
    Teuchos::ArrayRCP<Scalar>& CRS_vals,
    const Teuchos::ArrayView<const int>& SourcePids,
    Teuchos::Array<int>& TargetPids)
{
  using execution_space = typename Node::execution_space;
  using Tpetra::Details::PackTraits;

  using Kokkos::View;
  using Kokkos::deep_copy;

  using Teuchos::ArrayView;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MAX;
  using Teuchos::reduceAll;

  typedef typename Node::device_type DT;

  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef typename matrix_type::impl_scalar_type ST;

  const char prefix[] = "Tpetra::Details::unpackAndCombineIntoCrsArrays_new: ";
#  ifdef HAVE_TPETRA_MMM_TIMINGS
   using Teuchos::TimeMonitor;
   Teuchos::RCP<TimeMonitor> tm;
#  endif

  using Kokkos::MemoryUnmanaged;

  TEUCHOS_TEST_FOR_EXCEPTION
    (permute_to_lids_d.size () != permute_from_lids_d.size (), std::invalid_argument,
     prefix << "permute_to_lids_d.size() = " << permute_to_lids_d.size () << " != "
     "permute_from_lids_d.size() = " << permute_from_lids_d.size() << ".");
  // FIXME (mfh 26 Jan 2015) If there are no entries on the calling
  // process, then the matrix is neither locally nor globally indexed.
  const bool locallyIndexed = sourceMatrix.isLocallyIndexed ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (! locallyIndexed, std::invalid_argument, prefix << "The input "
    "CrsMatrix 'sourceMatrix' must be locally indexed.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (((size_t)import_lids_d.size ()) != num_packets_per_lid_d.size (), std::invalid_argument,
     prefix << "import_lids_d.size() = " << import_lids_d.size () << " != "
     "num_packets_per_lid_d.size() = " << num_packets_per_lid_d.size () << ".");

  auto local_matrix = sourceMatrix.getLocalMatrixDevice ();

  // TargetNumNonzeros is number of nonzeros in local matrix.
# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("unpackAndCombineWithOwningPIDsCount"))));
# endif
  size_t TargetNumNonzeros =
     UnpackAndCombineCrsMatrixImpl::unpackAndCombineWithOwningPIDsCount(
      local_matrix, permute_from_lids_d, imports_d,
      num_packets_per_lid_d, numSameIDs);
# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::null;
# endif

# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("resize CRS pointers"))));
# endif
  CRS_rowptr.resize (TargetNumRows+1);
  CRS_colind.resize(TargetNumNonzeros);
  CRS_vals.resize(TargetNumNonzeros);
  Teuchos::ArrayRCP<ST> const & CRS_vals_impl_scalar_type = Teuchos::arcp_reinterpret_cast<ST>(CRS_vals);
# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::null;
# endif

  TEUCHOS_TEST_FOR_EXCEPTION(
    permute_to_lids_d.size () != permute_from_lids_d.size (), std::invalid_argument,
    prefix << "permuteToLIDs.size() = " << permute_to_lids_d.size ()
    << "!= permute_from_lids_d.size() = " << permute_from_lids_d.size () << ".");

  // Preseed TargetPids with -1 for local
  if (static_cast<size_t> (TargetPids.size ()) != TargetNumNonzeros) {
    TargetPids.resize (TargetNumNonzeros);
  }
  TargetPids.assign (TargetNumNonzeros, -1);

  // Grab pointers for sourceMatrix
  auto local_col_map = sourceMatrix.getColMap()->getLocalMap();

# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("create mirror views from inputs"))));
# endif
  // Convert input arrays to Kokkos::Views
  DT outputDevice;

  auto crs_rowptr_d =
    create_mirror_view_from_raw_host_array(outputDevice, CRS_rowptr.getRawPtr(),
        CRS_rowptr.size(), true, "crs_rowptr");

  auto crs_colind_d =
    create_mirror_view_from_raw_host_array(outputDevice, CRS_colind.getRawPtr(),
        CRS_colind.size(), true, "crs_colidx");
#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  static_assert (! std::is_same<
      typename std::remove_const<
        typename std::decay<
          decltype (CRS_vals_impl_scalar_type)
        >::type::value_type
      >::type,
      std::complex<double> >::value,
    "CRS_vals::value_type is std::complex<double>; this should never happen"
    ", since std::complex does not work in Kokkos::View objects.");
#endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE

  auto crs_vals_d =
    create_mirror_view_from_raw_host_array(outputDevice, CRS_vals_impl_scalar_type.getRawPtr(),
        CRS_vals_impl_scalar_type.size(), true, "crs_vals");

#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  static_assert (! std::is_same<
      typename decltype (crs_vals_d)::non_const_value_type,
      std::complex<double> >::value,
    "crs_vals_d::non_const_value_type is std::complex<double>; this should "
    "never happen, since std::complex does not work in Kokkos::View objects.");
#endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE

  auto src_pids_d =
    create_mirror_view_from_raw_host_array(outputDevice, SourcePids.getRawPtr(),
        SourcePids.size(), true, "src_pids");

  auto tgt_pids_d =
    create_mirror_view_from_raw_host_array(outputDevice, TargetPids.getRawPtr(),
        TargetPids.size(), true, "tgt_pids");

# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::null;
# endif

  size_t bytes_per_value = 0;
  if (PackTraits<ST>::compileTimeSize) {
    // assume that ST is default constructible
    bytes_per_value = PackTraits<ST>::packValueCount(ST());
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
    size_t bytes_per_value_l = 0;
    if (local_matrix.values.extent(0) > 0) {
      const ST& val = local_matrix.values(0);
      bytes_per_value_l = PackTraits<ST>::packValueCount(val);
    } else {
      const ST& val = crs_vals_d(0);
      bytes_per_value_l = PackTraits<ST>::packValueCount(val);
    }
    Teuchos::reduceAll<int, size_t>(*(sourceMatrix.getComm()),
                                    Teuchos::REDUCE_MAX,
                                    bytes_per_value_l,
                                    outArg(bytes_per_value));
  }

#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  static_assert (! std::is_same<
      typename decltype (crs_vals_d)::non_const_value_type,
      std::complex<double> >::value,
    "crs_vals_d::non_const_value_type is std::complex<double>; this should "
    "never happen, since std::complex does not work in Kokkos::View objects.");
#endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE

# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("unpackAndCombineIntoCrsArrays"))));
# endif
  UnpackAndCombineCrsMatrixImpl::unpackAndCombineIntoCrsArrays(
      local_matrix, local_col_map, import_lids_d, imports_d,
      num_packets_per_lid_d, permute_to_lids_d, permute_from_lids_d,
      crs_rowptr_d, crs_colind_d, crs_vals_d, src_pids_d, tgt_pids_d,
      numSameIDs, TargetNumRows, TargetNumNonzeros, MyTargetPID,
      bytes_per_value);
# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::null;
# endif

  // Copy outputs back to host
# ifdef HAVE_TPETRA_MMM_TIMINGS
  tm = Teuchos::rcp(new TimeMonitor(*TimeMonitor::getNewTimer(prefix + std::string("copy back to host"))));
# endif
  typename decltype(crs_rowptr_d)::HostMirror crs_rowptr_h(
      CRS_rowptr.getRawPtr(), CRS_rowptr.size());
  // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
  deep_copy(execution_space(), crs_rowptr_h, crs_rowptr_d);

  typename decltype(crs_colind_d)::HostMirror crs_colind_h(
      CRS_colind.getRawPtr(), CRS_colind.size());
  // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
  deep_copy(execution_space(), crs_colind_h, crs_colind_d);

  typename decltype(crs_vals_d)::HostMirror crs_vals_h(
      CRS_vals_impl_scalar_type.getRawPtr(), CRS_vals_impl_scalar_type.size());
  // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
  deep_copy(execution_space(), crs_vals_h, crs_vals_d);

  typename decltype(tgt_pids_d)::HostMirror tgt_pids_h(
      TargetPids.getRawPtr(), TargetPids.size());
  // DEEP_COPY REVIEW - DEVICE-TO-HOSTMIRROR
  deep_copy(execution_space(), tgt_pids_h, tgt_pids_d);

} //unpackAndCombineIntoCrsArrays


} // namespace Details
} // namespace Tpetra

#define TPETRA_DETAILS_UNPACKCRSMATRIXANDCOMBINE_INSTANT( ST, LO, GO, NT ) \
  template void \
  Details::unpackCrsMatrixAndCombine<ST, LO, GO, NT> ( \
    const CrsMatrix<ST, LO, GO, NT>&, \
    const Teuchos::ArrayView<const char>&, \
    const Teuchos::ArrayView<const size_t>&, \
    const Teuchos::ArrayView<const LO>&, \
    size_t, \
    CombineMode); \
  template size_t \
  Details::unpackAndCombineWithOwningPIDsCount<ST, LO, GO, NT> ( \
    const CrsMatrix<ST, LO, GO, NT> &, \
    const Teuchos::ArrayView<const LO> &, \
    const Teuchos::ArrayView<const char> &, \
    const Teuchos::ArrayView<const size_t>&, \
    size_t, \
    CombineMode, \
    size_t, \
    const Teuchos::ArrayView<const LO>&, \
    const Teuchos::ArrayView<const LO>&); \
  template void \
  Details::unpackCrsMatrixAndCombineNew<ST, LO, GO, NT> ( \
    const CrsMatrix<ST, LO, GO, NT>&, \
    Kokkos::DualView<char*, typename DistObject<char, LO, GO, NT>::buffer_device_type>, \
    Kokkos::DualView<size_t*, typename DistObject<char, LO, GO, NT>::buffer_device_type>, \
    const Kokkos::DualView<const LO*, typename DistObject<char, LO, GO, NT>::buffer_device_type>&, \
    const size_t, \
    const CombineMode); \
  template void \
  Details::unpackAndCombineIntoCrsArrays<ST, LO, GO, NT> ( \
    const CrsMatrix<ST, LO, GO, NT> &, \
    const Kokkos::View<LO const *, \
                       Kokkos::Device<typename NT::device_type::execution_space, \
                                      Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>>,\
                       void, void >, \
    const Kokkos::View<const char*, \
                       Kokkos::Device<typename NT::device_type::execution_space, \
                                      Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>>, \
                       void, void >, \
    const Kokkos::View<const size_t*, \
                       Kokkos::Device<typename NT::device_type::execution_space, \
                                      Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>>, \
                       void, void >, \
    const size_t, \
    const Kokkos::View<LO const *, \
                       Kokkos::Device<typename NT::device_type::execution_space, \
                                      Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>>, \
                       void, void >, \
    const Kokkos::View<LO const *, \
                       Kokkos::Device<typename NT::device_type::execution_space, \
                                      Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>>, \
                       void, void >, \
    size_t, \
    const int, \
    Kokkos::View<size_t*,typename NT::device_type>&, \
    Kokkos::View<GO*,typename NT::device_type>&, \
    Kokkos::View<typename CrsMatrix<ST, LO, GO, NT>::impl_scalar_type*,typename NT::device_type>&, \
    const Teuchos::ArrayView<const int>&, \
    Kokkos::View<int*,typename NT::device_type>&); \
  template void \
  Details::unpackAndCombineIntoCrsArrays<ST, LO, GO, NT> ( \
    const CrsMatrix<ST, LO, GO, NT> &, \
    const Kokkos::View<LO const *, \
                       Kokkos::Device<typename NT::device_type::execution_space, \
                                      Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>>,\
                       void, void >, \
    const Kokkos::View<const char*, \
                       Kokkos::Device<typename NT::device_type::execution_space, \
                                      Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>>, \
                       void, void >, \
    const Kokkos::View<const size_t*, \
                       Kokkos::Device<typename NT::device_type::execution_space, \
                                      Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>>, \
                       void, void >, \
    const size_t, \
    const Kokkos::View<LO const *, \
                       Kokkos::Device<typename NT::device_type::execution_space, \
                                      Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>>, \
                       void, void >, \
    const Kokkos::View<LO const *, \
                       Kokkos::Device<typename NT::device_type::execution_space, \
                                      Tpetra::Details::DefaultTypes::comm_buffer_memory_space<typename NT::device_type>>, \
                       void, void >, \
    size_t, \
    const int, \
    Teuchos::ArrayRCP<size_t>&, \
    Teuchos::ArrayRCP<GO>&, \
    Teuchos::ArrayRCP<ST>&, \
    const Teuchos::ArrayView<const int>&, \
    Teuchos::Array<int>&);

#endif // TPETRA_DETAILS_UNPACKCRSMATRIXANDCOMBINE_DEF_HPP
