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

#ifndef TPETRA_DETAILS_UNPACKCRSMATRIXANDCOMBINE_DEF_HPP
#define TPETRA_DETAILS_UNPACKCRSMATRIXANDCOMBINE_DEF_HPP

#include "TpetraCore_config.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Tpetra_Details_castAwayConstDualView.hpp"
#include "Tpetra_Details_computeOffsets.hpp"
#include "Tpetra_Details_createMirrorView.hpp"
#include "Tpetra_Details_OrdinalTraits.hpp"
#include "Tpetra_Details_PackTraits.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_Details_getEntryOnHost.hpp"
#include "Kokkos_Core.hpp"
#include <memory>
#include <string>

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

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// Forward declaration of Distributor
class Distributor;
#endif // DOXYGEN_SHOULD_SKIP_THIS

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
/// \tparam DT The Kokkos device type.  See the documentation of Map
///   for requirements.
/// \tparam BDT The "buffer device type."
template<class ST, class LO, class GO, class DT, class BDT>
KOKKOS_FUNCTION int
unpackRow(typename PackTraits<GO, DT>::output_array_type& gids_out,
          typename PackTraits<int, DT>::output_array_type& pids_out,
          typename PackTraits<ST, DT>::output_array_type& vals_out,
          const Kokkos::View<const char*, BDT>& imports,
          const size_t offset,
          const size_t num_bytes,
          const size_t num_ent,
          const size_t num_bytes_per_value)
{
  if (num_ent == 0) {
    // Empty rows always take zero bytes, to ensure sparsity.
    return 0;
  }
  bool unpack_pids = pids_out.size() > 0;

  const size_t num_ent_beg = offset;
  const size_t num_ent_len = PackTraits<LO, BDT>::packValueCount (LO (0));

  const size_t gids_beg = num_ent_beg + num_ent_len;
  const size_t gids_len =
    num_ent * PackTraits<GO, BDT>::packValueCount (GO (0));

  const size_t pids_beg = gids_beg + gids_len;
  const size_t pids_len = unpack_pids ?
    size_t (num_ent * PackTraits<int, BDT>::packValueCount (int (0))) :
    size_t (0);

  const size_t vals_beg = gids_beg + gids_len + pids_len;
  const size_t vals_len = num_ent * num_bytes_per_value;

  const char* const num_ent_in = imports.data () + num_ent_beg;
  const char* const gids_in = imports.data () + gids_beg;
  const char* const pids_in = unpack_pids ? imports.data () + pids_beg : NULL;
  const char* const vals_in = imports.data () + vals_beg;

  size_t num_bytes_out = 0;
  LO num_ent_out;
  num_bytes_out += PackTraits<LO, BDT>::unpackValue (num_ent_out, num_ent_in);
  if (static_cast<size_t> (num_ent_out) != num_ent) {
    return 20; // error code
  }

  {
    Kokkos::pair<int, size_t> p;
    p = PackTraits<GO, BDT>::unpackArray (gids_out.data (), gids_in, num_ent);
    if (p.first != 0) {
      return 21; // error code
    }
    num_bytes_out += p.second;

    if (unpack_pids) {
      p = PackTraits<int, BDT>::unpackArray (pids_out.data (), pids_in, num_ent);
      if (p.first != 0) {
        return 22; // error code
      }
      num_bytes_out += p.second;
    }

    p = PackTraits<ST, BDT>::unpackArray (vals_out.data (), vals_in, num_ent);
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
}

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
  typedef Kokkos::View<const LO*, DT> import_lids_type;

  typedef Kokkos::View<LO*, DT> lids_scratch_type;
  typedef Kokkos::View<GO*, DT> gids_scratch_type;
  typedef Kokkos::View<int*,DT> pids_scratch_type;
  typedef Kokkos::View<ST*, DT> vals_scratch_type;

  typedef Kokkos::pair<int, LO> value_type;

  static_assert (std::is_same<LO, typename local_matrix_type::ordinal_type>::value,
                 "LocalMap::local_ordinal_type and "
                 "LocalMatrix::ordinal_type must be the same.");

  local_matrix_type local_matrix;
  local_map_type local_col_map;
  input_buffer_type imports;
  num_packets_per_lid_type num_packets_per_lid;
  import_lids_type import_lids;
  offsets_type offsets;
  Tpetra::CombineMode combine_mode;
  size_t max_num_ent;
  bool unpack_pids;
  size_t num_bytes_per_value;
  bool atomic;
  Kokkos::Experimental::UniqueToken<XS, Kokkos::Experimental::UniqueTokenScope::Global> tokens;
  lids_scratch_type lids_scratch;
  gids_scratch_type gids_scratch;
  pids_scratch_type pids_scratch;
  vals_scratch_type vals_scratch;

  UnpackCrsMatrixAndCombineFunctor(
      const local_matrix_type& local_matrix_in,
      const local_map_type& local_col_map_in,
      const input_buffer_type& imports_in,
      const num_packets_per_lid_type& num_packets_per_lid_in,
      const import_lids_type& import_lids_in,
      const offsets_type& offsets_in,
      const Tpetra::CombineMode combine_mode_in,
      const size_t max_num_ent_in,
      const bool unpack_pids_in,
      const size_t num_bytes_per_value_in,
      const bool atomic_in) :
    local_matrix (local_matrix_in),
    local_col_map (local_col_map_in),
    imports (imports_in),
    num_packets_per_lid (num_packets_per_lid_in),
    import_lids (import_lids_in),
    offsets (offsets_in),
    combine_mode (combine_mode_in),
    max_num_ent (max_num_ent_in),
    unpack_pids (unpack_pids_in),
    num_bytes_per_value (num_bytes_per_value_in),
    atomic (atomic_in),
    tokens (XS()),
    lids_scratch ("pids_scratch", tokens.size() * max_num_ent),
    gids_scratch ("gids_scratch", tokens.size() * max_num_ent),
    pids_scratch ("lids_scratch", tokens.size() * max_num_ent),
    vals_scratch ("vals_scratch", tokens.size() * max_num_ent)
  {}

  KOKKOS_INLINE_FUNCTION void init(value_type& dst) const
  {
    using Tpetra::Details::OrdinalTraits;
    dst = Kokkos::make_pair (0, OrdinalTraits<LO>::invalid ());
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst, const volatile value_type& src) const
  {
    // `dst` should reflect the first (least) bad index and
    // all other associated error codes and data.  Thus, we need only
    // check if the `src` object shows an error and if its associated
    // bad index is less than `dst`'s bad index.
    using Tpetra::Details::OrdinalTraits;
    if (src.second != OrdinalTraits<LO>::invalid ()) {
      // An error in the src; check if
      //   1. `dst` shows errors
      //   2. If `dst` does show errors, if src's bad index is less than
      //      *this' bad index
      if (dst.second == OrdinalTraits<LO>::invalid () ||
          src.second < dst.second) {
        dst = src;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const LO i, value_type& dst) const
  {
    using Kokkos::View;
    using Kokkos::subview;
    using Kokkos::MemoryUnmanaged;
    typedef typename XS::size_type size_type;
    typedef typename Kokkos::pair<size_type, size_type> slice;
    typedef BufferDeviceType BDT;

    typedef View<LO*, DT, MemoryUnmanaged> lids_out_type;
    typedef View<int*,DT, MemoryUnmanaged> pids_out_type;
    typedef View<GO*, DT, MemoryUnmanaged> gids_out_type;
    typedef View<ST*, DT, MemoryUnmanaged> vals_out_type;

    const size_t num_bytes = num_packets_per_lid(i);

    // Only unpack data if there is a nonzero number of bytes.
    if (num_bytes == 0) {
      return;
    }

    // there is actually something in the row
    const LO import_lid = import_lids[i];
    const size_t buf_size = imports.size();
    const size_t offset = offsets(i);

    // Get the number of entries to expect in the received data for this row.
    LO num_ent_LO = 0;
    const char* const in_buf = imports.data () + offset;
    (void) PackTraits<LO, BDT>::unpackValue (num_ent_LO, in_buf);
    const size_t num_ent = static_cast<size_t> (num_ent_LO);

    // Count the number of bytes expected to unpack
    size_t expected_num_bytes = 0;
    {
      expected_num_bytes += PackTraits<LO, BDT>::packValueCount (LO (0));
      expected_num_bytes += num_ent * PackTraits<GO, BDT>::packValueCount (GO (0));
      if (unpack_pids) {
        expected_num_bytes += num_ent * PackTraits<int, BDT>::packValueCount (int (0));
      }
      expected_num_bytes += num_ent * PackTraits<ST, BDT>::packValueCount (ST ());
    }

    if (expected_num_bytes > num_bytes) {
      dst = Kokkos::make_pair (1, i); // wrong number of bytes
      return;
    }

    if (offset > buf_size || offset + num_bytes > buf_size) {
      dst = Kokkos::make_pair (2, i); // out of bounds
      return;
    }

    // Get subviews in to the scratch arrays.  The token returned from acquire
    // is an integer in [0, tokens.size()).  It is used to grab a unique (to
    // this thread) subview of the scratch arrays.
    const size_type token = tokens.acquire();
    const size_t a = static_cast<size_t>(token) * max_num_ent;
    const size_t b = a + num_ent;
    lids_out_type lids_out = subview(lids_scratch, slice(a, b));
    gids_out_type gids_out = subview(gids_scratch, slice(a, b));
    pids_out_type pids_out = subview(pids_scratch, slice(a, (unpack_pids ? b : a)));
    vals_out_type vals_out = subview(vals_scratch, slice(a, b));

    // Unpack this row!
    int unpack_err =
      unpackRow<ST,LO,GO,DT,BDT>(gids_out, pids_out, vals_out,
                                 imports, offset, num_bytes,
                                 num_ent, num_bytes_per_value);
    if (unpack_err != 0) {
      dst = Kokkos::make_pair (unpack_err, i); // unpack error
      tokens.release (token);
      return;
    }

    // Column indices come in as global indices, in case the
    // source object's column Map differs from the target object's
    // (this's) column Map, and must be converted local index values
    for (size_t k = 0; k < num_ent; ++k) {
      lids_out(k) = local_col_map.getLocalElement (gids_out(k));
    }

    // Combine the values according to the combine_mode
    const LO* const lids_raw = const_cast<const LO*> (lids_out.data ());
    const ST* const vals_raw = const_cast<const ST*> (vals_out.data ());
    LO num_modified = 0;
    if (combine_mode == ADD) {
      num_modified +=
        local_matrix.sumIntoValues (import_lid, lids_raw, num_ent,
                                    vals_raw, false, atomic);
    }
    else if (combine_mode == REPLACE) {
      num_modified +=
        local_matrix.replaceValues (import_lid, lids_raw, num_ent,
                                    vals_raw, false, atomic);
    }
    else {
      dst = Kokkos::make_pair (4, i); // invalid combine mode
      tokens.release (token);
      return;
    }

    tokens.release (token);
  }
};

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
  typedef Kokkos::pair<size_t,size_t> slice;

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
      (void) PackTraits<LO, BDT>::unpackValue (num_ent_LO, in_buf);
      const size_t num_ent = static_cast<size_t> (num_ent_LO);

      update = (update < num_ent) ? num_ent : update;
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (const MaxNumEntTag,
        volatile value_type& dst,
        const volatile value_type& src) const
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
      (void) PackTraits<LO, BDT>::unpackValue (num_ent_LO, in_buf);
      tot_num_ent += static_cast<size_t> (num_ent_LO);
    }
  }
};

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
    const typename PackTraits<typename LocalMap::local_ordinal_type, typename LocalMap::device_type>::input_array_type import_lids,
    const Tpetra::CombineMode combine_mode,
    const bool unpack_pids,
    const bool atomic)
{
  typedef typename LocalMatrix::value_type ST;
  typedef typename LocalMap::local_ordinal_type LO;
  typedef typename LocalMap::device_type DT;
  typedef typename DT::execution_space XS;
  typedef Kokkos::RangePolicy<XS, Kokkos::IndexType<LO> > range_policy;
  typedef UnpackCrsMatrixAndCombineFunctor<LocalMatrix, LocalMap, BufferDeviceType> unpack_functor_type;

  const char prefix[] =
    "Tpetra::Details::UnpackAndCombineCrsMatrixImpl::unpackAndCombineIntoCrsMatrix: ";

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

  // Determine the maximum number of entries in any row in the matrix.  The
  // maximum number of entries is needed to allocate unpack buffers on the
  // device.
  size_t max_num_ent = compute_maximum_num_entries<LO,DT>(
      num_packets_per_lid, offsets, imports);

  // FIXME (TJF SEP 2017)
  // The scalar type is not necessarily default constructible
  size_t num_bytes_per_value = PackTraits<ST, DT>::packValueCount(ST());

  // Now do the actual unpack!
  unpack_functor_type f(local_matrix, local_map,
      imports, num_packets_per_lid, import_lids, offsets, combine_mode,
      max_num_ent, unpack_pids, num_bytes_per_value, atomic);

  typename unpack_functor_type::value_type x;
  Kokkos::parallel_reduce(range_policy(0, static_cast<LO>(num_import_lids)), f, x);
  auto x_h = x.to_std_pair();
  TEUCHOS_TEST_FOR_EXCEPTION(x_h.first != 0, std::runtime_error,
      prefix << "UnpackCrsMatrixAndCombineFunctor reported error code "
             << x_h.first << " for the first bad row " << x_h.second);

  return;
}

template<class LocalMatrix, class BufferDeviceType>
size_t
unpackAndCombineWithOwningPIDsCount(
  const LocalMatrix& local_matrix,
  const typename PackTraits<typename LocalMatrix::ordinal_type, typename LocalMatrix::device_type>::input_array_type permute_from_lids,
  const Kokkos::View<const char*, BufferDeviceType>& imports,
  const Kokkos::View<const size_t*, BufferDeviceType>& num_packets_per_lid,
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
}

template<class LO, class DT, class BDT>
KOKKOS_INLINE_FUNCTION
size_t
unpackRowCount(const Kokkos::View<const char*, BDT>& imports,
               const size_t offset,
               const size_t num_bytes)
{
  LO num_ent_LO = 0;
  if (num_bytes > 0) {
    const size_t p_num_bytes = PackTraits<LO, DT>::packValueCount(num_ent_LO);
    if (p_num_bytes > num_bytes) {
      return OrdinalTraits<size_t>::invalid();
    }
    const char* const in_buf = imports.data () + offset;
    (void) PackTraits<LO,DT>::unpackValue(num_ent_LO, in_buf);
  }
  return static_cast<size_t>(num_ent_LO);
}

/// \brief Setup row pointers for remotes
template<class LO, class DT, class BDT>
int
setupRowPointersForRemotes(
  const typename PackTraits<size_t, DT>::output_array_type& tgt_rowptr,
  const typename PackTraits<LO, DT>::input_array_type& import_lids,
  const Kokkos::View<const char*, BDT>& imports,
  const Kokkos::View<const size_t*, BDT>& num_packets_per_lid,
  const typename PackTraits<size_t, DT>::input_array_type& offsets)
{
  using Kokkos::parallel_reduce;
  typedef typename DT::execution_space XS;
  typedef typename PackTraits<size_t,DT>::input_array_type::size_type size_type;
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
      const size_t num_ent = unpackRowCount<LO, DT, BDT> (imports, offset, num_bytes);
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
    const typename PackTraits<size_t,DT>::output_array_type& tgt_rowptr,
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
    const typename PackTraits<typename LocalMap::global_ordinal_type, typename LocalMap::device_type>::output_array_type& tgt_colind,
    const typename PackTraits<int, typename LocalMap::device_type>::output_array_type& tgt_pids,
    const typename PackTraits<typename LocalMatrix::value_type, typename LocalMap::device_type>::output_array_type& tgt_vals,
    const Kokkos::View<size_t*,typename LocalMap::device_type>& new_start_row,
    const typename PackTraits<size_t, typename LocalMap::device_type>::output_array_type& tgt_rowptr,
    const typename PackTraits<int, typename LocalMap::device_type>::input_array_type& src_pids,
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
    const typename PackTraits<typename LocalMap::global_ordinal_type, typename LocalMap::device_type>::output_array_type& tgt_colind,
    const typename PackTraits<int, typename LocalMap::device_type>::output_array_type& tgt_pids,
    const typename PackTraits<typename LocalMatrix::value_type, typename LocalMap::device_type>::output_array_type& tgt_vals,
    const Kokkos::View<size_t*,typename LocalMap::device_type>& new_start_row,
    const typename PackTraits<size_t, typename LocalMap::device_type>::output_array_type& tgt_rowptr,
    const typename PackTraits<int, typename LocalMap::device_type>::input_array_type& src_pids,
    const typename PackTraits<typename LocalMap::local_ordinal_type, typename LocalMap::device_type>::input_array_type& permute_to_lids,
    const typename PackTraits<typename LocalMap::local_ordinal_type, typename LocalMap::device_type>::input_array_type& permute_from_lids,
    const LocalMatrix& local_matrix,
    const LocalMap& local_col_map,
    const int my_pid)
{
  using Kokkos::parallel_for;
  typedef typename LocalMap::device_type DT;
  typedef typename LocalMap::local_ordinal_type LO;
  typedef typename DT::execution_space XS;
  typedef typename PackTraits<LO,DT>::input_array_type::size_type size_type;
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
    const typename PackTraits<typename LocalMap::global_ordinal_type, typename LocalMap::device_type>::output_array_type& tgt_colind,
    const typename PackTraits<int, typename LocalMap::device_type>::output_array_type& tgt_pids,
    const typename PackTraits<typename LocalMatrix::value_type, typename LocalMap::device_type>::output_array_type& tgt_vals,
    const Kokkos::View<size_t*,typename LocalMap::device_type>& new_start_row,
    const typename PackTraits<size_t, typename LocalMap::device_type>::input_array_type& offsets,
    const typename PackTraits<typename LocalMap::local_ordinal_type, typename LocalMap::device_type>::input_array_type& import_lids,
    const Kokkos::View<const char*, BufferDeviceType>& imports,
    const Kokkos::View<const size_t*, BufferDeviceType>& num_packets_per_lid,
    const LocalMatrix& local_matrix,
    const LocalMap /*& local_col_map*/,
    const int my_pid,
    const size_t num_bytes_per_value)
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
  typedef BufferDeviceType BDT;

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
      size_t num_ent = unpackRowCount<LO,DT,BDT>(imports, offset, num_bytes);
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

      k_error += unpackRow<ST,LO,GO,DT,BDT>(gids_out, pids_out, vals_out,
                                            imports, offset, num_bytes,
                                            num_ent, num_bytes_per_value);

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
    const typename PackTraits<typename LocalMap::local_ordinal_type, typename LocalMap::device_type>::input_array_type& import_lids,
    const Kokkos::View<const char*, BufferDeviceType>& imports,
    const Kokkos::View<const size_t*, BufferDeviceType>& num_packets_per_lid,
    const typename PackTraits<typename LocalMap::local_ordinal_type, typename LocalMap::device_type>::input_array_type& permute_to_lids,
    const typename PackTraits<typename LocalMap::local_ordinal_type, typename LocalMap::device_type>::input_array_type& permute_from_lids,
    const typename PackTraits<size_t, typename LocalMap::device_type>::output_array_type& tgt_rowptr,
    const typename PackTraits<typename LocalMap::global_ordinal_type, typename LocalMap::device_type>::output_array_type& tgt_colind,
    const typename PackTraits<typename LocalMatrix::value_type, typename LocalMap::device_type>::output_array_type& tgt_vals,
    const typename PackTraits<int, typename LocalMap::device_type>::input_array_type& src_pids,
    const typename PackTraits<int, typename LocalMap::device_type>::output_array_type& tgt_pids,
    const size_t num_same_ids,
    const size_t tgt_num_rows,
    const size_t tgt_num_nonzeros,
    const int my_tgt_pid,
    const size_t num_bytes_per_value)
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
  const size_t mynnz = tgt_num_nonzeros;

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
  {
    auto nth_tgt_rowptr_h = getEntryOnHost(tgt_rowptr, N);
    bool condition = nth_tgt_rowptr_h != mynnz;
    TEUCHOS_TEST_FOR_EXCEPTION(condition, std::invalid_argument,
      prefix << "CRS_rowptr[last] = " <<
      nth_tgt_rowptr_h << "!= mynnz = " << mynnz << ".");
  }

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
      local_matrix, local_col_map, my_pid, num_bytes_per_value);
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
    size_t constantNumPackets,
    Distributor & distor,
    CombineMode combineMode,
    const bool atomic)
{
  using Kokkos::View;
  typedef typename Node::device_type device_type;
  typedef typename CrsMatrix<ST, LO, GO, Node>::local_matrix_type local_matrix_type;
  static_assert (std::is_same<device_type, typename local_matrix_type::device_type>::value,
                 "Node::device_type and LocalMatrix::device_type must be the same.");

  // Execution space.
  typedef typename device_type::execution_space XS;

  // Convert all Teuchos::Array to Kokkos::View.
  typename XS::device_type outputDevice;

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

  auto local_matrix = sourceMatrix.getLocalMatrix();
  auto local_col_map = sourceMatrix.getColMap()->getLocalMap();

  // Now do the actual unpack!
  UnpackAndCombineCrsMatrixImpl::unpackAndCombineIntoCrsMatrix(
      local_matrix, local_col_map, imports_d, num_packets_per_lid_d,
      import_lids_d, combineMode, false, atomic);

  return;
}

template<typename ST, typename LO, typename GO, typename NT>
void
unpackCrsMatrixAndCombineNew (const CrsMatrix<ST, LO, GO, NT>& sourceMatrix,
                              const Kokkos::DualView<const char*, typename DistObject<char, LO, GO, NT>::buffer_device_type>& imports,
                              const Kokkos::DualView<const size_t*, typename DistObject<char, LO, GO, NT>::buffer_device_type>& numPacketsPerLID,
                              const Kokkos::DualView<const LO*, typename NT::device_type>& importLIDs,
                              const size_t constantNumPackets,
                              Distributor& distor,
                              const CombineMode combineMode,
                              const bool atomic)
{
  using Tpetra::Details::castAwayConstDualView;
  using Kokkos::View;
  typedef typename NT::device_type device_type;
  typedef CrsMatrix<ST, LO, GO, NT> crs_matrix_type;
  typedef typename crs_matrix_type::local_matrix_type local_matrix_type;
  typedef DistObject<char, LO, GO, NT> dist_object_type;
  typedef typename dist_object_type::buffer_device_type buffer_device_type;
  typedef typename buffer_device_type::memory_space BMS;
  typedef typename device_type::memory_space MS;

  static_assert (std::is_same<device_type,
                   typename local_matrix_type::device_type>::value,
                 "NT::device_type and LocalMatrix::device_type must be "
                 "the same.");

  {
    auto numPacketsPerLID_nc = castAwayConstDualView (numPacketsPerLID);
    numPacketsPerLID_nc.template sync<BMS> ();
  }
  auto num_packets_per_lid_d = numPacketsPerLID.template view<BMS> ();

  {
    auto importLIDs_nc = castAwayConstDualView (importLIDs);
    importLIDs_nc.template sync<MS> ();
  }
  auto import_lids_d = importLIDs.template view<MS> ();

  {
    auto imports_nc = castAwayConstDualView (imports);
    imports_nc.template sync<BMS> ();
  }
  auto imports_d = imports.template view<BMS> ();

  auto local_matrix = sourceMatrix.getLocalMatrix ();
  auto local_col_map = sourceMatrix.getColMap ()->getLocalMap ();
  typedef decltype (local_col_map) local_map_type;

  // Now do the actual unpack!
  UnpackAndCombineCrsMatrixImpl::unpackAndCombineIntoCrsMatrix<
      local_matrix_type,
      local_map_type,
      buffer_device_type
    > (local_matrix, local_col_map, imports_d, num_packets_per_lid_d,
       import_lids_d, combineMode, false, atomic);
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
    size_t constantNumPackets,
    Distributor &distor,
    CombineMode combineMode,
    size_t numSameIDs,
    const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
    const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs)
{
  using Kokkos::MemoryUnmanaged;
  using Kokkos::View;
  typedef typename Node::device_type DT;
  typedef typename DistObject<char, LocalOrdinal, GlobalOrdinal, Node>::buffer_device_type BDT;
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

  auto local_matrix = sourceMatrix.getLocalMatrix ();
  auto permute_from_lids_d =
    create_mirror_view_from_raw_host_array (DT (),
                                            permuteFromLIDs.getRawPtr (),
                                            permuteFromLIDs.size (), true,
                                            "permute_from_lids");
  auto imports_d =
    create_mirror_view_from_raw_host_array (BDT (),
                                            imports.getRawPtr (),
                                            imports.size (), true,
                                            "imports");
  auto num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array (BDT (),
                                            numPacketsPerLID.getRawPtr (),
                                            numPacketsPerLID.size (), true,
                                            "num_packets_per_lid");

  return UnpackAndCombineCrsMatrixImpl::unpackAndCombineWithOwningPIDsCount(
      local_matrix, permute_from_lids_d, imports_d,
      num_packets_per_lid_d, numSameIDs);
}

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
    const Teuchos::ArrayView<const LocalOrdinal>& importLIDs,
    const Teuchos::ArrayView<const char>& imports,
    const Teuchos::ArrayView<const size_t>& numPacketsPerLID,
    const size_t constantNumPackets,
    Distributor& distor,
    const CombineMode combineMode,
    const size_t numSameIDs,
    const Teuchos::ArrayView<const LocalOrdinal>& permuteToLIDs,
    const Teuchos::ArrayView<const LocalOrdinal>& permuteFromLIDs,
    size_t TargetNumRows,
    size_t TargetNumNonzeros,
    const int MyTargetPID,
    const Teuchos::ArrayView<size_t>& CRS_rowptr,
    const Teuchos::ArrayView<GlobalOrdinal>& CRS_colind,
    const Teuchos::ArrayView<typename CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::impl_scalar_type>& CRS_vals,
    const Teuchos::ArrayView<const int>& SourcePids,
    Teuchos::Array<int>& TargetPids)
{
  using Tpetra::Details::PackTraits;

  using Kokkos::View;
  using Kokkos::deep_copy;

  using Teuchos::ArrayView;
  using Teuchos::outArg;
  using Teuchos::REDUCE_MAX;
  using Teuchos::reduceAll;

  typedef LocalOrdinal LO;

  typedef typename Node::device_type DT;
  typedef typename DT::execution_space XS;

  typedef CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef typename matrix_type::impl_scalar_type ST;
  typedef typename ArrayView<const LO>::size_type size_type;

  const char prefix[] = "Tpetra::Details::unpackAndCombineIntoCrsArrays: ";

  TEUCHOS_TEST_FOR_EXCEPTION(
    TargetNumRows + 1 != static_cast<size_t> (CRS_rowptr.size ()),
    std::invalid_argument, prefix << "CRS_rowptr.size() = " <<
    CRS_rowptr.size () << "!= TargetNumRows+1 = " << TargetNumRows+1 << ".");

  TEUCHOS_TEST_FOR_EXCEPTION(
    permuteToLIDs.size () != permuteFromLIDs.size (), std::invalid_argument,
    prefix << "permuteToLIDs.size() = " << permuteToLIDs.size ()
    << "!= permuteFromLIDs.size() = " << permuteFromLIDs.size () << ".");
  const size_type numImportLIDs = importLIDs.size ();

  TEUCHOS_TEST_FOR_EXCEPTION(
    numImportLIDs != numPacketsPerLID.size (), std::invalid_argument,
    prefix << "importLIDs.size() = " << numImportLIDs << " != "
    "numPacketsPerLID.size() = " << numPacketsPerLID.size() << ".");

  // Preseed TargetPids with -1 for local
  if (static_cast<size_t> (TargetPids.size ()) != TargetNumNonzeros) {
    TargetPids.resize (TargetNumNonzeros);
  }
  TargetPids.assign (TargetNumNonzeros, -1);

  // Grab pointers for sourceMatrix
  auto local_matrix = sourceMatrix.getLocalMatrix();
  auto local_col_map = sourceMatrix.getColMap()->getLocalMap();

  // Convert input arrays to Kokkos::View
  typename XS::device_type outputDevice;
  auto import_lids_d =
    create_mirror_view_from_raw_host_array(outputDevice, importLIDs.getRawPtr(),
        importLIDs.size(), true, "import_lids");

  auto imports_d =
    create_mirror_view_from_raw_host_array(outputDevice, imports.getRawPtr(),
        imports.size(), true, "imports");

  auto num_packets_per_lid_d =
    create_mirror_view_from_raw_host_array(outputDevice, numPacketsPerLID.getRawPtr(),
        numPacketsPerLID.size(), true, "num_packets_per_lid");

  auto permute_from_lids_d =
    create_mirror_view_from_raw_host_array(outputDevice, permuteFromLIDs.getRawPtr(),
        permuteFromLIDs.size(), true, "permute_from_lids");

  auto permute_to_lids_d =
    create_mirror_view_from_raw_host_array(outputDevice, permuteToLIDs.getRawPtr(),
        permuteToLIDs.size(), true, "permute_to_lids");

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
          decltype (CRS_vals)
        >::type::value_type
      >::type,
      std::complex<double> >::value,
    "CRS_vals::value_type is std::complex<double>; this should never happen"
    ", since std::complex does not work in Kokkos::View objects.");
#endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE

  auto crs_vals_d =
    create_mirror_view_from_raw_host_array(outputDevice, CRS_vals.getRawPtr(),
        CRS_vals.size(), true, "crs_vals");

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

  size_t num_bytes_per_value = 0;
  if (PackTraits<ST,DT>::compileTimeSize) {
    // assume that ST is default constructible
    num_bytes_per_value = PackTraits<ST,DT>::packValueCount(ST());
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
    if (local_matrix.values.extent(0) > 0) {
      const ST& val = local_matrix.values(0);
      num_bytes_per_value_l = PackTraits<ST,DT>::packValueCount(val);
    } else {
      const ST& val = crs_vals_d(0);
      num_bytes_per_value_l = PackTraits<ST,DT>::packValueCount(val);
    }
    Teuchos::reduceAll<int, size_t>(*(sourceMatrix.getComm()),
                                    Teuchos::REDUCE_MAX,
                                    num_bytes_per_value_l,
                                    outArg(num_bytes_per_value));
  }

#ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
  static_assert (! std::is_same<
      typename decltype (crs_vals_d)::non_const_value_type,
      std::complex<double> >::value,
    "crs_vals_d::non_const_value_type is std::complex<double>; this should "
    "never happen, since std::complex does not work in Kokkos::View objects.");
#endif // HAVE_TPETRA_INST_COMPLEX_DOUBLE

  UnpackAndCombineCrsMatrixImpl::unpackAndCombineIntoCrsArrays(
      local_matrix, local_col_map, import_lids_d, imports_d,
      num_packets_per_lid_d, permute_to_lids_d, permute_from_lids_d,
      crs_rowptr_d, crs_colind_d, crs_vals_d, src_pids_d, tgt_pids_d,
      numSameIDs, TargetNumRows, TargetNumNonzeros, MyTargetPID,
      num_bytes_per_value);

  // Copy outputs back to host
  typename decltype(crs_rowptr_d)::HostMirror crs_rowptr_h(
      CRS_rowptr.getRawPtr(), CRS_rowptr.size());
  deep_copy(crs_rowptr_h, crs_rowptr_d);

  typename decltype(crs_colind_d)::HostMirror crs_colind_h(
      CRS_colind.getRawPtr(), CRS_colind.size());
  deep_copy(crs_colind_h, crs_colind_d);

  typename decltype(crs_vals_d)::HostMirror crs_vals_h(
      CRS_vals.getRawPtr(), CRS_vals.size());
  deep_copy(crs_vals_h, crs_vals_d);

  typename decltype(tgt_pids_d)::HostMirror tgt_pids_h(
      TargetPids.getRawPtr(), TargetPids.size());
  deep_copy(tgt_pids_h, tgt_pids_d);

}

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
    Distributor&, \
    CombineMode, \
    const bool); \
  template void \
  Details::unpackCrsMatrixAndCombineNew<ST, LO, GO, NT> ( \
    const CrsMatrix<ST, LO, GO, NT>&, \
    const Kokkos::DualView<const char*, typename DistObject<char, LO, GO, NT>::buffer_device_type>&, \
    const Kokkos::DualView<const size_t*, typename DistObject<char, LO, GO, NT>::buffer_device_type>&, \
    const Kokkos::DualView<const LO*, NT::device_type>&, \
    const size_t, \
    Distributor&, \
    const CombineMode, \
    const bool); \
  template void \
  Details::unpackAndCombineIntoCrsArrays<ST, LO, GO, NT> ( \
    const CrsMatrix<ST, LO, GO, NT> &, \
    const Teuchos::ArrayView<const LO>&, \
    const Teuchos::ArrayView<const char>&, \
    const Teuchos::ArrayView<const size_t>&, \
    const size_t, \
    Distributor&, \
    const CombineMode, \
    const size_t, \
    const Teuchos::ArrayView<const LO>&, \
    const Teuchos::ArrayView<const LO>&, \
    size_t, \
    size_t, \
    const int, \
    const Teuchos::ArrayView<size_t>&, \
    const Teuchos::ArrayView<GO>&, \
    const Teuchos::ArrayView<CrsMatrix<ST, LO, GO, NT>::impl_scalar_type>&, \
    const Teuchos::ArrayView<const int>&, \
    Teuchos::Array<int>&); \
  template size_t \
  Details::unpackAndCombineWithOwningPIDsCount<ST, LO, GO, NT> ( \
    const CrsMatrix<ST, LO, GO, NT> &, \
    const Teuchos::ArrayView<const LO> &, \
    const Teuchos::ArrayView<const char> &, \
    const Teuchos::ArrayView<const size_t>&, \
    size_t, \
    Distributor &, \
    CombineMode, \
    size_t, \
    const Teuchos::ArrayView<const LO>&, \
    const Teuchos::ArrayView<const LO>&);

#endif // TPETRA_DETAILS_UNPACKCRSMATRIXANDCOMBINE_DEF_HPP
