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
#include "Tpetra_Details_PackTraits.hpp"
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

/// \brief Unpack a single row of a CrsMatrix
///
/// \tparam ST The type of the numerical entries of the matrix.
///   (You can use real-valued or complex-valued types here, unlike
///   in Epetra, where the scalar type is always \c double.)
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam GO The type of global indices.  See the
///   documentation of Map for requirements.
/// \tparam Device The Kokkos Device type.  See the documentation of Map
///   for requirements.
template<class ST, class LO, class GO, class Device>
KOKKOS_FUNCTION int
unpack_crs_matrix_row(typename PackTraits<GO, Device>::output_array_type& gids_out,
                      typename PackTraits<int, Device>::output_array_type& pids_out,
                      typename PackTraits<ST, Device>::output_array_type& vals_out,
                      const typename PackTraits<LO, Device>::input_buffer_type& imports,
                      const size_t offset, const size_t num_bytes, const size_t num_ent,
                      size_t num_bytes_per_value=0)
{
  using Kokkos::subview;

  // NOTE (mfh 02 Feb 2015) This assumes that output_buffer_type is
  // the same, no matter what type we're packing.  It's a reasonable
  // assumption, given that we go through the trouble of PackTraits
  // just so that we can pack data of different types in the same
  // buffer.
  typedef typename PackTraits<LO, Device>::input_buffer_type input_buffer_type;
  typedef typename input_buffer_type::size_type size_type;
  typedef typename Kokkos::pair<size_type, size_type> slice;

  if (num_ent == 0) {
    // Empty rows always take zero bytes, to ensure sparsity.
    return 0;
  }
  bool unpack_pids = pids_out.size() > 0;

  const LO num_ent_LO = static_cast<LO>(num_ent); // packValueCount wants this
  const size_t num_ent_beg = offset;
  const size_t num_ent_len = PackTraits<LO, Device>::packValueCount(num_ent_LO);

  const GO gid = 0; // packValueCount wants this
  const size_t gids_beg = num_ent_beg + num_ent_len;
  const size_t gids_len = num_ent * PackTraits<GO, Device>::packValueCount(gid);

  const int pid = 0; // packValueCount wants this
  const size_t pids_beg = gids_beg + gids_len;
  const size_t pids_len = (unpack_pids) ? num_ent * PackTraits<int, Device>::packValueCount(pid) : 0;

  if (num_bytes_per_value == 0) {
    ST val; // packValueCount wants this
    num_bytes_per_value = PackTraits<ST, Device>::packValueCount(val);
  }
  const size_t vals_beg = gids_beg + gids_len + pids_len;
  const size_t vals_len = num_ent * num_bytes_per_value;

  input_buffer_type num_ent_in =
    subview(imports, slice(num_ent_beg, num_ent_beg + num_ent_len));
  input_buffer_type gids_in =
    subview(imports, slice(gids_beg, gids_beg + gids_len));
  input_buffer_type pids_in;
  if (unpack_pids) {
    pids_in = subview(imports, slice(pids_beg, pids_beg + pids_len));
  }
  input_buffer_type vals_in =
    subview(imports, slice(vals_beg, vals_beg + vals_len));

  size_t num_bytes_out = 0;

  LO num_ent_out;
  num_bytes_out += PackTraits<LO, Device>::unpackValue(num_ent_out, num_ent_in);

  if (static_cast<size_t>(num_ent_out) != num_ent) {
    return 20;
  }

  {
    Kokkos::pair<int, size_t> p;
    p = PackTraits<GO, Device>::unpackArray(gids_out, gids_in, num_ent);
    if (p.first != 0) return 21;
    num_bytes_out += p.second;

    if (unpack_pids) {
      p = PackTraits<int, Device>::unpackArray(pids_out, pids_in, num_ent);
      if (p.first != 0) return 22;
      num_bytes_out += p.second;
    }

    p = PackTraits<ST, Device>::unpackArray(vals_out, vals_in, num_ent);
    if (p.first != 0) return 23;
    num_bytes_out += p.second;
  }

  const size_t expected_num_bytes = num_ent_len + gids_len + pids_len + vals_len;

  if (num_bytes_out != expected_num_bytes) return 24;

  return 0;
}

/// \brief Unpacks and combines a single row of the CrsMatrix.
///
/// \tparam LocalMatrix the specialization of the KokkosSparse::CrsMatrix
///   local matrix
/// \tparam LocalMap the type of the local column map
///
/// Data (bytes) describing the row of the CRS matrix are "unpacked"
/// from a single (concatenated) (view of) char* directly in to the row of the matrix
template<class LocalMatrix, class LocalMap>
struct UnpackCrsMatrixAndCombineFunctor {

  typedef LocalMatrix local_matrix_type;
  typedef LocalMap local_map_type;

  typedef typename local_matrix_type::value_type ST;
  typedef typename local_map_type::local_ordinal_type LO;
  typedef typename local_map_type::global_ordinal_type GO;
  typedef typename local_map_type::device_type Device;
  typedef typename Device::execution_space execution_space;

  typedef typename PackTraits<size_t,Device>::input_array_type num_packets_per_lid_type;
  typedef typename PackTraits<size_t,Device>::input_array_type offsets_type;
  typedef typename PackTraits<LO,Device>::input_buffer_type input_buffer_type;
  typedef typename PackTraits<LO,Device>::input_array_type import_lids_type;

  typedef Kokkos::View<LO*, Device> lids_scratch_type;
  typedef Kokkos::View<GO*, Device> gids_scratch_type;
  typedef Kokkos::View<int*,Device> pids_scratch_type;
  typedef Kokkos::View<ST*, Device> vals_scratch_type;

  typedef Kokkos::pair<int, LO> value_type;

  static_assert (std::is_same<LO, typename local_matrix_type::ordinal_type>::value,
                 "LocalMap::local_ordinal_type and "
                 "LocalMatrix::ordinal_type must be the same.");

  LO InvalidIndex = OrdinalTraits<LO>::invalid();

  local_matrix_type local_matrix;
  local_map_type local_col_map;
  input_buffer_type imports;
  num_packets_per_lid_type num_packets_per_lid;
  import_lids_type import_lids;
  offsets_type offsets;
  Tpetra::CombineMode combine_mode;
  size_t max_num_ent;
  bool unpack_pids;
  bool atomic;
  Kokkos::Experimental::UniqueToken<execution_space , Kokkos::Experimental::UniqueTokenScope::Global> tokens;
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
    atomic (atomic_in),
    tokens (execution_space()),
    lids_scratch ("pids_scratch", tokens.size() * max_num_ent),
    gids_scratch ("gids_scratch", tokens.size() * max_num_ent),
    pids_scratch ("lids_scratch", tokens.size() * max_num_ent),
    vals_scratch ("vals_scratch", tokens.size() * max_num_ent)
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
  void operator()(const LO i, value_type& dst) const
  {
    using Kokkos::View;
    using Kokkos::subview;
    using Kokkos::MemoryUnmanaged;
    typedef typename PackTraits<LO,Device>::input_buffer_type input_buffer_type;
    typedef typename execution_space::size_type size_type;
    typedef typename Kokkos::pair<size_type, size_type> slice;

    typedef View<LO*, Device, MemoryUnmanaged> lids_out_type;
    typedef View<int*,Device, MemoryUnmanaged> pids_out_type;
    typedef View<GO*, Device, MemoryUnmanaged> gids_out_type;
    typedef View<ST*, Device, MemoryUnmanaged> vals_out_type;

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
    const size_t num_ent_len = PackTraits<LO,Device>::packValueCount(num_ent_LO);
    input_buffer_type in_buf = subview(imports, slice(offset, offset+num_ent_len));
    (void) PackTraits<LO,Device>::unpackValue(num_ent_LO, in_buf);
    const size_t num_ent = static_cast<size_t>(num_ent_LO);

    // Count the number of bytes expected to unpack
    size_t expected_num_bytes = 0;
    {
      LO lid = 0; // packValueCount wants these
      expected_num_bytes += PackTraits<LO,Device>::packValueCount(lid);

      GO gid = 0; // packValueCount wants these
      expected_num_bytes += num_ent * PackTraits<GO,Device>::packValueCount(gid);

      if (unpack_pids) {
        int pid = 0;
        expected_num_bytes += num_ent * PackTraits<int,Device>::packValueCount(pid);
      }

      ST val; // packValueCount wants these
      expected_num_bytes += num_ent * PackTraits<ST,Device>::packValueCount(val);
    }

    if (expected_num_bytes > num_bytes) {
      // Wrong number of bytes
      dst = Kokkos::make_pair(1, i);
      return;
    }

    if (offset > buf_size || offset + num_bytes > buf_size) {
      // Out of bounds
      dst = Kokkos::make_pair(2, i);
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
    int unpack_err = unpack_crs_matrix_row<ST,LO,GO,Device>(gids_out, pids_out, vals_out,
        imports, offset, num_bytes, num_ent);

    if (unpack_err != 0) {
      // Unpack error
      dst = Kokkos::make_pair(unpack_err, i);
      tokens.release(token);
      return;
    }

    // Column indices come in as global indices, in case the
    // source object's column Map differs from the target object's
    // (this's) column Map, and must be converted local index values
    for (size_t k = 0; k < num_ent; ++k) {
      lids_out(k) = local_col_map.getLocalElement(gids_out(k));
    }

    // Combine the values according to the combine_mode
    const LO* const lids_raw = const_cast<const LO*> (lids_out.data());
    const ST* const vals_raw = const_cast<const ST*> (vals_out.data());
    LO num_modified = 0;
    if (combine_mode == ADD) {
      num_modified += local_matrix.sumIntoValues(import_lid,
          lids_raw, num_ent, vals_raw, false, atomic);
    }
    else if (combine_mode == REPLACE) {
      num_modified += local_matrix.replaceValues(import_lid,
          lids_raw, num_ent, vals_raw, false, atomic);
    }
    else {
      // Invalid combine mode
      dst = Kokkos::make_pair(4, i);
      tokens.release(token);
      return;
    }

    tokens.release(token);

    // Unpack successful.
  }

};

/// \brief Functor to determine the maximum number of entries in any row of the
///   matrix using Kokkos::parallel_reduce
///
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam Device The Kokkos Device type.  See the documentation of Map
///   for requirements.
///
/// The maximum number of entries is required for allocated device scratch space
template<class LO, class Device>
struct MaxNumEntriesFunctor {
  typedef size_t value_type;
  typedef typename PackTraits<LO,Device>::input_buffer_type input_buffer_type;
  typedef typename PackTraits<size_t,Device>::input_array_type offsets_type;
  typedef typename PackTraits<size_t,Device>::input_array_type num_packets_per_lid_type;
  typedef Kokkos::pair<size_t,size_t> slice;
  num_packets_per_lid_type num_packets_per_lid;
  offsets_type offsets;
  input_buffer_type imports;
  MaxNumEntriesFunctor(const num_packets_per_lid_type num_packets_per_lid_in,
                       const offsets_type& offsets_in,
                       const input_buffer_type& imports_in) :
    num_packets_per_lid(num_packets_per_lid_in),
    offsets(offsets_in),
    imports(imports_in)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const LO i, value_type& update) const {
    // Get how many entries to expect in the received data for this row.
    const size_t num_bytes = num_packets_per_lid(i);
    if (num_bytes > 0) {
      LO num_ent_LO = 0;
      const size_t offset = offsets(i);
      const size_t num_ent_len = PackTraits<LO,Device>::packValueCount(num_ent_LO);
      input_buffer_type in_buf = subview(imports, slice(offset, offset+num_ent_len));
      (void) PackTraits<LO,Device>::unpackValue(num_ent_LO, in_buf);
      size_t num_ent = static_cast<size_t>(num_ent_LO);
      if (update < num_ent) update = num_ent;
    }
  }

  KOKKOS_INLINE_FUNCTION void
  join (volatile value_type& dst, const volatile value_type& src) const
  {if (dst < src) dst = src; }
};

/// \brief Functor to determine the maximum number of entries in any row of the
///   matrix using Kokkos::parallel_reduce
///
/// \tparam LO The type of local indices.  See the
///   documentation of Map for requirements.
/// \tparam Device The Kokkos Device type.  See the documentation of Map
///   for requirements.
///
/// This procedure is a high level interface to MaxNumEntriesFunctor
template<class LO, class Device>
size_t compute_maximum_num_entries(
    const typename PackTraits<size_t,Device>::input_array_type num_packets_per_lid,
    const typename PackTraits<size_t,Device>::input_array_type offsets,
    const typename PackTraits<LO,Device>::input_buffer_type imports)
{
  typedef typename Device::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, LO> range_type;
  size_t n = 0;
  MaxNumEntriesFunctor<LO,Device> functor(num_packets_per_lid, offsets, imports);
  Kokkos::parallel_reduce(range_type(0, num_packets_per_lid.dimension_0()), functor, n);
  return n;
}

/// \brief Perform the unpack operation for the matrix
///
/// \tparam LocalMatrix the specialization of the KokkosSparse::CrsMatrix
///   local matrix
/// \tparam LocalMap the type of the local column map
///
/// This is a higher level interface to the UnpackCrsMatrixAndCombineFunctor
template<class LocalMatrix, class LocalMap>
void
do_unpack_and_combine(const LocalMatrix& local_matrix,
                      const LocalMap& local_map,
                      const typename PackTraits<typename LocalMap::local_ordinal_type, typename LocalMap::device_type>::input_buffer_type imports,
                      const typename PackTraits<size_t, typename LocalMap::device_type>::input_array_type num_packets_per_lid,
                      const typename PackTraits<typename LocalMap::local_ordinal_type, typename LocalMap::device_type>::input_array_type import_lids,
                      const Tpetra::CombineMode combine_mode,
                      const bool unpack_pids,
                      const bool atomic)
{
  typedef typename LocalMap::local_ordinal_type LO;
  typedef typename LocalMap::device_type Device;
  typedef typename Device::execution_space execution_space;
  typedef Kokkos::RangePolicy<execution_space, LO> range_type;
  typedef UnpackCrsMatrixAndCombineFunctor<LocalMatrix, LocalMap> unpack_functor_type;

  const char prefix[] = "Tpetra::Details::do_unpack_and_combine: ";

  const size_t num_import_lids = static_cast<size_t>(import_lids.dimension_0());
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
      num_import_lids != static_cast<size_t>(num_packets_per_lid.dimension_0());
    TEUCHOS_TEST_FOR_EXCEPTION(bad_num_import_lids,
        std::invalid_argument,
        prefix << "importLIDs.size() (" << num_import_lids << ") != "
        "numPacketsPerLID.size() (" << num_packets_per_lid.dimension_0() << ").");
  } // end QA error checking

  // Get the offsets
  Kokkos::View<size_t*, Device> offsets("offsets_d", num_import_lids+1);
  computeOffsetsFromCounts(offsets, num_packets_per_lid);

  // Determine the maximum number of entries in any row in the matrix.  The
  // maximum number of entries is needed to allocate unpack buffers on the
  // device.
  size_t max_num_ent = compute_maximum_num_entries<LO,Device>(
      num_packets_per_lid, offsets, imports);

  // Now do the actual unpack!
  unpack_functor_type f(local_matrix, local_map,
      imports, num_packets_per_lid, import_lids, offsets, combine_mode,
      max_num_ent, unpack_pids, atomic);

  typename unpack_functor_type::value_type x;
  Kokkos::parallel_reduce(range_type(0, num_import_lids), f, x);
  auto x_h = x.to_std_pair();
  TEUCHOS_TEST_FOR_EXCEPTION(x_h.first != 0, std::runtime_error,
      prefix << "UnpackCrsMatrixAndCombineFunctor reported error code "
             << x_h.first << " for the first bad row " << x_h.second);

  return;
}

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
/// \param ixports [in/out] Output pack buffer; resized if needed.
///
/// \param numPacketsPerLID [out] Entry k gives the number of bytes
///   packed for row exportLIDs[k] of the local matrix.
///
/// \param ixportLIDs [in] Local indices of the rows to pack.
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
unpackCrsMatrixAndCombine(const CrsMatrix<ST, LO, GO, Node>& sourceMatrix,
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
  typedef typename CrsMatrix<ST,LO,GO,Node>::local_matrix_type local_matrix_type;
  static_assert (std::is_same<device_type, typename local_matrix_type::device_type>::value,
                 "Node::device_type and LocalMatrix::device_type must be the same.");

  // Execution space.
  typedef typename device_type::execution_space execution_space;

  // Convert all Teuchos::Array to Kokkos::View.
  typename execution_space::device_type outputDevice;

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
  do_unpack_and_combine(local_matrix, local_col_map, imports_d,
      num_packets_per_lid_d, import_lids_d, combineMode, false, atomic);

  return;
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_UNPACKCRSMATRIX_HPP
