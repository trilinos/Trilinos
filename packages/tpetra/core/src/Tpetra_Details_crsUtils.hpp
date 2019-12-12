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

#ifndef TPETRA_DETAILS_CRSUTILS_HPP
#define TPETRA_DETAILS_CRSUTILS_HPP

#include <numeric>
#include <type_traits>

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_UnorderedMap.hpp"

/// \file Tpetra_Details_crsUtils.hpp
/// \brief Functions for manipulating CRS arrays
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {
namespace Details {

namespace impl {

template<class ViewType>
ViewType
uninitialized_view(const std::string& name, const size_t& size)
{
  return ViewType (Kokkos::view_alloc(name, Kokkos::WithoutInitializing), size);
}

/// \brief Implementation of padCrsArrays
template<class RowPtr, class Indices, class Values, class Padding>
void
pad_crs_arrays(
               const RowPtr& row_ptr_beg,
               const RowPtr& row_ptr_end,
    Indices& indices,
    Values& values,
    const Padding& padding)
{

  if (padding.size() == 0 || row_ptr_beg.size() == 0) {
    // Nothing to do
    return;
  }

  auto pad_values = values.extent(0) == indices.extent(0);

  // Determine if the indices array is large enough
  auto num_row = row_ptr_beg.size() - 1;
  auto policy = Kokkos::RangePolicy<typename Padding::execution_space>(0, num_row);
  RowPtr entries_this_row("entries_this_row", num_row);
  Kokkos::deep_copy(entries_this_row, 0);
  size_t additional_size_needed = 0;
  Kokkos::parallel_reduce("Determine additional size needed", policy,
    KOKKOS_LAMBDA(const int i, size_t& ladditional_size_needed) {

      auto allocated_this_row = row_ptr_beg(i+1) - row_ptr_beg(i);
      auto used_this_row = row_ptr_end(i) - row_ptr_beg(i);
      auto free_this_row = allocated_this_row - used_this_row;
      entries_this_row(i) = allocated_this_row;

      auto k = padding.find(static_cast<typename Padding::key_type>(i));
      if (padding.valid_at(k)) {
        // Additional padding was requested for this LID
        auto num_ent = padding.value_at(k);
        auto n = (num_ent > free_this_row) ? num_ent - free_this_row : 0;
        entries_this_row(i) += n;
        ladditional_size_needed += n;
      }
    }, additional_size_needed);

  if (additional_size_needed == 0)
    return;

  using ptrs_value_type = typename RowPtr::non_const_value_type;
  using inds_value_type = typename Indices::non_const_value_type;
  using vals_value_type = typename Values::non_const_value_type;

  // The indices array must be resized and the row_ptr arrays shuffled
  auto indices_new = uninitialized_view<Indices>("ind new", indices.size()+additional_size_needed);
  auto values_new = uninitialized_view<Values>("val new", pad_values ? values.size()+additional_size_needed : 0);
  Kokkos::deep_copy(values_new, vals_value_type(0.0));

  // mfh: Not so fussy about this not being a kernel initially,
  // since we're adding a new feature upon which existing code does not rely,
  // namely Export/Import to a StaticProfile CrsGraph.  However, watch out
  // for fence()ing relating to UVM.
  auto this_row_beg = row_ptr_beg(0);
  auto this_row_end = row_ptr_end(0);
  using range = Kokkos::pair<ptrs_value_type, ptrs_value_type>;
  for (typename RowPtr::size_type i=0; i<num_row-1; i++) {

    auto used_this_row = this_row_end - this_row_beg;

    // Copy over indices
    {
      auto indices_old_subview = subview(indices, range(this_row_beg, this_row_beg+used_this_row));
      auto indices_new_subview = subview(indices_new, range(row_ptr_beg(i), row_ptr_beg(i)+used_this_row));
      // just call memcpy; it works fine on device if this becomes a kernel
      memcpy(indices_new_subview.data(), indices_old_subview.data(), used_this_row * sizeof(inds_value_type));
    }

    // And then the values
    if (pad_values) {
      auto values_old_subview = subview(values, range(this_row_beg, this_row_beg+used_this_row));
      auto values_new_subview = subview(values_new, range(row_ptr_beg(i), row_ptr_beg(i)+used_this_row));
      memcpy(values_new_subview.data(), values_old_subview.data(), used_this_row * sizeof(vals_value_type));
    }

    // Before modifying the row_ptr arrays, save current beg, end for next iteration
    this_row_beg = row_ptr_beg(i+1);
    this_row_end = row_ptr_end(i+1);

    auto used_next_row = row_ptr_end(i+1) - row_ptr_beg(i+1);
    row_ptr_beg(i+1) = row_ptr_beg(i) + entries_this_row(i);
    row_ptr_end(i+1) = row_ptr_beg(i+1) + used_next_row;
  }

  {
    // Copy indices/values for last row
    row_ptr_beg(num_row) = indices_new.size();
    auto n = num_row - 1;
    auto used_this_row = row_ptr_end(n) - row_ptr_beg(n);

    {
      auto indices_old_subview = subview(indices, range(this_row_beg, this_row_beg+used_this_row));
      auto indices_new_subview = subview(indices_new, range(row_ptr_beg(n), row_ptr_beg(n)+used_this_row));
      memcpy(indices_new_subview.data(), indices_old_subview.data(), used_this_row * sizeof(inds_value_type));
    }

    if (pad_values) {
      auto values_old_subview = subview(values, range(this_row_beg, this_row_beg+used_this_row));
      auto values_new_subview = subview(values_new, range(row_ptr_beg(n), row_ptr_beg(n)+used_this_row));
      memcpy(values_new_subview.data(), values_old_subview.data(), used_this_row * sizeof(vals_value_type));
    }
  }

  indices = indices_new;
  values = values_new;
}

/// \brief Implementation of insertCrsIndices
template <class Pointers, class InOutIndices, class InIndices, class IndexMap>
size_t
insert_crs_indices(
    typename Pointers::value_type const row,
    Pointers const& row_ptrs,
    InOutIndices& cur_indices,
    size_t& num_assigned,
    InIndices const& new_indices,
    IndexMap&& map,
    std::function<void(size_t const, size_t const, size_t const)> cb)
{
  if (new_indices.size() == 0) {
    return 0;
  }

  if (cur_indices.size() == 0) {
    // No room to insert new indices
    return Teuchos::OrdinalTraits<size_t>::invalid();
  }

  using offset_type = typename std::decay<decltype (row_ptrs[0])>::type;
  using ordinal_type = typename std::decay<decltype (cur_indices[0])>::type;

  const offset_type start = row_ptrs[row];
  offset_type end = start + static_cast<offset_type> (num_assigned);
  const size_t num_avail = (row_ptrs[row + 1] < end) ? size_t (0) :
    row_ptrs[row + 1] - end;
  const size_t num_new_indices = static_cast<size_t> (new_indices.size ());
  size_t num_inserted = 0;

  for (size_t k = 0; k < num_new_indices; ++k) {
    const ordinal_type idx = std::forward<IndexMap>(map)(new_indices[k]);
    offset_type row_offset = start;
    for (; row_offset < end; ++row_offset) {
      if (idx == cur_indices[row_offset]) {
        break;
      }
    }

    if (row_offset == end) {
      if (num_inserted >= num_avail) { // not enough room
        return Teuchos::OrdinalTraits<size_t>::invalid();
      }
      // This index is not yet in indices
      cur_indices[end++] = idx;
      num_inserted++;
    }
    if (cb) {
      cb(k, start, row_offset - start);
    }
  }
  num_assigned += num_inserted;
  return num_inserted;
}

/// \brief Implementation of findCrsIndices
template <class Pointers, class Indices1, class Indices2, class IndexMap, class Callback>
size_t
find_crs_indices(
    typename Pointers::value_type const row,
    Pointers const& row_ptrs,
    const size_t curNumEntries,
    Indices1 const& cur_indices,
    Indices2 const& new_indices,
    IndexMap&& map,
    Callback&& cb)
{
  if (new_indices.size() == 0)
    return 0;

  using ordinal = typename Indices1::value_type;
  auto invalid_ordinal = Teuchos::OrdinalTraits<ordinal>::invalid();

  const size_t start = static_cast<size_t> (row_ptrs[row]);
  const size_t end = start + curNumEntries;
  size_t num_found = 0;
  for (size_t k = 0; k < new_indices.size(); k++)
  {
    auto row_offset = start;
    auto idx = std::forward<IndexMap>(map)(new_indices[k]);
    if (idx == invalid_ordinal)
      continue;
    for (; row_offset < end; row_offset++)
    {
      if (idx == cur_indices[row_offset])
      {
        std::forward<Callback>(cb)(k, start, row_offset - start);
        num_found++;
      }
    }
  }
  return num_found;
}

} // namespace impl


/// \brief Determine if the row pointers and indices arrays need to be resized
///   to accommodate new entries. If they do need to be resized, resize the
///   indices arrays and shift the existing contents to accommodate new entries.
///   Modify values in the row pointers array to point to the newly shifted
///   locations in the indices arrays.
///
///   This routine is called to resize/shift the CRS arrays before attempting to
///   insert new values if the number of new values exceeds the amount of free
///   space in the CRS arrays.
///
/// \param [in/out] rowPtrBeg - rowPtrBeg[i] points to the first
///        column index (in the indices array) of row i.
/// \param [in/out] rowPtrEnd - rowPtrEnd[i] points to the last
///        column index (in the indices array) of row i.
/// \param [in/out] indices - array containing columns indices of nonzeros in
///        CRS representation.
/// \param [in] padding - Kookos::UnorderedMap. padding[i] is the amount of free
///        space required for row i. If the distance between row_ptr_end[i] and
///        row_ptr_beg[i] does not accommodate padding[i], we resize and shift
///        indices to accommodate.
///
template<class RowPtr, class Indices, class Padding>
void
padCrsArrays(
             const RowPtr& rowPtrBeg,
             const RowPtr& rowPtrEnd,
    Indices& indices,
    const Padding& padding)
{
  using impl::pad_crs_arrays;
  // send empty values array
  Indices values;
  pad_crs_arrays<RowPtr, Indices, Indices, Padding>(rowPtrBeg, rowPtrEnd, indices, values, padding);
}

template<class RowPtr, class Indices, class Values, class Padding>
void
padCrsArrays(
             const RowPtr& rowPtrBeg,
             const RowPtr& rowPtrEnd,
    Indices& indices,
    Values& values,
    const Padding& padding)
{
  using impl::pad_crs_arrays;
  pad_crs_arrays<RowPtr, Indices, Values, Padding>(rowPtrBeg, rowPtrEnd, indices, values, padding);
}

/// \brief Insert new indices in to current list of indices
///
/// \param row [in] The row in which to insert
/// \param rowPtrs [in] "Pointers" to beginning of each row
/// \param curIndices [in/out] The current indices
/// \param numAssigned [in/out] The number of currently assigned indices in row \c row
/// \param newIndices [in] The indices to insert
/// \param map [in] An optional function mapping newIndices[k] to its actual index
/// \param cb [in] An optional callback function called on every insertion at the local
///     index and the offset in to the inserted location
/// \return numInserted The number of indices inserted. If there is not
///     capacity in curIndices for newIndices, return -1;
///
/// \bf Notes
/// \c curIndices is the current list of CRS indices. it is not assumed to be sorted, but
/// entries are unique. For each \c newIndices[k], we look to see if the index exists in
/// \c cur_indices. If it does, we do not insert it (no repeats). If it does not exist, we
/// first check to make sure there is capacity in \c curIndices and if there is we insert
/// it at the end.
///
/// The actual value of \c newIndices[k] that is inserted is the value returned from \c
/// map(newIndices[k]). If an identity map is provided, \c newIndices[k] is directly
/// inserted. However, any other map can be provided. For instance, for a locally indexed
/// graph on which \c insertGlobalIndices is called, the \c curIndices array can be a
/// view of the graph's local indices, the \c newIndices array are the new *global*
/// indices, and \c map is the graph's column map to convert global indices to local.
/// If this function is called through the overload below without the \c map
/// argument, the identity map is provided.
///
/// The optional function \c cb is called on every valid index. \c cb is sent the
/// current loop index \c k, \c rowPtrs[k] (the start of the row), and the relative
/// offset from \c start in to the \c curIndices array for \c newIndices[k]. This
/// function could, for example, be used by \c CrsMatrix to fill the values array during
/// \c sumInto*Values or \c replace*Values; Eg, \c CrsMatrix::sumIntoLocalValues
/// might have the following:
///
/// <code>
/// CrsMatrix::sumIntoLocalValues(LO row, array<LO> cols, array<S> vals)
/// {
///   this->graph_->insertLocalValues(row, cols,
///       [&](size_t const k, size_t const start, size_t const offset){
///           this->values_[start+offset] += vals[k]; });
/// }
/// </code>
///
template <class Pointers, class InOutIndices, class InIndices>
size_t
insertCrsIndices(
    typename Pointers::value_type const row,
    Pointers const& rowPtrs,
    InOutIndices& curIndices,
    size_t& numAssigned,
    InIndices const& newIndices,
    std::function<void(const size_t, const size_t, const size_t)> cb =
        std::function<void(const size_t, const size_t, const size_t)>())
{
  static_assert(std::is_same<typename std::remove_const<typename InOutIndices::value_type>::type,
                             typename std::remove_const<typename InIndices::value_type>::type>::value,
    "Expected views to have same value type");

  // Provide a unit map for the more general insert_indices
  using ordinal = typename InOutIndices::value_type;
  auto numInserted = impl::insert_crs_indices(row, rowPtrs, curIndices,
    numAssigned, newIndices, [](ordinal const idx) { return idx; }, cb);
  return numInserted;
}

template <class Pointers, class InOutIndices, class InIndices>
size_t
insertCrsIndices(
    typename Pointers::value_type const row,
    Pointers const& rowPtrs,
    InOutIndices& curIndices,
    size_t& numAssigned,
    InIndices const& newIndices,
    std::function<typename InOutIndices::value_type(const typename InIndices::value_type)> map,
    std::function<void(const size_t, const size_t, const size_t)> cb =
      std::function<void(const size_t, const size_t, const size_t)>())
{
  auto numInserted = impl::insert_crs_indices(row, rowPtrs, curIndices,
    numAssigned, newIndices, map, cb);
  return numInserted;
}


/// \brief Finds offsets in to current list of indices
///
/// \param row [in] The row in which to insert
/// \param rowPtrs [in] "Pointers" to beginning of each row
/// \param curIndices [in] The current indices
/// \param numAssigned [in] The number of currently assigned indices in row \c row
/// \param newIndices [in] The indices to insert
/// \param cb [in] An optional function called on every insertion at the local
///     index and the offset in to the inserted location
/// \return numFound The number of indices found.
///
/// \bf Notes
/// \c curIndices is the current list of CRS indices. it is not assumed to be sorted, but
/// entries are unique. For each \c newIndices[k], we look to see if the index exists in
/// \c curIndices. If it does, we do not insert it (no repeats). If it does not exist, we
/// first check to make sure there is capacity in \c curIndices and if there is we insert
/// it at the end.
///
/// The actual value of \c newIndices[k] that is inserted is the value returned from \c
/// map(newIndices[k]). If an identity map is provided, \c newIndices[k] is directly
/// inserted. However, any other map can be provided. For instance, for a locally indexed
/// graph on which \c insertGlobalIndices is called, the \c curIndices array can be a
/// view of the graph's local indices, the \c newIndices array are the new *global*
/// indices, and \c map is the graph's column map to convert global indices to local.
/// If this function is called through the overload below without the \c map
/// argument, the identity map is provided.
///
/// The function \c cb is called on every valid index.
///
template <class Pointers, class Indices1, class Indices2, class Callback>
size_t
findCrsIndices(
    typename Pointers::value_type const row,
    Pointers const& rowPtrs,
    const size_t curNumEntries,
    Indices1 const& curIndices,
    Indices2 const& newIndices,
    Callback&& cb)
{
  static_assert(std::is_same<typename std::remove_const<typename Indices1::value_type>::type,
                             typename std::remove_const<typename Indices2::value_type>::type>::value,
    "Expected views to have same value type");
  // Provide a unit map for the more general find_crs_indices
  using ordinal = typename Indices2::value_type;
  auto numFound = impl::find_crs_indices(row, rowPtrs, curNumEntries, curIndices, newIndices,
    [=](ordinal ind){ return ind; }, cb);
  return numFound;
}

template <class Pointers, class Indices1, class Indices2, class IndexMap, class Callback>
size_t
findCrsIndices(
    typename Pointers::value_type const row,
    Pointers const& rowPtrs,
    const size_t curNumEntries,
    Indices1 const& curIndices,
    Indices2 const& newIndices,
    IndexMap&& map,
    Callback&& cb)
{
  return impl::find_crs_indices(row, rowPtrs, curNumEntries, curIndices, newIndices, map, cb);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CRSUTILS_HPP
