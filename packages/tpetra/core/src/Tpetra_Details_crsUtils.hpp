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

/// \brief Implementation of numpy.argsort
template <class Ordinal>
Teuchos::Array<size_t>
argsort(Ordinal const * const a, size_t const n)
{
  Teuchos::Array<size_t> ix(n);
  std::iota(ix.begin(), ix.end(), 0);
  std::sort(ix.begin(), ix.end(), [&](size_t i, size_t j) {return a[i] < a[j];});
  return ix;
}

/// \brief Implementation of padCrsArrays
template<class RowPtr, class Indices, class Padding>
void
pad_crs_arrays(
    RowPtr& row_ptr_beg,
    RowPtr& row_ptr_end,
    Indices& indices,
    const Padding& padding)
{

  if (padding.size() == 0 || row_ptr_beg.size() == 0) {
    // Nothing to do
    return;
  }

  // Determine if the indices array is large enough
  auto num_row = row_ptr_beg.size() - 1;
  auto policy = Kokkos::RangePolicy<typename Padding::execution_space>(0, num_row);
  RowPtr entries_this_row("entries_this_row", num_row);
  Kokkos::deep_copy(entries_this_row, 0);
  size_t additional_size_needed = 0;
  Kokkos::parallel_reduce("Determine additional size needed", policy,
    KOKKOS_LAMBDA(const int& i, size_t& ladditional_size_needed) {

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

  if (additional_size_needed == 0) return;

  // The indices array must be resized and the row_ptr arrays shuffled
  auto indices_new = uninitialized_view<Indices>("ind new", indices.size()+additional_size_needed);

  // mfh: Not so fussy about this not being a kernel initially,
  // since we're adding a new feature upon which existing code does not rely,
  // namely Export/Import to a StaticProfile CrsGraph.  However, watch out
  // for fence()ing relating to UVM.
  auto this_row_beg = row_ptr_beg(0);
  auto this_row_end = row_ptr_end(0);
  using range = Kokkos::pair<typename RowPtr::value_type, typename RowPtr::value_type>;
  for (typename RowPtr::size_type i=0; i<num_row-1; i++) {

    // First, copy over indices for this row
    auto used_this_row = this_row_end - this_row_beg;
    auto indices_old_subview =
      subview(indices, range(this_row_beg, this_row_beg+used_this_row));
    auto indices_new_subview =
      subview(indices_new, range(row_ptr_beg(i), row_ptr_beg(i)+used_this_row));

    // just call memcpy; it works fine on device if this becomes a kernel
    using value_type = typename Indices::non_const_value_type;
    memcpy(indices_new_subview.data(), indices_old_subview.data(),
           used_this_row * sizeof(value_type));

    // TODO could this actually have extra entries at the end of each row?
    // If yes, then fill those entries with an "invalid" value.

    // Before modifying the row_ptr arrays, save current beg, end for next iteration
    this_row_beg = row_ptr_beg(i+1);
    this_row_end = row_ptr_end(i+1);

    auto used_next_row = row_ptr_end(i+1) - row_ptr_beg(i+1);
    row_ptr_beg(i+1) = row_ptr_beg(i) + entries_this_row(i);
    row_ptr_end(i+1) = row_ptr_beg(i+1) + used_next_row;
  }
  {
    row_ptr_beg(num_row) = indices_new.size();
    // Copy indices for last row
    auto n = num_row - 1;
    auto used_this_row = row_ptr_end(n) - row_ptr_beg(n);
    auto indices_old_subview =
      subview(indices, range(this_row_beg, this_row_beg+used_this_row));
    auto indices_new_subview =
      subview(indices_new, range(row_ptr_beg(n), row_ptr_beg(n)+used_this_row));
    using value_type = typename Indices::non_const_value_type;
    memcpy(indices_new_subview.data(), indices_old_subview.data(),
           used_this_row * sizeof(value_type));
  }
  indices = indices_new;
}


/// \brief Merge two arrays \c indices1 and \c indices2 representing CRS index arrays.
///   The result is returned in one sorted range with unique elements. Elements
///   are compared using operator<.
///
///   On input, indices1 need not be sorted, but indices2 must be sorted. The
///   intended use case is indices1 contains new elements to merge in to the
///   existing CRS indices contained in indices2 - which should already be
///   sorted.
///
///   This function is modeled after std::merge but takes arrays and lengths
///   instead of iterators. It also does not require that both input arrays be
///   sorted (just the second).
///
/// \param [in] indices1 - the first array of indices
/// \param [in] n1 - the length of indices1
/// \param [in] indices2 - the second array of indices
/// \param [in] n2 - the length of indices2
///
/// \return merged - The sorted and unique elements
///
template <class Ordinal>
Teuchos::Array<Ordinal>
ind_merge(
    Ordinal const * const indices1,
    size_t const n1,
    Ordinal const * const indices2,
    size_t const n2)
{

  // We assume that indices1 is not sorted. It is the index array sent to us by
  // users. indices2, on the other hand, should be sorted. It is the existing
  // CRS index array.
  auto ix = argsort(indices1, n1);

  Teuchos::Array<Ordinal> merged;
  merged.reserve(n1 + n2);
  if (n2 == 0)
  {
    for (size_t i = 0; i < n1; i++)
      merged.push_back(indices1[ix[i]]);
  }
  else
  {
    size_t i1 = 0;
    size_t i2 = 0;
    while (i1 < n1)
    {
      if (i2 == n2)
      {
        for (size_t i = i1; i < n1; i++)
          merged.push_back(indices1[ix[i]]);
        break;
      }
      if (indices2[i2] < indices1[ix[i1]])
        merged.push_back(indices2[i2++]);
      else
        merged.push_back(indices1[ix[i1++]]);
    }
    std::copy(indices2+i2, indices2+n2, std::back_inserter(merged));
  }

  // Return only the unique elements
  auto last = std::unique(merged.begin(), merged.end());
  merged.erase(last, merged.end());
  return merged;
}

/// \brief Copies the elements from the array \c indices1 which are
///   not found in the array \c indices2 to the array \c diff.
///
///   On input, indices1 need not be sorted, but indices2 must be sorted. The
///   intended use case is indices1 contains new elements to merge in to the
///   existing CRS indices contained in indices2 - which should already be
///   sorted.
///
///   This function is modeled after std::set_difference but takes arrays and
///   lengths instead of iterators. It also does not require that both input
///   arrays be sorted (just the second).
///
/// \param [in] indices1 - the first array of indices
/// \param [in] n1 - the length of indices1
/// \param [in] indices2 - the second array of indices
/// \param [in] n2 - the length of indices2
///
/// \return diff - The difference between \c indices1 and \c indices2
///
template <class Ordinal>
Teuchos::Array<Ordinal>
ind_difference(
    Ordinal const * const indices1,
    size_t const n1,
    Ordinal const * const indices2,
    size_t const n2)
{
  Teuchos::Array<Ordinal> diff;

  // We assume that indices1 is not sorted. It is the index array sent to us by
  // users. indices2, on the other hand, should be sorted. It is the existing
  // CRS index array.
  auto ix = argsort(indices1, n1);

  size_t i1 = 0;
  size_t i2 = 0;
  while (i1 < n1)
  {
    if (i2 == n2)
    {
      for (size_t i = i1; i < n1; i++)
        diff.push_back(indices1[ix[i]]);
      break;
    }
    auto ind1 = indices1[ix[i1]];
    auto ind2 = indices2[i2];
    if (ind1 < ind2)
    {
      diff.push_back(ind1);
      i1++;
    }
    else
    {
      if (! (ind2 < ind1))
        i1++;
      i2++;
    }
  }

  // Return only unique elements
  auto last = std::unique(diff.begin(), diff.end());
  diff.erase(last, diff.end());
  return diff;
}

/// \brief Implementation of insertCrsIndices
template<class T>
int
insert_crs_indices(
    T& indices,
    size_t const num_assigned,
    typename T::value_type const in_indices[],
    size_t const num_in)
{
  if (num_in == 0)
    return 0;

  // Rather than loop through the two arrays of indices and looking for their
  // difference and then checking if their is enough capacity in indices to
  // accommodate the new elements (not already in indices, excluding
  // duplicates), we go ahead and merge the two arrays in to one. After the
  // arrays are merged, we can check if the length of the merged array exceeds
  // the capacity of the indices array. If it does, return -1 to indicate the
  // error. If it does not, however, we can just copy the merged elements back
  // in to \c indices.
  auto merged = ind_merge(in_indices, num_in, indices.data(), num_assigned);

  // Check if we have the capacity for the incoming indices
  if (static_cast<size_t>(merged.size()) > indices.size())
    return -1;
  std::copy(merged.begin(), merged.end(), indices.data());
  return merged.size() - num_assigned;
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
    RowPtr& rowPtrBeg,
    RowPtr& rowPtrEnd,
    Indices& indices,
    const Padding& padding)
{
  using impl::pad_crs_arrays;
  pad_crs_arrays<RowPtr, Indices, Padding>(rowPtrBeg, rowPtrEnd, indices, padding);
}

/// \brief Insert new indices in \c inIndices in to \c indices
///
///   On input, inIndices need not be sorted, but indices must be sorted. The
///   intended use case is inIndices contains new elements to insert in to the
///   existing CRS indices contained in indices - which should already be
///   sorted.
///
///   Only elements in inIndices that are not already in indices are inserted
///   (duplicates are skipped).
///
/// \param indices [in/out] CRS indices. The CRS indices must be sorted on input
/// \param numAssigned [in].The number of indices that have already been assigned.
/// \param inIndices [in] The indices to insert. Only those entries not already
///    present in \c indices will be inserted.
/// \param numIn [in] The length of in_indices
///
/// \return numInserted [out] The number of indices inserted.
///    If numInserted == -1, there was not enough capacity for all of the incoming
///    indices and none were inserted.
template <class T>
int
insertCrsIndices(
    T& indices,
    size_t const numAssigned,
    typename T::value_type const inIndices[],
    size_t const numIn)
{
  return impl::insert_crs_indices(indices, numAssigned, inIndices, numIn);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_CRSUTILS_HPP
