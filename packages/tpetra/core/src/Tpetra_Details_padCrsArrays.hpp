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

#ifndef TPETRA_DETAILS_PADCRSARRAYS_HPP
#define TPETRA_DETAILS_PADCRSARRAYS_HPP

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"
#include "Kokkos_UnorderedMap.hpp"

/// \file Tpetra_Details_padCrsArrays.hpp
/// \brief Functions that pad CRS arrays by resizing indices and inserting empty
///   space and adjusting the row pointers as needed.
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {
namespace Details {

namespace pad_crs_impl {

template<class ViewType>
ViewType uninitialized_view(const std::string& name, const size_t& size) {
  return ViewType (Kokkos::view_alloc(name, Kokkos::WithoutInitializing), size);
}

// Determine if row_ptr and indices arrays need to be resized to accommodate
// new entries
template<class RowPtr, class Indices, class Padding>
void
pad_crs_arrays(RowPtr& row_ptr_beg, RowPtr& row_ptr_end, Indices& indices, const Padding& padding)
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
} // namespace pad_crs_impl

// Determine if row_ptr and indices arrays need to be resized to accommodate
// new entries
template<class RowPtr, class Indices, class Padding>
void
padCrsArrays(RowPtr& row_ptr_beg, RowPtr& row_ptr_end, Indices& indices,
             const Padding& padding)
{
  using pad_crs_impl::pad_crs_arrays;
  pad_crs_arrays<RowPtr, Indices, Padding>(row_ptr_beg, row_ptr_end, indices, padding);
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_PADCRSARRAYS_HPP
