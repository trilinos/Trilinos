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

#ifndef TPETRA_DETAILS_RESIZEROWPTR_HPP
#define TPETRA_DETAILS_RESIZEROWPTR_HPP

#include "TpetraCore_config.h"
#include "Kokkos_Core.hpp"

/// \file Tpetra_Details_resizeRowPtr.hpp
/// \brief Functions that resizes CSR row pointers and indices
/// \warning This file, and its contents, are implementation details
///   of Tpetra.  The file itself or its contents may disappear or
///   change at any time.

namespace Tpetra {
namespace Details {


// Determine if row_ptrs and indices arrays need to be resized to accommodate
// new entries
template<class RowView, class EntriesView, class NumPacketsView>
KOKKOS_FUNCTION void
resizeRowPtrsAndIndices(RowView& row_ptrs_beg, RowView& row_ptrs_end, EntriesView& indices,
                        const NumPacketsView& num_packets_per_lid, const bool unpack_pids)
{

  // Determine if the indices array is large enough
  RowView entries_this_row("entries_this_row", num_packets_per_lid.size());
  Kokkos::deep_copy(entries_this_row, 0);
  size_t additional_size_needed = 0;
  Kokkos::parallel_reduce("Determine additional size needed", num_packets_per_lid.size(),
    KOKKOS_LAMBDA(const int& i, size_t& ladditional_size_needed) {
      size_t num_packets_this_lid = num_packets_per_lid(i);
      size_t num_ent = (unpack_pids) ? num_packets_this_lid/2
                                     : num_packets_this_lid;
      auto allocated_this_row = row_ptrs_beg(i+1) - row_ptrs_beg(i);
      auto used_this_row = row_ptrs_end(i) - row_ptrs_beg(i);
      auto free_this_row = allocated_this_row - used_this_row;
      auto n = (num_ent > free_this_row) ? num_ent - free_this_row : 0;
      entries_this_row(i) = allocated_this_row + n;
      ladditional_size_needed += n;
    }, additional_size_needed);

  if (additional_size_needed == 0) return;

  // The row_ptrs and indices arrays must be resized
  EntriesView indices_copy("", indices.size());
  Kokkos::deep_copy(indices_copy, indices);

  auto this_row_beg = row_ptrs_beg(0);
  auto next_row_beg = row_ptrs_beg(1);
  auto this_row_end = row_ptrs_end(0);
  Kokkos::resize(indices, indices.size()+additional_size_needed);
  Kokkos::deep_copy(indices, Teuchos::OrdinalTraits<typename EntriesView::value_type>::invalid());
  for (size_t i=0; i<num_packets_per_lid.size()-1; i++) {

    auto allocated_this_row = next_row_beg - this_row_beg;
    auto used_this_row = this_row_end - this_row_beg;
    auto additional_entries_needed = entries_this_row(i) - allocated_this_row;

    // First, copy over indices for this row
    const auto slice1 = std::pair<size_t,size_t>(this_row_beg, this_row_beg+used_this_row);
    auto indices_copy_subview = subview(EntriesView(indices_copy), slice1);

    // subview of new indices array
    const auto slice2 = std::pair<size_t,size_t>(row_ptrs_beg(i), row_ptrs_beg(i)+used_this_row);
    auto indices_subview = subview(EntriesView(indices), slice2);

    Kokkos::deep_copy(indices_subview, indices_copy_subview);

    // Before modifying the row_ptrs arrays, save current beg, end for next iteration
    this_row_beg = row_ptrs_beg(i+1);
    next_row_beg = row_ptrs_beg(i+2);
    this_row_end = row_ptrs_end(i+1);

    if (additional_entries_needed <= 0) {
      // Nothing more to do
      continue;
    }

    // Shift the row_ptrs array to accommodate the extra space needed
    auto used_next_row = row_ptrs_end(i+1) - row_ptrs_beg(i+1);
    row_ptrs_beg(i+1) = row_ptrs_beg(i) + entries_this_row(i);
    row_ptrs_end(i+1) = row_ptrs_beg(i+1) + used_next_row;
  }
  // Copy last row
  {
    auto used_this_row = this_row_end - this_row_beg;
    auto n = row_ptrs_beg.extent(0);
    const auto slice1 = std::pair<size_t,size_t>(this_row_beg, this_row_beg+used_this_row);
    auto indices_copy_subview = subview(EntriesView(indices_copy), slice1);
    const auto slice2 = std::pair<size_t,size_t>(row_ptrs_beg(n-2), row_ptrs_beg(n-2)+used_this_row);
    auto indices_subview = subview(EntriesView(indices), slice2);
    Kokkos::deep_copy(indices_subview, indices_copy_subview);
    row_ptrs_beg(n-1) = indices.extent(0);
  }

}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_RESIZEROWPTR_HPP
