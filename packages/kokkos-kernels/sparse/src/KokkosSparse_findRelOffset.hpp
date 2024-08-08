//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_SPARSE_FINDRELOFFSET_HPP
#define KOKKOS_SPARSE_FINDRELOFFSET_HPP

/// \file KokkosSparse_findRelOffset.hpp
/// \brief Find the relative offset of a column index in a sparse
///   graph's or sparse matrix's row.

#include "Kokkos_Macros.hpp"  // KOKKOS_FUNCTION
#include <type_traits>

namespace KokkosSparse {
/// \brief Search <tt>indsToSearch[0 .. numEnt-1]</tt> for
///   \c indToFind, using equality comparison.
///
/// \return If found, return index of \c indToFind in \c indsToSearch;
///   else, return \c numEnt (by analogy with C++ Standard Library
///   functions like std::find, that return "the end of the sequence"
///   in this case).
///
/// \tparam OffsetType Integer type that can be used to represent any
///   valid index in \c indsToSearch, up to and including \c numEnt.
/// \tparam IndexViewType 1-D array of equality-comparable entries
///   (generally intended to be column indices).  This may a 1-D
///   Kokkos::View, a raw 1-D array, or any type that implements
///   operator[](OffsetType).
///
/// \param indsToSearch [in] Array of indices to search.  For a
///   sparse graph or matrix, this is the array of all the column
///   indices for some row of the graph / matrix.
/// \param numEnt [in] Number of entries in \c indsToSearch to
///   search.  This is a separate argument, first so that this
///   function works with raw arrays as well as Kokkos::View, and
///   second so that users don't have to incur the overhead of
///   calling Kokkos::subview to limit the length of a View.  The
///   latter may be particularly helpful for the case of the
///   begin/end-pointer variant of CSR graph/matrix storage.
/// \param indToFind [in] (Local) column index for which to find the
///   offset.  This has the same type as that of each entry in
///   \c indsToSearch.
/// \param hint [in] Hint for where to find \c indToFind in the array.
///   If <tt>indsToSearch[hint] == indToFind</tt>, then the hint is
///   correct.  The hint is ignored if it is out of range (that is,
///   greater than or equal to the number of entries in the given
///   row).
/// \param isSorted [in] Whether the input array of indices to search
///   is sorted in increasing order.
///
/// The hint optimizes for the case of calling this method several
/// times with the same sparse graph / matrix row, when several
/// index inputs occur in consecutive sequence.  This may occur (for
/// example) when there are multiple degrees of freedom per mesh
/// point, and users are handling the assignment of degrees of
/// freedom to global indices manually (rather than letting some
/// other class take care of it).  In that case, users might choose
/// to assign the degrees of freedom for a mesh point to consecutive
/// global indices.  Epetra implements the hint for this reason.
///
/// The hint only costs two comparisons (one to check range, and the
/// other to see if the hint was correct), and it can save searching
/// for the indices (which may take a lot more than two
/// comparisons).
///
/// \note To implementers: We put <tt>indsToSearch</tt> before
///   <tt>indToFind</tt> so that we can derive the type of
///   <tt>indToFind</tt> directly from that of each entry of
///   <tt>indsToSearch</tt>, without needing <tt>IndexViewType</tt>
///   to be a Kokkos::View.  Thankfully, arguments to a C++ function
///   behave more like LET* than LET (in ANSI Common Lisp terms).
template <class OffsetType, class IndexViewType>
KOKKOS_FUNCTION OffsetType findRelOffset(const IndexViewType& indsToSearch, const OffsetType numEnt,
                                         /* typename IndexViewType::const_value_type */
                                         const typename std::decay<decltype(indsToSearch[0])>::type indToFind,
                                         const OffsetType hint, const bool isSorted) {
  // IndexViewType doesn't have to be a Kokkos::View; it just has to
  // implement operator[] like a 1-D array.
  //
  // static_assert (Kokkos::is_view<IndexViewType>::value,
  //                 "IndexViewType must be a Kokkos::View");
  // static_assert (static_cast<int> (IndexViewType::rank) == 1,
  //                 "IndexViewType must be a rank-1 Kokkos::View");
  static_assert(std::is_integral<OffsetType>::value, "OffsetType must be an integer.");

  if (hint < numEnt && indsToSearch[hint] == indToFind) {
    return hint;  // hint was correct
  }

  // Even if the array is sorted, use linear search if the number of
  // entries is small ("small" is a tuning parameter; feel free to
  // tune for your architecture).  'constexpr' promises the compiler
  // that it can bake this constant as a literal into the code.
  constexpr OffsetType linearSearchThreshold = 16;

  if (!isSorted || numEnt < linearSearchThreshold) {
    for (OffsetType k = 0; k < numEnt; ++k) {
      if (indsToSearch[k] == indToFind) {
        return k;
      }
    }
  } else {  // use binary search
    OffsetType start = 0;
    OffsetType end   = numEnt;
    // Compare epetra/src/Epetra_Util.cpp, Epetra_Util_binary_search.
    // Unlike that function, I don't use end = numEnt-1, because I
    // want this code to work also for unsigned OffsetType (signed is
    // preferred, though).  Thus, in my code, end is always "one past
    // the last valid index."
    while (end > start) {
      // Invariants: 0 <= start < end, thus start + end > 0.
      const OffsetType mid = (start + end - 1) / 2;
      // Invariants: 0 <= start <= mid < end.
      if (indsToSearch[mid] < indToFind) {
        // Invariant: start < mid+1 (thus, recursion terminates),
        // and for all k <= mid, indsToSearch[k] < indToFind.
        start = mid + 1;  // Invariant: 0 < mid < start <= end.
      } else {            // indsToSearch[mid] >= indToFind
        // Invariant: mid < end (thus, recursion terminates),
        // and for all k <= mid, indsToSearch[k] >= indToFind.
        end = mid;  // Invariant: 0 <= start <= mid <= end.
      }
    }
    // Invariant: 0 <= start == end.

    // Don't check if we've already passed the end.
    if (start < numEnt && indsToSearch[start] == indToFind) {
      return start;
    }
  }

  return numEnt;  // "end of sequence"
}

}  // namespace KokkosSparse

#endif  // KOKKOS_SPARSE_FINDRELOFFSET_HPP
