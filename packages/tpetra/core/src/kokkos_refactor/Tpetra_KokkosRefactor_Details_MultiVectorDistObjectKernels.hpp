/*
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
*/

// mfh 13/14 Sep 2013 The "should use as<size_t>" comments are both
// incorrect (as() is not a device function) and usually irrelevant
// (it would only matter if LocalOrdinal were bigger than size_t on a
// particular platform, which is unlikely).

#ifndef TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_HPP
#define TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_HPP

#include "Kokkos_Core.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "impl/Kokkos_Atomic_Generic.hpp"
#include <sstream>
#include <stdexcept>

namespace Tpetra {
namespace KokkosRefactor {
namespace Details {

/// \brief Implementation details <i>of</i> implementation details
///
/// \warning DO NOT USE ANYTHING HERE.  WE MAKE NO PROMISES OF
///   BACKWARDS COMPATIBILITY FOR ANYTHING IN THIS NAMESPACE.  IT MAY
///   DISAPPEAR OR CHANGE AT ANY TIME.  THIS IS <i>NOT</i> FOR USERS.
namespace Impl {

/// \brief Is x out of bounds?  That is, is x less than zero, or
///   greater than or equal to the given exclusive upper bound?
///
/// We go through all this trouble to avoid the compiler warnings that
/// may result from asking whether x is less than zero, when x is an
/// unsigned integer.
template<class IntegerType,
         const bool isSigned = std::numeric_limits<IntegerType>::is_signed>
struct OutOfBounds {
  static KOKKOS_INLINE_FUNCTION bool
  test (const IntegerType x,
        const IntegerType exclusiveUpperBound);
};

// Partial specialization for the case where IntegerType IS signed.
template<class IntegerType>
struct OutOfBounds<IntegerType, true> {
  static KOKKOS_INLINE_FUNCTION bool
  test (const IntegerType x,
        const IntegerType exclusiveUpperBound)
  {
    return x < static_cast<IntegerType> (0) || x >= exclusiveUpperBound;
  }
};

// Partial specialization for the case where IntegerType is NOT signed.
template<class IntegerType>
struct OutOfBounds<IntegerType, false> {
  static KOKKOS_INLINE_FUNCTION bool
  test (const IntegerType x,
        const IntegerType exclusiveUpperBound)
  {
    return x >= exclusiveUpperBound;
  }
};

/// \brief Is x out of bounds?  That is, is x less than zero, or
///   greater than or equal to the given exclusive upper bound?
template<class IntegerType>
KOKKOS_INLINE_FUNCTION bool
outOfBounds (const IntegerType x, const IntegerType exclusiveUpperBound)
{
  return OutOfBounds<IntegerType>::test (x, exclusiveUpperBound);
}

} // namespace Impl

  // Functors for implementing packAndPrepare and unpackAndCombine
  // through parallel_for

  template <typename DstView, typename SrcView, typename IdxView>
  struct PackArraySingleColumn {
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    DstView dst;
    SrcView src;
    IdxView idx;
    size_t col;

    PackArraySingleColumn (const DstView& dst_,
                           const SrcView& src_,
                           const IdxView& idx_,
                           const size_t col_) :
      dst(dst_), src(src_), idx(idx_), col(col_) {}

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k) const {
      dst(k) = src(idx(k), col);
    }

    static void
    pack (const DstView& dst,
          const SrcView& src,
          const IdxView& idx,
          const size_t col)
    {
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      Kokkos::parallel_for
        ("Tpetra::MultiVector pack one col",
         range_type (0, idx.size ()),
         PackArraySingleColumn (dst, src, idx, col));
    }
  };

  template <typename DstView,
            typename SrcView,
            typename IdxView,
            typename SizeType = typename DstView::execution_space::size_type>
  class PackArraySingleColumnWithBoundsCheck {
  private:
    static_assert (Kokkos::Impl::is_view<DstView>::value,
                   "DstView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<SrcView>::value,
                   "SrcView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<IdxView>::value,
                   "IdxView must be a Kokkos::View.");
    static_assert (static_cast<int> (DstView::rank) == 1,
                   "DstView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 2,
                   "SrcView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");
    static_assert (std::is_integral<SizeType>::value,
                   "SizeType must be a built-in integer type.");
  public:
    typedef SizeType size_type;
    using value_type = size_t;

  private:
    DstView dst;
    SrcView src;
    IdxView idx;
    size_type col;

  public:
    PackArraySingleColumnWithBoundsCheck (const DstView& dst_,
                                          const SrcView& src_,
                                          const IdxView& idx_,
                                          const size_type col_) :
      dst (dst_), src (src_), idx (idx_), col (col_) {}

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k, value_type& lclErrCount) const {
      using index_type = typename IdxView::non_const_value_type;

      const index_type lclRow = idx(k);
      if (lclRow < static_cast<index_type> (0) ||
          lclRow >= static_cast<index_type> (src.extent (0))) {
        ++lclErrCount;
      }
      else {
        dst(k) = src(lclRow, col);
      }
    }

    KOKKOS_INLINE_FUNCTION
    void init (value_type& initialErrorCount) const {
      initialErrorCount = 0;
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& dstErrorCount,
          const volatile value_type& srcErrorCount) const
    {
      dstErrorCount += srcErrorCount;
    }

    static void
    pack (const DstView& dst,
          const SrcView& src,
          const IdxView& idx,
          const size_type col)
    {
      typedef typename DstView::execution_space execution_space;
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      typedef typename IdxView::non_const_value_type index_type;

      size_t errorCount = 0;
      Kokkos::parallel_reduce
        ("Tpetra::MultiVector pack one col debug only",
         range_type (0, idx.size ()),
         PackArraySingleColumnWithBoundsCheck (dst, src, idx, col),
         errorCount);

      if (errorCount != 0) {
        // Go back and find the out-of-bounds entries in the index
        // array.  Performance doesn't matter since we are already in
        // an error state, so we can do this sequentially, on host.
        auto idx_h = Kokkos::create_mirror_view (idx);
        Kokkos::deep_copy (idx_h, idx);

        std::vector<index_type> badIndices;
        const size_type numInds = idx_h.extent (0);
        for (size_type k = 0; k < numInds; ++k) {
          if (idx_h(k) < static_cast<index_type> (0) ||
              idx_h(k) >= static_cast<index_type> (src.extent (0))) {
            badIndices.push_back (idx_h(k));
          }
        }

        TEUCHOS_TEST_FOR_EXCEPTION
          (errorCount != badIndices.size (), std::logic_error,
           "PackArraySingleColumnWithBoundsCheck: errorCount = " << errorCount
           << " != badIndices.size() = " << badIndices.size () << ".  This sho"
           "uld never happen.  Please report this to the Tpetra developers.");

        std::ostringstream os;
        os << "MultiVector single-column pack kernel had "
           << badIndices.size () << " out-of bounds index/ices.  "
          "Here they are: [";
        for (size_t k = 0; k < badIndices.size (); ++k) {
          os << badIndices[k];
          if (k + 1 < badIndices.size ()) {
            os << ", ";
          }
        }
        os << "].";
        throw std::runtime_error (os.str ());
      }
    }
  };


  template <typename DstView, typename SrcView, typename IdxView>
  void
  pack_array_single_column (const DstView& dst,
                            const SrcView& src,
                            const IdxView& idx,
                            const size_t col,
                            const bool debug = true)
  {
    static_assert (Kokkos::Impl::is_view<DstView>::value,
                   "DstView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<SrcView>::value,
                   "SrcView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<IdxView>::value,
                   "IdxView must be a Kokkos::View.");
    static_assert (static_cast<int> (DstView::rank) == 1,
                   "DstView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 2,
                   "SrcView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");

    if (debug) {
      typedef PackArraySingleColumnWithBoundsCheck<DstView,SrcView,IdxView> impl_type;
      impl_type::pack (dst, src, idx, col);
    }
    else {
      typedef PackArraySingleColumn<DstView,SrcView,IdxView> impl_type;
      impl_type::pack (dst, src, idx, col);
    }
  }

  template <typename DstView, typename SrcView, typename IdxView>
  struct PackArrayMultiColumn {
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    DstView dst;
    SrcView src;
    IdxView idx;
    size_t numCols;

    PackArrayMultiColumn (const DstView& dst_,
                          const SrcView& src_,
                          const IdxView& idx_,
                          const size_t numCols_) :
      dst(dst_), src(src_), idx(idx_), numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        dst(offset + j) = src(localRow, j);
      }
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     size_t numCols) {
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      Kokkos::parallel_for
        ("Tpetra::MultiVector pack multicol const stride",
         range_type (0, idx.size ()),
         PackArrayMultiColumn (dst, src, idx, numCols));
    }
  };

  template <typename DstView,
            typename SrcView,
            typename IdxView,
            typename SizeType = typename DstView::execution_space::size_type>
  class PackArrayMultiColumnWithBoundsCheck {
  public:
    using size_type = SizeType;
    using value_type = size_t;

  private:
    DstView dst;
    SrcView src;
    IdxView idx;
    size_type numCols;

  public:
    PackArrayMultiColumnWithBoundsCheck (const DstView& dst_,
                                         const SrcView& src_,
                                         const IdxView& idx_,
                                         const size_type numCols_) :
      dst (dst_), src (src_), idx (idx_), numCols (numCols_) {}

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k, value_type& lclErrorCount) const {
      typedef typename IdxView::non_const_value_type index_type;

      const index_type lclRow = idx(k);
      if (lclRow < static_cast<index_type> (0) ||
          lclRow >= static_cast<index_type> (src.extent (0))) {
        ++lclErrorCount; // failed
      }
      else {
        const size_type offset = k*numCols;
        for (size_type j = 0; j < numCols; ++j) {
          dst(offset + j) = src(lclRow, j);
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void init (value_type& initialErrorCount) const {
      initialErrorCount = 0;
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& dstErrorCount,
          const volatile value_type& srcErrorCount) const
    {
      dstErrorCount += srcErrorCount;
    }

    static void
    pack (const DstView& dst,
          const SrcView& src,
          const IdxView& idx,
          const size_type numCols)
    {
      typedef typename DstView::execution_space execution_space;
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      typedef typename IdxView::non_const_value_type index_type;

      size_t errorCount = 0;
      Kokkos::parallel_reduce
        ("Tpetra::MultiVector pack multicol const stride debug only",
         range_type (0, idx.size ()),
         PackArrayMultiColumnWithBoundsCheck (dst, src, idx, numCols),
         errorCount);
      if (errorCount != 0) {
        // Go back and find the out-of-bounds entries in the index
        // array.  Performance doesn't matter since we are already in
        // an error state, so we can do this sequentially, on host.
        auto idx_h = Kokkos::create_mirror_view (idx);
        Kokkos::deep_copy (idx_h, idx);

        std::vector<index_type> badIndices;
        const size_type numInds = idx_h.extent (0);
        for (size_type k = 0; k < numInds; ++k) {
          if (idx_h(k) < static_cast<index_type> (0) ||
              idx_h(k) >= static_cast<index_type> (src.extent (0))) {
            badIndices.push_back (idx_h(k));
          }
        }

        TEUCHOS_TEST_FOR_EXCEPTION
          (errorCount != badIndices.size (), std::logic_error,
           "PackArraySingleColumnWithBoundsCheck: errorCount = " << errorCount
           << " != badIndices.size() = " << badIndices.size () << ".  This sho"
           "uld never happen.  Please report this to the Tpetra developers.");

        std::ostringstream os;
        os << "Tpetra::MultiVector multiple-column pack kernel had "
           << badIndices.size () << " out-of bounds index/ices (errorCount = "
           << errorCount << "): [";
        for (size_t k = 0; k < badIndices.size (); ++k) {
          os << badIndices[k];
          if (k + 1 < badIndices.size ()) {
            os << ", ";
          }
        }
        os << "].";
        throw std::runtime_error (os.str ());
      }
    }
  };


  template <typename DstView,
            typename SrcView,
            typename IdxView>
  void
  pack_array_multi_column (const DstView& dst,
                           const SrcView& src,
                           const IdxView& idx,
                           const size_t numCols,
                           const bool debug = true)
  {
    static_assert (Kokkos::Impl::is_view<DstView>::value,
                   "DstView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<SrcView>::value,
                   "SrcView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<IdxView>::value,
                   "IdxView must be a Kokkos::View.");
    static_assert (static_cast<int> (DstView::rank) == 1,
                   "DstView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 2,
                   "SrcView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");

    if (debug) {
      typedef PackArrayMultiColumnWithBoundsCheck<DstView,
        SrcView, IdxView> impl_type;
      impl_type::pack (dst, src, idx, numCols);
    }
    else {
      typedef PackArrayMultiColumn<DstView, SrcView, IdxView> impl_type;
      impl_type::pack (dst, src, idx, numCols);
    }
  }

  template <typename DstView, typename SrcView, typename IdxView,
            typename ColView>
  struct PackArrayMultiColumnVariableStride {
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    DstView dst;
    SrcView src;
    IdxView idx;
    ColView col;
    size_t numCols;

    PackArrayMultiColumnVariableStride (const DstView& dst_,
                                        const SrcView& src_,
                                        const IdxView& idx_,
                                        const ColView& col_,
                                        const size_t numCols_) :
      dst(dst_), src(src_), idx(idx_), col(col_), numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_type k) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        dst(offset + j) = src(localRow, col(j));
      }
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     const ColView& col,
                     size_t numCols) {
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      Kokkos::parallel_for
        ("Tpetra::MultiVector pack multicol var stride",
         range_type (0, idx.size ()),
         PackArrayMultiColumnVariableStride (dst, src, idx, col, numCols));
    }
  };

  template <typename DstView,
            typename SrcView,
            typename IdxView,
            typename ColView,
            typename SizeType = typename DstView::execution_space::size_type>
  class PackArrayMultiColumnVariableStrideWithBoundsCheck {
  public:
    using size_type = SizeType;
    using value_type = size_t;

  private:
    DstView dst;
    SrcView src;
    IdxView idx;
    ColView col;
    size_type numCols;

  public:
    PackArrayMultiColumnVariableStrideWithBoundsCheck (const DstView& dst_,
                                                       const SrcView& src_,
                                                       const IdxView& idx_,
                                                       const ColView& col_,
                                                       const size_type numCols_) :
      dst (dst_), src (src_), idx (idx_), col (col_), numCols (numCols_) {}

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k, value_type& lclErrorCount) const {
      typedef typename IdxView::non_const_value_type row_index_type;
      typedef typename ColView::non_const_value_type col_index_type;

      const row_index_type lclRow = idx(k);
      if (lclRow < static_cast<row_index_type> (0) ||
          lclRow >= static_cast<row_index_type> (src.extent (0))) {
        ++lclErrorCount = 0;
      }
      else {
        const size_type offset = k*numCols;
        for (size_type j = 0; j < numCols; ++j) {
          const col_index_type lclCol = col(j);
          if (Impl::outOfBounds<col_index_type> (lclCol, src.extent (1))) {
            ++lclErrorCount = 0;
          }
          else { // all indices are valid; do the assignment
            dst(offset + j) = src(lclRow, lclCol);
          }
        }
      }
    }

    KOKKOS_INLINE_FUNCTION void
    init (value_type& initialErrorCount) const {
      initialErrorCount = 0;
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& dstErrorCount,
          const volatile value_type& srcErrorCount) const
    {
      dstErrorCount += srcErrorCount;
    }

    static void
    pack (const DstView& dst,
          const SrcView& src,
          const IdxView& idx,
          const ColView& col,
          const size_type numCols)
    {
      using execution_space = typename DstView::execution_space;
      using range_type = Kokkos::RangePolicy<execution_space, size_type>;
      using row_index_type = typename IdxView::non_const_value_type;
      using col_index_type = typename ColView::non_const_value_type;

      size_t errorCount = 0;
      Kokkos::parallel_reduce
        ("Tpetra::MultiVector pack multicol var stride debug only",
         range_type (0, idx.size ()),
         PackArrayMultiColumnVariableStrideWithBoundsCheck (dst, src, idx,
                                                            col, numCols),
         errorCount);
      if (errorCount != 0) {
        constexpr size_t maxNumBadIndicesToPrint = 100;

        std::ostringstream os; // for error reporting
        os << "Tpetra::MultiVector multicolumn variable stride pack kernel "
          "found " << errorCount
          << " error" << (errorCount != size_t (1) ? "s" : "") << ".  ";

        // Go back and find any out-of-bounds entries in the array of
        // row indices.  Performance doesn't matter since we are already
        // in an error state, so we can do this sequentially, on host.
        auto idx_h = Kokkos::create_mirror_view (idx);
        Kokkos::deep_copy (idx_h, idx);

        std::vector<row_index_type> badRows;
        const size_type numRowInds = idx_h.extent (0);
        for (size_type k = 0; k < numRowInds; ++k) {
          if (Impl::outOfBounds<row_index_type> (idx_h(k), src.extent (0))) {
            badRows.push_back (idx_h(k));
          }
        }

        if (badRows.size () != 0) {
          os << badRows.size () << " out-of-bounds row ind"
             << (badRows.size () != size_t (1) ? "ices" : "ex");
          if (badRows.size () <= maxNumBadIndicesToPrint) {
            os << ": [";
            for (size_t k = 0; k < badRows.size (); ++k) {
              os << badRows[k];
              if (k + 1 < badRows.size ()) {
                os << ", ";
              }
            }
            os << "].  ";
          }
          else {
            os << ".  ";
          }
        }
        else {
          os << "No out-of-bounds row indices.  ";
        }

        // Go back and find any out-of-bounds entries in the array
        // of column indices.
        auto col_h = Kokkos::create_mirror_view (col);
        Kokkos::deep_copy (col_h, col);

        std::vector<col_index_type> badCols;
        const size_type numColInds = col_h.extent (0);
        for (size_type k = 0; k < numColInds; ++k) {
          if (Impl::outOfBounds<col_index_type> (col_h(k), src.extent (1))) {
            badCols.push_back (col_h(k));
          }
        }

        if (badCols.size () != 0) {
          os << badCols.size () << " out-of-bounds column ind"
             << (badCols.size () != size_t (1) ? "ices" : "ex");
          if (badCols.size () <= maxNumBadIndicesToPrint) {
            os << ": [";
            for (size_t k = 0; k < badCols.size (); ++k) {
              os << badCols[k];
              if (k + 1 < badCols.size ()) {
                os << ", ";
              }
            }
            os << "].  ";
          }
          else {
            os << ".  ";
          }
        }
        else {
          os << "No out-of-bounds column indices.  ";
        }

        TEUCHOS_TEST_FOR_EXCEPTION
          (errorCount != 0 && badRows.size () == 0 && badCols.size () == 0,
           std::logic_error, "Tpetra::MultiVector variable stride pack "
           "kernel reports errorCount=" << errorCount << ", but we failed "
           "to find any bad rows or columns.  This should never happen.  "
           "Please report this to the Tpetra developers.");

        throw std::runtime_error (os.str ());
      } // hasErr
    }
  };

  template <typename DstView,
            typename SrcView,
            typename IdxView,
            typename ColView>
  void
  pack_array_multi_column_variable_stride (const DstView& dst,
                                           const SrcView& src,
                                           const IdxView& idx,
                                           const ColView& col,
                                           const size_t numCols,
                                           const bool debug = true)
  {
    static_assert (Kokkos::Impl::is_view<DstView>::value,
                   "DstView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<SrcView>::value,
                   "SrcView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<IdxView>::value,
                   "IdxView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ColView>::value,
                   "ColView must be a Kokkos::View.");
    static_assert (static_cast<int> (DstView::rank) == 1,
                   "DstView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 2,
                   "SrcView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (ColView::rank) == 1,
                   "ColView must be a rank-1 Kokkos::View.");

    if (debug) {
      typedef PackArrayMultiColumnVariableStrideWithBoundsCheck<DstView,
        SrcView, IdxView, ColView> impl_type;
      impl_type::pack (dst, src, idx, col, numCols);
    }
    else {
      typedef PackArrayMultiColumnVariableStride<DstView,
        SrcView, IdxView, ColView> impl_type;
      impl_type::pack (dst, src, idx, col, numCols);
    }
  }

  // Tag types to indicate whether to use atomic updates in the
  // various CombineMode "Op"s.
  struct atomic_tag {};
  struct nonatomic_tag {};

  template<class SC>
  struct AddOp {
    KOKKOS_INLINE_FUNCTION
    void operator() (atomic_tag, SC& dest, const SC& src) const {
      Kokkos::atomic_add (&dest, src);
    }

    KOKKOS_INLINE_FUNCTION
    void operator() (nonatomic_tag, SC& dest, const SC& src) const {
      dest += src;
    }
  };

  template<class SC>
  struct InsertOp {
    // There's no point to using Kokkos::atomic_assign for the REPLACE
    // or INSERT CombineModes, since this is not a well-defined
    // reduction for MultiVector anyway.  See GitHub Issue #4417
    // (which this fixes).
    KOKKOS_INLINE_FUNCTION
    void operator() (atomic_tag, SC& dest, const SC& src) const {
      dest = src;
    }

    KOKKOS_INLINE_FUNCTION
    void operator() (nonatomic_tag, SC& dest, const SC& src) const {
      dest = src;
    }
  };

  // Kokkos::Impl::atomic_fetch_oper wants a class like this.
  template<class Scalar1, class Scalar2>
  struct AbsMaxOper {
    KOKKOS_INLINE_FUNCTION
    static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
      const auto val1_abs = Kokkos::ArithTraits<Scalar1>::abs(val1);
      const auto val2_abs = Kokkos::ArithTraits<Scalar2>::abs(val2);
      return val1_abs > val2_abs ? Scalar1(val1_abs) : Scalar1(val2_abs);
    }
  };

  template <typename SC>
  struct AbsMaxOp {
    KOKKOS_INLINE_FUNCTION
    void operator() (atomic_tag, SC& dest, const SC& src) const {
      Kokkos::Impl::atomic_fetch_oper (AbsMaxOper<SC, SC> (), &dest, src);
    }

    KOKKOS_INLINE_FUNCTION
    void operator() (nonatomic_tag, SC& dest, const SC& src) const {
      dest = AbsMaxOper<SC, SC> ().apply (dest, src);
    }
  };

  template <typename ExecutionSpace,
            typename DstView,
            typename SrcView,
            typename IdxView,
            typename Op>
  class UnpackArrayMultiColumn {
  private:
    static_assert (Kokkos::Impl::is_view<DstView>::value,
                   "DstView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<SrcView>::value,
                   "SrcView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<IdxView>::value,
                   "IdxView must be a Kokkos::View.");
    static_assert (static_cast<int> (DstView::rank) == 2,
                   "DstView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 1,
                   "SrcView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");

  public:
    typedef typename ExecutionSpace::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

  private:
    DstView dst;
    SrcView src;
    IdxView idx;
    Op op;
    size_t numCols;

  public:
    UnpackArrayMultiColumn (const ExecutionSpace& /* execSpace */,
                            const DstView& dst_,
                            const SrcView& src_,
                            const IdxView& idx_,
                            const Op& op_,
                            const size_t numCols_) :
      dst (dst_),
      src (src_),
      idx (idx_),
      op (op_),
      numCols (numCols_)
    {}

    template<class TagType>
    KOKKOS_INLINE_FUNCTION void
    operator() (TagType tag, const size_type k) const
    {
      static_assert
        (std::is_same<TagType, atomic_tag>::value ||
         std::is_same<TagType, nonatomic_tag>::value,
         "TagType must be atomic_tag or nonatomic_tag.");

      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (tag, dst(localRow, j), src(offset+j));
      }
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const Op& op,
            const size_t numCols,
            const bool use_atomic_updates)
    {
      if (use_atomic_updates) {
        using range_type =
          Kokkos::RangePolicy<atomic_tag, execution_space, size_type>;
        Kokkos::parallel_for
          ("Tpetra::MultiVector unpack const stride atomic",
           range_type (0, idx.size ()),
           UnpackArrayMultiColumn (execSpace, dst, src, idx, op, numCols));
      }
      else {
        using range_type =
          Kokkos::RangePolicy<nonatomic_tag, execution_space, size_type>;
        Kokkos::parallel_for
          ("Tpetra::MultiVector unpack const stride nonatomic",
           range_type (0, idx.size ()),
           UnpackArrayMultiColumn (execSpace, dst, src, idx, op, numCols));
      }
    }
  };

  template <typename ExecutionSpace,
            typename DstView,
            typename SrcView,
            typename IdxView,
            typename Op,
            typename SizeType = typename ExecutionSpace::execution_space::size_type>
  class UnpackArrayMultiColumnWithBoundsCheck {
  private:
    static_assert (Kokkos::Impl::is_view<DstView>::value,
                   "DstView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<SrcView>::value,
                   "SrcView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<IdxView>::value,
                   "IdxView must be a Kokkos::View.");
    static_assert (static_cast<int> (DstView::rank) == 2,
                   "DstView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 1,
                   "SrcView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");
    static_assert (std::is_integral<SizeType>::value,
                   "SizeType must be a built-in integer type.");

  public:
    using execution_space = typename ExecutionSpace::execution_space;
    using size_type = SizeType;
    using value_type = size_t;

  private:
    DstView dst;
    SrcView src;
    IdxView idx;
    Op op;
    size_type numCols;

  public:
    UnpackArrayMultiColumnWithBoundsCheck (const ExecutionSpace& /* execSpace */,
                                           const DstView& dst_,
                                           const SrcView& src_,
                                           const IdxView& idx_,
                                           const Op& op_,
                                           const size_type numCols_) :
      dst (dst_),
      src (src_),
      idx (idx_),
      op (op_),
      numCols (numCols_)
    {}

    template<class TagType>
    KOKKOS_INLINE_FUNCTION void
    operator() (TagType tag,
                const size_type k,
                size_t& lclErrCount) const
    {
      static_assert
        (std::is_same<TagType, atomic_tag>::value ||
         std::is_same<TagType, nonatomic_tag>::value,
         "TagType must be atomic_tag or nonatomic_tag.");
      using index_type = typename IdxView::non_const_value_type;

      const index_type lclRow = idx(k);
      if (lclRow < static_cast<index_type> (0) ||
          lclRow >= static_cast<index_type> (dst.extent (0))) {
        ++lclErrCount;
      }
      else {
        const size_type offset = k*numCols;
        for (size_type j = 0; j < numCols; ++j) {
          op (tag, dst(lclRow,j), src(offset+j));
        }
      }
    }

    template<class TagType>
    KOKKOS_INLINE_FUNCTION void
    init (TagType, size_t& initialErrorCount) const {
      initialErrorCount = 0;
    }

    template<class TagType>
    KOKKOS_INLINE_FUNCTION void
    join (TagType,
          volatile size_t& dstErrorCount,
          const volatile size_t& srcErrorCount) const
    {
      dstErrorCount += srcErrorCount;
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const Op& op,
            const size_type numCols,
            const bool use_atomic_updates)
    {
      using index_type = typename IdxView::non_const_value_type;

      size_t errorCount = 0;
      if (use_atomic_updates) {
        using range_type =
          Kokkos::RangePolicy<atomic_tag, execution_space, size_type>;
        Kokkos::parallel_reduce
          ("Tpetra::MultiVector unpack multicol const stride atomic debug only",
           range_type (0, idx.size ()),
           UnpackArrayMultiColumnWithBoundsCheck (execSpace, dst, src,
                                                  idx, op, numCols),
           errorCount);
      }
      else {
        using range_type =
          Kokkos::RangePolicy<nonatomic_tag, execution_space, size_type>;
        Kokkos::parallel_reduce
          ("Tpetra::MultiVector unpack multicol const stride nonatomic debug only",
           range_type (0, idx.size ()),
           UnpackArrayMultiColumnWithBoundsCheck (execSpace, dst, src,
                                                  idx, op, numCols),
           errorCount);
      }

      if (errorCount != 0) {
        // Go back and find the out-of-bounds entries in the index
        // array.  Performance doesn't matter since we are already in
        // an error state, so we can do this sequentially, on host.
        auto idx_h = Kokkos::create_mirror_view (idx);
        Kokkos::deep_copy (idx_h, idx);

        std::vector<index_type> badIndices;
        const size_type numInds = idx_h.extent (0);
        for (size_type k = 0; k < numInds; ++k) {
          if (idx_h(k) < static_cast<index_type> (0) ||
              idx_h(k) >= static_cast<index_type> (dst.extent (0))) {
            badIndices.push_back (idx_h(k));
          }
        }

        if (errorCount != badIndices.size ()) {
          std::ostringstream os;
          os << "MultiVector unpack kernel: errorCount = " << errorCount
             << " != badIndices.size() = " << badIndices.size ()
             << ".  This should never happen.  "
            "Please report this to the Tpetra developers.";
          throw std::logic_error (os.str ());
        }

        std::ostringstream os;
        os << "MultiVector unpack kernel had " << badIndices.size ()
           << " out-of bounds index/ices.  Here they are: [";
        for (size_t k = 0; k < badIndices.size (); ++k) {
          os << badIndices[k];
          if (k + 1 < badIndices.size ()) {
            os << ", ";
          }
        }
        os << "].";
        throw std::runtime_error (os.str ());
      }
    }
  };

  template <typename ExecutionSpace,
            typename DstView,
            typename SrcView,
            typename IdxView,
            typename Op>
  void
  unpack_array_multi_column (const ExecutionSpace& execSpace,
                             const DstView& dst,
                             const SrcView& src,
                             const IdxView& idx,
                             const Op& op,
                             const size_t numCols,
                             const bool use_atomic_updates,
                             const bool debug)
  {
    static_assert (Kokkos::Impl::is_view<DstView>::value,
                   "DstView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<SrcView>::value,
                   "SrcView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<IdxView>::value,
                   "IdxView must be a Kokkos::View.");
    static_assert (static_cast<int> (DstView::rank) == 2,
                   "DstView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 1,
                   "SrcView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");

    if (debug) {
      typedef UnpackArrayMultiColumnWithBoundsCheck<ExecutionSpace,
        DstView, SrcView, IdxView, Op> impl_type;
      impl_type::unpack (execSpace, dst, src, idx, op, numCols,
                         use_atomic_updates);
    }
    else {
      typedef UnpackArrayMultiColumn<ExecutionSpace,
        DstView, SrcView, IdxView, Op> impl_type;
      impl_type::unpack (execSpace, dst, src, idx, op, numCols,
                         use_atomic_updates);
    }
  }

  template <typename ExecutionSpace,
            typename DstView,
            typename SrcView,
            typename IdxView,
            typename ColView,
            typename Op>
  class UnpackArrayMultiColumnVariableStride {
  private:
    static_assert (Kokkos::Impl::is_view<DstView>::value,
                   "DstView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<SrcView>::value,
                   "SrcView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<IdxView>::value,
                   "IdxView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ColView>::value,
                   "ColView must be a Kokkos::View.");
    static_assert (static_cast<int> (DstView::rank) == 2,
                   "DstView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 1,
                   "SrcView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (ColView::rank) == 1,
                   "ColView must be a rank-1 Kokkos::View.");

  public:
    using execution_space = typename ExecutionSpace::execution_space;
    using size_type = typename execution_space::size_type;

  private:
    DstView dst;
    SrcView src;
    IdxView idx;
    ColView col;
    Op op;
    size_t numCols;

  public:
    UnpackArrayMultiColumnVariableStride (const ExecutionSpace& /* execSpace */,
                                          const DstView& dst_,
                                          const SrcView& src_,
                                          const IdxView& idx_,
                                          const ColView& col_,
                                          const Op& op_,
                                          const size_t numCols_) :
      dst (dst_),
      src (src_),
      idx (idx_),
      col (col_),
      op (op_),
      numCols (numCols_)
    {}

    template<class TagType>
    KOKKOS_INLINE_FUNCTION void
    operator() (TagType tag, const size_type k) const
    {
      static_assert
        (std::is_same<TagType, atomic_tag>::value ||
         std::is_same<TagType, nonatomic_tag>::value,
         "TagType must be atomic_tag or nonatomic_tag.");

      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (tag, dst(localRow, col(j)), src(offset+j));
      }
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const ColView& col,
            const Op& op,
            const size_t numCols,
            const bool use_atomic_updates)
    {
      if (use_atomic_updates) {
        using range_type =
          Kokkos::RangePolicy<atomic_tag, execution_space, size_type>;
        Kokkos::parallel_for
          ("Tpetra::MultiVector unpack var stride atomic",
           range_type (0, idx.size ()),
           UnpackArrayMultiColumnVariableStride (execSpace, dst, src,
                                                 idx, col, op, numCols));
      }
      else {
        using range_type =
          Kokkos::RangePolicy<nonatomic_tag, execution_space, size_type>;
        Kokkos::parallel_for
          ("Tpetra::MultiVector unpack var stride nonatomic",
           range_type (0, idx.size ()),
           UnpackArrayMultiColumnVariableStride (execSpace, dst, src,
                                                 idx, col, op, numCols));
      }
    }
  };

  template <typename ExecutionSpace,
            typename DstView,
            typename SrcView,
            typename IdxView,
            typename ColView,
            typename Op,
            typename SizeType = typename ExecutionSpace::execution_space::size_type>
  class UnpackArrayMultiColumnVariableStrideWithBoundsCheck {
  private:
    static_assert (Kokkos::Impl::is_view<DstView>::value,
                   "DstView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<SrcView>::value,
                   "SrcView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<IdxView>::value,
                   "IdxView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ColView>::value,
                   "ColView must be a Kokkos::View.");
    static_assert (static_cast<int> (DstView::rank) == 2,
                   "DstView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 1,
                   "SrcView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (ColView::rank) == 1,
                   "ColView must be a rank-1 Kokkos::View.");
    static_assert (std::is_integral<SizeType>::value,
                   "SizeType must be a built-in integer type.");

  public:
    using execution_space = typename ExecutionSpace::execution_space;
    using size_type = SizeType;
    using value_type = size_t;

  private:
    DstView dst;
    SrcView src;
    IdxView idx;
    ColView col;
    Op op;
    size_type numCols;

  public:
    UnpackArrayMultiColumnVariableStrideWithBoundsCheck
      (const ExecutionSpace& /* execSpace */,
       const DstView& dst_,
       const SrcView& src_,
       const IdxView& idx_,
       const ColView& col_,
       const Op& op_,
       const size_t numCols_) :
        dst (dst_),
        src (src_),
        idx (idx_),
        col (col_),
        op (op_),
        numCols (numCols_)
    {}

    template<class TagType>
    KOKKOS_INLINE_FUNCTION void
    operator() (TagType tag,
                const size_type k,
                value_type& lclErrorCount) const
    {
      static_assert
        (std::is_same<TagType, atomic_tag>::value ||
         std::is_same<TagType, nonatomic_tag>::value,
         "TagType must be atomic_tag or nonatomic_tag.");
      using row_index_type = typename IdxView::non_const_value_type;
      using col_index_type = typename ColView::non_const_value_type;

      const row_index_type lclRow = idx(k);
      if (lclRow < static_cast<row_index_type> (0) ||
          lclRow >= static_cast<row_index_type> (dst.extent (0))) {
        ++lclErrorCount;
      }
      else {
        const size_type offset = k * numCols;
        for (size_type j = 0; j < numCols; ++j) {
          const col_index_type lclCol = col(j);
          if (Impl::outOfBounds<col_index_type> (lclCol, dst.extent (1))) {
            ++lclErrorCount;
          }
          else { // all indices are valid; apply the op
            op (tag, dst(lclRow, col(j)), src(offset+j));
          }
        }
      }
    }

    KOKKOS_INLINE_FUNCTION void
    init (value_type& initialErrorCount) const {
      initialErrorCount = 0;
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& dstErrorCount,
          const volatile value_type& srcErrorCount) const
    {
      dstErrorCount += srcErrorCount;
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const ColView& col,
            const Op& op,
            const size_type numCols,
            const bool use_atomic_updates)
    {
      using row_index_type = typename IdxView::non_const_value_type;
      using col_index_type = typename ColView::non_const_value_type;

      size_t errorCount = 0;
      if (use_atomic_updates) {
        using range_type =
          Kokkos::RangePolicy<atomic_tag, execution_space, size_type>;
        Kokkos::parallel_reduce
          ("Tpetra::MultiVector unpack var stride atomic debug only",
           range_type (0, idx.size ()),
           UnpackArrayMultiColumnVariableStrideWithBoundsCheck
             (execSpace, dst, src, idx, col, op, numCols),
           errorCount);
      }
      else {
        using range_type =
          Kokkos::RangePolicy<nonatomic_tag, execution_space, size_type>;
        Kokkos::parallel_reduce
          ("Tpetra::MultiVector unpack var stride nonatomic debug only",
           range_type (0, idx.size ()),
           UnpackArrayMultiColumnVariableStrideWithBoundsCheck
             (execSpace, dst, src, idx, col, op, numCols),
           errorCount);
      }

      if (errorCount != 0) {
        constexpr size_t maxNumBadIndicesToPrint = 100;

        std::ostringstream os; // for error reporting
        os << "Tpetra::MultiVector multicolumn variable stride unpack kernel "
          "found " << errorCount
          << " error" << (errorCount != size_t (1) ? "s" : "") << ".  ";

        // Go back and find any out-of-bounds entries in the array of
        // row indices.  Performance doesn't matter since we are
        // already in an error state, so we can do this sequentially,
        // on host.
        auto idx_h = Kokkos::create_mirror_view (idx);
        Kokkos::deep_copy (idx_h, idx);

        std::vector<row_index_type> badRows;
        const size_type numRowInds = idx_h.extent (0);
        for (size_type k = 0; k < numRowInds; ++k) {
          if (idx_h(k) < static_cast<row_index_type> (0) ||
              idx_h(k) >= static_cast<row_index_type> (dst.extent (0))) {
            badRows.push_back (idx_h(k));
          }
        }

        if (badRows.size () != 0) {
          os << badRows.size () << " out-of-bounds row ind"
             << (badRows.size () != size_t (1) ? "ices" : "ex");
          if (badRows.size () <= maxNumBadIndicesToPrint) {
            os << ": [";
            for (size_t k = 0; k < badRows.size (); ++k) {
              os << badRows[k];
              if (k + 1 < badRows.size ()) {
                os << ", ";
              }
            }
            os << "].  ";
          }
          else {
            os << ".  ";
          }
        }
        else {
          os << "No out-of-bounds row indices.  ";
        }

        // Go back and find any out-of-bounds entries in the array
        // of column indices.
        auto col_h = Kokkos::create_mirror_view (col);
        Kokkos::deep_copy (col_h, col);

        std::vector<col_index_type> badCols;
        const size_type numColInds = col_h.extent (0);
        for (size_type k = 0; k < numColInds; ++k) {
          if (Impl::outOfBounds<col_index_type> (col_h(k), dst.extent (1))) {
            badCols.push_back (col_h(k));
          }
        }

        if (badCols.size () != 0) {
          os << badCols.size () << " out-of-bounds column ind"
             << (badCols.size () != size_t (1) ? "ices" : "ex");
          if (badCols.size () <= maxNumBadIndicesToPrint) {
            for (size_t k = 0; k < badCols.size (); ++k) {
              os << ": [";
              os << badCols[k];
              if (k + 1 < badCols.size ()) {
                os << ", ";
              }
            }
            os << "].  ";
          }
          else {
            os << ".  ";
          }
        }
        else {
          os << "No out-of-bounds column indices.  ";
        }

        TEUCHOS_TEST_FOR_EXCEPTION
          (errorCount != 0 && badRows.size () == 0 && badCols.size () == 0,
           std::logic_error, "Tpetra::MultiVector variable stride unpack "
           "kernel reports errorCount=" << errorCount << ", but we failed "
           "to find any bad rows or columns.  This should never happen.  "
           "Please report this to the Tpetra developers.");

        throw std::runtime_error (os.str ());
      } // hasErr
    }
  };

  template <typename ExecutionSpace,
            typename DstView,
            typename SrcView,
            typename IdxView,
            typename ColView,
            typename Op>
  void
  unpack_array_multi_column_variable_stride (const ExecutionSpace& execSpace,
                                             const DstView& dst,
                                             const SrcView& src,
                                             const IdxView& idx,
                                             const ColView& col,
                                             const Op& op,
                                             const size_t numCols,
                                             const bool use_atomic_updates,
                                             const bool debug)
  {
    static_assert (Kokkos::Impl::is_view<DstView>::value,
                   "DstView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<SrcView>::value,
                   "SrcView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<IdxView>::value,
                   "IdxView must be a Kokkos::View.");
    static_assert (Kokkos::Impl::is_view<ColView>::value,
                   "ColView must be a Kokkos::View.");
    static_assert (static_cast<int> (DstView::rank) == 2,
                   "DstView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 1,
                   "SrcView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (ColView::rank) == 1,
                   "ColView must be a rank-1 Kokkos::View.");

    if (debug) {
      using impl_type =
        UnpackArrayMultiColumnVariableStrideWithBoundsCheck<ExecutionSpace,
          DstView, SrcView, IdxView, ColView, Op>;
      impl_type::unpack (execSpace, dst, src, idx, col, op, numCols,
                         use_atomic_updates);
    }
    else {
      using impl_type = UnpackArrayMultiColumnVariableStride<ExecutionSpace,
        DstView, SrcView, IdxView, ColView, Op>;
      impl_type::unpack (execSpace, dst, src, idx, col, op, numCols,
                         use_atomic_updates);
    }
  }

  template <typename DstView, typename SrcView,
            typename DstIdxView, typename SrcIdxView>
  struct PermuteArrayMultiColumn {
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    DstView dst;
    SrcView src;
    DstIdxView dst_idx;
    SrcIdxView src_idx;
    size_t numCols;

    PermuteArrayMultiColumn (const DstView& dst_,
                             const SrcView& src_,
                             const DstIdxView& dst_idx_,
                             const SrcIdxView& src_idx_,
                             const size_t numCols_) :
      dst(dst_), src(src_), dst_idx(dst_idx_), src_idx(src_idx_),
      numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      for (size_t j = 0; j < numCols; ++j) {
        dst(toRow, j) = src(fromRow, j);
      }
    }

    static void
    permute (const DstView& dst,
	     const SrcView& src,
	     const DstIdxView& dst_idx,
	     const SrcIdxView& src_idx,
	     const size_t numCols)
    {
      using range_type = Kokkos::RangePolicy<execution_space, size_type>;      
      const size_type n = std::min (dst_idx.size (), src_idx.size ());
      Kokkos::parallel_for
	("Tpetra::MultiVector permute multicol const stride",
	 range_type (0, n),
	 PermuteArrayMultiColumn (dst, src, dst_idx, src_idx, numCols));
    }
  };

  // To do:  Add enable_if<> restrictions on DstView::Rank == 1,
  // SrcView::Rank == 2
  template <typename DstView, typename SrcView,
            typename DstIdxView, typename SrcIdxView>
  void permute_array_multi_column(const DstView& dst,
                                  const SrcView& src,
                                  const DstIdxView& dst_idx,
                                  const SrcIdxView& src_idx,
                                  size_t numCols) {
    PermuteArrayMultiColumn<DstView,SrcView,DstIdxView,SrcIdxView>::permute(
      dst, src, dst_idx, src_idx, numCols);
  }

  template <typename DstView, typename SrcView,
            typename DstIdxView, typename SrcIdxView,
            typename DstColView, typename SrcColView>
  struct PermuteArrayMultiColumnVariableStride {
    typedef typename DstView::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

    DstView dst;
    SrcView src;
    DstIdxView dst_idx;
    SrcIdxView src_idx;
    DstColView dst_col;
    SrcColView src_col;
    size_t numCols;

    PermuteArrayMultiColumnVariableStride(const DstView& dst_,
                                          const SrcView& src_,
                                          const DstIdxView& dst_idx_,
                                          const SrcIdxView& src_idx_,
                                          const DstColView& dst_col_,
                                          const SrcColView& src_col_,
                                          const size_t numCols_) :
      dst(dst_), src(src_), dst_idx(dst_idx_), src_idx(src_idx_),
      dst_col(dst_col_), src_col(src_col_),
      numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      for (size_t j = 0; j < numCols; ++j) {
        dst(toRow, dst_col(j)) = src(fromRow, src_col(j));
      }
    }

    static void
    permute (const DstView& dst,
	     const SrcView& src,
	     const DstIdxView& dst_idx,
	     const SrcIdxView& src_idx,
	     const DstColView& dst_col,
	     const SrcColView& src_col,
	     const size_t numCols)
    {
      using range_type = Kokkos::RangePolicy<execution_space, size_type>;      
      const size_type n = std::min (dst_idx.size (), src_idx.size ());
      Kokkos::parallel_for
	("Tpetra::MultiVector permute multicol var stride",
	 range_type (0, n),
	 PermuteArrayMultiColumnVariableStride (dst, src, dst_idx, src_idx,
						dst_col, src_col, numCols));
    }
  };

  // To do:  Add enable_if<> restrictions on DstView::Rank == 1,
  // SrcView::Rank == 2
  template <typename DstView, typename SrcView,
            typename DstIdxView, typename SrcIdxView,
            typename DstColView, typename SrcColView>
  void permute_array_multi_column_variable_stride(const DstView& dst,
                                                  const SrcView& src,
                                                  const DstIdxView& dst_idx,
                                                  const SrcIdxView& src_idx,
                                                  const DstColView& dst_col,
                                                  const SrcColView& src_col,
                                                  size_t numCols) {
    PermuteArrayMultiColumnVariableStride<DstView,SrcView,
      DstIdxView,SrcIdxView,DstColView,SrcColView>::permute(
      dst, src, dst_idx, src_idx, dst_col, src_col, numCols);
  }

} // Details namespace
} // KokkosRefactor namespace
} // Tpetra namespace

#endif // TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_DIST_OBJECT_KERNELS_HPP
