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

    PackArraySingleColumn(const DstView& dst_,
                          const SrcView& src_,
                          const IdxView& idx_,
                          size_t col_) :
      dst(dst_), src(src_), idx(idx_), col(col_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      dst(k) = src(idx(k), col);
    }

    static void
    pack (const DstView& dst,
          const SrcView& src,
          const IdxView& idx,
          size_t col)
    {
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      Kokkos::parallel_for (range_type (0, idx.size ()),
                            PackArraySingleColumn (dst,src,idx,col));
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
    //! We use int as a Boolean (Kokkos doesn't allow bool reduction types).
    typedef int value_type;

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
    operator() (const size_type& k, value_type& result) const {
      typedef typename IdxView::non_const_value_type index_type;

      const index_type lclRow = idx(k);
      if (lclRow < static_cast<index_type> (0) ||
          lclRow >= static_cast<index_type> (src.extent (0))) {
        result = 0; // failed!
      }
      else {
        dst(k) = src(lclRow, col);
      }
    }

    KOKKOS_INLINE_FUNCTION
    void init (value_type& initialResult) const {
      initialResult = 1; // success
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& dstResult,
          const volatile value_type& srcResult) const
    {
      dstResult = (dstResult == 0 || srcResult == 0) ? 0 : 1;
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

      int result = 1;
      Kokkos::parallel_reduce (range_type (0, idx.size ()),
                               PackArraySingleColumnWithBoundsCheck (dst, src,
                                                                     idx, col),
                               result);
      if (result != 1) {
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

    PackArrayMultiColumn(const DstView& dst_,
                         const SrcView& src_,
                         const IdxView& idx_,
                         size_t numCols_) :
      dst(dst_), src(src_), idx(idx_), numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        dst(offset + j) = src(localRow, j);
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     size_t numCols) {
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      Kokkos::parallel_for (range_type (0, idx.size ()),
                            PackArrayMultiColumn (dst,src,idx,numCols));
    }
  };

  template <typename DstView,
            typename SrcView,
            typename IdxView,
            typename SizeType = typename DstView::execution_space::size_type>
  class PackArrayMultiColumnWithBoundsCheck {
  public:
    typedef SizeType size_type;
    //! We use int as a Boolean (Kokkos doesn't allow bool reduction types).
    typedef int value_type;

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
    operator() (const size_type& k, value_type& result) const {
      typedef typename IdxView::non_const_value_type index_type;

      const index_type lclRow = idx(k);
      if (lclRow < static_cast<index_type> (0) ||
          lclRow >= static_cast<index_type> (src.extent (0))) {
        result = 0; // failed!
      }
      else {
        const size_type offset = k*numCols;
        for (size_type j = 0; j < numCols; ++j) {
          dst(offset + j) = src(lclRow, j);
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void init (value_type& initialResult) const {
      initialResult = 1; // success
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& dstResult,
          const volatile value_type& srcResult) const
    {
      dstResult = (dstResult == 0 || srcResult == 0) ? 0 : 1;
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

      int result = 1;
      Kokkos::parallel_reduce (range_type (0, idx.size ()),
                               PackArrayMultiColumnWithBoundsCheck (dst, src,
                                                                    idx, numCols),
                               result);
      if (result != 1) {
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

        std::ostringstream os;
        os << "MultiVector multiple-column pack kernel had "
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

    PackArrayMultiColumnVariableStride(const DstView& dst_,
                                       const SrcView& src_,
                                       const IdxView& idx_,
                                       const ColView& col_,
                                       size_t numCols_) :
      dst(dst_), src(src_), idx(idx_), col(col_), numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j)
        dst(offset + j) = src(localRow, col(j));
    }

    static void pack(const DstView& dst,
                     const SrcView& src,
                     const IdxView& idx,
                     const ColView& col,
                     size_t numCols) {
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      Kokkos::parallel_for (range_type (0, idx.size ()),
                            PackArrayMultiColumnVariableStride(
                              dst,src,idx,col,numCols) );
    }
  };

  template <typename DstView,
            typename SrcView,
            typename IdxView,
            typename ColView,
            typename SizeType = typename DstView::execution_space::size_type>
  class PackArrayMultiColumnVariableStrideWithBoundsCheck {
  public:
    typedef SizeType size_type;
    //! We use int as a Boolean (Kokkos doesn't allow bool reduction types).
    typedef Kokkos::pair<int, int> value_type;

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
    operator() (const size_type& k, value_type& result) const {
      typedef typename IdxView::non_const_value_type row_index_type;
      typedef typename ColView::non_const_value_type col_index_type;

      const row_index_type lclRow = idx(k);
      if (lclRow < static_cast<row_index_type> (0) ||
          lclRow >= static_cast<row_index_type> (src.extent (0))) {
        result.first = 0; // failed!
      }
      else {
        const size_type offset = k*numCols;
        for (size_type j = 0; j < numCols; ++j) {
          const col_index_type lclCol = col(j);
          if (Impl::outOfBounds<col_index_type> (lclCol, src.extent (1))) {
            result.second = 0; // failed!
          }
          else { // all indices are valid; do the assignment
            dst(offset + j) = src(lclRow, lclCol);
          }
        }
      }
    }

    KOKKOS_INLINE_FUNCTION void
    init (value_type& initialResult) const {
      initialResult.first = 1; // success
      initialResult.second = 1; // success
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& dstResult,
          const volatile value_type& srcResult) const
    {
      dstResult.first = (dstResult.first == 0 || srcResult.first == 0) ? 0 : 1;
      dstResult.second = (dstResult.second == 0 || srcResult.second == 0) ? 0 : 1;
    }

    static void
    pack (const DstView& dst,
          const SrcView& src,
          const IdxView& idx,
          const ColView& col,
          const size_type numCols)
    {
      using Kokkos::parallel_reduce;
      typedef typename DstView::execution_space execution_space;
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      typedef typename IdxView::non_const_value_type row_index_type;
      typedef typename ColView::non_const_value_type col_index_type;

      Kokkos::pair<int, int> result (1, 1);
      parallel_reduce (range_type (0, idx.size ()),
                       PackArrayMultiColumnVariableStrideWithBoundsCheck (dst, src,
                                                                          idx, col,
                                                                          numCols),
                       result);
      const bool hasBadRows = (result.first != 1);
      const bool hasBadCols = (result.second != 1);
      const bool hasErr = hasBadRows || hasBadCols;
      if (hasErr) {
        std::ostringstream os; // for error reporting

        if (hasBadRows) {
          // Go back and find the out-of-bounds entries in the array of
          // row indices.  Performance doesn't matter since we are already
          // in an error state, so we can do this sequentially, on host.
          auto idx_h = Kokkos::create_mirror_view (idx);
          Kokkos::deep_copy (idx_h, idx);

          std::vector<row_index_type> badRows;
          const size_type numInds = idx_h.extent (0);
          for (size_type k = 0; k < numInds; ++k) {
            if (Impl::outOfBounds<row_index_type> (idx_h(k), src.extent (0))) {
              badRows.push_back (idx_h(k));
            }
          }
          os << "MultiVector multiple-column pack kernel had "
             << badRows.size () << " out-of bounds row index/ices: [";
          for (size_t k = 0; k < badRows.size (); ++k) {
            os << badRows[k];
            if (k + 1 < badRows.size ()) {
              os << ", ";
            }
          }
          os << "].";
        } // hasBadRows

        if (hasBadCols) {
          // Go back and find the out-of-bounds entries in the array
          // of column indices.  Performance doesn't matter since we
          // are already in an error state, so we can do this
          // sequentially, on host.
          auto col_h = Kokkos::create_mirror_view (col);
          Kokkos::deep_copy (col_h, col);

          std::vector<col_index_type> badCols;
          const size_type numInds = col_h.extent (0);
          for (size_type k = 0; k < numInds; ++k) {
            if (Impl::outOfBounds<col_index_type> (col_h(k), src.extent (1))) {
              badCols.push_back (col_h(k));
            }
          }

          if (hasBadRows) {
            os << "  ";
          }
          os << "MultiVector multiple-column pack kernel had "
             << badCols.size () << " out-of bounds column index/ices: [";
          for (size_t k = 0; k < badCols.size (); ++k) {
            os << badCols[k];
            if (k + 1 < badCols.size ()) {
              os << ", ";
            }
          }
          os << "].";
        } // hasBadCols

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

  // FIXME (mfh 16 Dec 2016) It looks like these are totally generic
  // and don't get specialized in Stokhos.  This suggests we can
  // template them on the execution space, and specialize them for
  // Kokkos::Serial so they don't do atomic update operations.  It
  // might be better to make the unpack kernels (the pack kernels
  // don't use atomic updates) functions that we prebuild, in the
  // manner of KokkosKernels, or at least handle via ETI.

  template<class ExecutionSpace>
  struct InsertOp {
    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      Kokkos::atomic_assign(&dest, src);
    }
  };

#ifdef KOKKOS_ENABLE_SERIAL
  template<>
  struct InsertOp< ::Kokkos::Serial > {
    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      dest = src; // no need for an atomic operation here
    }
  };
#endif // KOKKOS_ENABLE_SERIAL

  template<class ExecutionSpace>
  struct AddOp {
    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      Kokkos::atomic_add(&dest, src);
    }
  };

#ifdef KOKKOS_ENABLE_SERIAL
  template<>
  struct AddOp< ::Kokkos::Serial > {
    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      dest += src; // no need for an atomic operation here
    }
  };
#endif // KOKKOS_ENABLE_SERIAL

  template<class Scalar1, class Scalar2>
  struct AbsMaxOper {
    KOKKOS_FORCEINLINE_FUNCTION
    static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
      const auto val1_abs = Kokkos::Details::ArithTraits<Scalar1>::abs(val1);
      const auto val2_abs = Kokkos::Details::ArithTraits<Scalar2>::abs(val2);
      return val1_abs > val2_abs ? Scalar1(val1_abs) : Scalar1(val2_abs);
    }
  };

  template<class ExecutionSpace>
  struct AbsMaxOp {
    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      Kokkos::Impl::atomic_fetch_oper (AbsMaxOper<Scalar,Scalar>(), &dest, src);
    }
  };

#ifdef KOKKOS_ENABLE_SERIAL
  template<>
  struct AbsMaxOp< ::Kokkos::Serial > {
    // ETP:  Is this really what we want?  This seems very odd if
    // Scalar != SCT::mag_type (e.g., Scalar == std::complex<T>)
    //
    // mfh: I didn't write this code, but note that we don't use T =
    // Scalar here, we use T = ArithTraits<Scalar>::mag_type.  That
    // makes this code reasonable.
    template <typename T>
    KOKKOS_INLINE_FUNCTION
    T max(const T& a, const T& b) const { return a > b ? a : b; }

    template <typename Scalar>
    KOKKOS_INLINE_FUNCTION
    void operator() (Scalar& dest, const Scalar& src) const {
      typedef Kokkos::Details::ArithTraits<Scalar> SCT;
      // no need for an atomic operation here
      dest = static_cast<Scalar> (max (SCT::abs (dest), SCT::abs (src)));
    }
  };
#endif // KOKKOS_ENABLE_SERIAL

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

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k) const
    {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (dst(localRow, j), src(offset+j));
      }
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const Op& op,
            const size_t numCols)
    {
      Kokkos::parallel_for
        ("Tpetra::MultiVector unpack (constant stride)",
         Kokkos::RangePolicy<execution_space, size_type> (0, idx.size ()),
         UnpackArrayMultiColumn (execSpace, dst, src, idx, op, numCols));
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
    typedef typename ExecutionSpace::execution_space execution_space;
    typedef SizeType size_type;
    //! We use int as a Boolean (Kokkos doesn't allow bool reduction types).
    typedef int value_type;

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

    KOKKOS_INLINE_FUNCTION
    void operator() (const size_type& k, value_type& result) const
    {
      typedef typename IdxView::non_const_value_type index_type;

      const index_type lclRow = idx(k);
      if (lclRow < static_cast<index_type> (0) ||
          lclRow >= static_cast<index_type> (dst.extent (0))) {
        result = 0; // failed!
      }
      else {
        const size_type offset = k*numCols;
        for (size_type j = 0; j < numCols; ++j) {
          op (dst(lclRow,j), src(offset+j));
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void init (value_type& initialResult) const {
      initialResult = 1; // success
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& dstResult,
          const volatile value_type& srcResult) const
    {
      dstResult = (dstResult == 0 || srcResult == 0) ? 0 : 1;
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const Op& op,
            const size_type numCols)
    {
      typedef typename IdxView::non_const_value_type index_type;

      int result = 1;
      Kokkos::parallel_reduce
        ("Tpetra::MultiVector unpack (constant stride) (with bounds check)",
         Kokkos::RangePolicy<execution_space, size_type> (0, idx.size ()),
         UnpackArrayMultiColumnWithBoundsCheck (execSpace, dst, src,
                                                idx, op, numCols),
         result);

      if (result != 1) {
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
                             const bool debug = true)
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
      impl_type::unpack (execSpace, dst, src, idx, op, numCols);
    }
    else {
      typedef UnpackArrayMultiColumn<ExecutionSpace,
        DstView, SrcView, IdxView, Op> impl_type;
      impl_type::unpack (execSpace, dst, src, idx, op, numCols);
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
    typedef typename ExecutionSpace::execution_space execution_space;
    typedef typename execution_space::size_type size_type;

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
                                          size_t numCols_) :
      dst (dst_),
      src (src_),
      idx (idx_),
      col (col_),
      op (op_),
      numCols (numCols_)
    {}

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type k) const
    {
      const typename IdxView::value_type localRow = idx(k);
      const size_t offset = k*numCols;
      for (size_t j = 0; j < numCols; ++j) {
        op (dst(localRow, col(j)), src(offset+j));
      }
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const ColView& col,
            const Op& op,
            const size_t numCols)
    {
      Kokkos::parallel_for
        ("Tpetra::MultiVector unpack (nonconstant stride)",
         Kokkos::RangePolicy<execution_space, size_type> (0, idx.size ()),
         UnpackArrayMultiColumnVariableStride (execSpace, dst, src,
                                               idx, col, op, numCols));
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
    typedef typename ExecutionSpace::execution_space execution_space;
    typedef SizeType size_type;
    //! We use int as a Boolean (Kokkos doesn't allow bool reduction types).
    typedef Kokkos::pair<int, int> value_type;

  private:
    DstView dst;
    SrcView src;
    IdxView idx;
    ColView col;
    Op op;
    size_type numCols;

  public:
    UnpackArrayMultiColumnVariableStrideWithBoundsCheck (const ExecutionSpace& /* execSpace */,
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

    KOKKOS_INLINE_FUNCTION void
    operator() (const size_type& k, value_type& result) const
    {
      typedef typename IdxView::non_const_value_type row_index_type;
      typedef typename ColView::non_const_value_type col_index_type;

      const row_index_type lclRow = idx(k);
      if (lclRow < static_cast<row_index_type> (0) ||
          lclRow >= static_cast<row_index_type> (dst.extent (0))) {
        result.first = 0; // failed!
      }
      else {
        const size_type offset = k*numCols;
        for (size_type j = 0; j < numCols; ++j) {
          const col_index_type lclCol = col(j);

          if (Impl::outOfBounds<col_index_type> (lclCol, dst.extent (1))) {
            result.second = 0; // failed!
          }
          else { // all indices are valid; apply the op
            op (dst(lclRow, col(j)), src(offset+j));
          }
        }
      }
    }

    KOKKOS_INLINE_FUNCTION void
    init (value_type& initialResult) const {
      initialResult.first = 1; // success
      initialResult.second = 1; // success
    }

    KOKKOS_INLINE_FUNCTION void
    join (volatile value_type& dstResult,
          const volatile value_type& srcResult) const
    {
      dstResult.first = (dstResult.first == 0 || srcResult.first == 0) ? 0 : 1;
      dstResult.second = (dstResult.second == 0 || srcResult.second == 0) ? 0 : 1;
    }

    static void
    unpack (const ExecutionSpace& execSpace,
            const DstView& dst,
            const SrcView& src,
            const IdxView& idx,
            const ColView& col,
            const Op& op,
            const size_type numCols)
    {
      typedef typename IdxView::non_const_value_type row_index_type;
      typedef typename ColView::non_const_value_type col_index_type;

      Kokkos::pair<int, int> result (1, 1);
      Kokkos::parallel_reduce
        ("Tpetra::MultiVector unpack (nonconstant stride) (with bounds check)",
         Kokkos::RangePolicy<execution_space, size_type> (0, idx.size ()),
         UnpackArrayMultiColumnVariableStrideWithBoundsCheck (execSpace, dst,
                                                              src, idx, col,
                                                              op, numCols),
         result);

      const bool hasBadRows = (result.first != 1);
      const bool hasBadCols = (result.second != 1);
      const bool hasErr = hasBadRows || hasBadCols;
      if (hasErr) {
        std::ostringstream os; // for error reporting

        if (hasBadRows) {
          // Go back and find the out-of-bounds entries in the array
          // of row indices.  Performance doesn't matter since we are
          // already in an error state, so we can do this
          // sequentially, on host.
          auto idx_h = Kokkos::create_mirror_view (idx);
          Kokkos::deep_copy (idx_h, idx);

          std::vector<row_index_type> badRows;
          const size_type numInds = idx_h.extent (0);
          for (size_type k = 0; k < numInds; ++k) {
            if (idx_h(k) < static_cast<row_index_type> (0) ||
                idx_h(k) >= static_cast<row_index_type> (dst.extent (0))) {
              badRows.push_back (idx_h(k));
            }
          }
          os << "MultiVector multiple-column unpack kernel had "
             << badRows.size () << " out-of bounds row index/ices: [";
          for (size_t k = 0; k < badRows.size (); ++k) {
            os << badRows[k];
            if (k + 1 < badRows.size ()) {
              os << ", ";
            }
          }
          os << "].";
        } // hasBadRows

        if (hasBadCols) {
          // Go back and find the out-of-bounds entries in the array
          // of column indices.  Performance doesn't matter since we
          // are already in an error state, so we can do this
          // sequentially, on host.
          auto col_h = Kokkos::create_mirror_view (col);
          Kokkos::deep_copy (col_h, col);

          std::vector<col_index_type> badCols;
          const size_type numInds = col_h.extent (0);
          for (size_type k = 0; k < numInds; ++k) {
            if (Impl::outOfBounds<col_index_type> (col_h(k), dst.extent (1))) {
              badCols.push_back (col_h(k));
            }
          }

          if (hasBadRows) {
            os << "  ";
          }
          os << "MultiVector multiple-column unpack kernel had "
             << badCols.size () << " out-of bounds column index/ices: [";
          for (size_t k = 0; k < badCols.size (); ++k) {
            os << badCols[k];
            if (k + 1 < badCols.size ()) {
              os << ", ";
            }
          }
          os << "].";
        } // hasBadCols

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
    static_assert (static_cast<int> (DstView::rank) == 2,
                   "DstView must be a rank-2 Kokkos::View.");
    static_assert (static_cast<int> (SrcView::rank) == 1,
                   "SrcView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (IdxView::rank) == 1,
                   "IdxView must be a rank-1 Kokkos::View.");
    static_assert (static_cast<int> (ColView::rank) == 1,
                   "ColView must be a rank-1 Kokkos::View.");

    if (debug) {
      typedef UnpackArrayMultiColumnVariableStrideWithBoundsCheck<ExecutionSpace,
        DstView, SrcView, IdxView, ColView, Op> impl_type;
      impl_type::unpack (execSpace, dst, src, idx, col, op, numCols);
    }
    else {
      typedef UnpackArrayMultiColumnVariableStride<ExecutionSpace,
        DstView, SrcView, IdxView, ColView, Op> impl_type;
      impl_type::unpack (execSpace, dst, src, idx, col, op, numCols);
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

    PermuteArrayMultiColumn(const DstView& dst_,
                            const SrcView& src_,
                            const DstIdxView& dst_idx_,
                            const SrcIdxView& src_idx_,
                            size_t numCols_) :
      dst(dst_), src(src_), dst_idx(dst_idx_), src_idx(src_idx_),
      numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      for (size_t j = 0; j < numCols; ++j)
        dst(toRow, j) = src(fromRow, j);
    }

    static void permute(const DstView& dst,
                        const SrcView& src,
                        const DstIdxView& dst_idx,
                        const SrcIdxView& src_idx,
                        size_t numCols) {
      const size_type n = std::min( dst_idx.size(), src_idx.size() );
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      Kokkos::parallel_for (range_type (0, n),
                            PermuteArrayMultiColumn (dst,src,dst_idx,src_idx,numCols));
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
                                          size_t numCols_) :
      dst(dst_), src(src_), dst_idx(dst_idx_), src_idx(src_idx_),
      dst_col(dst_col_), src_col(src_col_),
      numCols(numCols_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()( const size_type k ) const {
      const typename DstIdxView::value_type toRow = dst_idx(k);
      const typename SrcIdxView::value_type fromRow = src_idx(k);
      for (size_t j = 0; j < numCols; ++j)
        dst(toRow, dst_col(j)) = src(fromRow, src_col(j));
    }

    static void permute(const DstView& dst,
                        const SrcView& src,
                        const DstIdxView& dst_idx,
                        const SrcIdxView& src_idx,
                        const DstColView& dst_col,
                        const SrcColView& src_col,
                        size_t numCols) {
      const size_type n = std::min( dst_idx.size(), src_idx.size() );
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;
      Kokkos::parallel_for (range_type (0, n),
                            PermuteArrayMultiColumnVariableStride (dst, src,
                                                                   dst_idx,
                                                                   src_idx,
                                                                   dst_col,
                                                                   src_col,
                                                                   numCols));
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
