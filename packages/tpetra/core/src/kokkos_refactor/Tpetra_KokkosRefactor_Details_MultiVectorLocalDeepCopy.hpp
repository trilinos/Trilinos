// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_HPP
#define TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_HPP

#include "Kokkos_Core.hpp"
#include "Tpetra_Details_copyConvert.hpp"

namespace Tpetra {
namespace Details {

namespace { // (anonymous)

template<class DstViewType, class SrcViewType>
void
copyConvertResolvingPossibleAliasing (const DstViewType& dst,
                                      const SrcViewType& src)
{
  // NOTE: It's important to do the addition _inside_ the
  // reinterpret-cast.  If you reinterpret_cast the separate results,
  // you may get the wrong answer (e.g., because ptrdiff_t is signed,
  // and pointers may have arbitrary 64-bit virtual addresses).  I'm
  // speaking from experience here.
  const ptrdiff_t dst_beg = reinterpret_cast<ptrdiff_t> (dst.data ());
  const ptrdiff_t dst_end =
    reinterpret_cast<ptrdiff_t> (dst.data () + dst.span ());
  const ptrdiff_t src_beg = reinterpret_cast<ptrdiff_t> (src.data ());
  const ptrdiff_t src_end =
    reinterpret_cast<ptrdiff_t> (src.data () + src.span ());

  if (src_beg == dst_beg && src_end == dst_end) {
    // Do nothing; there's no need to copy
  }
  else if (dst_end <= src_beg || src_end <= dst_beg) { // no aliasing
    ::Tpetra::Details::copyConvert (dst, src);
  }
  else {
    // dst and src alias each other, so we can't call
    // Kokkos::deep_copy(dst,src) directly (Kokkos detects this and
    // throws, at least in debug mode).  Instead, we make temporary
    // host storage (create_mirror always makes a new allocation,
    // unlike create_mirror_view).  Use host because it's cheaper to
    // allocate.  Hopefully users aren't doing aliased copies in a
    // tight loop.
    auto src_copy = Kokkos::create_mirror (Kokkos::HostSpace (), src);

    // DEEP_COPY REVIEW - NOT TESTED
    Kokkos::deep_copy (src_copy, src);
    ::Tpetra::Details::copyConvert (dst, src_copy);
  }
}

} // namespace (anonymous)

/// \brief Implementation of Tpetra::MultiVector deep copy of local data.
///
/// This implements <tt>Tpetra::MultiVector</tt> deep copy, as in
/// <ul>
/// <li> <tt>Tpetra::deep_copy</tt> </li>
/// <li> <tt>Tpetra::MultiVector::assign</tt> </li>
/// <li> <tt>Tpetra::MultiVector::createCopy</tt> </li>
/// <li> <tt>The two-argument MultiVector copy constructor with
///          <tt>Teuchos::Copy</tt> as the second argument </li>
/// </ul>
///
/// \param dst [in/out] Rank-2 <tt>Kokkos::View</tt>; destination of
///   the copy
///
/// \param src [in] Rank-2 <tt>Kokkos::View</tt>; source of the copy
///
/// \param dstConstStride [in] Whether <tt>dst</tt> is "constant
///   stride."  If so, then the <tt>j</tt>-th column of <tt>dst</tt>
///   has index <tt>j</tt>.  If not, then it has index
///   <tt>dstWhichVecs[j]</tt>.
///
/// \param srcConstStride [in] Whether <tt>src</tt> is "constant
///   stride."  If so, then the <tt>j</tt>-th column of <tt>src</tt>
///   has index <tt>j</tt>.  If not, then it has index
///   <tt>srcWhichVecs[j]</tt>.
///
/// \param dstWhichVecs [in] Host-readable Rank-1 array of some kind,
///   corresponding to <tt>dst.whichVectors_</tt>.  Need only be
///   readable (from host) if <tt>dstConstStride</tt> is true.
///
/// \param srcWhichVecs [in] Host-readable Rank-1 array of some kind,
///   corresponding to <tt>src.whichVectors_</tt>.  Need only be
///   readable (from host) if <tt>srcConstStride</tt> is true.
template<class DstViewType,
         class SrcViewType,
         class DstWhichVecsType,
         class SrcWhichVecsType>
void
localDeepCopy (const DstViewType& dst,
               const SrcViewType& src,
               const bool dstConstStride,
               const bool srcConstStride,
               const DstWhichVecsType& dstWhichVecs,
               const SrcWhichVecsType& srcWhichVecs)
{
  using Kokkos::ALL;
  using Kokkos::subview;
  using size_type = typename DstViewType::size_type;

  if (dstConstStride && srcConstStride) {
    copyConvertResolvingPossibleAliasing (dst, src);
  }
  else {
    const size_type numCols = dstConstStride ?
      static_cast<size_type> (srcWhichVecs.size ()) :
      static_cast<size_type> (dstWhichVecs.size ());
    for (size_type j = 0; j < numCols; ++j) {
      const size_type dst_col = dstConstStride ? j :
        static_cast<size_type> (dstWhichVecs[j]);
      const auto dst_j = subview (dst, ALL (), dst_col);
      const size_type src_col = srcConstStride ? j :
        static_cast<size_type> (srcWhichVecs[j]);
      const auto src_j = subview (src, ALL (), src_col);

      copyConvertResolvingPossibleAliasing (dst_j, src_j);
    }
  }
}

/// \brief Implementation of Tpetra::MultiVector deep copy of local
///   data, for when both the source and destination MultiVector
///   objects have constant stride (isConstantStride() is true).
template<class DstViewType,
         class SrcViewType>
void
localDeepCopyConstStride (const DstViewType& dst,
                          const SrcViewType& src)
{
  return copyConvertResolvingPossibleAliasing (dst, src);
}

} // Details namespace
} // Tpetra namespace

#endif // TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_HPP
