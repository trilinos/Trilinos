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

#ifndef TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_HPP
#define TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_HPP

#include "Kokkos_Core.hpp"
#include "Tpetra_Details_copyConvert.hpp"

namespace Tpetra {
namespace Details {

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
    ::Tpetra::Details::copyConvert (dst, src);
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

      ::Tpetra::Details::copyConvert (dst_j, src_j);
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
  return ::Tpetra::Details::copyConvert (dst, src);
}

} // Details namespace
} // Tpetra namespace

#endif // TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_HPP
