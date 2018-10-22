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
#include "Kokkos_ArithTraits.hpp"

namespace Tpetra {
namespace Details {

  // Functor for deep-copying between two 2-D Kokkos Views.
  // This implements Tpetra::MultiVector deep copy, as in
  //   - Tpetra::deep_copy
  //   - Tpetra::MultiVector::assign
  //   - Tpetra::MultiVector::createCopy
  //   - The two-argument MultiVector copy constructor with
  //     Teuchos::Copy as the second argument
  //
  // DstViewType and SrcViewType must be 2-D Kokkos Views.
  // DstWhichVecsType and SrcWhichVecsType must be 1-D Kokkos Views
  // whose value type is an integer.  They correspond to the "which
  // vectors?" arrays for the destination and source Views,
  // respectively.
  //
  // If DstConstStride is true, dstWhichVecs is ignored.  If
  // SrcConstStride is true, srcWhichVecs is ignored.  These bool
  // template parameters take advantage of compilers' ability to
  // disable code inside an "if (false) { ... }" scope.
  template<class DstViewType,
           class SrcViewType,
           class DstWhichVecsType,
           class SrcWhichVecsType,
           const bool DstConstStride,
           const bool SrcConstStride>
  struct LocalDeepCopyFunctor {
    typedef typename DstViewType::execution_space execution_space;
    typedef typename DstViewType::size_type index_type;

    DstViewType dst_;
    SrcViewType src_;
    DstWhichVecsType dstWhichVecs_;
    SrcWhichVecsType srcWhichVecs_;
    const index_type numVecs_;

    // You may always call the 4-argument constructor.
    // dstWhichVecs is ignored if DstConstStride is true.
    // srcWhichVecs is ignored if SrcConstStride is true.
    LocalDeepCopyFunctor (const DstViewType& dst,
                          const SrcViewType& src,
                          const DstWhichVecsType& dstWhichVecs,
                          const SrcWhichVecsType& srcWhichVecs) :
      dst_ (dst),
      src_ (src),
      dstWhichVecs_ (dstWhichVecs),
      srcWhichVecs_ (srcWhichVecs),
      numVecs_ (DstConstStride ? dst.extent (1) : dstWhichVecs.extent (0))
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! DstConstStride && ! SrcConstStride &&
        dstWhichVecs.extent (0) != srcWhichVecs.extent (0),
        std::invalid_argument, "LocalDeepCopyFunctor (4-arg constructor): "
        "Neither src nor dst have constant stride, but "
        "dstWhichVecs.extent(0) = " << dstWhichVecs.extent (0)
        << " != srcWhichVecs.extent(0) = " << srcWhichVecs.extent (0)
        << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        DstConstStride && ! SrcConstStride &&
        srcWhichVecs.extent (0) != dst.extent (1),
        std::invalid_argument, "LocalDeepCopyFunctor (4-arg constructor): "
        "src does not have constant stride, but srcWhichVecs.extent(0) = "
        << srcWhichVecs.extent (0) << " != dst.extent(1) = "
        << dst.extent (1) << ".");
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! DstConstStride && SrcConstStride &&
        dstWhichVecs.extent (0) != src.extent (1),
        std::invalid_argument, "LocalDeepCopyFunctor (4-arg constructor): "
        "dst does not have constant stride, but dstWhichVecs.extent(0) = "
        << dstWhichVecs.extent (0) << " != src.extent(1) = "
        << src.extent (1) << ".");
    }

    // You may only call the 2-argument constructor if DstConstStride
    // and SrcConstStride are both true.
    LocalDeepCopyFunctor (const DstViewType& dst, const SrcViewType& src) :
      dst_ (dst),
      src_ (src),
      numVecs_ (dst.extent (1))
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! DstConstStride || ! SrcConstStride, std::logic_error,
        "Tpetra::LocalDeepCopyFunctor: You may not use the constant-stride "
        "constructor if either of the Boolean template parameters is false.");
    }

    void KOKKOS_INLINE_FUNCTION operator () (const index_type i) const {
      if (DstConstStride) {
        if (SrcConstStride) {
          for (index_type j = 0; j < numVecs_; ++j) {
            dst_(i,j) = src_(i,j);
          }
        } else {
          for (index_type j = 0; j < numVecs_; ++j) {
            dst_(i,j) = src_(i,srcWhichVecs_(j));
          }
        }
      } else {
        if (SrcConstStride) {
          for (index_type j = 0; j < numVecs_; ++j) {
            dst_(i,dstWhichVecs_(j)) = src_(i,j);
          }
        } else {
          for (index_type j = 0; j < numVecs_; ++j) {
            dst_(i,dstWhichVecs_(j)) = src_(i,srcWhichVecs_(j));
          }
        }
      }
    }
  };

  template<class DstViewType,
           class SrcViewType = DstViewType,
           class DstWhichVecsType = Kokkos::View<const typename DstViewType::size_type*, typename DstViewType::execution_space>,
           class SrcWhichVecsType = Kokkos::View<const typename SrcViewType::size_type*, typename SrcViewType::execution_space> >
  struct LocalDeepCopy {
    typedef typename DstViewType::size_type index_type;

    static void
    run (const DstViewType& dst, const SrcViewType& src,
         const bool dstConstStride, const bool srcConstStride)
    {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! dstConstStride || ! srcConstStride, std::invalid_argument,
        "LocalDeepCopy::run: You may only call the 4-argument version of this "
        "function if dstConstStride and srcConstStride are both true.");

      // FIXME (mfh 22 Jul 2014, 10 Dec 2014) Currently, it doesn't
      // work to do a 2-D copy, even if both MultiVectors have
      // constant stride.  This is because Kokkos can't currently tell
      // the difference between padding (which permits a single
      // deep_copy for the whole 2-D View) and stride > numRows (which
      // does NOT permit a single deep_copy for the whole 2-D View).
      // Carter is working on this, but for now, the temporary fix is
      // to copy one column at a time.
      //
      // FIXME (mfh 10 Dec 2014) Call Kokkos::deep_copy if appropriate
      // for dst and src.  See note above.
      typedef LocalDeepCopyFunctor<DstViewType, SrcViewType,
        DstWhichVecsType, SrcWhichVecsType, true, true> functor_type;
      functor_type f (dst, src);
      typedef typename DstViewType::execution_space execution_space;
      typedef decltype (dst.extent (0)) size_type;
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;

      Kokkos::parallel_for ("Tpetra::Details::LocalDeepCopy(x,y,stride)",range_type (0, dst.extent (0)), f);
    }

    static void
    run (const DstViewType& dst, const SrcViewType& src,
         const bool dstConstStride, const bool srcConstStride,
         const DstWhichVecsType& dstWhichVecs,
         const SrcWhichVecsType& srcWhichVecs)
    {
      // FIXME (mfh 22 Jul 2014, 10 Dec 2014) Currently, it doesn't
      // work to do a 2-D copy, even if both MultiVectors have
      // constant stride.  This is because Kokkos can't currently tell
      // the difference between padding (which permits a single
      // deep_copy for the whole 2-D View) and stride > numRows (which
      // does NOT permit a single deep_copy for the whole 2-D View).
      // Carter is working on this, but for now, the temporary fix is
      // to copy one column at a time.

      typedef typename DstViewType::execution_space execution_space;
      typedef decltype (dst.extent (0)) size_type;
      typedef Kokkos::RangePolicy<execution_space, size_type> range_type;

      if (dstConstStride) {
        if (srcConstStride) {
          // FIXME (mfh 10 Dec 2014) Do a Kokkos::deep_copy if
          // appropriate for dst and src.  See note above.
          typedef LocalDeepCopyFunctor<DstViewType, SrcViewType,
            DstWhichVecsType, SrcWhichVecsType, true, true> functor_type;
          functor_type f (dst, src);
          Kokkos::parallel_for ("Tpetra::Details::LocalDeepCopy(x,y,stride,whichVecs,0)",range_type (0, dst.extent (0)), f);
        }
        else { // ! srcConstStride
          typedef LocalDeepCopyFunctor<DstViewType, SrcViewType,
            DstWhichVecsType, SrcWhichVecsType, true, false> functor_type;
          functor_type f (dst, src, srcWhichVecs, srcWhichVecs);
          Kokkos::parallel_for ("Tpetra::Details::LocalDeepCopy(x,y,stride,whichVecs,1)",range_type (0, dst.extent (0)), f);
        }
      }
      else { // ! dstConstStride
        if (srcConstStride) {
          typedef LocalDeepCopyFunctor<DstViewType, SrcViewType,
            DstWhichVecsType, SrcWhichVecsType, false, true> functor_type;
          functor_type f (dst, src, dstWhichVecs, dstWhichVecs);
          Kokkos::parallel_for ("Tpetra::Details::LocalDeepCopy(x,y,stride,whichVecs,2)",range_type (0, dst.extent (0)), f);
        }
        else { // ! srcConstStride
          typedef LocalDeepCopyFunctor<DstViewType, SrcViewType,
            DstWhichVecsType, SrcWhichVecsType, false, false> functor_type;
          functor_type f (dst, src, dstWhichVecs, srcWhichVecs);
          Kokkos::parallel_for ("Tpetra::Details::LocalDeepCopy(x,y,stride,whichVecs,3)",range_type (0, dst.extent (0)), f);
        }
      }
    }
  };


  template<class DstViewType,
           class SrcViewType,
           class DstWhichVecsType,
           class SrcWhichVecsType>
  void
  localDeepCopy (const DstViewType& dst, const SrcViewType& src,
                 const bool dstConstStride, const bool srcConstStride,
                 const DstWhichVecsType& dstWhichVecs,
                 const SrcWhichVecsType& srcWhichVecs)
  {
    typedef LocalDeepCopy<DstViewType, SrcViewType,
      DstWhichVecsType, SrcWhichVecsType> impl_type;
    impl_type::run (dst, src, dstConstStride, srcConstStride,
                    dstWhichVecs, srcWhichVecs);
  }

  template<class DstViewType, class SrcViewType>
  void
  localDeepCopyConstStride (const DstViewType& dst, const SrcViewType& src)
  {
    typedef LocalDeepCopy<DstViewType, SrcViewType> impl_type;
    impl_type::run (dst, src, true, true);
  }

} // Details namespace
} // Tpetra namespace

#endif // TPETRA_KOKKOS_REFACTOR_DETAILS_MULTI_VECTOR_LOCAL_DEEP_COPY_HPP
