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

#ifndef TPETRA_EXPERIMENTAL_BLOCKVIEW_HPP
#define TPETRA_EXPERIMENTAL_BLOCKVIEW_HPP

/// \file Tpetra_Experimental_BlockView.hpp
/// \brief Linear algebra kernels for small dense matrices and vectors
///
/// This file declares and defines generic computational kernels for
/// small dense linear algebra operations, with matrices and vectors
/// stored as Kokkos::View.  The operations are meant as helpers for
/// Tpetra::Experimental::BlockCrsMatrix.

#include "TpetraCore_config.h"
#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_Complex.hpp"

namespace Tpetra {

/// \brief Namespace for new Tpetra features that are not ready for
///   public release, but are ready for evaluation by friendly expert
///   users.
///
/// \warning Expect header files, classes, functions, and other
///   interfaces to change or disappear.  Anything in this namespace
///   is under active development and evaluation.  Documentation may
///   be sparse or not exist yet.  Generally, unit tests will exist,
///   but coverage may be lacking.  If you understand these caveats
///   and accept them, please feel free to take a look inside and try
///   things out.
namespace Experimental {

namespace Impl {

/// \brief Implementation of Tpetra's ABSMAX CombineMode for the small
///   dense blocks in BlockCrsMatrix, or the small dense vectors in
///   BlockMultiVector and BlockVector.
///
/// This is the "generic" version that we don't implement.
/// We actually implement versions for ViewType rank 1 or rank 2.
template<class ViewType1,
         class ViewType2,
         const int rank1 = ViewType1::rank>
struct AbsMax {
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType2& Y, const ViewType1& X);
};

/// \brief Implementation of Tpetra's ABSMAX CombineMode for the small
///   dense blocks in BlockCrsMatrix.
///
/// Tpetra uses this operation to implement the ABSMAX CombineMode.
template<class ViewType1,
         class ViewType2>
struct AbsMax<ViewType1, ViewType2, 2> {
  /// \brief <tt>(*this)(i,j) := max(abs((*this)(i,j)), abs(X(i,j)))</tt>
  ///   for all (i,j).
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType2& Y, const ViewType1& X)
  {
    static_assert (static_cast<int> (ViewType1::rank) == static_cast<int> (ViewType2::rank),
                   "AbsMax: ViewType1 and ViewType2 must have the same rank.");
    typedef typename std::remove_reference<decltype (Y(0,0)) >::type STY;
    static_assert(! std::is_const<STY>::value,
      "AbsMax: The type of each entry of Y must be nonconst.");
    typedef typename std::decay<decltype (X(0,0)) >::type STX;
    static_assert(  std::is_same<STX, STY>::value,
      "AbsMax: The type of each entry of X and Y must be the same.");
    typedef Kokkos::Details::ArithTraits<STY> KAT;

    const int numCols = Y.extent (1);
    const int numRows = Y.extent (0);
    for (int j = 0; j < numCols; ++j) {
      for (int i = 0; i < numRows; ++i) {
        STY& Y_ij = Y(i,j); // use ref here to avoid 2nd op() call on Y
        const STX X_ij = X(i,j);
        // NOTE: no std::max (not a CUDA __device__ function); must
        // cast back up to complex.
        const auto Y_ij_abs = KAT::abs (Y_ij);
        const auto X_ij_abs = KAT::abs (X_ij);
        Y_ij = (Y_ij_abs >= X_ij_abs) ?
          static_cast<STY> (Y_ij_abs) :
          static_cast<STY> (X_ij_abs);
      }
    }
  }
};

/// \brief Implementation of Tpetra's ABSMAX CombineMode for the small
///   dense vectors in BlockMultiVector and BlockVector.
///
/// Tpetra uses this operation to implement the ABSMAX CombineMode.
template<class ViewType1,
         class ViewType2>
struct AbsMax<ViewType1, ViewType2, 1> {
  /// \brief <tt>(*this)(i) := max(abs((*this)(i)), abs(X(i)))</tt>
  ///   for all i.
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType2& Y, const ViewType1& X)
  {
    static_assert (static_cast<int> (ViewType1::rank) == static_cast<int> (ViewType2::rank),
                   "AbsMax: ViewType1 and ViewType2 must have the same rank.");

    typedef typename std::remove_reference<decltype (Y(0)) >::type STY;
    static_assert(! std::is_const<STY>::value,
      "AbsMax: The type of each entry of Y must be nonconst.");
    typedef typename std::remove_const<typename std::remove_reference<decltype (X(0)) >::type>::type STX;
    static_assert(  std::is_same<STX, STY>::value,
      "AbsMax: The type of each entry of X and Y must be the same.");
    typedef Kokkos::Details::ArithTraits<STY> KAT;

    const int numRows = Y.extent (0);
    for (int i = 0; i < numRows; ++i) {
      STY& Y_i = Y(i); // use ref here to avoid 2nd op() call on Y
      const STX X_i = X(i);
      // NOTE: no std::max (not a CUDA __device__ function); must
      // cast back up to complex.
      const auto Y_i_abs = KAT::abs (Y_i);
      const auto X_i_abs = KAT::abs (X_i);
      Y_i = (Y_i_abs >= X_i_abs) ?
        static_cast<STY> (Y_i_abs) :
        static_cast<STY> (X_i_abs);
    }
  }
};

/// \brief Implementation of Tpetra's ABSMAX CombineMode for the small
///   dense blocks in BlockCrsMatrix, and the small dense vectors in
///   BlockMultiVector and BlockVector.
///
/// This is the function that Tpetra actually uses to implement the
/// ABSMAX CombineMode.
template<class ViewType1, class ViewType2, const int rank = ViewType1::rank>
KOKKOS_INLINE_FUNCTION void
absMax (const ViewType2& Y, const ViewType1& X)
{
  static_assert (static_cast<int> (ViewType1::rank) == static_cast<int> (ViewType2::rank),
                 "absMax: ViewType1 and ViewType2 must have the same rank.");
  AbsMax<ViewType1, ViewType2, rank>::run (Y, X);
}

/// \brief Implementation of Tpetra::Experimental::SCAL function.
///
/// This is the "generic" version that we don't implement.
/// We actually implement versions for ViewType rank 1 or rank 2.
template<class ViewType,
         class CoefficientType,
         class LayoutType = typename ViewType::array_layout,
         class IndexType = int,
         const int rank = ViewType::rank>
struct SCAL {
  static KOKKOS_INLINE_FUNCTION void
  run (const CoefficientType& alpha, const ViewType& x);
};

/// \brief Implementation of Tpetra::Experimental::SCAL function, for
///   ViewType rank 1 (i.e., a vector).
template<class ViewType,
         class CoefficientType,
         class LayoutType,
         class IndexType>
struct SCAL<ViewType, CoefficientType, LayoutType, IndexType, 1> {
  /// \brief x := alpha*x (rank-1 x, i.e., a vector)
  static KOKKOS_INLINE_FUNCTION void
  run (const CoefficientType& alpha, const ViewType& x)
  {
    const IndexType numRows = static_cast<IndexType> (x.extent (0));
    // BLAS _SCAL doesn't check whether alpha is 0.
    for (IndexType i = 0; i < numRows; ++i) {
      x(i) = alpha * x(i);
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::SCAL function, for
///   ViewType rank 2 (i.e., a matrix).
template<class ViewType,
         class CoefficientType,
         class LayoutType,
         class IndexType>
struct SCAL<ViewType, CoefficientType, LayoutType, IndexType, 2> {
  /// \brief A := alpha*A (rank-2 A, i.e., a matrix)
  static KOKKOS_INLINE_FUNCTION void
  run (const CoefficientType& alpha, const ViewType& A)
  {
    const IndexType numRows = static_cast<IndexType> (A.extent (0));
    const IndexType numCols = static_cast<IndexType> (A.extent (1));

    // BLAS _SCAL doesn't check whether alpha is 0.
    for (IndexType i = 0; i < numRows; ++i) {
      for (IndexType j = 0; j < numCols; ++j) {
        A(i,j) = alpha * A(i,j);
      }
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::SCAL function, for
///   ViewType rank 2 (i.e., a matrix), and LayoutType = LayoutRight.
///
/// For LayoutRight (or LayoutLeft) input, we can flatten indexing
/// from 2-D to 1-D.
template<class ViewType,
         class CoefficientType,
         class IndexType>
struct SCAL<ViewType, CoefficientType, Kokkos::LayoutRight, IndexType, 2> {
  /// \brief A := alpha*A (rank-2 A, i.e., a matrix)
  static KOKKOS_INLINE_FUNCTION void
  run (const CoefficientType& alpha, const ViewType& A)
  {
    const IndexType N = A.size ();
    typedef typename std::decay<decltype (A(0,0)) >::type scalar_type;
    scalar_type* const A_raw = A.data ();

    for (IndexType i = 0; i < N; ++i) {
      A_raw[i] = alpha * A_raw[i];
    }
  }
};


/// \brief Implementation of Tpetra::Experimental::FILL function.
///
/// This is the "generic" version that we don't implement.
/// We actually implement versions for ViewType rank 1 or rank 2.
template<class ViewType,
         class InputType,
         class LayoutType = typename ViewType::array_layout,
         class IndexType = int,
         const int rank = ViewType::rank>
struct FILL {
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType& x, const InputType& val);
};

/// \brief Implementation of Tpetra::Experimental::FILL function, for
///   ViewType rank 1 (i.e., a vector).
template<class ViewType,
         class InputType,
         class LayoutType,
         class IndexType>
struct FILL<ViewType, InputType, LayoutType, IndexType, 1> {
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType& x, const InputType& val)
  {
    const IndexType numRows = static_cast<IndexType> (x.extent (0));
    for (IndexType i = 0; i < numRows; ++i) {
      x(i) = val;
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::FILL function, for
///   ViewType rank 2 (i.e., a matrix).
template<class ViewType,
         class InputType,
         class LayoutType,
         class IndexType>
struct FILL<ViewType, InputType, LayoutType, IndexType, 2> {
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType& X, const InputType& val)
  {
    const IndexType numRows = static_cast<IndexType> (X.extent (0));
    const IndexType numCols = static_cast<IndexType> (X.extent (1));
    for (IndexType j = 0; j < numCols; ++j) {
      for (IndexType i = 0; i < numRows; ++i) {
        X(i,j) = val;
      }
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::AXPY function.
///
/// This is the "generic" version that we don't implement.
/// We actually implement versions for ViewType rank 1 or rank 2.
template<class CoefficientType,
         class ViewType1,
         class ViewType2,
         class LayoutType1 = typename ViewType1::array_layout,
         class LayoutType2 = typename ViewType2::array_layout,
         class IndexType = int,
         const int rank = ViewType1::rank>
struct AXPY {
  static KOKKOS_INLINE_FUNCTION void
  run (const CoefficientType& alpha,
       const ViewType1& x,
       const ViewType2& y);
};

/// \brief Implementation of Tpetra::Experimental::AXPY function, for
///   ViewType1 and ViewType2 rank 1 (i.e., vectors).
template<class CoefficientType,
         class ViewType1,
         class ViewType2,
         class LayoutType1,
         class LayoutType2,
         class IndexType>
struct AXPY<CoefficientType, ViewType1, ViewType2, LayoutType1, LayoutType2, IndexType, 1> {
  /// \brief y := y + alpha*x (rank-1 x and y, i.e., vectors)
  static KOKKOS_INLINE_FUNCTION void
  run (const CoefficientType& alpha,
       const ViewType1& x,
       const ViewType2& y)
  {
    using Kokkos::Details::ArithTraits;
    static_assert (static_cast<int> (ViewType1::rank) == static_cast<int> (ViewType2::rank),
                   "AXPY: x and y must have the same rank.");

    const IndexType numRows = static_cast<IndexType> (y.extent (0));
    if (alpha != ArithTraits<CoefficientType>::zero ()) {
      for (IndexType i = 0; i < numRows; ++i) {
        y(i) += alpha * x(i);
      }
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::AXPY function, for
///   ViewType1 and ViewType2 rank 2 (i.e., matrices).
template<class CoefficientType,
         class ViewType1,
         class ViewType2,
         class LayoutType1,
         class LayoutType2,
         class IndexType>
struct AXPY<CoefficientType, ViewType1, ViewType2, LayoutType1, LayoutType2, IndexType, 2> {
  /// \brief Y := Y + alpha*X (rank-2 X and Y, i.e., matrices)
  static KOKKOS_INLINE_FUNCTION void
  run (const CoefficientType& alpha,
       const ViewType1& X,
       const ViewType2& Y)
  {
    using Kokkos::Details::ArithTraits;
    static_assert (ViewType1::rank == ViewType2::rank,
                   "AXPY: X and Y must have the same rank.");
    const IndexType numRows = static_cast<IndexType> (Y.extent (0));
    const IndexType numCols = static_cast<IndexType> (Y.extent (1));

    if (alpha != ArithTraits<CoefficientType>::zero ()) {
      for (IndexType i = 0; i < numRows; ++i) {
        for (IndexType j = 0; j < numCols; ++j) {
          Y(i,j) += alpha * X(i,j);
        }
      }
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::AXPY function, for
///   ViewType1 and ViewType2 rank 2 (i.e., matrices), when both
///   ViewType1 and ViewType2 have LayoutRight.
template<class CoefficientType,
         class ViewType1,
         class ViewType2,
         class IndexType>
struct AXPY<CoefficientType, ViewType1, ViewType2, Kokkos::LayoutRight, Kokkos::LayoutRight, IndexType, 2> {
  /// \brief Y := Y + alpha*X (rank-2 X and Y, i.e., matrices)
  static KOKKOS_INLINE_FUNCTION void
  run (const CoefficientType& alpha,
       const ViewType1& X,
       const ViewType2& Y)
  {
    using Kokkos::Details::ArithTraits;
    static_assert (static_cast<int> (ViewType1::rank) ==
                   static_cast<int> (ViewType2::rank),
                   "AXPY: X and Y must have the same rank.");
    typedef typename std::decay<decltype (X(0,0)) >::type SX;
    typedef typename std::decay<decltype (Y(0,0)) >::type SY;

    const IndexType N = static_cast<IndexType> (Y.size ());
    const SX* const X_raw = X.data ();
    SY* const Y_raw = Y.data ();

    if (alpha != ArithTraits<CoefficientType>::zero ()) {
      for (IndexType i = 0; i < N; ++i) {
        Y_raw[i] += alpha * X_raw[i];
      }
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::AXPY function, for
///   ViewType1 and ViewType2 rank 2 (i.e., matrices), when both
///   ViewType1 and ViewType2 have LayoutLeft.
template<class CoefficientType,
         class ViewType1,
         class ViewType2,
         class IndexType>
struct AXPY<CoefficientType, ViewType1, ViewType2, Kokkos::LayoutLeft, Kokkos::LayoutLeft, IndexType, 2> {
  /// \brief Y := Y + alpha*X (rank-2 X and Y, i.e., matrices)
  static KOKKOS_INLINE_FUNCTION void
  run (const CoefficientType& alpha,
       const ViewType1& X,
       const ViewType2& Y)
  {
    using Kokkos::Details::ArithTraits;
    static_assert (ViewType1::rank == ViewType2::rank,
                   "AXPY: X and Y must have the same rank.");
    typedef typename std::decay<decltype (X(0,0)) >::type SX;
    typedef typename std::decay<decltype (Y(0,0)) >::type SY;

    const IndexType N = static_cast<IndexType> (Y.size ());
    const SX* const X_raw = X.data ();
    SY* const Y_raw = Y.data ();

    if (alpha != ArithTraits<CoefficientType>::zero ()) {
      for (IndexType i = 0; i < N; ++i) {
        Y_raw[i] += alpha * X_raw[i];
      }
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::COPY function.
///
/// This is the "generic" version that we don't implement.
/// We actually implement versions for ViewType rank 1 or rank 2.
template<class ViewType1,
         class ViewType2,
         class LayoutType1 = typename ViewType1::array_layout,
         class LayoutType2 = typename ViewType2::array_layout,
         class IndexType = int,
         const int rank = ViewType1::rank>
struct COPY {
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType1& x, const ViewType2& y);
};

/// \brief Implementation of Tpetra::Experimental::COPY function, for
///   ViewType1 and ViewType2 rank 1 (i.e., vectors).
template<class ViewType1,
         class ViewType2,
         class LayoutType1,
         class LayoutType2,
         class IndexType>
struct COPY<ViewType1, ViewType2, LayoutType1, LayoutType2, IndexType, 1> {
  /// \brief y := x (rank-1 x and y, i.e., vectors)
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType1& x, const ViewType2& y)
  {
    const IndexType numRows = static_cast<IndexType> (x.extent (0));
    for (IndexType i = 0; i < numRows; ++i) {
      y(i) = x(i);
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::COPY function, for
///   ViewType1 and ViewType2 rank 2 (i.e., matrices).
template<class ViewType1,
         class ViewType2,
         class LayoutType1,
         class LayoutType2,
         class IndexType>
struct COPY<ViewType1, ViewType2, LayoutType1, LayoutType2, IndexType, 2> {
  /// \brief Y := X (rank-2 X and Y, i.e., matrices)
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType1& X, const ViewType2& Y)
  {
    const IndexType numRows = static_cast<IndexType> (Y.extent (0));
    const IndexType numCols = static_cast<IndexType> (Y.extent (1));

    // BLAS _SCAL doesn't check whether alpha is 0.
    for (IndexType i = 0; i < numRows; ++i) {
      for (IndexType j = 0; j < numCols; ++j) {
        Y(i,j) = X(i,j);
      }
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::COPY function, for
///   ViewType1 and ViewType2 rank 2 (i.e., matrices), where both have
///   LayoutRight (row-major order, with contiguous storage).
template<class ViewType1,
         class ViewType2,
         class IndexType>
struct COPY<ViewType1, ViewType2, Kokkos::LayoutRight, Kokkos::LayoutRight, IndexType, 2> {
  /// \brief Y := X (rank-2 X and Y, i.e., matrices)
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType1& X, const ViewType2& Y)
  {
    typedef typename std::decay<decltype (X(0,0)) >::type SX;
    typedef typename std::decay<decltype (Y(0,0)) >::type SY;

    const IndexType N = static_cast<IndexType> (Y.size ());
    const SX* const X_raw = X.data ();
    SY* const Y_raw = Y.data ();

    // BLAS _SCAL doesn't check whether alpha is 0.
    for (IndexType i = 0; i < N; ++i) {
      Y_raw[i] = X_raw[i];
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::COPY function, for
///   ViewType1 and ViewType2 rank 2 (i.e., matrices), where both have
///   LayoutLeft (row-major order, with contiguous storage).
template<class ViewType1,
         class ViewType2,
         class IndexType>
struct COPY<ViewType1, ViewType2, Kokkos::LayoutLeft, Kokkos::LayoutLeft, IndexType, 2> {
  /// \brief Y := X (rank-2 X and Y, i.e., matrices)
  static KOKKOS_INLINE_FUNCTION void
  run (const ViewType1& X, const ViewType2& Y)
  {
    typedef typename std::decay<decltype (X(0,0)) >::type SX;
    typedef typename std::decay<decltype (Y(0,0)) >::type SY;

    const IndexType N = static_cast<IndexType> (Y.size ());
    const SX* const X_raw = X.data ();
    SY* const Y_raw = Y.data ();

    // BLAS _SCAL doesn't check whether alpha is 0.
    for (IndexType i = 0; i < N; ++i) {
      Y_raw[i] = X_raw[i];
    }
  }
};


template<class VecType1,
         class BlkType,
         class VecType2,
         class CoeffType,
         class IndexType = int,
         class VecLayoutType1 = typename VecType1::array_layout,
         class BlkLayoutType = typename BlkType::array_layout,
         class VecLayoutType2 = typename VecType2::array_layout>
struct GEMV {
  static KOKKOS_INLINE_FUNCTION void
  run (const CoeffType& alpha,
       const BlkType& A,
       const VecType1& x,
       const VecType2& y)
  {
    static_assert (VecType1::rank == 1, "GEMV: VecType1 must have rank 1.");
    static_assert (BlkType::rank == 2, "GEMV: BlkType must have rank 2.");
    static_assert (VecType2::rank == 1, "GEMV: VecType2 must have rank 1.");
    typedef typename std::decay<decltype (y(0)) >::type y_elt_type;

    const IndexType numRows = static_cast<IndexType> (A.extent (0));
    const IndexType numCols = static_cast<IndexType> (A.extent (1));

    for (IndexType i = 0; i < numRows; ++i) {
      y_elt_type y_i = y(i);
      for (IndexType j = 0; j < numCols; ++j) {
        y_i += alpha * A(i,j) * x(j);
      }
      y(i) = y_i;
    }
  }
};

template<class VecType1,
         class BlkType,
         class VecType2,
         class CoeffType,
         class IndexType>
struct GEMV<VecType1, BlkType, VecType2, CoeffType, IndexType,
            Kokkos::LayoutRight, Kokkos::LayoutRight, Kokkos::LayoutRight>
{
  static KOKKOS_INLINE_FUNCTION void
  run (const CoeffType& alpha,
       const BlkType& A,
       const VecType1& x,
       const VecType2& y)
  {
    static_assert (VecType1::rank == 1, "GEMV: VecType1 must have rank 1.");
    static_assert (BlkType::rank == 2, "GEMV: BlkType must have rank 2.");
    static_assert (VecType2::rank == 1, "GEMV: VecType2 must have rank 1.");
    typedef typename std::decay<decltype (y(0)) >::type y_elt_type;
    typedef typename std::decay<decltype (A(0,0)) >::type A_elt_type;

    const IndexType numRows = static_cast<IndexType> (A.extent (0));
    const IndexType numCols = static_cast<IndexType> (A.extent (1));
    const A_elt_type* const A_raw = A.data ();

    for (IndexType i = 0; i < numRows; ++i) {
      y_elt_type y_i = y(i);
      const A_elt_type* const A_i = A_raw + i*numCols;
      for (IndexType j = 0; j < numCols; ++j) {
        y_i += alpha * A_i[j] * x(j);
      }
      y(i) = y_i;
    }
  }
};

template<class VecType1,
         class BlkType,
         class VecType2,
         class CoeffType,
         class IndexType>
struct GEMV<VecType1, BlkType, VecType2, CoeffType, IndexType,
            Kokkos::LayoutLeft, Kokkos::LayoutLeft, Kokkos::LayoutLeft>
{
  static KOKKOS_INLINE_FUNCTION void
  run (const CoeffType& alpha,
       const BlkType& A,
       const VecType1& x,
       const VecType2& y)
  {
    static_assert (VecType1::rank == 1, "GEMV: VecType1 must have rank 1.");
    static_assert (BlkType::rank == 2, "GEMV: BlkType must have rank 2.");
    static_assert (VecType2::rank == 1, "GEMV: VecType2 must have rank 1.");
    typedef typename std::decay<decltype (A(0,0)) >::type A_elt_type;

    const A_elt_type* const A_raw = A.data ();
    const IndexType numRows = static_cast<IndexType> (A.extent (0));
    const IndexType numCols = static_cast<IndexType> (A.extent (1));
    for (IndexType j = 0; j < numCols; ++j) {
      const A_elt_type* const A_j = A_raw + j*numRows;
      for (IndexType i = 0; i < numRows; ++i) {
        y(i) += alpha * A_j[i] * x(i);
      }
    }
  }
};

} // namespace Impl

/// \brief x := alpha*x, where x is either rank 1 (a vector) or rank 2
///   (a matrix).
template<class ViewType,
         class CoefficientType,
         class LayoutType = typename ViewType::array_layout,
         class IndexType = int,
         const int rank = ViewType::rank>
KOKKOS_INLINE_FUNCTION void
SCAL (const CoefficientType& alpha, const ViewType& x)
{
  Impl::SCAL<ViewType, CoefficientType, LayoutType, IndexType, rank>::run (alpha, x);
}

/// \brief Set every entry of x to val.
template<class ViewType,
         class InputType,
         class LayoutType = typename ViewType::array_layout,
         class IndexType = int,
         const int rank = ViewType::rank>
KOKKOS_INLINE_FUNCTION void
FILL (const ViewType& x, const InputType& val)
{
  Impl::FILL<ViewType, InputType, LayoutType, IndexType, rank>::run (x, val);
}

/// \brief <tt>y := y + alpha * x</tt> (dense vector or matrix update)
///
/// This function follows the BLAS convention that if alpha == 0, then
/// it does nothing.  (This matters only if x contains Inf or NaN
/// values.)
template<class CoefficientType,
         class ViewType1,
         class ViewType2,
         class LayoutType1 = typename ViewType1::array_layout,
         class LayoutType2 = typename ViewType2::array_layout,
         class IndexType = int,
         const int rank = ViewType1::rank>
KOKKOS_INLINE_FUNCTION void
AXPY (const CoefficientType& alpha,
      const ViewType1& x,
      const ViewType2& y)
{
  static_assert (static_cast<int> (ViewType1::rank) ==
                 static_cast<int> (ViewType2::rank),
                 "AXPY: x and y must have the same rank.");
  Impl::AXPY<CoefficientType, ViewType1, ViewType2, LayoutType1, LayoutType2,
    IndexType, rank>::run (alpha, x, y);
}

/// \brief Deep copy x into y, where x and y are either rank 1
///   (vectors) or rank 2 (matrices) with the same dimension(s).
///
/// \param x [in] The input vector / matrix.
/// \param y [out] The output vector / matrix.
///
/// We put the output argument last, because that's what the BLAS
/// functions _COPY (replace _ with "S", "D", "C", or "Z") do.
template<class ViewType1,
         class ViewType2,
         class LayoutType1 = typename ViewType1::array_layout,
         class LayoutType2 = typename ViewType2::array_layout,
         class IndexType = int,
         const int rank = ViewType1::rank>
KOKKOS_INLINE_FUNCTION void
COPY (const ViewType1& x, const ViewType2& y)
{
  static_assert (static_cast<int> (ViewType1::rank) ==
                 static_cast<int> (ViewType2::rank),
                 "COPY: x and y must have the same rank.");
  Impl::COPY<ViewType1, ViewType2, LayoutType1, LayoutType2, IndexType,
    rank>::run (x, y);
}

/// \brief <tt>y := y + alpha * A * x</tt> (dense matrix-vector multiply)
///
/// \param alpha [in] Coefficient by which to multiply A*x (this does
///   NOT necessarily follow BLAS rules; the caller is responsible for
///   checking whether alpha == 0 and implementing BLAS rules in that
///   case).
/// \param A [in] Small dense matrix (must have rank 2)
/// \param x [in] Small dense vector input (must have rank 1 and at
///   least as many rows as A has columns)
/// \param y [in/out] Small dense vector output (must have rank 1 and
///   at least as many rows as A has rows)
template<class VecType1,
         class BlkType,
         class VecType2,
         class CoeffType,
         class IndexType = int>
KOKKOS_INLINE_FUNCTION void
GEMV (const CoeffType& alpha,
      const BlkType& A,
      const VecType1& x,
      const VecType2& y)
{
  Impl::GEMV<VecType1, BlkType, VecType2, CoeffType,
    IndexType>::run (alpha, A, x, y);
}

/// \brief Small dense matrix-matrix multiply: <tt>C := alpha*A*B + beta*C</tt>
///
/// \tparam ViewType1 Type of the first matrix input A.
/// \tparam ViewType2 Type of the second matrix input B.
/// \tparam ViewType3 Type of the third matrix input/output C.
/// \tparam CoefficientType Type of the scalar coefficients alpha and beta.
/// \tparam IndexType Type of the index used in for loops; defaults to \c int.
template<class ViewType1,
         class ViewType2,
         class ViewType3,
         class CoefficientType,
         class IndexType = int>
KOKKOS_INLINE_FUNCTION void
GEMM (const char transA[],
      const char transB[],
      const CoefficientType& alpha,
      const ViewType1& A,
      const ViewType2& B,
      const CoefficientType& beta,
      const ViewType3& C)
{
  // Assert that A, B, and C are in fact matrices
  static_assert (ViewType1::rank == 2, "GEMM: A must have rank 2 (be a matrix).");
  static_assert (ViewType2::rank == 2, "GEMM: B must have rank 2 (be a matrix).");
  static_assert (ViewType3::rank == 2, "GEMM: C must have rank 2 (be a matrix).");

  typedef typename std::remove_reference<decltype (A(0,0))>::type Scalar;
  typedef Kokkos::Details::ArithTraits<Scalar> STS;
  const Scalar ZERO = STS::zero();
  const Scalar ONE = STS::one();

  // Get the dimensions
  IndexType m, n, k;
  if(transA[0] == 'N' || transA[0] == 'n') {
    m = static_cast<IndexType> (A.extent (0));
    n = static_cast<IndexType> (A.extent (1));
  }
  else {
    m = static_cast<IndexType> (A.extent (1));
    n = static_cast<IndexType> (A.extent (0));
  }
  k = static_cast<IndexType> (C.extent (1));

  // quick return if possible
  if(alpha == ZERO && beta == ONE)
    return;

  // And if alpha equals zero...
  if(alpha == ZERO) {
    if(beta == ZERO) {
      for(IndexType i=0; i<m; i++) {
        for(IndexType j=0; j<k; j++) {
          C(i,j) = ZERO;
        }
      }
    }
    else {
      for(IndexType i=0; i<m; i++) {
        for(IndexType j=0; j<k; j++) {
          C(i,j) = beta*C(i,j);
        }
      }
    }
  }

  // Start the operations
  if(transB[0] == 'n' || transB[0] == 'N') {
    if(transA[0] == 'n' || transA[0] == 'N') {
      // Form C = alpha*A*B + beta*C
      for(IndexType j=0; j<n; j++) {
        if(beta == ZERO) {
          for(IndexType i=0; i<m; i++) {
            C(i,j) = ZERO;
          }
        }
        else if(beta != ONE) {
          for(IndexType i=0; i<m; i++) {
            C(i,j) = beta*C(i,j);
          }
        }
        for(IndexType l=0; l<k; l++) {
          Scalar temp = alpha*B(l,j);
          for(IndexType i=0; i<m; i++) {
            C(i,j) = C(i,j) + temp*A(i,l);
          }
        }
      }
    }
    else {
      // Form C = alpha*A**T*B + beta*C
      for(IndexType j=0; j<n; j++) {
        for(IndexType i=0; i<m; i++) {
          Scalar temp = ZERO;
          for(IndexType l=0; l<k; l++) {
            temp = temp + A(l,i)*B(l,j);
          }
          if(beta == ZERO) {
            C(i,j) = alpha*temp;
          }
          else {
            C(i,j) = alpha*temp + beta*C(i,j);
          }
        }
      }
    }
  }
  else {
    if(transA[0] == 'n' || transA[0] == 'N') {
      // Form C = alpha*A*B**T + beta*C
      for(IndexType j=0; j<n; j++) {
        if(beta == ZERO) {
          for(IndexType i=0; i<m; i++) {
            C(i,j) = ZERO;
          }
        }
        else if(beta != ONE) {
          for(IndexType i=0; i<m; i++) {
            C(i,j) = beta*C(i,j);
          }
        }
        for(IndexType l=0; l<k; l++) {
          Scalar temp = alpha*B(j,l);
          for(IndexType i=0; i<m; i++) {
            C(i,j) = C(i,j) + temp*A(i,l);
          }
        }
      }
    }
    else {
      // Form C = alpha*A**T*B**T + beta*C
      for(IndexType j=0; j<n; j++) {
        for(IndexType i=0; i<m; i++) {
          Scalar temp = ZERO;
          for(IndexType l=0; l<k; l++) {
            temp = temp + A(l,i)*B(j,l);
          }
          if(beta == ZERO) {
            C(i,j) = alpha*temp;
          }
          else {
            C(i,j) = alpha*temp + beta*C(i,j);
          }
        }
      }
    }
  }
}

/// \brief Computes A = P*L*U
template<class LittleBlockType,
         class LittleVectorType>
KOKKOS_INLINE_FUNCTION void
GETF2 (const LittleBlockType& A, const LittleVectorType& ipiv, int& info)
{
  // The type of an entry of ipiv is the index type.
  typedef typename std::decay<decltype (ipiv(0)) >::type IndexType;
  static_assert (std::is_integral<IndexType>::value,
                 "GETF2: The type of each entry of ipiv must be an integer type.");
  typedef typename std::remove_reference<decltype (A(0,0))>::type Scalar;
  static_assert (! std::is_const<Scalar>::value,
                 "GETF2: A must not be a const View (or LittleBlock).");
  static_assert (! std::is_const<std::remove_reference<decltype (ipiv(0))>>::value,
                 "GETF2: ipiv must not be a const View (or LittleBlock).");
  static_assert (LittleBlockType::rank == 2, "GETF2: A must have rank 2 (be a matrix).");
  typedef Kokkos::Details::ArithTraits<Scalar> STS;
  const Scalar ZERO = STS::zero();

  const IndexType numRows = static_cast<IndexType> (A.extent (0));
  const IndexType numCols = static_cast<IndexType> (A.extent (1));
  const IndexType pivDim = static_cast<IndexType> (ipiv.extent (0));

  // std::min is not a CUDA device function
  const IndexType minPivDim = (numRows < numCols) ? numRows : numCols;
  if (pivDim < minPivDim) {
    info = -2;
    return;
  }

  // Initialize info
  info = 0;

  for(IndexType j=0; j < pivDim; j++)
  {
    // Find pivot and test for singularity
    IndexType jp = j;
    for(IndexType i=j+1; i<numRows; i++)
    {
      if(STS::abs(A(i,j)) > STS::abs(A(jp,j))) {
        jp = i;
      }
    }
    ipiv(j) = jp+1;

    if(A(jp,j) != ZERO)
    {
      // Apply the interchange to columns 1:N
      if(jp != j)
      {
        for(IndexType i=0; i < numCols; i++)
        {
          Scalar temp = A(jp,i);
          A(jp,i) = A(j,i);
          A(j,i) = temp;
        }
      }

      // Compute elements J+1:M of J-th column
      for(IndexType i=j+1; i<numRows; i++) {
        A(i,j) = A(i,j) / A(j,j);
      }
    }
    else if(info == 0) {
      info = j;
    }

    // Update trailing submatrix
    for(IndexType r=j+1; r < numRows; r++)
    {
      for(IndexType c=j+1; c < numCols; c++) {
        A(r,c) = A(r,c) - A(r,j) * A(j,c);
      }
    }
  }
}

namespace Impl {

/// \brief Computes the solution to Ax=b
///
/// We have not implemented transpose yet, or multiple RHS
template<class LittleBlockType,
         class LittleIntVectorType,
         class LittleScalarVectorType,
         const int rank = LittleScalarVectorType::rank>
struct GETRS {
  static KOKKOS_INLINE_FUNCTION void
  run (const char mode[],
       const LittleBlockType& A,
       const LittleIntVectorType& ipiv,
       const LittleScalarVectorType& B,
       int& info);
};

//! Special case of GETRS for a single right-hand side.
template<class LittleBlockType,
         class LittleIntVectorType,
         class LittleScalarVectorType>
struct GETRS<LittleBlockType, LittleIntVectorType, LittleScalarVectorType, 1> {
  static KOKKOS_INLINE_FUNCTION void
  run (const char mode[],
       const LittleBlockType& A,
       const LittleIntVectorType& ipiv,
       const LittleScalarVectorType& B,
       int& info)
  {
    // The type of an entry of ipiv is the index type.
    typedef typename std::decay<decltype (ipiv(0))>::type IndexType;
    // IndexType must be signed, because this code does a countdown loop
    // to zero.  Unsigned integers are always >= 0, even on underflow.
    static_assert (std::is_integral<IndexType>::value &&
                   std::is_signed<IndexType>::value,
                   "GETRS: The type of each entry of ipiv must be a signed integer.");
    typedef typename std::decay<decltype (A(0,0))>::type Scalar;
    static_assert (! std::is_const<std::remove_reference<decltype (B(0))>>::value,
                   "GETRS: B must not be a const View (or LittleBlock).");
    static_assert (LittleBlockType::rank == 2, "GETRS: A must have rank 2 (be a matrix).");
    static_assert (LittleIntVectorType::rank == 1, "GETRS: ipiv must have rank 1.");
    static_assert (LittleScalarVectorType::rank == 1, "GETRS: For this specialization, B must have rank 1.");

    typedef Kokkos::Details::ArithTraits<Scalar> STS;
    const Scalar ZERO = STS::zero();
    const IndexType numRows = static_cast<IndexType> (A.extent (0));
    const IndexType numCols = static_cast<IndexType> (A.extent (1));
    const IndexType pivDim = static_cast<IndexType> (ipiv.extent (0));

    info = 0;

    // Ensure that the matrix is square
    if (numRows != numCols) {
      info = -2;
      return;
    }

    // Ensure that the pivot array is sufficiently large
    if (pivDim < numRows) {
      info = -3;
      return;
    }

    // No transpose case
    if(mode[0] == 'n' || mode[0] == 'N') {
      // Apply row interchanges to the RHS
      for(IndexType i=0; i<numRows; i++) {
        if(ipiv(i) != i+1) {
          Scalar temp = B(i);
          B(i) = B(ipiv(i)-1);
          B(ipiv(i)-1) = temp;
        }
      }

      // Solve Lx=b, overwriting b with x
      for(IndexType r=1; r < numRows; r++) {
        for(IndexType c=0; c < r; c++) {
          B(r) = B(r) - A(r,c)*B(c);
        }
      }

      // Solve Ux=b, overwriting b with x
      for(IndexType r=numRows-1; r >= 0; r--) {
        // Check whether U is singular
        if(A(r,r) == ZERO) {
          info = r+1;
          return;
        }

        for(IndexType c=r+1; c < numCols; c++) {
          B(r) = B(r) - A(r,c)*B(c);
        }
        B(r) = B(r) / A(r,r);
      }
    }
    // Transpose case
    else if(mode[0] == 't' || mode[0] == 'T') {
      info = -1; // NOT YET IMPLEMENTED
      return;
    }
    // Conjugate transpose case
    else if (mode[0] == 'c' || mode[0] == 'C') {
      info = -1; // NOT YET IMPLEMENTED
      return;
    }
    else { // invalid mode
      info = -1;
      return;
    }
  }
};


//! Special case of GETRS for multiple right-hand sides.
template<class LittleBlockType,
         class LittleIntVectorType,
         class LittleScalarVectorType>
struct GETRS<LittleBlockType, LittleIntVectorType, LittleScalarVectorType, 2> {
  static KOKKOS_INLINE_FUNCTION void
  run (const char mode[],
       const LittleBlockType& A,
       const LittleIntVectorType& ipiv,
       const LittleScalarVectorType& B,
       int& info)
  {
    // The type of an entry of ipiv is the index type.
    typedef typename std::decay<decltype (ipiv(0)) >::type IndexType;
    static_assert (std::is_integral<IndexType>::value,
                   "GETRS: The type of each entry of ipiv must be an integer type.");
    static_assert (! std::is_const<std::remove_reference<decltype (B(0)) > >::value,
                   "GETRS: B must not be a const View (or LittleBlock).");
    static_assert (LittleBlockType::rank == 2, "GETRS: A must have rank 2 (be a matrix).");
    static_assert (LittleIntVectorType::rank == 1, "GETRS: ipiv must have rank 1.");
    static_assert (LittleScalarVectorType::rank == 2, "GETRS: For this specialization, B must have rank 2.");

    // The current implementation iterates over one right-hand side at
    // a time.  It might be faster to do this differently, but this
    // should work for now.
    const IndexType numRhs = B.extent (1);
    info = 0;

    for (IndexType rhs = 0; rhs < numRhs; ++rhs) {
      auto B_cur = Kokkos::subview (B, Kokkos::ALL (), rhs);
      GETRS<LittleBlockType, LittleIntVectorType, decltype (B_cur), 1>::run (mode, A, ipiv, B_cur, info);
      if (info != 0) {
        return;
      }
    }
  }
};

} // namespace Impl

/// \brief Solve the linear system(s) AX=B, using the result of GETRF or GETF2.
///
/// \warning We have not implemented transpose yet, or multiple right-hand sides.
template<class LittleBlockType,
         class LittleIntVectorType,
         class LittleScalarVectorType>
KOKKOS_INLINE_FUNCTION void
GETRS (const char mode[],
       const LittleBlockType& A,
       const LittleIntVectorType& ipiv,
       const LittleScalarVectorType& B,
       int& info)
{
  Impl::GETRS<LittleBlockType, LittleIntVectorType, LittleScalarVectorType,
    LittleScalarVectorType::rank>::run (mode, A, ipiv, B, info);
}


/// \brief Compute inverse of A, using result of GETRF or GETF2.
///
/// \tparam LittleBlockType Type of dense matrix \c A
/// \tparam LittleBlockType Type of 1-D pivot array \c ipiv
/// \tparam LittleScalarVectorType Type of 1-D work array \c work
///
/// \param A [in/out] On input: output matrix resulting from running
///   GETRF or GETF2 on a square matrix A.  On output: inverse of the
///   original matrix A.
/// \param ipiv [in] Pivot array from the LU factorization.
/// \param work [out] Temporary workspace; must be at least as long as
///   the number of rows in A.
/// \param info [out] On output, 0 if the routine was successful, else
///   nonzero.
template<class LittleBlockType,
         class LittleIntVectorType,
         class LittleScalarVectorType>
KOKKOS_INLINE_FUNCTION void
GETRI (const LittleBlockType& A,
       const LittleIntVectorType& ipiv,
       const LittleScalarVectorType& work,
       int& info)
{
  // The type of an entry of ipiv is the index type.
  typedef typename std::decay<decltype (ipiv(0))>::type IndexType;
  // IndexType must be signed, because this code does a countdown loop
  // to zero.  Unsigned integers are always >= 0, even on underflow.
  static_assert (std::is_integral<IndexType>::value &&
                 std::is_signed<IndexType>::value,
                 "GETRI: The type of each entry of ipiv must be a signed integer.");
  typedef typename std::remove_reference<decltype (A(0,0))>::type Scalar;
  static_assert (! std::is_const<std::remove_reference<decltype (A(0,0))>>::value,
                 "GETRI: A must not be a const View (or LittleBlock).");
  static_assert (! std::is_const<std::remove_reference<decltype (work(0))>>::value,
                 "GETRI: work must not be a const View (or LittleBlock).");
  static_assert (LittleBlockType::rank == 2, "GETRI: A must have rank 2 (be a matrix).");
  typedef Kokkos::Details::ArithTraits<Scalar> STS;
  const Scalar ZERO = STS::zero();
  const Scalar ONE = STS::one();

  const IndexType numRows = static_cast<IndexType> (A.extent (0));
  const IndexType numCols = static_cast<IndexType> (A.extent (1));
  const IndexType pivDim = static_cast<IndexType> (ipiv.extent (0));
  const IndexType workDim = static_cast<IndexType> (work.extent (0));

  info = 0;

  // Ensure that the matrix is square
  if (numRows != numCols) {
    info = -1;
    return;
  }

  // Ensure that the pivot array is sufficiently large
  if (pivDim < numRows) {
    info = -2;
    return;
  }

  // Ensure that the work array is sufficiently large
  if (workDim < numRows) {
    info = -3;
    return;
  }

  // Form Uinv in place
  for(IndexType j=0; j < numRows; j++) {
    if(A(j,j) == ZERO) {
      info = j+1;
      return;
    }

    A(j,j) = ONE / A(j,j);

    // Compute elements 1:j-1 of j-th column
    for(IndexType r=0; r < j; r++) {
      A(r,j) = A(r,r)*A(r,j);
      for(IndexType c=r+1; c < j; c++) {
        A(r,j) = A(r,j) + A(r,c)*A(c,j);
      }
    }
    for(IndexType r=0; r < j; r++) {
      A(r,j) = -A(j,j)*A(r,j);
    }
  }

  // Compute Ainv by solving A\L = Uinv
  for(IndexType j = numCols-2; j >= 0; j--) {
    // Copy lower triangular data to work array and replace with 0
    for(IndexType r=j+1; r < numRows; r++) {
      work(r) = A(r,j);
      A(r,j) = 0;
    }

    for(IndexType r=0; r < numRows; r++) {
      for(IndexType i=j+1; i < numRows; i++) {
        A(r,j) = A(r,j) - work(i)*A(r,i);
      }
    }
  }

  // Apply column interchanges
  for(IndexType j=numRows-1; j >= 0; j--) {
    IndexType jp = ipiv(j)-1;
    if(j != jp) {
      for(IndexType r=0; r < numRows; r++) {
        Scalar temp = A(r,j);
        A(r,j) = A(r,jp);
        A(r,jp) = temp;
      }
    }
  }
}


// mfh 08 Nov 2015: I haven't tested this overload yet.  It also needs
// an implementation for trans != 'N' (the transpose and conjugate
// transpose cases).
#if 0
template<class LittleBlockType,
         class LittleVectorTypeX,
         class LittleVectorTypeY,
         class CoefficientType,
         class IndexType = int>
KOKKOS_INLINE_FUNCTION void
GEMV (const char trans,
      const CoefficientType& alpha,
      const LittleBlockType& A,
      const LittleVectorTypeX& x,
      const CoefficientType& beta,
      const LittleVectorTypeY& y)
{
  // y(0) returns a reference to the 0-th entry of y.  Remove that
  // reference to get the type of each entry of y.  It's OK if y has
  // zero entries -- this doesn't actually do y(i), it just returns
  // the type of that expression.
  typedef typename std::remove_reference<decltype (y(0)) >::type y_value_type;
  const IndexType numRows = static_cast<IndexType> (A.extent (0));
  const IndexType numCols = static_cast<IndexType> (A.extent (1));

  if (beta == 0.0) {
    if (alpha == 0.0) {
      for (IndexType i = 0; i < numRows; ++i) {
        y(i) = 0.0;
      }
    }
    else {
      for (IndexType i = 0; i < numRows; ++i) {
        y_value_type y_i = 0.0;
        for (IndexType j = 0; j < numCols; ++j) {
          y_i += A(i,j) * x(j);
        }
        y(i) = y_i;
      }
    }
  }
  else { // beta != 0
    if (alpha == 0.0) {
      if (beta == 0.0) {
        for (IndexType i = 0; i < numRows; ++i) {
          y(i) = 0.0;
        }
      }
      else {
        for (IndexType i = 0; i < numRows; ++i) {
          y(i) *= beta;
        }
      }
    }
    else {
      for (IndexType i = 0; i < numRows; ++i) {
        y_value_type y_i = beta * y(i);
        for (IndexType j = 0; j < numCols; ++j) {
          y_i += alpha * A(i,j) * x(j);
        }
        y(i) = y_i;
      }
    }
  }
}

#endif // 0

} // namespace Experimental
} // namespace Tpetra

#endif // TPETRA_EXPERIMENTAL_BLOCKVIEW_HPP
