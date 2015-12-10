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
/// \brief Declaration and definition of LittleBlock and LittleVector

#include "Tpetra_ConfigDefs.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_LAPACK.hpp"
#ifdef HAVE_TPETRA_INST_FLOAT128
#  include "Teuchos_BLAS.hpp"
#endif // HAVE_TPETRA_INST_FLOAT128

#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_Complex.hpp"

#ifdef HAVE_TPETRA_INST_FLOAT128
#  include "Teuchos_Details_Lapack128.hpp"
#endif // HAVE_TPETRA_INST_FLOAT128

namespace Tpetra {
namespace Details {

  /// \brief Return the Teuchos::LAPACK specialization corresponding
  ///   to the given Scalar type.
  ///
  /// The reason this exists is the same reason why the
  /// impl_scalar_type typedef in Tpetra::MultiVector may differ from
  /// its Scalar template parameter.  For example, Scalar =
  /// std::complex<T> corresponds to impl_scalar_type =
  /// Kokkos::complex<T>.  The latter has no Teuchos::LAPACK
  /// specialization, so we have to map it back to std::complex<T>.
  template<class Scalar>
  struct GetLapackType {
    typedef Scalar lapack_scalar_type;
    typedef Teuchos::LAPACK<int, Scalar> lapack_type;
  };

  template<class T>
  struct GetLapackType<Kokkos::complex<T> > {
    typedef std::complex<T> lapack_scalar_type;
    typedef Teuchos::LAPACK<int, std::complex<T> > lapack_type;
  };

#ifdef HAVE_TPETRA_INST_FLOAT128
  template<>
  struct GetLapackType<__float128> {
    typedef __float128 lapack_scalar_type;
    // Use the Lapack128 class we declared above to implement the
    // linear algebra operations needed for small dense blocks and
    // vectors.
    typedef Teuchos::Details::Lapack128 lapack_type;
  };
#endif // HAVE_TPETRA_INST_FLOAT128

} // namespace Details
} // namespace Tpetra


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
         const int rank = ViewType1::rank>
struct AbsMax {
  static void run (const ViewType2& Y, const ViewType1& X);
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
  static void run (const ViewType2& Y, const ViewType1& X)
  {
    typedef typename std::remove_reference<decltype (Y(0,0)) >::type STY;
    static_assert(! std::is_const<STY>::value,
      "AbsMax: The type of each entry of Y must be nonconst.");
    typedef typename std::remove_const<typename std::remove_reference<decltype (X(0,0)) >::type>::type STX;
    static_assert(  std::is_same<STX, STY>::value,
      "AbsMax: The type of each entry of X and Y must be the same.");
    typedef Kokkos::Details::ArithTraits<STY> KAT;

    for (int j = 0; j < Y.dimension_1 (); ++j) {
      for (int i = 0; i < Y.dimension_0 (); ++i) {
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
  static void run (const ViewType2& Y, const ViewType1& X)
  {
    typedef typename std::remove_reference<decltype (Y(0)) >::type STY;
    static_assert(! std::is_const<STY>::value,
      "AbsMax: The type of each entry of Y must be nonconst.");
    typedef typename std::remove_const<typename std::remove_reference<decltype (X(0)) >::type>::type STX;
    static_assert(  std::is_same<STX, STY>::value,
      "AbsMax: The type of each entry of X and Y must be the same.");
    typedef Kokkos::Details::ArithTraits<STY> KAT;

    for (int i = 0; i < Y.dimension_0 (); ++i) {
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
void absMax (const ViewType2& Y, const ViewType1& X) {
  AbsMax<ViewType1, ViewType2, rank>::run (Y, X);
}

/// \brief Implementation of Tpetra::Experimental::SCAL function.
///
/// This is the "generic" version that we don't implement.
/// We actually implement versions for ViewType rank 1 or rank 2.
template<class ViewType,
         class CoefficientType,
         class IndexType = int,
         const int rank = ViewType::rank>
struct SCAL {
  static void run (const CoefficientType& alpha, const ViewType& x);
};

/// \brief Implementation of Tpetra::Experimental::SCAL function, for
///   ViewType rank 1 (i.e., a vector).
template<class ViewType,
         class CoefficientType,
         class IndexType>
struct SCAL<ViewType, CoefficientType, IndexType, 1> {
  /// \brief x := alpha*x (rank-1 x, i.e., a vector)
  static void run (const CoefficientType& alpha, const ViewType& x)
  {
    const IndexType numRows = static_cast<IndexType> (x.dimension_0 ());
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
         class IndexType>
struct SCAL<ViewType, CoefficientType, IndexType, 2> {
  /// \brief A := alpha*A (rank-2 A, i.e., a matrix)
  static void run (const CoefficientType& alpha, const ViewType& A)
  {
    const IndexType numRows = static_cast<IndexType> (A.dimension_0 ());
    const IndexType numCols = static_cast<IndexType> (A.dimension_1 ());

    // BLAS _SCAL doesn't check whether alpha is 0.
    for (IndexType i = 0; i < numRows; ++i) {
      for (IndexType j = 0; j < numCols; ++j) {
        A(i,j) = alpha * A(i,j);
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
         class IndexType = int,
         const int rank = ViewType1::rank>
struct AXPY {
  static void
  run (const CoefficientType& alpha,
       const ViewType1& x,
       const ViewType2& y);
};

/// \brief Implementation of Tpetra::Experimental::AXPY function, for
///   ViewType1 and ViewType2 rank 1 (i.e., vectors).
template<class CoefficientType,
         class ViewType1,
         class ViewType2,
         class IndexType>
struct AXPY<CoefficientType, ViewType1, ViewType2, IndexType, 1> {
  /// \brief y := y + alpha*x (rank-1 x and y, i.e., vectors)
  static void
  run (const CoefficientType& alpha,
       const ViewType1& x,
       const ViewType2& y)
  {
    static_assert (ViewType1::rank == ViewType2::rank,
                   "AXPY: x and y must have the same rank.");
    const IndexType numRows = static_cast<IndexType> (y.dimension_0 ());
    if (alpha != 0.0) {
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
         class IndexType>
struct AXPY<CoefficientType, ViewType1, ViewType2, IndexType, 2> {
  /// \brief Y := Y + alpha*X (rank-2 X and Y, i.e., matrices)
  static void
  run (const CoefficientType& alpha,
       const ViewType1& X,
       const ViewType2& Y)
  {
    static_assert (ViewType1::rank == ViewType2::rank,
                   "AXPY: X and Y must have the same rank.");
    const IndexType numRows = static_cast<IndexType> (Y.dimension_0 ());
    const IndexType numCols = static_cast<IndexType> (Y.dimension_1 ());

    if (alpha != 0.0) {
      for (IndexType i = 0; i < numRows; ++i) {
        for (IndexType j = 0; j < numCols; ++j) {
          Y(i,j) += alpha * X(i,j);
        }
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
         class IndexType,
         const int rank = ViewType1::rank>
struct COPY {
  static void run (const ViewType1& x, const ViewType2& y);
};

/// \brief Implementation of Tpetra::Experimental::COPY function, for
///   ViewType1 and ViewType2 rank 1 (i.e., vectors).
template<class ViewType1,
         class ViewType2,
         class IndexType>
struct COPY<ViewType1, ViewType2, IndexType, 1> {
  /// \brief y := x (rank-1 x and y, i.e., vectors)
  static void run (const ViewType1& x, const ViewType2& y)
  {
    const IndexType numRows = static_cast<IndexType> (x.dimension_0 ());
    for (IndexType i = 0; i < numRows; ++i) {
      y(i) = x(i);
    }
  }
};

/// \brief Implementation of Tpetra::Experimental::COPY function, for
///   ViewType1 and ViewType2 rank 2 (i.e., matrices).
template<class ViewType1,
         class ViewType2,
         class IndexType>
struct COPY<ViewType1, ViewType2, IndexType, 2> {
  /// \brief Y := X (rank-2 X and Y, i.e., matrices)
  static void run (const ViewType1& X, const ViewType2& Y)
  {
    const IndexType numRows = static_cast<IndexType> (Y.dimension_0 ());
    const IndexType numCols = static_cast<IndexType> (Y.dimension_1 ());

    // BLAS _SCAL doesn't check whether alpha is 0.
    for (IndexType i = 0; i < numRows; ++i) {
      for (IndexType j = 0; j < numCols; ++j) {
        Y(i,j) = X(i,j);
      }
    }
  }
};

} // namespace Impl

/// \brief x := alpha*x, where x is either rank 1 (a vector) or rank 2
///   (a matrix).
template<class ViewType,
         class CoefficientType,
         class IndexType = int,
         const int rank = ViewType::rank>
void SCAL (const CoefficientType& alpha, const ViewType& x) {
  Impl::SCAL<ViewType, CoefficientType, IndexType, rank>::run (alpha, x);
}

/// \brief y := y + alpha * x
///
/// This function follows the BLAS convention that if alpha == 0, then
/// it does nothing.  (This matters only if x contains Inf or NaN
/// values.)
template<class CoefficientType,
         class ViewType1,
         class ViewType2,
         class IndexType = int,
         const int rank = ViewType1::rank>
void
AXPY (const CoefficientType& alpha,
      const ViewType1& x,
      const ViewType2& y)
{
  static_assert (ViewType1::rank == ViewType1::rank,
                 "AXPY: x and y must have the same rank.");
  Impl::AXPY<CoefficientType, ViewType1, ViewType2, IndexType, rank>::run (alpha, x, y);
}

/// \brief x := alpha*x, where x is either rank 1 (a vector) or rank 2
///   (a matrix).
template<class ViewType1,
         class ViewType2,
         class IndexType = int,
         const int rank = ViewType1::rank>
void COPY (const ViewType1& x, const ViewType2& y) {
  static_assert (ViewType1::rank == ViewType1::rank,
                 "COPY: x and y must have the same rank.");
  Impl::COPY<ViewType1, ViewType2, IndexType, rank>::run (x, y);
}

/// \brief y := y + alpha * A * x
///
/// This function follows the BLAS convention that if alpha == 0, then
/// it does nothing.  (This matters only if x contains Inf or NaN
/// values.)
template<class LittleVectorType1,
         class LittleBlockType,
         class LittleVectorType2,
         class CoefficientType,
         class IndexType = int>
void
GEMV (const CoefficientType& alpha,
      const LittleBlockType& A,
      const LittleVectorType1& x,
      const LittleVectorType2& y)
{
  const IndexType numRows = static_cast<IndexType> (A.dimension_0 ());
  const IndexType numCols = static_cast<IndexType> (A.dimension_1 ());

  if (alpha != 0.0) {
    for (IndexType i = 0; i < numRows; ++i) {
      typename std::remove_reference<decltype (y(i)) >::type y_i = y(i);
      for (IndexType j = 0; j < numCols; ++j) {
        y_i += alpha * A(i,j) * x(j);
      }
      y(i) = y_i;
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
void
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
  const IndexType numRows = static_cast<IndexType> (A.dimension_0 ());
  const IndexType numCols = static_cast<IndexType> (A.dimension_1 ());

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

/// \class LittleBlock
/// \brief Nonowning view of a square dense block in a block matrix.
/// \tparam Scalar The type of entries in the block.
/// \tparam LO The type of local indices.  See the documentation of
///   the first template parameter of Map for requirements.
///
/// "Little" means local (not distributed over multiple MPI processes;
/// stored to maximize locality) and small (think 3x3, not 1000x1000).
///
/// The \c Scalar template parameter may be const or nonconst.  This
/// is one reason why instance methods below that take a LittleBlock
/// accept it as a template parameter: that lets you add a const
/// LittleBlock (e.g., LittleBlock<const double, int>) to a nonconst
/// LittleBlock (e.g., LittleBlock<double, int>).
template<class Scalar, class LO>
class LittleBlock {
public:
  typedef Scalar scalar_type;
  typedef typename Kokkos::Details::ArithTraits<Scalar>::val_type impl_scalar_type;

private:
  typedef Kokkos::Details::ArithTraits<impl_scalar_type> STS;

public:
  //! Number of dimensions
  static const int rank = 2;

  /// \brief Constructor
  /// \param A [in] Pointer to the block's entries
  /// \param blockSize [in] Dimension of the block (all blocks are square)
  /// \param strideX [in] Stride between consecutive entries in a column
  /// \param strideY [in] Stride between consecutive entries in a row
  LittleBlock (Scalar* const A,
               const LO blockSize,
               const LO strideX,
               const LO strideY) :
    A_ (reinterpret_cast<impl_scalar_type*> (A)),
    blockSize_ (blockSize),
    strideX_ (strideX),
    strideY_ (strideY)
  {}

  /// \brief Constructor that takes an \c impl_scalar_type pointer.
  ///
  /// \param A [in] Pointer to the block's entries, as
  ///   <tt>impl_scalar_type*</tt> rather than <tt>Scalar*</tt>
  /// \param blockSize [in] Dimension of the block (all blocks are square)
  /// \param strideX [in] Stride between consecutive entries in a column
  /// \param strideY [in] Stride between consecutive entries in a row
  ///
  /// While this constructor is templated on a type \c T, the intent
  /// is that <tt>T == impl_scalar_type</tt>.  (We must template on T
  /// rather than using <tt>impl_scalar_type</tt> directly, because of
  /// how std::enable_if works.)  The long, complicated std::enable_if
  /// expression ensures that this constructor only exists if
  /// <tt>Scalar</tt> differs from <tt>impl_scalar_type</tt>, but the
  /// two types are mutually compatible and have the same size.  (They
  /// must be bitwise compatible, so that \c reinterpret_cast makes
  /// sense between them.)
  template<class T>
  LittleBlock (T* const A,
               const LO blockSize,
               const LO strideX,
               const LO strideY,
               typename std::enable_if<
                 ! std::is_same<Scalar, T>::value &&
                 std::is_convertible<Scalar, T>::value &&
                 sizeof (Scalar) == sizeof (T),
               int*>::type ignoreMe = NULL) :
    A_ (reinterpret_cast<impl_scalar_type*> (A)),
    blockSize_ (blockSize),
    strideX_ (strideX),
    strideY_ (strideY)
  {}

  //! The block size (number of rows, and number of columns).
  LO getBlockSize () const {
    return blockSize_;
  }

  //! Number of rows in the block.
  LO dimension_0 () const {
    return blockSize_;
  }

  //! Number of columns in the block.
  LO dimension_1 () const {
    return blockSize_;
  }

  template<class IntegerType>
  void stride (IntegerType* const s) const {
    s[0] = strideX_;
    s[1] = strideY_;
  }

  //! Pointer to the block's entries, as <tt>Scalar*</tt>.
  Scalar* ptr_on_device () const {
    return reinterpret_cast<Scalar*> (A_);
  }

  //! Pointer to the block's entries, as <tt>Scalar*</tt>.
  Scalar* getRawPtr () const {
    return reinterpret_cast<Scalar*> (A_);
  }

  /// \brief Reference to entry (i,j) of the block.
  ///
  /// \note To Tpetra developers: This is returned as
  ///   <tt>impl_scalar_type</tt> and not as \c Scalar, in order to
  ///   avoid a lot of reinterpret_cast calls in the inner loop of the
  ///   sparse matrix-vector multiply kernel of
  ///   Tpetra::Experimental::BlockCrsMatrix.  Any pair of types
  ///   <tt>impl_scalar_type</tt>, \c Scalar used here should always
  ///   be convertible in either direction, so the return type should
  ///   not pose any issues in practice.
  impl_scalar_type& operator() (const LO i, const LO j) const {
    return A_[i * strideX_ + j * strideY_];
  }

  //! <tt>*this := *this + alpha * X</tt>.
  template<class LittleBlockType>
  void update (const Scalar& alpha, const LittleBlockType& X) const {
    AXPY (static_cast<impl_scalar_type> (alpha), X, *this);
  }

  //! <tt>*this := X</tt>.
  template<class LittleBlockType>
  void assign (const LittleBlockType& X) const {
    COPY (X, *this);
  }

  //! <tt>(*this)(i,j) := alpha * (*this)(i,j)</tt> for all (i,j).
  void scale (const Scalar& alpha) const {
    SCAL (static_cast<impl_scalar_type> (alpha), *this);
  }

  //! <tt>(*this)(i,j) := alpha</tt> for all (i,j).
  void fill (const Scalar& alpha) const {
    const impl_scalar_type theAlpha = static_cast<Scalar> (alpha);
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) = theAlpha;
      }
    }
  }

  void factorize (int* ipiv, int & info)
  {
    typedef typename Tpetra::Details::GetLapackType<Scalar>::lapack_scalar_type LST;
    typedef typename Tpetra::Details::GetLapackType<Scalar>::lapack_type lapack_type;

    LST* const A_raw = reinterpret_cast<LST*> (A_);
    lapack_type lapack;
    // NOTE (mfh 03 Jan 2015) This method doesn't check the 'info'
    // output argument, but it returns info, so the user is
    // responsible for checking.
    lapack.GETRF(blockSize_, blockSize_, A_raw, blockSize_, ipiv, &info);
  }

  template<class LittleVectorType>
  void solve (LittleVectorType & X, const int* ipiv) const
  {
    typedef typename Tpetra::Details::GetLapackType<Scalar>::lapack_scalar_type LST;
    typedef typename Tpetra::Details::GetLapackType<Scalar>::lapack_type lapack_type;

    // FIXME (mfh 03 Jan 2015) Check using enable_if that Scalar can
    // be safely converted to LST.

    lapack_type lapack;
    LST* const A_raw = reinterpret_cast<LST*> (A_);
    LST* const X_raw = reinterpret_cast<LST*> (X.getRawPtr ());
    int info = 0;
    char trans = 'T';
    // FIXME (mfh 03 Jan 2015) Either check the 'info' output
    // argument, or return it.
    lapack.GETRS(trans, blockSize_, 1, A_raw, blockSize_, ipiv, X_raw, blockSize_, &info);
  }

private:
  impl_scalar_type* const A_;
  const LO blockSize_;
  const LO strideX_;
  const LO strideY_;
};


/// \class LittleVector
/// \brief Nonowning view of a set of degrees of freedom corresponding
///   to a mesh point in a block vector or multivector.
/// \tparam Scalar The type of entries.
/// \tparam LO The type of local indices.  See the documentation of
///   the first template parameter of Map for requirements.
///
/// "Little" means local (not distributed over multiple MPI processes;
/// stored to maximize locality) and small (think length 3, not length
/// 1000).
///
/// The \c Scalar template parameter may be const or nonconst.  This
/// is one reason why instance methods below that take a LittleVector
/// accept it as a template parameter: that lets you add a const
/// LittleVector (e.g., LittleVector<const double, int>) to a nonconst
/// LittleVector (e.g., LittleVector<double, int>).
template<class Scalar, class LO>
class LittleVector {
public:
  typedef Scalar scalar_type;
  typedef typename Kokkos::Details::ArithTraits<Scalar>::val_type impl_scalar_type;

private:
  typedef Kokkos::Details::ArithTraits<impl_scalar_type> STS;

public:
  //! Number of dimensions
  static const int rank = 1;

  /// \brief Constructor
  /// \param A [in] Pointer to the vector's entries
  /// \param blockSize [in] Dimension of the vector
  /// \param strideX [in] Stride between consecutive entries
  LittleVector (Scalar* const A, const LO blockSize, const LO strideX) :
    A_ (reinterpret_cast<impl_scalar_type*> (A)),
    blockSize_ (blockSize),
    strideX_ (strideX)
  {}

  /// \brief Constructor that takes an \c impl_scalar_type pointer.
  ///
  /// \param A [in] Pointer to the vector's entries, as
  ///   <tt>impl_scalar_type*</tt> rather than <tt>Scalar*</tt>
  /// \param blockSize [in] Dimension of the vector
  /// \param strideX [in] Stride between consecutive entries
  ///
  /// While this constructor is templated on a type \c T, the intent
  /// is that <tt>T == impl_scalar_type</tt>.  (We must template on T
  /// rather than using <tt>impl_scalar_type</tt> directly, because of
  /// how std::enable_if works.)  The long, complicated std::enable_if
  /// expression ensures that this constructor only exists if
  /// <tt>Scalar</tt> differs from <tt>impl_scalar_type</tt>, but the
  /// two types are mutually compatible and have the same size.  (They
  /// must be bitwise compatible, so that \c reinterpret_cast makes
  /// sense between them.)
  template<class T>
  LittleVector (T* const A,
                const LO blockSize,
                const LO strideX,
                typename std::enable_if<
                  ! std::is_same<Scalar, T>::value &&
                  std::is_convertible<Scalar, T>::value &&
                  sizeof (Scalar) == sizeof (T),
                int*>::type ignoreMe = NULL) :
    A_ (reinterpret_cast<impl_scalar_type*> (A)),
    blockSize_ (blockSize),
    strideX_ (strideX)
  {}

  //! Pointer to the vector's entries.
  Scalar* getRawPtr () const {
    return reinterpret_cast<Scalar*> (A_);
  }

  //! Pointer to the vector's entries.
  Scalar* ptr_on_device () const {
    return reinterpret_cast<Scalar*> (A_);
  }

  //! The block size (number of degrees of freedom per mesh point).
  LO getBlockSize () const {
    return blockSize_;
  }

  //! Number of entries in the vector.
  LO dimension_0 () const {
    return blockSize_;
  }

  //! Stride between consecutive entries.
  LO getStride () const {
    return strideX_;
  }

  /// \brief Stride between consecutive entries.
  ///
  /// This exists for compatibility with Kokkos::View.
  template<class IntegerType>
  void stride (IntegerType* const s) const {
    s[0] = strideX_;
  }

  /// \brief Reference to entry (i) of the vector.
  ///
  /// \note To Tpetra developers: This is returned as
  ///   <tt>impl_scalar_type</tt> and not as \c Scalar, in order to
  ///   avoid a lot of reinterpret_cast calls in the inner loop of the
  ///   sparse matrix-vector multiply kernel of
  ///   Tpetra::Experimental::BlockCrsMatrix.  Any pair of types
  ///   <tt>impl_scalar_type</tt>, \c Scalar used here should always
  ///   be convertible in either direction, so the return type should
  ///   not pose any issues in practice.
  impl_scalar_type& operator() (const LO i) const {
    return A_[i * strideX_];
  }

  //! <tt>*this := *this + alpha * X</tt>.
  template<class LittleVectorType>
  void update (const Scalar& alpha, const LittleVectorType& X) const {
    AXPY (static_cast<impl_scalar_type> (alpha), X, *this);
  }

  //! <tt>*this := X</tt>.
  template<class LittleVectorType>
  void assign (const LittleVectorType& X) const {
    COPY (X, *this);
  }

  //! <tt>(*this)(i) := alpha * (*this)(i)</tt> for all (i,j).
  void scale (const Scalar& alpha) const {
    SCAL (static_cast<impl_scalar_type> (alpha), *this);
  }

  //! <tt>(*this)(i,j) := alpha</tt> for all (i,j).
  void fill (const Scalar& alpha) const {
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) = theAlpha;
    }
  }

  //! true if and only if all entries of this equal all entries of X.
  template<class LittleVectorType>
  bool equal (const LittleVectorType& X) const {
    if (getBlockSize () != X.getBlockSize ()) {
      return false;
    }
    for (LO i = 0; i < blockSize_; ++i) {
      if ((*this)(i) != X(i)) {
        return false;
      }
    }
    return true;
  }

  //! <tt>(*this) := (*this) + alpha * A * X</tt> (matrix-vector multiply).
  template<class LittleBlockType, class LittleVectorType>
  void
  matvecUpdate (const Scalar& alpha,
                const LittleBlockType& A,
                const LittleVectorType& X) const
  {
    GEMV (static_cast<impl_scalar_type> (alpha), A, X, *this);
  }

private:
  impl_scalar_type* const A_;
  const LO blockSize_;
  const LO strideX_;
};




} // namespace Experimental
} // namespace Tpetra

#endif // TPETRA_EXPERIMENTAL_BLOCKVIEW_HPP
