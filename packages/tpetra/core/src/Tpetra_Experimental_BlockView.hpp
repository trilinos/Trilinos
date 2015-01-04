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

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_LAPACK.hpp>

#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
#  include <Kokkos_ArithTraits.hpp>
#  include <Kokkos_complex.hpp>
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

namespace { // anonymous

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

#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
  template<class T>
  struct GetLapackType<Kokkos::complex<T> > {
    typedef std::complex<T> lapack_scalar_type;
    typedef Teuchos::LAPACK<int, std::complex<T> > lapack_type;
  };
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

} // namespace (anonymous)

//#include "Teuchos_LAPACK_wrappers.hpp"
//extern "C" {int DGETRF_F77(const int *, const int *, double *, const int*, int *, int*);}

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
#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
  typedef typename Kokkos::Details::ArithTraits<Scalar>::val_type impl_scalar_type;
#else
  typedef Scalar impl_scalar_type;
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

private:
#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
  typedef Kokkos::Details::ArithTraits<impl_scalar_type> STS;
#else
  typedef Teuchos::ScalarTraits<impl_scalar_type> STS;
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

public:
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

#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
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
#  ifdef KOKKOS_HAVE_CXX11
               typename std::enable_if<
                 ! std::is_same<Scalar, T>::value &&
                 std::is_convertible<Scalar, T>::value &&
                 sizeof (Scalar) == sizeof (T),
#  else
               typename Kokkos::Impl::enable_if<
                 ! Kokkos::Impl::is_same<Scalar, T>::value &&
                 sizeof (Scalar) == sizeof (T),
#  endif // KOKKOS_HAVE_CXX11
               int*>::type ignoreMe = NULL) :
    A_ (reinterpret_cast<impl_scalar_type*> (A)),
    blockSize_ (blockSize),
    strideX_ (strideX),
    strideY_ (strideY)
  {}
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

  //! The block size (number of rows, and number of columns).
  LO getBlockSize () const {
    return blockSize_;
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
    const impl_scalar_type theAlpha = static_cast<Scalar> (alpha);
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) += theAlpha * X(i,j);
      }
    }
  }

  //! <tt>*this := X</tt>.
  template<class LittleBlockType>
  void assign (const LittleBlockType& X) const {
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) = X(i,j);
      }
    }
  }

  //! <tt>(*this)(i,j) := alpha * (*this)(i,j)</tt> for all (i,j).
  void scale (const Scalar& alpha) const {
    const impl_scalar_type theAlpha = static_cast<Scalar> (alpha);
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) *= theAlpha;
      }
    }
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

  /// \brief <tt>(*this)(i,j) := max(abs((*this)(i,j)), abs(X(i,j)))</tt>
  ///   for all (i,j).
  ///
  /// Tpetra uses this operation to implement the ABSMAX CombineMode.
  template<class LittleBlockType>
  void absmax (const LittleBlockType& X) const {
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        impl_scalar_type& Y_ij = (*this)(i,j);
        const impl_scalar_type X_ij = X(i,j);
        Y_ij = std::max (STS::magnitude (Y_ij), STS::magnitude (X_ij));
      }
    }
  }

  void factorize (int* ipiv, int & info)
  {
    typedef typename GetLapackType<Scalar>::lapack_scalar_type LST;
    typedef typename GetLapackType<Scalar>::lapack_type lapack_type;

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
    typedef typename GetLapackType<Scalar>::lapack_scalar_type LST;
    typedef typename GetLapackType<Scalar>::lapack_type lapack_type;

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
  // FIXME (mfh 04 Jan 2015) I strongly object to putting the LU
  // factorization's pivot array here.  Pivot arrays should be stored
  // separately in the preconditioner.  std::vector adds a void* and
  // two size_t values (size and capacity) to the struct.
  std::vector<int> ipiv_;
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
#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
  typedef typename Kokkos::Details::ArithTraits<Scalar>::val_type impl_scalar_type;
#else
  typedef Scalar impl_scalar_type;
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

private:
#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
  typedef Kokkos::Details::ArithTraits<impl_scalar_type> STS;
#else
  typedef Teuchos::ScalarTraits<impl_scalar_type> STS;
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

public:
  /// \brief Constructor
  /// \param A [in] Pointer to the vector's entries
  /// \param blockSize [in] Dimension of the vector
  /// \param stride [in] Stride between consecutive entries
  LittleVector (Scalar* const A, const LO blockSize, const LO stride) :
    A_ (reinterpret_cast<impl_scalar_type*> (A)),
    blockSize_ (blockSize),
    strideX_ (stride)
  {}

#ifdef TPETRA_HAVE_KOKKOS_REFACTOR
  /// \brief Constructor that takes an \c impl_scalar_type pointer.
  ///
  /// \param A [in] Pointer to the vector's entries, as
  ///   <tt>impl_scalar_type*</tt> rather than <tt>Scalar*</tt>
  /// \param blockSize [in] Dimension of the vector
  /// \param stride [in] Stride between consecutive entries
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
                const LO stride,
#  ifdef KOKKOS_HAVE_CXX11
                typename std::enable_if<
                  ! std::is_same<Scalar, T>::value &&
                  std::is_convertible<Scalar, T>::value &&
                  sizeof (Scalar) == sizeof (T),
#  else
                typename Kokkos::Impl::enable_if<
                  ! Kokkos::Impl::is_same<Scalar, T>::value &&
                  sizeof (Scalar) == sizeof (T),
#  endif // KOKKOS_HAVE_CXX11
                int*>::type ignoreMe = NULL) :
    A_ (reinterpret_cast<impl_scalar_type*> (A)),
    blockSize_ (blockSize),
    strideX_ (stride)
  {}
#endif // TPETRA_HAVE_KOKKOS_REFACTOR

  //! Pointer to the block's entries.
  Scalar* getRawPtr () const {
    return reinterpret_cast<Scalar*> (A_);
  }

  //! The block size (number of degrees of freedom per mesh point).
  LO getBlockSize () const {
    return blockSize_;
  }

  //! Stride between consecutive entries.
  LO getStride () const {
    return strideX_;
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
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) += theAlpha * X(i);
    }
  }

  //! <tt>*this := X</tt>.
  template<class LittleVectorType>
  void assign (const LittleVectorType& X) const {
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) = X(i);
    }
  }

  //! <tt>(*this)(i,j) := alpha * (*this)(i,j)</tt> for all (i,j).
  void scale (const Scalar& alpha) const {
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) *= theAlpha;
    }
  }

  //! <tt>(*this)(i,j) := alpha</tt> for all (i,j).
  void fill (const Scalar& alpha) const {
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) = theAlpha;
    }
  }

  /// \brief <tt>(*this)(i,j) := max(abs((*this)(i,j)), abs(X(i,j)))</tt>
  ///   for all (i,j).
  ///
  /// Tpetra uses this operation to implement the ABSMAX CombineMode.
  template<class LittleVectorType>
  void absmax (const LittleVectorType& X) const {
    for (LO i = 0; i < blockSize_; ++i) {
      impl_scalar_type& Y_i = (*this)(i);
      Y_i = std::max (STS::magnitude (Y_i), STS::magnitude (X (i)));
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
    const impl_scalar_type theAlpha = static_cast<impl_scalar_type> (alpha);
    // FIXME (mfh 07 May 2014) This is suitable for column major, not
    // for row major.  Of course, we'll have to change other loops
    // above as well to make row major faster.
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i) += theAlpha * A(i,j) * X(j);
      }
    }
  }

private:
  impl_scalar_type* const A_;
  const LO blockSize_;
  const LO strideX_;
};

} // namespace Experimental
} // namespace Tpetra

#endif // TPETRA_EXPERIMENTAL_BLOCKVIEW_HPP
