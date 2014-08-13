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

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_LAPACK.hpp>


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
private:
  typedef Teuchos::ScalarTraits<Scalar> STS;

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
    A_ (A), blockSize_ (blockSize), strideX_ (strideX), strideY_ (strideY)
  {
  }

  //! The block size (number of rows, and number of columns).
  LO getBlockSize () const {
    return blockSize_;
  }

  //! Pointer to the block's entries.
  Scalar* getRawPtr () const {
    return A_;
  }

  //! Reference to entry (i,j) of the block.
  Scalar& operator() (const LO i, const LO j) const {
    return A_[i * strideX_ + j * strideY_];
  }

  //! <tt>*this := *this + alpha * X</tt>.
  template<class LittleBlockType>
  void update (const Scalar& alpha, const LittleBlockType& X) const {
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) += alpha * X(i,j);
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
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) *= alpha;
      }
    }
  }

  //! <tt>(*this)(i,j) := alpha</tt> for all (i,j).
  void fill (const Scalar& alpha) const {
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i,j) = alpha;
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
        Scalar& Y_ij = (*this)(i,j);
        const Scalar X_ij = X(i,j);
        Y_ij = std::max (STS::magnitude (Y_ij), STS::magnitude (X_ij));
      }
    }
  }



  void factorize (int* ipiv, int & info)
  {
    Teuchos::LAPACK<int, Scalar> lapack;
    lapack.GETRF(blockSize_, blockSize_, A_, blockSize_, ipiv, &info);
  }

  template<class LittleVectorType>
  void solve (LittleVectorType & X, const int* ipiv) const
  {
    int info;
    Teuchos::LAPACK<int, Scalar> lapack;
    char trans = 'T';
    lapack.GETRS(trans, blockSize_, 1, A_, blockSize_, ipiv, X.getRawPtr(), blockSize_, &info);
  }

private:
  Scalar* const A_;
  const LO blockSize_;
  const LO strideX_;
  const LO strideY_;
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
private:
  typedef Teuchos::ScalarTraits<Scalar> STS;

public:
  /// \brief Constructor
  /// \param A [in] Pointer to the vector's entries
  /// \param blockSize [in] Dimension of the vector
  /// \param stride [in] Stride between consecutive entries
  LittleVector (Scalar* const A, const LO blockSize, const LO stride) :
    A_ (A), blockSize_ (blockSize), strideX_ (stride)
  {}

  //! Pointer to the block's entries.
  Scalar* getRawPtr () const {
    return A_;
  }

  //! The block size (number of degrees of freedom per mesh point).
  LO getBlockSize () const {
    return blockSize_;
  }

  //! Stride between consecutive entries.
  LO getStride () const {
    return strideX_;
  }

  //! Reference to entry (i) of the vector.
  Scalar& operator() (const LO i) const {
    return A_[i * strideX_];
  }

  //! <tt>*this := *this + alpha * X</tt>.
  template<class LittleVectorType>
  void update (const Scalar& alpha, const LittleVectorType& X) const {
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) += alpha * X(i);
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
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) *= alpha;
    }
  }

  //! <tt>(*this)(i,j) := alpha</tt> for all (i,j).
  void fill (const Scalar& alpha) const {
    for (LO i = 0; i < blockSize_; ++i) {
      (*this)(i) = alpha;
    }
  }

  /// \brief <tt>(*this)(i,j) := max(abs((*this)(i,j)), abs(X(i,j)))</tt>
  ///   for all (i,j).
  ///
  /// Tpetra uses this operation to implement the ABSMAX CombineMode.
  template<class LittleVectorType>
  void absmax (const LittleVectorType& X) const {
    for (LO i = 0; i < blockSize_; ++i) {
      Scalar& Y_i = (*this)(i);
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
    // FIXME (mfh 07 May 2014) This is suitable for column major, not
    // for row major.  Of course, we'll have to change other loops
    // above as well to make row major faster.
    for (LO j = 0; j < blockSize_; ++j) {
      for (LO i = 0; i < blockSize_; ++i) {
        (*this)(i) += alpha * A(i,j) * X(j);
      }
    }
  }

private:
  Scalar* const A_;
  const LO blockSize_;
  const LO strideX_;
};

} // namespace Experimental
} // namespace Tpetra

#endif // TPETRA_EXPERIMENTAL_BLOCKVIEW_HPP
