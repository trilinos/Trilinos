/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_DIAGONAL_DEF_HPP
#define IFPACK2_DIAGONAL_DEF_HPP

#include "Ifpack2_Diagonal_decl.hpp"
#include "Tpetra_CrsMatrix_def.hpp"
#include "Ifpack2_Condest.hpp"

namespace Ifpack2 {

template<class MatrixType>
Diagonal<MatrixType>::Diagonal (const Teuchos::RCP<const row_matrix_type>& A) :
  matrix_ (A),
  initializeTime_ (0.0),
  computeTime_ (0.0),
  applyTime_ (0.0),
  numInitialize_ (0),
  numCompute_ (0),
  numApply_ (0),
  condEst_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  isInitialized_ (false),
  isComputed_ (false)
{}

template<class MatrixType>
Diagonal<MatrixType>::Diagonal (const Teuchos::RCP<const crs_matrix_type>& A) :
  matrix_ (A),
  initializeTime_ (0.0),
  computeTime_ (0.0),
  applyTime_ (0.0),
  numInitialize_ (0),
  numCompute_ (0),
  numApply_ (0),
  condEst_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  isInitialized_ (false),
  isComputed_ (false)
{}

template<class MatrixType>
Diagonal<MatrixType>::Diagonal (const Teuchos::RCP<const vector_type>& diag) :
  userInverseDiag_ (diag),
  inverseDiag_ (diag),
  initializeTime_ (0.0),
  computeTime_ (0.0),
  applyTime_ (0.0),
  numInitialize_ (0),
  numCompute_ (0),
  numApply_ (0),
  condEst_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  isInitialized_ (false),
  isComputed_ (false)
{}

template<class MatrixType>
Diagonal<MatrixType>::~Diagonal ()
{}

template<class MatrixType>
Teuchos::RCP<const typename Diagonal<MatrixType>::map_type>
Diagonal<MatrixType>::getDomainMap () const
{
  if (matrix_.is_null ()) {
    if (userInverseDiag_.is_null ()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Ifpack2::Diagonal::getDomainMap: "
        "The input matrix A is null, and you did not provide a vector of "
        "inverse diagonal entries.  Please call setMatrix() with a nonnull "
        "input matrix before calling this method.");
    } else {
      return userInverseDiag_->getMap ();
    }
  } else {
    return matrix_->getDomainMap ();
  }
}

template<class MatrixType>
Teuchos::RCP<const typename Diagonal<MatrixType>::map_type>
Diagonal<MatrixType>::getRangeMap () const
{
  if (matrix_.is_null ()) {
    if (userInverseDiag_.is_null ()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::runtime_error, "Ifpack2::Diagonal::getRangeMap: "
        "The input matrix A is null, and you did not provide a vector of "
        "inverse diagonal entries.  Please call setMatrix() with a nonnull "
        "input matrix before calling this method.");
    } else {
      return userInverseDiag_->getMap ();
    }
  } else {
    return matrix_->getRangeMap ();
  }
}

template<class MatrixType>
void Diagonal<MatrixType>::
setParameters (const Teuchos::ParameterList& /*params*/)
{}

template<class MatrixType>
void Diagonal<MatrixType>::reset ()
{
  inverseDiag_ = Teuchos::null;
  offsets_ = Teuchos::null;
  condEst_ = -Teuchos::ScalarTraits<magnitude_type>::one ();
  isInitialized_ = false;
  isComputed_ = false;
}

template<class MatrixType>
void Diagonal<MatrixType>::
setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A.getRawPtr () != matrix_.getRawPtr ()) { // it's a different matrix
    reset ();
    matrix_ = A;
  }
}

template<class MatrixType>
void Diagonal<MatrixType>::initialize ()
{
  // Either the matrix to precondition must be nonnull, or the user
  // must have provided a Vector of diagonal entries.
  TEUCHOS_TEST_FOR_EXCEPTION(
    matrix_.is_null () && userInverseDiag_.is_null (), std::runtime_error,
    "Ifpack2::Diagonal::initialize: The matrix to precondition is null, "
    "and you did not provide a Tpetra::Vector of diagonal entries.  "
    "Please call setMatrix() with a nonnull input before calling this method.");

  // If the user did provide an input matrix, then that takes
  // precedence over the vector of inverse diagonal entries, if they
  // provided one earlier.  This is only possible if they created this
  // Diagonal instance using the constructor that takes a
  // Tpetra::Vector pointer, and then called setMatrix() with a
  // nonnull input matrix.
  if (! matrix_.is_null ()) {
    // If you call initialize(), it means that you are asserting that
    // the structure of the input sparse matrix may have changed.
    // This means we should always recompute the diagonal offsets, if
    // the input matrix is a Tpetra::CrsMatrix.
    Teuchos::RCP<const crs_matrix_type> A_crs =
      Teuchos::rcp_dynamic_cast<const crs_matrix_type> (matrix_);
    if (A_crs.is_null ()) {
      offsets_ = Teuchos::null; // offsets are no longer valid
    }
    else {
      A_crs->getLocalDiagOffsets (offsets_);
    }
  }

  isInitialized_ = true;
  ++numInitialize_;
}

template<class MatrixType>
void Diagonal<MatrixType>::compute ()
{
  // Either the matrix to precondition must be nonnull, or the user
  // must have provided a Vector of diagonal entries.
  TEUCHOS_TEST_FOR_EXCEPTION(
    matrix_.is_null () && userInverseDiag_.is_null (), std::runtime_error,
    "Ifpack2::Diagonal::compute: The matrix to precondition is null, "
    "and you did not provide a Tpetra::Vector of diagonal entries.  "
    "Please call setMatrix() with a nonnull input before calling this method.");

  if (! isInitialized_) {
    initialize ();
  }

  // If the user did provide an input matrix, then that takes
  // precedence over the vector of inverse diagonal entries, if they
  // provided one earlier.  This is only possible if they created this
  // Diagonal instance using the constructor that takes a
  // Tpetra::Vector pointer, and then called setMatrix() with a
  // nonnull input matrix.
  if (matrix_.is_null ()) { // accept the user's diagonal
    inverseDiag_ = userInverseDiag_;
  }
  else {
    Teuchos::RCP<vector_type> tmpVec (new vector_type (matrix_->getRowMap ()));
    Teuchos::RCP<const crs_matrix_type> A_crs =
      Teuchos::rcp_dynamic_cast<const crs_matrix_type> (matrix_);
    if (A_crs.is_null ()) {
      // Get the diagonal entries from the Tpetra::RowMatrix.
      matrix_->getLocalDiagCopy (*tmpVec);
    }
    else {
      // Get the diagonal entries from the Tpetra::CrsMatrix using the
      // precomputed offsets.
      A_crs->getLocalDiagCopy (*tmpVec, offsets_ ());
    }
    tmpVec->reciprocal (*tmpVec); // invert the diagonal entries
    inverseDiag_ = tmpVec;
  }

  isComputed_ = true;
  ++numCompute_;
}

template<class MatrixType>
void Diagonal<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp /*mode*/,
       scalar_type alpha,
       scalar_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (), std::runtime_error, "Ifpack2::Diagonal::apply: You "
    "must first call compute() before you may call apply().  Once you have "
    "called compute(), you need not call it again unless the values in the "
    "matrix have changed, or unless you have called setMatrix().");

  // FIXME (mfh 12 Sep 2014) This assumes that row Map == range Map ==
  // domain Map.  If the preconditioner has a matrix, we should ask
  // the matrix whether we need to do an Import before and/or an
  // Export after.

  Y.elementWiseMultiply (alpha, *inverseDiag_, X, beta);
  ++numApply_;
}

template<class MatrixType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
Diagonal<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const row_matrix_type>& matrix)
{
  const magnitude_type minusOne = -Teuchos::ScalarTraits<magnitude_type>::one ();

  if (! isComputed ()) { // cannot compute right now
    return minusOne;
  }
  // NOTE: this is computing the *local* condest
  if (condEst_ == minusOne) {
    condEst_ = Ifpack2::Condest (*this, CT, MaxIters, Tol, matrix);
  }
  return condEst_;
}

template <class MatrixType>
int Diagonal<MatrixType>::getNumInitialize() const {
  return numInitialize_;
}

template <class MatrixType>
int Diagonal<MatrixType>::getNumCompute() const {
  return numCompute_;
}

template <class MatrixType>
int Diagonal<MatrixType>::getNumApply() const {
  return numApply_;
}

template <class MatrixType>
double Diagonal<MatrixType>::getInitializeTime() const {
  return initializeTime_;
}

template<class MatrixType>
double Diagonal<MatrixType>::getComputeTime() const {
  return computeTime_;
}

template<class MatrixType>
double Diagonal<MatrixType>::getApplyTime() const {
  return applyTime_;
}

template <class MatrixType>
std::string Diagonal<MatrixType>::description () const
{
  std::ostringstream out;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  out << "\"Ifpack2::Diagonal\": "
      << "{";
  if (this->getObjectLabel () != "") {
    out << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  if (matrix_.is_null ()) {
    out << "Matrix: null";
  }
  else {
    out << "Matrix: not null"
        << ", Global matrix dimensions: ["
        << matrix_->getGlobalNumRows () << ", "
        << matrix_->getGlobalNumCols () << "]";
  }

  out << "}";
  return out.str ();
}

template <class MatrixType>
void Diagonal<MatrixType>::
describe (Teuchos::FancyOStream &out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;
  if (vl != Teuchos::VERB_NONE) {
    Teuchos::OSTab tab0 (out);
    out << "\"Ifpack2::Diagonal\":";
    Teuchos::OSTab tab1 (out);
    out << "Template parameter: "
        << Teuchos::TypeNameTraits<MatrixType>::name () << endl;
    if (this->getObjectLabel () != "") {
      out << "Label: \"" << this->getObjectLabel () << "\", ";
    }
    out << "Number of initialize calls: " << numInitialize_ << endl
        << "Number of compute calls: " << numCompute_ << endl
        << "Number of apply calls: " << numApply_ << endl;
  }
}

} // namespace Ifpack2

#define IFPACK2_DIAGONAL_INSTANT(S,LO,GO,N)                            \
  template class Ifpack2::Diagonal< Tpetra::CrsMatrix<S, LO, GO, N> >; \
  template class Ifpack2::Diagonal< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif
