// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DIAGONAL_DEF_HPP
#define IFPACK2_DIAGONAL_DEF_HPP

#include "Ifpack2_Diagonal_decl.hpp"
#include "Tpetra_CrsMatrix.hpp"

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
  offsets_ = offsets_type ();
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
      offsets_ = offsets_type (); // offsets are no longer valid
    }
    else {
      const size_t lclNumRows = A_crs->getLocalNumRows ();
      if (offsets_.extent (0) < lclNumRows) {
        offsets_ = offsets_type (); // clear first to save memory
        offsets_ = offsets_type ("offsets", lclNumRows);
      }
      A_crs->getCrsGraph ()->getLocalDiagOffsets (offsets_);
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
      A_crs->getLocalDiagCopy (*tmpVec, offsets_);
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
  template class Ifpack2::Diagonal< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif
