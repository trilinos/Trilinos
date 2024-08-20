// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_IDENTITY_SOLVER_DEF_HPP
#define IFPACK2_IDENTITY_SOLVER_DEF_HPP

#include "Ifpack2_IdentitySolver_decl.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Export.hpp"

namespace Ifpack2 {

template<class MatrixType>
IdentitySolver<MatrixType>::
IdentitySolver (const Teuchos::RCP<const row_matrix_type>& A)
  : matrix_ (A),
    isInitialized_ (false),
    isComputed_ (false),
    numInitialize_ (0),
    numCompute_ (0),
    numApply_ (0),
    initializeTime_(0.0),
    computeTime_(0.0),
    applyTime_(0.0)
{
}

template<class MatrixType>
IdentitySolver<MatrixType>::~IdentitySolver ()
{
}

template<class MatrixType>
void IdentitySolver<MatrixType>::setParameters (const Teuchos::ParameterList& /*params*/)
{
}

template<class MatrixType>
void IdentitySolver<MatrixType>::initialize ()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    matrix_.is_null (), std::runtime_error, "Ifpack2::IdentitySolver: "
    "You must call setMatrix() with a nonnull input matrix "
    "before you may call initialize() or compute().");

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! matrix_->getDomainMap ()->isCompatible (* (matrix_->getRangeMap ())),
    std::invalid_argument,
    "Ifpack2::IdentitySolver: The domain and range Maps "
    "of the input matrix must be compatible.");

  // If the domain and range Maps are not the same, then we need to
  // construct an Export from the domain Map to the range Map, so that
  // this operator is really the identity and not a permutation.
  if (! matrix_->getDomainMap ()->isSameAs (* (matrix_->getRangeMap ()))) {
    export_ = Teuchos::rcp (new export_type (matrix_->getDomainMap (),
                                             matrix_->getRangeMap ()));
  }
  else {
    // If the Export is null, we won't do the Export in apply().
    // Thus, we need to set it to null here as a flag.
    export_ = Teuchos::null;
  }

  isInitialized_ = true;
  ++numInitialize_;
}

template<class MatrixType>
void IdentitySolver<MatrixType>::compute ()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    matrix_.is_null (), std::runtime_error, "Ifpack2::IdentitySolver: "
    "You must call setMatrix() with a nonnull input matrix "
    "before you may call initialize() or compute().");

  if (! isInitialized_) {
    initialize ();
  }

  isComputed_ = true;
  ++numCompute_;
}

template<class MatrixType>
void IdentitySolver<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp /*mode*/,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::RCP;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> MV;

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (), std::runtime_error,
    "Ifpack2::IdentitySolver::apply: If compute() has not yet been called, "
    "or if you have changed the matrix via setMatrix(), "
    "you must call compute() before you may call this method.");

  // "Identity solver" does what it says: it's the identity operator.
  // We have to Export if the domain and range Maps are not the same.
  // Otherwise, this operator would be a permutation, not the identity.
  if (export_.is_null ()) {
    Y.update (alpha, X, beta);
  }
  else {
    if (alpha == STS::one () && beta == STS::zero ()) { // the common case
      Y.doExport (X, *export_, Tpetra::REPLACE);
    }
    else {
      // We know that the domain and range Maps are compatible.  First
      // bring X into the range Map via Export.  Then compute in place
      // in Y.
      MV X_tmp (Y.getMap (), Y.getNumVectors ());
      X_tmp.doExport (X, *export_, Tpetra::REPLACE);
      Y.update (alpha, X_tmp, beta);
    }
  }
  ++numApply_;
}

template <class MatrixType>
int IdentitySolver<MatrixType>::getNumInitialize() const {
  return(numInitialize_);
}

template <class MatrixType>
int IdentitySolver<MatrixType>::getNumCompute() const {
  return(numCompute_);
}

template <class MatrixType>
int IdentitySolver<MatrixType>::getNumApply() const {
  return(numApply_);
}

template <class MatrixType>
double IdentitySolver<MatrixType>::getInitializeTime() const {
  return(initializeTime_);
}

template<class MatrixType>
double IdentitySolver<MatrixType>::getComputeTime() const {
  return(computeTime_);
}

template<class MatrixType>
double IdentitySolver<MatrixType>::getApplyTime() const {
  return(applyTime_);
}

template <class MatrixType>
std::string IdentitySolver<MatrixType>::description () const
{
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::IdentitySolver\": {";
  if (this->getObjectLabel () != "") {
    os << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  os << "Initialized: " << (isInitialized () ? "true" : "false") << ", "
     << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  if (matrix_.is_null ()) {
    os << "Matrix: null";
  }
  else {
    os << "Matrix: not null"
       << ", Global matrix dimensions: ["
       << matrix_->getGlobalNumRows () << ", "
       << matrix_->getGlobalNumCols () << "]";
  }

  os << "}";
  return os.str ();
}

template <class MatrixType>
void IdentitySolver<MatrixType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  const Teuchos::EVerbosityLevel vl
    = (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;

  if (vl != Teuchos::VERB_NONE) {
    // By convention, describe() should always begin with a tab.
    Teuchos::OSTab tab0 (out);
    out << "\"Ifpack2::IdentitySolver\":" << endl;
    Teuchos::OSTab tab1 (out);
    out << "MatrixType: " << Teuchos::TypeNameTraits<MatrixType>::name () << endl;
    out << "numInitialize: " << numInitialize_ << endl;
    out << "numCompute: " << numCompute_ << endl;
    out << "numApply: " << numApply_ << endl;
  }
}

template <class MatrixType>
Teuchos::RCP<const typename IdentitySolver<MatrixType>::map_type> IdentitySolver<MatrixType>::getDomainMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    matrix_.is_null (), std::runtime_error, "Ifpack2::IdentitySolver::getDomainMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return matrix_->getDomainMap ();
}

template <class MatrixType>
Teuchos::RCP<const typename IdentitySolver<MatrixType>::map_type> IdentitySolver<MatrixType>::getRangeMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    matrix_.is_null (), std::runtime_error, "Ifpack2::IdentitySolver::getRangeMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return matrix_->getRangeMap ();
}

template<class MatrixType>
void IdentitySolver<MatrixType>::
setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // Check in serial or one-process mode if the matrix is square.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A.is_null () && A->getComm ()->getSize () == 1 &&
    A->getLocalNumRows () != A->getLocalNumCols (),
    std::runtime_error, "Ifpack2::IdentitySolver::setMatrix: If A's communicator only "
    "contains one process, then A must be square.  Instead, you provided a "
    "matrix A with " << A->getLocalNumRows () << " rows and "
    << A->getLocalNumCols () << " columns.");

  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates the preconditioner.
  isInitialized_ = false;
  isComputed_ = false;
  export_ = Teuchos::null;

  matrix_ = A;
}

} // namespace Ifpack2

#define IFPACK2_IDENTITYSOLVER_INSTANT(S,LO,GO,N)                            \
  template class Ifpack2::IdentitySolver< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_IDENTITY_SOLVER_DEF_HPP
