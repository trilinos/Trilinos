/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_LOCALSPARSETRIANGULARSOLVER_DEF_HPP
#define IFPACK2_LOCALSPARSETRIANGULARSOLVER_DEF_HPP

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"

namespace Ifpack2 {

template<class MatrixType>
LocalSparseTriangularSolver<MatrixType>::
LocalSparseTriangularSolver (const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A),
  isInitialized_ (false),
  isComputed_ (false),
  numInitialize_ (0),
  numCompute_ (0),
  numApply_ (0),
  initializeTime_ (0.0),
  computeTime_ (0.0),
  applyTime_ (0.0)
{
  typedef typename Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> crs_matrix_type;
  if (! A.is_null ()) {
    Teuchos::RCP<const crs_matrix_type> A_crs =
      Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A);
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_crs.is_null (), std::invalid_argument,
       "Ifpack2::LocalSparseTriangularSolver constructor: "
       "The input matrix A is not a Tpetra::CrsMatrix.");
    A_crs_ = A_crs;
  }
}

template<class MatrixType>
LocalSparseTriangularSolver<MatrixType>::
~LocalSparseTriangularSolver ()
{}

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
setParameters (const Teuchos::ParameterList& /*params*/)
{}

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
initialize ()
{
  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::initialize: ";

  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "You must call "
     "setMatrix() with a nonnull input matrix before you may call "
     "initialize() or compute().");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_crs_.is_null (), std::logic_error, prefix << "A_ is nonnull, but "
     "A_crs_ is null.  Please report this bug to the Ifpack2 developers.");
  auto G = A_->getGraph ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (G.is_null (), std::logic_error, prefix << "A_ and A_crs_ are nonnull, "
     "but A_'s RowGraph G is null.  "
     "Please report this bug to the Ifpack2 developers.");

  isInitialized_ = true;
  ++numInitialize_;
}

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
compute ()
{
  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::compute: ";

  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "You must call "
     "setMatrix() with a nonnull input matrix before you may call "
     "initialize() or compute().");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_crs_.is_null (), std::logic_error, prefix << "A_ is nonnull, but "
     "A_crs_ is null.  Please report this bug to the Ifpack2 developers.");

  if (! isInitialized_) {
    initialize ();
  }
  TEUCHOS_TEST_FOR_EXCEPTION
    (! isInitialized_, std::logic_error, prefix << "initialize() should have "
     "been called by this point, but isInitialized_ is false.  "
     "Please report this bug to the Ifpack2 developers.");

  isComputed_ = true;
  ++numCompute_;
}

template<class MatrixType>
void LocalSparseTriangularSolver<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type,
         global_ordinal_type, node_type>& X,
       Tpetra::MultiVector<scalar_type, local_ordinal_type,
         global_ordinal_type, node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef scalar_type ST;
  typedef Teuchos::ScalarTraits<ST> STS;
  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::apply: ";

  TEUCHOS_TEST_FOR_EXCEPTION
    (! isComputed (), std::runtime_error, prefix << "If compute() has not yet "
     "been called, or if you have changed the matrix via setMatrix(), you must "
     "call compute() before you may call this method.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::logic_error, prefix << "A_ is null.  "
     "Please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_crs_.is_null (), std::logic_error, prefix << "A_crs_ is null.  "
     "Please report this bug to the Ifpack2 developers.");
  auto G = A_->getGraph ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (G.is_null (), std::logic_error, prefix << "A_ and A_crs_ are nonnull, "
     "but A_'s RowGraph G is null.  "
     "Please report this bug to the Ifpack2 developers.");
  auto importer = G->getImporter ();
  auto exporter = G->getExporter ();

  if (! importer.is_null ()) {
    if (X_colMap_.is_null () || X_colMap_->getNumVectors () != X.getNumVectors ()) {
      X_colMap_ = rcp (new MV (importer->getTargetMap (), X.getNumVectors ()));
    }
    else {
      X_colMap_->putScalar (STS::zero ());
    }
    X_colMap_->doImport (X, *importer, Tpetra::ADD);
  }
  RCP<const MV> X_cur = importer.is_null () ? rcpFromRef (X) :
    Teuchos::rcp_const_cast<const MV> (X_colMap_);

  if (! exporter.is_null ()) {
    if (Y_rowMap_.is_null () || Y_rowMap_->getNumVectors () != Y.getNumVectors ()) {
      Y_rowMap_ = rcp (new MV (exporter->getSourceMap (), Y.getNumVectors ()));
    }
    else {
      Y_rowMap_->putScalar (STS::zero ());
    }
    Y_rowMap_->doExport (Y, *importer, Tpetra::ADD);
  }
  RCP<MV> Y_cur = exporter.is_null () ? rcpFromRef (Y) : Y_rowMap_;

  localApply (*X_cur, *Y_cur, mode, alpha, beta);

  if (! exporter.is_null ()) {
    Y.putScalar (STS::zero ());
    Y.doExport (*Y_cur, *exporter, Tpetra::ADD);
  }

  ++numApply_;
}

template<class MatrixType>
void
LocalSparseTriangularSolver<MatrixType>::
localApply (const MV& X,
            MV& Y,
            const Teuchos::ETransp mode,
            const scalar_type& alpha,
            const scalar_type& beta) const
{
  using Teuchos::RCP;
  typedef scalar_type ST;
  typedef Teuchos::ScalarTraits<ST> STS;

  if (beta == STS::zero ()) {
    if (alpha == STS::zero ()) {
      Y.putScalar (STS::zero ()); // Y := 0 * Y (ignore contents of Y)
    }
    else { // alpha != 0
      A_crs_->template localSolve<ST, ST> (X, Y, mode);
      if (alpha != STS::one ()) {
        Y.scale (alpha);
      }
    }
  }
  else { // beta != 0
    if (alpha == STS::zero ()) {
      Y.scale (beta); // Y := beta * Y
    }
    else { // alpha != 0
      MV Y_tmp (Y, Teuchos::Copy);
      A_crs_->template localSolve<ST, ST> (X, Y_tmp, mode); // Y_tmp := M * X
      Y.update (alpha, Y_tmp, beta); // Y := beta * Y + alpha * Y_tmp
    }
  }
}


template <class MatrixType>
int
LocalSparseTriangularSolver<MatrixType>::
getNumInitialize () const {
  return numInitialize_;
}

template <class MatrixType>
int
LocalSparseTriangularSolver<MatrixType>::
getNumCompute () const {
  return numCompute_;
}

template <class MatrixType>
int
LocalSparseTriangularSolver<MatrixType>::
getNumApply () const {
  return numApply_;
}

template <class MatrixType>
double
LocalSparseTriangularSolver<MatrixType>::
getInitializeTime () const {
  return initializeTime_;
}

template<class MatrixType>
double
LocalSparseTriangularSolver<MatrixType>::
getComputeTime () const {
  return computeTime_;
}

template<class MatrixType>
double
LocalSparseTriangularSolver<MatrixType>::
getApplyTime() const {
  return applyTime_;
}

template <class MatrixType>
std::string
LocalSparseTriangularSolver<MatrixType>::
description () const
{
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::LocalSparseTriangularSolver\": {";
  if (this->getObjectLabel () != "") {
    os << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  os << "Initialized: " << (isInitialized () ? "true" : "false") << ", "
     << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  if (A_.is_null ()) {
    os << "Matrix: null";
  }
  else {
    os << "Matrix: not null"
       << ", Global matrix dimensions: ["
       << A_->getGlobalNumRows () << ", "
       << A_->getGlobalNumCols () << "]";
  }

  os << "}";
  return os.str ();
}

template <class MatrixType>
void LocalSparseTriangularSolver<MatrixType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  // Default verbosity level is VERB_LOW, which prints only on Process
  // 0 of the matrix's communicator.
  const Teuchos::EVerbosityLevel vl
    = (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;

  if (vl != Teuchos::VERB_NONE) {
    // Print only on Process 0 in the matrix's communicator.  If the
    // matrix is null, though, we have to get the communicator from
    // somewhere, so we ask Tpetra for its default communicator.  If
    // MPI is enabled, this wraps MPI_COMM_WORLD or a clone thereof.
    auto comm = A_.is_null () ?
      Tpetra::DefaultPlatform::getDefaultPlatform ().getComm () :
      A_->getComm ();

    // Users aren't supposed to do anything with the matrix on
    // processes where its communicator is null.
    if (! comm.is_null () && comm->getRank () == 0) {
      // By convention, describe() should always begin with a tab.
      Teuchos::OSTab tab0 (out);
      // Output is in YAML format.  We have to escape the class name,
      // because it has a colon.
      out << "\"Ifpack2::LocalSparseTriangularSolver\":" << endl;
      Teuchos::OSTab tab1 (out);
      out << "Scalar: " << Teuchos::TypeNameTraits<scalar_type>::name () << endl
          << "LocalOrdinal: " << Teuchos::TypeNameTraits<local_ordinal_type>::name () << endl
          << "GlobalOrdinal: " << Teuchos::TypeNameTraits<global_ordinal_type>::name () << endl
          << "Node: " << Teuchos::TypeNameTraits<node_type>::name () << endl;
    }
  }
}

template <class MatrixType>
Teuchos::RCP<const typename LocalSparseTriangularSolver<MatrixType>::map_type>
LocalSparseTriangularSolver<MatrixType>::
getDomainMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error,
     "Ifpack2::LocalSparseTriangularSolver::getDomainMap: "
     "The matrix is null.  Please call setMatrix() with a nonnull input "
     "before calling this method.");
  return A_->getDomainMap ();
}

template <class MatrixType>
Teuchos::RCP<const typename LocalSparseTriangularSolver<MatrixType>::map_type>
LocalSparseTriangularSolver<MatrixType>::
getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error,
     "Ifpack2::LocalSparseTriangularSolver::getRangeMap: "
     "The matrix is null.  Please call setMatrix() with a nonnull input "
     "before calling this method.");
  return A_->getRangeMap ();
}

template<class MatrixType>
void LocalSparseTriangularSolver<MatrixType>::
setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  const char prefix[] = "Ifpack2::LocalSparseTriangularSolver::setMatrix: ";

  // Check in serial or one-process mode if the matrix is square.
  TEUCHOS_TEST_FOR_EXCEPTION
    (! A.is_null () && A->getComm ()->getSize () == 1 &&
     A->getNodeNumRows () != A->getNodeNumCols (),
     std::runtime_error, prefix << "If A's communicator only contains one "
     "process, then A must be square.  Instead, you provided a matrix A with "
     << A->getNodeNumRows () << " rows and " << A->getNodeNumCols ()
     << " columns.");

  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates the preconditioner.
  isInitialized_ = false;
  isComputed_ = false;

  typedef typename Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> crs_matrix_type;
  if (A.is_null ()) {
    A_crs_ = Teuchos::null;
  }
  else {
    Teuchos::RCP<const crs_matrix_type> A_crs =
      Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A);
    TEUCHOS_TEST_FOR_EXCEPTION
      (A_crs.is_null (), std::invalid_argument, prefix <<
       "The input matrix A is not a Tpetra::CrsMatrix.");
    A_crs_ = A_crs;
  }
  A_ = A;

  TEUCHOS_TEST_FOR_EXCEPTION
    ((A.is_null () && A_.is_null () && A_crs_.is_null ()) ||
     (A_.getRawPtr () == A.getRawPtr () && ! A_crs_.is_null ()),
     std::logic_error, prefix << "This class' matrix pointers were set "
     "incorrectly.  Please report this bug to the Ifpack2 developers.");
}

} // namespace Ifpack2

#define IFPACK2_LOCALSPARSETRIANGULARSOLVER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::LocalSparseTriangularSolver< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_LOCALSPARSETRIANGULARSOLVER_DEF_HPP
