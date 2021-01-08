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

#ifndef IFPACK2_DETAILS_TRIDISOLVER_DEF_HPP
#define IFPACK2_DETAILS_TRIDISOLVER_DEF_HPP

#include "Ifpack2_LocalFilter.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"

#ifdef HAVE_MPI
#  include <mpi.h>
#  include "Teuchos_DefaultMpiComm.hpp"
#else
#  include "Teuchos_DefaultSerialComm.hpp"
#endif // HAVE_MPI


namespace Ifpack2 {
namespace Details {

//////////////////////////////////////////////////////////////////////
// Non-stub (full) implementation
//////////////////////////////////////////////////////////////////////

template<class MatrixType>
TriDiSolver<MatrixType, false>::
TriDiSolver (const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A),
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
Teuchos::RCP<const typename TriDiSolver<MatrixType, false>::map_type>
TriDiSolver<MatrixType, false>::getDomainMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Details::TriDiSolver::"
    "getDomainMap: The input matrix A is null.  Please call setMatrix() with a "
    "nonnull input matrix before calling this method.");
  // For an input matrix A, TriDiSolver solves Ax=b for x.
  // Thus, its Maps are reversed from those of the input matrix.
  return A_->getRangeMap ();
}


template<class MatrixType>
Teuchos::RCP<const typename TriDiSolver<MatrixType, false>::map_type>
TriDiSolver<MatrixType, false>::getRangeMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Details::TriDiSolver::"
    "getRangeMap: The input matrix A is null.  Please call setMatrix() with a "
    "nonnull input matrix before calling this method.");
  // For an input matrix A, TriDiSolver solves Ax=b for x.
  // Thus, its Maps are reversed from those of the input matrix.
  return A_->getDomainMap ();
}


template<class MatrixType>
void
TriDiSolver<MatrixType, false>::
setParameters (const Teuchos::ParameterList& params) {
  (void) params; // this preconditioner doesn't currently take any parameters
}


template<class MatrixType>
bool
TriDiSolver<MatrixType, false>::isInitialized () const {
  return isInitialized_;
}


template<class MatrixType>
bool
TriDiSolver<MatrixType, false>::isComputed () const {
  return isComputed_;
}


template<class MatrixType>
int
TriDiSolver<MatrixType, false>::getNumInitialize () const {
  return numInitialize_;
}


template<class MatrixType>
int
TriDiSolver<MatrixType, false>::getNumCompute () const {
  return numCompute_;
}


template<class MatrixType>
int
TriDiSolver<MatrixType, false>::getNumApply () const {
  return numApply_;
}


template<class MatrixType>
double
TriDiSolver<MatrixType, false>::getInitializeTime () const {
  return initializeTime_;
}


template<class MatrixType>
double
TriDiSolver<MatrixType, false>::getComputeTime () const {
  return computeTime_;
}


template<class MatrixType>
double
TriDiSolver<MatrixType, false>::getApplyTime () const {
  return applyTime_;
}


template<class MatrixType>
Teuchos::RCP<const typename TriDiSolver<MatrixType, false>::row_matrix_type>
TriDiSolver<MatrixType, false>::getMatrix () const {
  return A_;
}


template<class MatrixType>
void TriDiSolver<MatrixType, false>::
reset ()
{
  isInitialized_ = false;
  isComputed_ = false;
  A_local_ = Teuchos::null;
  A_local_tridi_.reshape (0);
  ipiv_.resize (0);
}


template<class MatrixType>
void TriDiSolver<MatrixType, false>::
setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // It's legitimate to call setMatrix() with a null input.  This has
  // the effect of resetting the preconditioner's internal state.
  if (! A_.is_null ()) {
    const global_size_t numRows = A->getRangeMap ()->getGlobalNumElements ();
    const global_size_t numCols = A->getDomainMap ()->getGlobalNumElements ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      numRows != numCols, std::invalid_argument, "Ifpack2::Details::TriDiSolver::"
      "setMatrix: Input matrix must be (globally) square.  "
      "The matrix you provided is " << numRows << " by " << numCols << ".");
  }
  // Clear any previously computed objects.
  reset ();

  // Now that we've cleared the state, we can keep the matrix.
  A_ = A;
}


template<class MatrixType>
void TriDiSolver<MatrixType, false>::initialize ()
{
  using Teuchos::Comm;
  using Teuchos::null;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  const std::string timerName ("Ifpack2::Details::TriDiSolver::initialize");

  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  double startTime = timer->wallTime();

  { // Begin timing here.
    Teuchos::TimeMonitor timeMon (*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null (), std::runtime_error, "Ifpack2::Details::TriDiSolver::"
      "initialize: The input matrix A is null.  Please call setMatrix() "
      "with a nonnull input before calling this method.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! A_->hasColMap (), std::invalid_argument, "Ifpack2::Details::TriDiSolver: "
      "The constructor's input matrix must have a column Map, "
      "so that it has local indices.");

    // Clear any previously computed objects.
    reset ();

    // Make the local filter of the input matrix A.
    if (A_->getComm ()->getSize () > 1) {
      A_local_ = rcp (new LocalFilter<row_matrix_type> (A_));
    } else {
      A_local_ = A_;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      A_local_.is_null (), std::logic_error, "Ifpack2::Details::TriDiSolver::"
      "initialize: A_local_ is null after it was supposed to have been "
      "initialized.  Please report this bug to the Ifpack2 developers.");

    // Allocate the TriDi local matrix and the pivot array.
    const size_t numRows = A_local_->getNodeNumRows ();
    const size_t numCols = A_local_->getNodeNumCols ();
    TEUCHOS_TEST_FOR_EXCEPTION(
      numRows != numCols, std::logic_error, "Ifpack2::Details::TriDiSolver::"
      "initialize: Local filter matrix is not square.  This should never happen.  "
      "Please report this bug to the Ifpack2 developers.");
    A_local_tridi_.reshape (numRows);
    ipiv_.resize (numRows);
    std::fill (ipiv_.begin (), ipiv_.end (), 0);

    isInitialized_ = true;
    ++numInitialize_;
  }

  initializeTime_ += (timer->wallTime() - startTime);
}


template<class MatrixType>
TriDiSolver<MatrixType, false>::~TriDiSolver()
{}


template<class MatrixType>
void TriDiSolver<MatrixType, false>::compute ()
{
  using Teuchos::RCP;
  const std::string timerName ("Ifpack2::Details::TriDiSolver::compute");

  RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }

  double startTime = timer->wallTime();

  // Begin timing here.
  {
    Teuchos::TimeMonitor timeMon (*timer);
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null (), std::runtime_error, "Ifpack2::Details::TriDiSolver::"
      "compute: The input matrix A is null.  Please call setMatrix() with a "
      "nonnull input, then call initialize(), before calling this method.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      A_local_.is_null (), std::logic_error, "Ifpack2::Details::TriDiSolver::"
      "compute: A_local_ is null.  Please report this bug to the Ifpack2 "
      "developers.");

    isComputed_ = false;
    if (! this->isInitialized ()) {
      this->initialize ();
    }
    extract (A_local_tridi_, *A_local_); // extract the tridi local matrix

    factor (A_local_tridi_, ipiv_ ()); // factor the tridi local matrix

    isComputed_ = true;
    ++numCompute_;
  }
  computeTime_ += (timer->wallTime() - startTime);
}

template<class MatrixType>
void TriDiSolver<MatrixType, false>::factor (Teuchos::SerialTriDiMatrix<int, scalar_type>& A,
        const Teuchos::ArrayView<int>& ipiv)
{
  // Fill the LU permutation array with zeros.
  std::fill (ipiv.begin (), ipiv.end (), 0);

  Teuchos::LAPACK<int, scalar_type> lapack;
  int INFO = 0;
  lapack.GTTRF (A.numRowsCols (),
                A.DL(),
                A.D(),
                A.DU(),
                A.DU2(),
                ipiv.getRawPtr (), &INFO);
  // INFO < 0 is a bug.
  TEUCHOS_TEST_FOR_EXCEPTION(
    INFO < 0, std::logic_error, "Ifpack2::Details::TriDiSolver::factor: "
    "LAPACK's _GTTRF (tridiagonal LU factorization with partial pivoting) "
    "was called incorrectly.  INFO = " << INFO << " < 0.  "
    "Please report this bug to the Ifpack2 developers.");
  // INFO > 0 means the matrix is singular.  This is probably an issue
  // either with the choice of rows the rows we extracted, or with the
  // input matrix itself.
  TEUCHOS_TEST_FOR_EXCEPTION(
    INFO > 0, std::runtime_error, "Ifpack2::Details::TriDiSolver::factor: "
    "LAPACK's _GTTRF (tridiagonal LU factorization with partial pivoting) "
    "reports that the computed U factor is exactly singular.  U(" << INFO <<
    "," << INFO << ") (one-based index i) is exactly zero.  This probably "
    "means that the input matrix has a singular diagonal block.");
}


template<class MatrixType>
void TriDiSolver<MatrixType, false>::
applyImpl (const MV& X,
           MV& Y,
           const Teuchos::ETransp mode,
           const scalar_type alpha,
           const scalar_type beta) const
{
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::CONJ_TRANS;
  using Teuchos::TRANS;

  const int numVecs = static_cast<int> (X.getNumVectors ());
  if (alpha == STS::zero ()) { // don't need to solve the linear system
    if (beta == STS::zero ()) {
      // Use BLAS AXPY semantics for beta == 0: overwrite, clobbering
      // any Inf or NaN values in Y (rather than multiplying them by
      // zero, resulting in NaN values).
      Y.putScalar (STS::zero ());
    }
    else { // beta != 0
      Y.scale (STS::zero ());
    }
  }
  else { // alpha != 0; must solve the linear system
    Teuchos::LAPACK<int, scalar_type> lapack;
    // If beta is nonzero, Y is not constant stride, or alpha != 1, we
    // have to use a temporary output multivector Y_tmp.  It gets a
    // copy of alpha*X, since GETRS overwrites its (multi)vector input
    // with its output.
    RCP<MV> Y_tmp;
    if (beta == STS::zero () && Y.isConstantStride () && alpha == STS::one ()) {
      deep_copy(Y, X);
      Y_tmp = rcpFromRef (Y);
    }
    else {
      Y_tmp = rcp (new MV (createCopy(X))); // constructor copies X
      if (alpha != STS::one ()) {
        Y_tmp->scale (alpha);
      }
    }
    const int Y_stride = static_cast<int> (Y_tmp->getStride ());
    ArrayRCP<scalar_type> Y_view = Y_tmp->get1dViewNonConst ();
    scalar_type* const Y_ptr = Y_view.getRawPtr ();
    int INFO = 0;
    const char trans =
      (mode == CONJ_TRANS ? 'C' : (mode == TRANS ? 'T' : 'N'));
    lapack.GTTRS (trans, A_local_tridi_.numRowsCols(), numVecs,
                  A_local_tridi_.DL(),
                  A_local_tridi_.D(),
                  A_local_tridi_.DU(),
                  A_local_tridi_.DU2(),
                  ipiv_.getRawPtr (), Y_ptr, Y_stride, &INFO);
    TEUCHOS_TEST_FOR_EXCEPTION(
      INFO != 0, std::runtime_error, "Ifpack2::Details::TriDiSolver::"
      "applyImpl: LAPACK's _GTTRS (tridiagonal solve using LU factorization "
      "with partial pivoting) failed with INFO = " << INFO << " != 0.");

    if (beta != STS::zero ()) {
      Y.update (alpha, *Y_tmp, beta);
    }
    else if (! Y.isConstantStride ()) {
      deep_copy(Y, *Y_tmp);
    }
  }
}


template<class MatrixType>
void TriDiSolver<MatrixType, false>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

  const std::string timerName ("Ifpack2::Details::TriDiSolver::apply");
  RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }

  double startTime = timer->wallTime();

  // Begin timing here.
  {
    Teuchos::TimeMonitor timeMon (*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isComputed_, std::runtime_error, "Ifpack2::Details::TriDiSolver::apply: "
      "You must have called the compute() method before you may call apply().  "
      "You may call the apply() method as many times as you want after calling "
      "compute() once, but you must have called compute() at least once.");

    const size_t numVecs = X.getNumVectors ();

    TEUCHOS_TEST_FOR_EXCEPTION(
      numVecs != Y.getNumVectors (), std::runtime_error,
      "Ifpack2::Details::TriDiSolver::apply: X and Y have different numbers "
      "of vectors.  X has " << X.getNumVectors () << ", but Y has "
      << X.getNumVectors () << ".");

    if (numVecs == 0) {
      return; // done! nothing to do
    }

    // Set up "local" views of X and Y.
    RCP<const MV> X_local;
    RCP<MV> Y_local;
    const bool multipleProcs = (A_->getRowMap ()->getComm ()->getSize () >= 1);
    if (multipleProcs) {
      // Interpret X and Y as "local" multivectors, that is, in the
      // local filter's domain resp. range Maps.  "Interpret" means that
      // we create views with different Maps; we don't have to copy.
      X_local = X.offsetView (A_local_->getDomainMap (), 0);
      Y_local = Y.offsetViewNonConst (A_local_->getRangeMap (), 0);
    }
    else { // only one process in A_'s communicator
      // X and Y are already "local"; no need to set up local views.
      X_local = rcpFromRef (X);
      Y_local = rcpFromRef (Y);
    }

    // Apply the local operator:
    // Y_local := beta*Y_local + alpha*M^{-1}*X_local
    this->applyImpl (*X_local, *Y_local, mode, alpha, beta);

    ++numApply_; // We've successfully finished the work of apply().
  }

  applyTime_ += (timer->wallTime() - startTime);
}


template<class MatrixType>
std::string TriDiSolver<MatrixType, false>::description () const
{
  std::ostringstream out;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  out << "\"Ifpack2::Details::TriDiSolver\": ";
  out << "{";
  if (this->getObjectLabel () != "") {
    out << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  out << "Initialized: " << (isInitialized () ? "true" : "false") << ", "
      << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  if (A_.is_null ()) {
    out << "Matrix: null";
  }
  else {
    out << "Matrix: not null"
        << ", Global matrix dimensions: ["
        << A_->getGlobalNumRows () << ", " << A_->getGlobalNumCols () << "]";
  }

  out << "}";
  return out.str ();
}


template<class MatrixType>
void TriDiSolver<MatrixType, false>::describeLocal (Teuchos::FancyOStream& out,
                                                    const Teuchos::EVerbosityLevel verbLevel) const {
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using std::endl;

  if (verbLevel == Teuchos::VERB_NONE) {
    return;
  }
  else {
    RCP<FancyOStream> ptrOut = rcpFromRef (out);
    OSTab tab1 (ptrOut);
    if (this->getObjectLabel () != "") {
      out << "label: " << this->getObjectLabel () << endl;
    }
    out << "initialized: " << (isInitialized_ ? "true" : "false") << endl
        << "computed: " << (isComputed_ ? "true" : "false") << endl
        << "number of initialize calls: " << numInitialize_ << endl
        << "number of compute calls: " << numCompute_ << endl
        << "number of apply calls: " << numApply_ << endl
        << "total time in seconds in initialize: " << initializeTime_ << endl
        << "total time in seconds in compute: " << computeTime_ << endl
        << "total time in seconds in apply: " << applyTime_ << endl;
    if (verbLevel >= Teuchos::VERB_EXTREME) {
      out << "A_local_tridi_:" << endl;
      A_local_tridi_.print(out);
      }
      out << "ipiv_: " << Teuchos::toString (ipiv_) << endl;
    }
}

template<class MatrixType>
void TriDiSolver<MatrixType, false>::describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using std::endl;

  RCP<FancyOStream> ptrOut = rcpFromRef (out);
  OSTab tab0 (ptrOut);
  if (A_.is_null ()) {
    // If A_ is null, we don't have a communicator, so we can't
    // safely print local data on all processes.  Just print the
    // local data without arbitration between processes, and hope
    // for the best.
    if (verbLevel > Teuchos::VERB_NONE) {
      out << "Ifpack2::Details::TriDiSolver:" << endl;
    }
    describeLocal (out, verbLevel);
  }
  else {
    // If A_ is not null, we have a communicator, so we can
    // arbitrate among all processes to print local data.
    const Teuchos::Comm<int>& comm = * (A_->getRowMap ()->getComm ());
    const int myRank = comm.getRank ();
    const int numProcs = comm.getSize ();
    if (verbLevel > Teuchos::VERB_NONE && myRank == 0) {
      out << "Ifpack2::Details::TriDiSolver:" << endl;
    }
    OSTab tab1 (ptrOut);
    for (int p = 0; p < numProcs; ++p) {
      if (myRank == p) {
        out << "Process " << myRank << ":" << endl;
        describeLocal (out, verbLevel);
      }
      comm.barrier ();
      comm.barrier ();
      comm.barrier ();
    } // for p = 0 .. numProcs-1
  }
}

template<class MatrixType>
void TriDiSolver<MatrixType, false>::extract (Teuchos::SerialTriDiMatrix<int, scalar_type>& A_local_tridi,
         const row_matrix_type& A_local)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  typedef local_ordinal_type LO;
  typedef typename Teuchos::ArrayView<LO>::size_type size_type;

  // Fill the local tridi matrix with zeros.
  A_local_tridi.putScalar (STS::zero ());

  //
  // Map both row and column indices to local indices.  We can use the
  // row Map's local indices for row indices, and the column Map's
  // local indices for column indices.  It doesn't really matter;
  // either way is just a permutation of rows and columns.
  //
  const map_type& rowMap = * (A_local.getRowMap ());

  // Temporary arrays to hold the indices and values of the entries in
  // each row of A_local.
  const size_type maxNumRowEntries =
    static_cast<size_type> (A_local.getNodeMaxNumRowEntries ());
  Array<LO> localIndices (maxNumRowEntries);
  Array<scalar_type> values (maxNumRowEntries);

  const LO numLocalRows = static_cast<LO> (rowMap.getNodeNumElements ());
  const LO minLocalRow = rowMap.getMinLocalIndex ();
  // This slight complication of computing the upper loop bound avoids
  // issues if the row Map has zero entries on the calling process.
  const LO maxLocalRow = minLocalRow + numLocalRows; // exclusive bound
  for (LO localRow = minLocalRow; localRow < maxLocalRow; ++localRow) {
    // The LocalFilter automatically excludes "off-process" entries.
    // That means all the column indices in this row belong to the
    // domain Map.  We can, therefore, just use the local row and
    // column indices to put each entry directly in the tridi matrix.
    // It's OK if the column Map puts the local indices in a different
    // order; the Import will bring them into the correct order.
    const size_type numEntriesInRow =
      static_cast<size_type> (A_local.getNumEntriesInLocalRow (localRow));
    size_t numEntriesOut = 0; // ignored
    A_local.getLocalRowCopy (localRow,
                             localIndices (0, numEntriesInRow),
                             values (0, numEntriesInRow),
                             numEntriesOut);
    for (LO k = 0; k < numEntriesInRow; ++k) {
      const LO localCol = localIndices[k];
      const scalar_type val = values[k];
      // We use += instead of =, in case there are duplicate entries
      // in the row.  There should not be, but why not be general?
      // NOTE: we only extract the TriDi part of the row matrix. Do not extract DU2
      if( localCol >= localRow-1 && localCol <= localRow+1 )
        A_local_tridi(localRow, localCol) += val;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// Stub implementation
//////////////////////////////////////////////////////////////////////

template<class MatrixType>
TriDiSolver<MatrixType, true>::TriDiSolver (const Teuchos::RCP<const row_matrix_type>& A) {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
Teuchos::RCP<const typename TriDiSolver<MatrixType, true>::map_type>
TriDiSolver<MatrixType, true>::getDomainMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
Teuchos::RCP<const typename TriDiSolver<MatrixType, true>::map_type>
TriDiSolver<MatrixType, true>::getRangeMap () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
void
TriDiSolver<MatrixType, true>::setParameters (const Teuchos::ParameterList& params) {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
bool
TriDiSolver<MatrixType, true>::isInitialized () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
bool
TriDiSolver<MatrixType, true>::isComputed () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
int
TriDiSolver<MatrixType, true>::getNumInitialize () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
int
TriDiSolver<MatrixType, true>::getNumCompute () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
int
TriDiSolver<MatrixType, true>::getNumApply () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
double
TriDiSolver<MatrixType, true>::getInitializeTime () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
double
TriDiSolver<MatrixType, true>::getComputeTime () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
double
TriDiSolver<MatrixType, true>::getApplyTime () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
Teuchos::RCP<const typename TriDiSolver<MatrixType, true>::row_matrix_type>
TriDiSolver<MatrixType, true>::getMatrix () const {
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
void TriDiSolver<MatrixType, true>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
void TriDiSolver<MatrixType, true>::initialize ()
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
TriDiSolver<MatrixType, true>::~TriDiSolver ()
{
  // Destructors should never throw exceptions.
}


template<class MatrixType>
void TriDiSolver<MatrixType, true>::compute ()
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
void TriDiSolver<MatrixType, true>::apply (
       const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
std::string
TriDiSolver<MatrixType, true>::description () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}


template<class MatrixType>
void TriDiSolver<MatrixType, true>::describe(Teuchos::FancyOStream& out,
                                              const Teuchos::EVerbosityLevel verbLevel) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Not implemented");
}

}// namespace Details
} // namespace Ifpack2

#define IFPACK2_DETAILS_TRIDISOLVER_INSTANT(S,LO,GO,N)                  \
  template class Ifpack2::Details::TriDiSolver< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_DETAILS_TRIDISOLVER_HPP
