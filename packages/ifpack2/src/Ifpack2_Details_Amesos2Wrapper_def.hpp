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

#ifndef IFPACK2_DETAILS_AMESOS2WRAPPER_DEF_HPP
#define IFPACK2_DETAILS_AMESOS2WRAPPER_DEF_HPP

// // disable clang warnings
// #ifdef __clang__
// #pragma clang system_header
// #endif

#include <Ifpack2_Heap.hpp>
#include <Ifpack2_Condest.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_AMESOS2)
#include <Amesos2.hpp>

namespace Ifpack2 {
namespace Details {

template <class MatrixType>
Amesos2Wrapper<MatrixType>::
Amesos2Wrapper (const Teuchos::RCP<const row_matrix_type>& A) :
  A_ (A),
  Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  IsInitialized_ (false),
  IsComputed_ (false)
{}

template <class MatrixType>
Amesos2Wrapper<MatrixType>::~Amesos2Wrapper()
{}

template <class MatrixType>
void Amesos2Wrapper<MatrixType>::setParameters (const Teuchos::ParameterList& params)
{
  // FIXME (mfh 10 Dec 2013) This class does not currently set parameters.
  (void) params;
}


template <class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
Amesos2Wrapper<MatrixType>::getComm () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Amesos2Wrapper::getComm: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getComm ();
}


template <class MatrixType>
Teuchos::RCP<const typename Amesos2Wrapper<MatrixType>::row_matrix_type>
Amesos2Wrapper<MatrixType>::getMatrix () const {
  return A_;
}


template <class MatrixType>
Teuchos::RCP<const typename Amesos2Wrapper<MatrixType>::map_type>
Amesos2Wrapper<MatrixType>::getDomainMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Amesos2Wrapper::getDomainMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getDomainMap ();
}


template <class MatrixType>
Teuchos::RCP<const typename Amesos2Wrapper<MatrixType>::map_type>
Amesos2Wrapper<MatrixType>::getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Amesos2Wrapper::getRangeMap: "
    "The matrix is null.  Please call setMatrix() with a nonnull input "
    "before calling this method.");
  return A_->getRangeMap ();
}


template <class MatrixType>
bool Amesos2Wrapper<MatrixType>::hasTransposeApply () const {
  return true;
}


template <class MatrixType>
int Amesos2Wrapper<MatrixType>::getNumInitialize () const {
  return NumInitialize_;
}


template <class MatrixType>
int Amesos2Wrapper<MatrixType>::getNumCompute () const {
  return NumCompute_;
}


template <class MatrixType>
int Amesos2Wrapper<MatrixType>::getNumApply () const {
  return NumApply_;
}


template <class MatrixType>
double Amesos2Wrapper<MatrixType>::getInitializeTime () const {
  return InitializeTime_;
}


template<class MatrixType>
double Amesos2Wrapper<MatrixType>::getComputeTime () const {
  return ComputeTime_;
}


template<class MatrixType>
double Amesos2Wrapper<MatrixType>::getApplyTime () const {
  return ApplyTime_;
}

template<class MatrixType>
typename Amesos2Wrapper<MatrixType>::magnitude_type
Amesos2Wrapper<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const row_matrix_type>& matrix)
{
  if (! isComputed ()) {
    return -STM::one ();
  }
  // NOTE: this is computing the *local* condest
  if (Condest_ == -STM::one ()) {
    Condest_ = Ifpack2::Condest (*this, CT, MaxIters, Tol, matrix);
  }
  return Condest_;
}


template<class MatrixType>
void Amesos2Wrapper<MatrixType>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  // Check in serial or one-process mode if the matrix is square.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A.is_null () && A->getComm ()->getSize () == 1 &&
    A->getNodeNumRows () != A->getNodeNumCols (),
    std::runtime_error, "Ifpack2::Amesos2Wrapper::setMatrix: If A's communicator only "
    "contains one process, then A must be square.  Instead, you provided a "
    "matrix A with " << A->getNodeNumRows () << " rows and "
    << A->getNodeNumCols () << " columns.");

  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous
  // factorization.
  IsInitialized_ = false;
  IsComputed_ = false;
  A_local_ = Teuchos::null;
  A_ = A;
}


template<class MatrixType>
void Amesos2Wrapper<MatrixType>::initialize ()
{
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  const std::string timerName ("Ifpack2::Amesos2Wrapper::initialize");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    // Check that the matrix is nonnull.
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_.is_null (), std::runtime_error, "Ifpack2::Amesos2Wrapper::initialize: "
      "The matrix to precondition is null.  Please call setMatrix() with a "
      "nonnull input before calling this method.");

    // Clear any previous computations.
    IsInitialized_ = false;
    IsComputed_ = false;
    Condest_ = -STM::one ();
    A_local_ = Teuchos::null;

    // Construct the local matrix.
    A_local_ = makeLocalMatrix (*A_);

    // FIXME (10 Dec 2013) This (the Amesos2 solver type) should be a
    // run-time parameter through the input ParameterList.

    std::string solverType;

#if defined(HAVE_AMESOS2_SUPERLU)
    solverType = "superlu";
#elif defined(HAVE_AMESOS2_KLU2)
    solverType = "klu";
#elif defined(HAVE_AMESOS2_SUPERLUDIST)
    solverType = "superludist";
#elif defined(HAVE_AMESOS2_LAPACK)
    solverType = "lapack";
#endif

    amesos2solver_ = Amesos2::create<MatrixType, MV> (solverType, A_local_);
    amesos2solver_->preOrdering ();

    // The symbolic factorization properly belongs to initialize(),
    // since initialize() is concerned with the matrix's structure
    // (and compute() with the matrix's values).
    amesos2solver_->symbolicFactorization ();
  } // Stop timing here.

  IsInitialized_ = true;
  ++NumInitialize_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  InitializeTime_ = timer->totalElapsedTime ();
}

template<class MatrixType>
void Amesos2Wrapper<MatrixType>::compute ()
{
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  // Don't count initialization in the compute() time.
  if (! isInitialized ()) {
    initialize ();
  }

  const std::string timerName ("Ifpack2::AdditiveSchwarz::compute");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    TimeMonitor timeMon (*timer);
    amesos2solver_->numericFactorization ();
  } // Stop timing here.

  IsComputed_ = true;
  ++NumCompute_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ComputeTime_ = timer->totalElapsedTime ();
}


template <class MatrixType>
void Amesos2Wrapper<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
       Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  const std::string timerName ("Ifpack2::Amesos2Wrapper::apply");
  RCP<Time> timer = TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = TimeMonitor::getNewCounter (timerName);
  }

  { // Start timing here.
    TimeMonitor timeMon (*timer);

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isComputed (), std::runtime_error,
      "Ifpack2::Amesos2Wrapper::apply: You must call compute() to compute the "
      "incomplete factorization, before calling apply().");

    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
      "Ifpack2::Amesos2Wrapper::apply: X and Y must have the same number of columns.  "
      "X has " << X.getNumVectors () << " columns, but Y has "
      << Y.getNumVectors () << " columns.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      ! mode != Teuchos::NO_TRANS, std::logic_error,
      "Ifpack2::Amesos2Wrapper::apply: Solving with the transpose (mode == "
      "Teuchos::TRANS) or conjugate transpose (Teuchos::CONJ_TRANS) is not "
      "implemented.");

    // If beta != 0, create a temporary multivector Y_temp to hold the
    // contents of alpha*M^{-1}*X.  Otherwise, alias Y_temp to Y.
    RCP<MV> Y_temp;
    Y_temp = rcp (new MV (Y.getMap (), Y.getNumVectors ()));

    // If X and Y are pointing to the same memory location, create an
    // auxiliary vector, X_temp, so that we don't clobber the input
    // when computing the output.  Otherwise, alias X_temp to X.
    RCP<const MV> X_temp;
    if (X.getLocalMV ().getValues () == Y.getLocalMV ().getValues ()) {
      X_temp = rcp (new MV (X));
    } else {
      X_temp = rcpFromRef (X);
    }

    // FIXME (mfh 10 Dec 2013) LocalFilter would take care of all this.
    // Replace this extraction stuff with use of LocalFilter.

    // construct local vectors
    size_t numvecs = X_temp->getNumVectors();
    RCP<const map_type> globalRowMap = A_->getRowMap ();
    const local_ordinal_type numRows = globalRowMap->getNodeNumElements ();
    RCP<MV> localY = rcp (new MV (A_local_->getRowMap (), numvecs));
    RCP<MV> localX = rcp (new MV (A_local_->getRowMap (), numvecs));
    // extract values
    for (size_t j = 0; j < numvecs; ++j) {
      Teuchos::ArrayRCP<const scalar_type> vecj = X_temp->getData (j);
      for(local_ordinal_type i = 0; i < numRows; ++i) {
        localX->replaceLocalValue (i, j, vecj[i]);
      }
    }

    // solve
    amesos2solver_->setX (localY);
    amesos2solver_->setB (localX);
    amesos2solver_->solve ();

    // FIXME (mfh 10 Dec 2013) LocalFilter would take care of all this.
    // Replace this extraction stuff with use of LocalFilter.

    // extract to global vector
    for (size_t j = 0; j < localY->getNumVectors (); ++j) {
      Teuchos::ArrayRCP<const scalar_type> localview = localY->getData (j);
      for (unsigned int i = 0; i<globalRowMap->getNodeNumElements (); ++i) {
        Y_temp->replaceLocalValue (i, j, localview[i]);
      }
    }

    Y.update (alpha, *Y_temp, beta);
  } // Stop timing here.

  ++NumApply_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ApplyTime_ = timer->totalElapsedTime ();
}


template <class MatrixType>
std::string Amesos2Wrapper<MatrixType>::description() const {
  using Teuchos::TypeNameTraits;
  std::ostringstream os;

  os << "Ifpack2::Amesos2Wrapper: {"
     << "MatrixType: \"" << TypeNameTraits<MatrixType>::name ()
     << "\", ";
  if (this->getObjectLabel () != "") {
    os << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  os << "Initialized: " << (isInitialized () ? "true" : "false")
     << ", "
     << "Computed: " << (isComputed () ? "true" : "false")
     << ", "
     << "Number of rows: " << A_->getGlobalNumRows ()
     << ", "
     << "Number of columns: " << A_->getGlobalNumCols ()
     << "}";
  return os.str();
}


template <class MatrixType>
void
Amesos2Wrapper<MatrixType>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::Comm;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::TypeNameTraits;
  using std::endl;

  const Teuchos::EVerbosityLevel vl = (verbLevel == Teuchos::VERB_DEFAULT) ?
    Teuchos::VERB_LOW : verbLevel;

  // describe() starts, by convention, with a tab before it prints anything.
  OSTab tab0 (out);
  if (vl > Teuchos::VERB_NONE) {
    out << "Ifpack2::Amesos2Wrapper:" << endl;
    OSTab tab1 (out);
    out << "MatrixType: \"" << TypeNameTraits<MatrixType>::name ()
        << "\"" << endl;

    if (this->getObjectLabel () != "") {
      out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
    }

    out << "Initialized: " << (isInitialized () ? "true" : "false") << endl;
    out << "Computed: " << (isComputed () ? "true" : "false") << endl;
    out << "Number of initialize calls: " << getNumInitialize () << endl;
    out << "Number of compute calls: " << getNumCompute () << endl;
    out << "Number of apply calls: " << getNumApply () << endl;
    out << "Total time in seconds for initialize: " << getInitializeTime () << endl;
    out << "Total time in seconds for compute: " << getComputeTime () << endl;
    out << "Total time in seconds for apply: " << getApplyTime () << endl;

    if (vl > Teuchos::VERB_LOW) {
      out << "Local matrix:" << endl;
      A_local_->describe (out, vl);
    }
  }
}

template <class MatrixType>
Teuchos::RCP<MatrixType>
Amesos2Wrapper<MatrixType>::makeLocalMatrix (const row_matrix_type& A)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // FIXME (mfh 10 Dec 2013) Why aren't you using LocalFilter here?

  // local communicator
  RCP<const Teuchos::Comm<int> > localComm;
#ifdef HAVE_MPI
  localComm = rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF));
#else
  localComm = rcp (new Teuchos::SerialComm<int> ());
#endif

  // get row map and setup local matrix
  RCP<const map_type> globalRowMap = A.getRowMap ();
  RCP<const map_type> globalColMap = A.getColMap ();
  const local_ordinal_type numRows = globalRowMap->getNodeNumElements ();

  RCP<const map_type> localRowMap =
    rcp (new map_type (numRows, 0, localComm, Tpetra::GloballyDistributed,
                       A.getNode ()));
  RCP<MatrixType> Alocal = rcp (new MatrixType (localRowMap, localRowMap, 100));

  // extract rows
  for (local_ordinal_type i = 0; i < numRows; ++i) {
    ArrayView<const local_ordinal_type> indices;
    ArrayView<const scalar_type> values;
    Array<local_ordinal_type> indices_vec;
    Array<scalar_type> values_vec;

    A.getLocalRowView (i, indices, values);
    indices_vec.resize (0);
    values_vec.resize (0);
    for (unsigned int j = 0; j < indices.size (); ++j) {
      const local_ordinal_type local_col = indices[j];
      const global_ordinal_type global_col = globalColMap->getGlobalElement (local_col);
      if (globalRowMap->isNodeGlobalElement (global_col)) {
        indices_vec.push_back (globalRowMap->getLocalElement (global_col));
        values_vec.push_back (values[j]);
      }
    }
    Alocal->insertLocalValues (i, indices_vec (), values_vec ());
  }
  Alocal->fillComplete ();
  return Alocal;
}

} // namespace Details
} // namespace Ifpack2

#endif // HAVE_IFPACK2_AMESOS2
#endif /* IFPACK2_DETAILS_AMESOS2WRAPPER_DEF_HPP */

