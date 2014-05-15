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

#ifdef HAVE_IFPACK2_AMESOS2
#include <Amesos2.hpp>

namespace Ifpack2 {
namespace Details {

template <class MatrixType>
Amesos2Wrapper<MatrixType>::
Amesos2Wrapper (const Teuchos::RCP<const row_matrix_type>& A) :
  Condest_ (-STM::one ()),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  IsInitialized_ (false),
  IsComputed_ (false)
{
  // The input matrix A (an instance of Tpetra::RowMatrix) must have
  // type MatrixType (a specialization of Tpetra::CrsMatrix).
  Teuchos::RCP<const MatrixType> A_crs =
    Teuchos::rcp_dynamic_cast<const MatrixType> (A);
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_crs.is_null (), std::invalid_argument, "Ifpack2::Details::Amesos2Wrapper "
    "constructor: The input matrix A is not a Tpetra::CrsMatrix instance, "
    "but it must be in order for Amesos2 to work correctly.");
  A_ = A_crs;
}

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
  // It's legal for A to be null; in that case, you may not call
  // initialize() until calling setMatrix() with a nonnull input.
  // Regardless, setting the matrix invalidates any previous
  // factorization.
  IsInitialized_ = false;
  IsComputed_ = false;

  if (A.is_null ()) {
    A_ = Teuchos::null;
  }
  else {
    // The input matrix A (an instance of Tpetra::RowMatrix) must have
    // type MatrixType (a specialization of Tpetra::CrsMatrix).
    Teuchos::RCP<const MatrixType> A_crs =
      Teuchos::rcp_dynamic_cast<const MatrixType> (A);
    TEUCHOS_TEST_FOR_EXCEPTION(
      A_crs.is_null (), std::invalid_argument, "Ifpack2::Details::Amesos2Wrapper "
      "constructor: The input matrix A is not a Tpetra::CrsMatrix instance, "
      "but it must be in order for Amesos2 to work correctly.");
    A_ = A_crs;
  }

  // FIXME (mfh 10 Dec 2013) Currently, initialize() recreates
  // amesos2solver_ unconditionally, so this code won't have any
  // effect.  Once we fix initialize() so that it keeps
  // amesos2solver_, the code below will be effective.
  if (! amesos2solver_.is_null ()) {
    amesos2solver_->setA (A_);
  }
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

    // FIXME (10 Dec 2013) This (the Amesos2 solver type) should be a
    // run-time parameter through the input ParameterList.
    // (9 May 2014) JJH Ifpack2 also shouldn't be checking the availability direct solvers.
    // It's up to Amesos2 to test for this and throw an exception
    // (which it does in Amesos2::Factory::create).

    std::string solverType;

#if defined(HAVE_AMESOS2_SUPERLU)
    solverType = "superlu";
#elif defined(HAVE_AMESOS2_KLU2)
    solverType = "klu";
#elif defined(HAVE_AMESOS2_SUPERLUDIST)
    solverType = "superludist";
#elif defined(HAVE_AMESOS2_CHOLMOD)
    solverType = "cholmod";
#elif defined(HAVE_AMESOS2_LAPACK)
    solverType = "lapack";
#else
    // FIXME (9 May 2014) JJH Amesos2 does not yet expose KLU2, its internal direct solver.
    // This means there's no fallback option, thus we throw an exception here.
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Amesos2 has not been configured with any direct solver support.");
#endif

    // FIXME (10 Dec 2013) It shouldn't be necessary to recreate the
    // solver each time, since Amesos2::Solver has a setA() method.
    // See the implementation of setMatrix().

    amesos2solver_ = Amesos2::create<MatrixType, MV> (solverType, A_);
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
      mode != Teuchos::NO_TRANS, std::logic_error,
      "Ifpack2::Amesos2Wrapper::apply: Solving with the transpose (mode == "
      "Teuchos::TRANS) or conjugate transpose (Teuchos::CONJ_TRANS) is not "
      "implemented.");

    // If alpha != 1 or beta != 0, create a temporary multivector
    // Y_temp to hold the contents of alpha*M^{-1}*X.  Otherwise,
    // alias Y_temp to Y.
    RCP<MV> Y_temp = (alpha != STS::one () || beta != STS::zero ()) ?
      rcp (new MV (Y.getMap (), Y.getNumVectors ())) :
      rcpFromRef (Y);

    // If X and Y are pointing to the same memory location, create an
    // auxiliary vector, X_temp, so that we don't clobber the input
    // when computing the output.  Otherwise, alias X_temp to X.
    RCP<const MV> X_temp;
    if (X.getLocalMV ().getValues () == Y.getLocalMV ().getValues ()) {
      X_temp = rcp (new MV (createCopy(X)));
    } else {
      X_temp = rcpFromRef (X);
    }

    // Use the precomputed factorization to solve.
    amesos2solver_->setX (Y_temp);
    amesos2solver_->setB (X_temp);
    amesos2solver_->solve ();

    if (alpha != STS::one () || beta != STS::zero ()) {
      Y.update (alpha, *Y_temp, beta);
    }
  } // Stop timing here.

  ++NumApply_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ApplyTime_ = timer->totalElapsedTime ();
}


template <class MatrixType>
std::string Amesos2Wrapper<MatrixType>::description () const {
  using Teuchos::TypeNameTraits;
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::Amesos2Wrapper\": {";
  if (this->getObjectLabel () != "") {
    os << "Label: \"" << this->getObjectLabel () << "\", ";
  }
  os << "Initialized: " << (isInitialized () ? "true" : "false")
     << ", Computed: " << (isComputed () ? "true" : "false");

  if (A_.is_null ()) {
    os << ", Matrix: null";
  }
  else {
    os << ", Matrix: not null"
       << ", Global matrix dimensions: ["
       << A_->getGlobalNumRows () << ", " << A_->getGlobalNumCols () << "]";
  }

  os << "}";
  return os.str ();
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
    out << "\"Ifpack2::Amesos2Wrapper\":" << endl;
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
      out << "Matrix:" << endl;
      A_->describe (out, vl);
    }
  }
}

} // namespace Details
} // namespace Ifpack2

#define IFPACK2_DETAILS_AMESOS2WRAPPER_INSTANT(S,LO,GO,N) \
  template class Ifpack2::Details::Amesos2Wrapper< Tpetra::CrsMatrix<S, LO, GO, N> >;

#endif // HAVE_IFPACK2_AMESOS2
#endif // IFPACK2_DETAILS_AMESOS2WRAPPER_DEF_HPP
