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

#ifndef IFPACK2_CHEBYSHEV_DEF_HPP
#define IFPACK2_CHEBYSHEV_DEF_HPP

#include <Ifpack2_Parameters.hpp>
#include <Ifpack2_Chebyshev.hpp>
#include <Teuchos_TimeMonitor.hpp>

namespace Ifpack2 {

template<class MatrixType>
Chebyshev<MatrixType>::
Chebyshev (const Teuchos::RCP<const row_matrix_type>& A)
  : impl_ (A),
    IsInitialized_ (false),
    IsComputed_ (false),
    NumInitialize_ (0),
    NumCompute_ (0),
    NumApply_ (0),
    InitializeTime_ (0.0),
    ComputeTime_ (0.0),
    ApplyTime_ (0.0),
    ComputeFlops_ (0.0),
    ApplyFlops_ (0.0)
{
  this->setObjectLabel ("Ifpack2::Chebyshev");
}


template<class MatrixType>
Chebyshev<MatrixType>::~Chebyshev() {
}


template<class MatrixType>
void Chebyshev<MatrixType>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A.getRawPtr () != impl_.getMatrix ().getRawPtr ()) {
    IsInitialized_ = false;
    IsComputed_ = false;
    impl_.setMatrix (A);
  }
}


template<class MatrixType>
void
Chebyshev<MatrixType>::setParameters (const Teuchos::ParameterList& List)
{
  // FIXME (mfh 25 Jan 2013) Casting away const is bad here.
  impl_.setParameters (const_cast<Teuchos::ParameterList&> (List));
}


template<class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
Chebyshev<MatrixType>::getComm () const
{
  Teuchos::RCP<const row_matrix_type> A = impl_.getMatrix ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null (), std::runtime_error, "Ifpack2::Chebyshev::getComm: The input "
    "matrix A is null.  Please call setMatrix() with a nonnull input matrix "
    "before calling this method.");
  return A->getRowMap ()->getComm ();
}


template<class MatrixType>
Teuchos::RCP<const typename Chebyshev<MatrixType>::row_matrix_type>
Chebyshev<MatrixType>::
getMatrix() const {
  return impl_.getMatrix ();
}


template<class MatrixType>
Teuchos::RCP<const MatrixType>
Chebyshev<MatrixType>::
getCrsMatrix() const {
  return Teuchos::rcp_dynamic_cast<const MatrixType> (impl_.getMatrix ());
}


template<class MatrixType>
Teuchos::RCP<const typename Chebyshev<MatrixType>::map_type>
Chebyshev<MatrixType>::
getDomainMap () const
{
  Teuchos::RCP<const row_matrix_type> A = impl_.getMatrix ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null (), std::runtime_error, "Ifpack2::Chebyshev::getDomainMap: The "
    "input matrix A is null.  Please call setMatrix() with a nonnull input "
    "matrix before calling this method.");
  return A->getDomainMap ();
}


template<class MatrixType>
Teuchos::RCP<const typename Chebyshev<MatrixType>::map_type>
Chebyshev<MatrixType>::
getRangeMap () const
{
  Teuchos::RCP<const row_matrix_type> A = impl_.getMatrix ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null (), std::runtime_error, "Ifpack2::Chebyshev::getRangeMap: The "
    "input matrix A is null.  Please call setMatrix() with a nonnull input "
    "matrix before calling this method.");
  return A->getRangeMap ();
}


template<class MatrixType>
bool Chebyshev<MatrixType>::hasTransposeApply() const {
  return impl_.hasTransposeApply ();
}


template<class MatrixType>
int Chebyshev<MatrixType>::getNumInitialize() const {
  return NumInitialize_;
}


template<class MatrixType>
int Chebyshev<MatrixType>::getNumCompute() const {
  return NumCompute_;
}


template<class MatrixType>
int Chebyshev<MatrixType>::getNumApply() const {
  return NumApply_;
}


template<class MatrixType>
double Chebyshev<MatrixType>::getInitializeTime() const {
  return InitializeTime_;
}


template<class MatrixType>
double Chebyshev<MatrixType>::getComputeTime() const {
  return ComputeTime_;
}


template<class MatrixType>
double Chebyshev<MatrixType>::getApplyTime() const {
  return ApplyTime_;
}


template<class MatrixType>
double Chebyshev<MatrixType>::getComputeFlops () const {
  return ComputeFlops_;
}


template<class MatrixType>
double Chebyshev<MatrixType>::getApplyFlops () const {
  return ApplyFlops_;
}


template<class MatrixType>
void
Chebyshev<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
       Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  const std::string timerName ("Ifpack2::Chebyshev::apply");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }

  // Start timing here.
  {
    Teuchos::TimeMonitor timeMon (*timer);

    // compute() calls initialize() if it hasn't already been called.
    // Thus, we only need to check isComputed().
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! isComputed (), std::runtime_error,
      "Ifpack2::Chebyshev::apply(): You must call the compute() method before "
      "you may call apply().");
    TEUCHOS_TEST_FOR_EXCEPTION(
      X.getNumVectors () != Y.getNumVectors (), std::runtime_error,
      "Ifpack2::Chebyshev::apply(): X and Y must have the same number of "
      "columns.  X.getNumVectors() = " << X.getNumVectors() << " != "
      << "Y.getNumVectors() = " << Y.getNumVectors() << ".");
    applyImpl (X, Y, mode, alpha, beta);
  }
  ++NumApply_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ApplyTime_ = timer->totalElapsedTime ();
}


template<class MatrixType>
void
Chebyshev<MatrixType>::
applyMat (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
          Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
          Teuchos::ETransp mode) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors () != Y.getNumVectors (), std::invalid_argument,
    "Ifpack2::Chebyshev::applyMat: X.getNumVectors() != Y.getNumVectors().");

  Teuchos::RCP<const row_matrix_type> A = impl_.getMatrix ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null (), std::runtime_error, "Ifpack2::Chebyshev::applyMat: The input "
    "matrix A is null.  Please call setMatrix() with a nonnull input matrix "
    "before calling this method.");

  A->apply (X, Y, mode);
}


template<class MatrixType>
void Chebyshev<MatrixType>::initialize () {
  // We create the timer, but this method doesn't do anything, so
  // there is no need to start the timer.  The resulting total time
  // will always be zero.
  const std::string timerName ("Ifpack2::Chebyshev::initialize");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }
  IsInitialized_ = true;
  ++NumInitialize_;
}


template<class MatrixType>
void Chebyshev<MatrixType>::compute ()
{
  const std::string timerName ("Ifpack2::Chebyshev::compute");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }

  // Start timing here.
  {
    Teuchos::TimeMonitor timeMon (*timer);
    if (! isInitialized ()) {
      initialize ();
    }
    IsComputed_ = false;
    impl_.compute ();
  }
  IsComputed_ = true;
  ++NumCompute_;

  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ComputeTime_ = timer->totalElapsedTime ();
}


template <class MatrixType>
std::string Chebyshev<MatrixType>::description () const {
  std::ostringstream out;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  out << "\"Ifpack2::Chebyshev\": {";
  out << "Initialized: " << (isInitialized () ? "true" : "false") << ", "
      << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  out << impl_.description() << ", ";

  if (impl_.getMatrix ().is_null ()) {
    out << "Matrix: null";
  }
  else {
    out << "Global matrix dimensions: ["
        << impl_.getMatrix ()->getGlobalNumRows () << ", "
        << impl_.getMatrix ()->getGlobalNumCols () << "]"
        << ", Global nnz: " << impl_.getMatrix ()->getGlobalNumEntries();
  }

  out << "}";
  return out.str ();
}


template <class MatrixType>
void Chebyshev<MatrixType>::
describe (Teuchos::FancyOStream &out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  const Teuchos::EVerbosityLevel vl =
    (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;
  const int myRank = this->getComm ()->getRank ();

  if (vl != Teuchos::VERB_NONE && myRank == 0) {
    // By convention, describe() starts with a tab.
    Teuchos::OSTab tab0 (out);
    out << description ();
  }

#if 0
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using std::endl;
  using std::setw;

  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) {
    vl = VERB_LOW;
  }
  RCP<const Comm<int> > comm = A_->getRowMap ()->getComm ();

  const int myImageID = comm->getRank();
  Teuchos::OSTab tab(out);

  scalar_type MinVal, MaxVal;
  if (IsComputed_) {
    Teuchos::ArrayRCP<const scalar_type> DiagView = InvDiagonal_->get1dView();
    scalar_type myMinVal = DiagView[0];
    scalar_type myMaxVal = DiagView[0];
    for(typename Teuchos::ArrayRCP<scalar_type>::size_type i=1; i<DiagView.size(); ++i) {
      if (STS::magnitude(myMinVal) > STS::magnitude(DiagView[i])) myMinVal = DiagView[i];
      if (STS::magnitude(myMaxVal) < STS::magnitude(DiagView[i])) myMaxVal = DiagView[i];
    }
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &myMinVal, &MinVal);
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, 1, &myMaxVal, &MaxVal);
  }

  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium:
  //    high:
  // extreme:
  if (vl != VERB_NONE && myImageID == 0) {
    out << this->description() << endl;
    out << endl;
    out << "===============================================================================" << std::endl;
    out << "Degree of polynomial      = " << PolyDegree_ << std::endl;
    if   (ZeroStartingSolution_) { out << "Using zero starting solution" << endl; }
    else                         { out << "Using input starting solution" << endl; }
    if (IsComputed_) {
      out << "Minimum value on stored inverse diagonal = " << MinVal << std::endl;
      out << "Maximum value on stored inverse diagonal = " << MaxVal << std::endl;
    }
    out << std::endl;
    out << "Phase           # calls    Total Time (s)     Total MFlops      MFlops/s       " << endl;
    out << "------------    -------    ---------------    ---------------   ---------------" << endl;
    out << setw(12) << "initialize()" << setw(5) << getNumInitialize() << "    " << setw(15) << getInitializeTime() << endl;
    out << setw(12) << "compute()" << setw(5) << getNumCompute()    << "    " << setw(15) << getComputeTime() << "    "
        << setw(15) << getComputeFlops() << "    "
        << setw(15) << (getComputeTime() != 0.0 ? getComputeFlops() / getComputeTime() * 1.0e-6 : 0.0) << endl;
    out << setw(12) << "apply()" << setw(5) << getNumApply()    << "    " << setw(15) << getApplyTime() << "    "
        << setw(15) << getApplyFlops() << "    "
        << setw(15) << (getApplyTime() != 0.0 ? getApplyFlops() / getApplyTime() * 1.0e-6 : 0.0) << endl;
    out << "===============================================================================" << std::endl;
    out << endl;
  }
#endif // 0
}

template<class MatrixType>
void
Chebyshev<MatrixType>::
applyImpl (const MV& X,
           MV& Y,
           Teuchos::ETransp mode,
           scalar_type alpha,
           scalar_type beta) const
{
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcpFromRef;

  const scalar_type zero = STS::zero();
  const scalar_type one = STS::one();

  // Y = beta*Y + alpha*M*X.

  // If alpha == 0, then we don't need to do Chebyshev at all.
  if (alpha == zero) {
    if (beta == zero) { // Obey Sparse BLAS rules; avoid 0*NaN.
      Y.putScalar (zero);
    }
    else {
      Y.scale (beta);
    }
    return;
  }

  // If beta != 0, then we need to keep a (deep) copy of the initial
  // value of Y, so that we can add beta*it to the Chebyshev result at
  // the end.  Usually this method is called with beta == 0, so we
  // don't have to worry about caching Y_org.
  RCP<MV> Y_orig;
  if (beta != zero) {
    Y_orig = rcp (new MV (Y, Teuchos::Copy));
  }

  // If X and Y point to the same memory location, we need to use a
  // (deep) copy of X (X_copy) as the input MV.  Otherwise, just let
  // X_copy point to X.
  //
  // This is hopefully an uncommon use case, so we don't bother to
  // optimize for it by caching X_copy.
  RCP<const MV> X_copy;
  bool copiedInput = false;
  if (X.getLocalMV ().getValues () == Y.getLocalMV ().getValues ()) {
    X_copy = rcp (new MV (X, Teuchos::Copy));
    copiedInput = true;
  }
  else {
    X_copy = rcpFromRef (X);
  }

  // If alpha != 1, fold alpha into (a deep copy of) X.
  //
  // This is an uncommon use case, so we don't bother to optimize for
  // it by caching X_copy.  However, we do check whether we've already
  // copied X above, to avoid a second copy.
  if (alpha != one) {
    RCP<MV> X_copy_nonConst = rcp_const_cast<MV> (X_copy);
    if (! copiedInput) {
      X_copy_nonConst = rcp (new MV (X, Teuchos::Copy));
      copiedInput = true;
    }
    X_copy_nonConst->scale (alpha);
    X_copy = rcp_const_cast<const MV> (X_copy_nonConst);
  }

  impl_.apply (*X_copy, Y);

  if (beta != zero) {
    Y.update (beta, *Y_orig, one); // Y = beta * Y_orig + 1 * Y
  }
}


template<class MatrixType>
typename MatrixType::scalar_type Chebyshev<MatrixType>::getLambdaMaxForApply () const {
  return impl_.getLambdaMaxForApply ();
}



}//namespace Ifpack2

#define IFPACK2_CHEBYSHEV_INSTANT(S,LO,GO,N)                            \
  template class Ifpack2::Chebyshev< Tpetra::CrsMatrix<S, LO, GO, N> >; \
  template class Ifpack2::Chebyshev< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_CHEBYSHEV_DEF_HPP
