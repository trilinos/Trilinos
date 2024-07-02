// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_CHEBYSHEV_DEF_HPP
#define IFPACK2_CHEBYSHEV_DEF_HPP

#include "Ifpack2_Parameters.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include <iostream>
#include <sstream>


namespace Ifpack2 {

template<class MatrixType>
Chebyshev<MatrixType>::
Chebyshev (const Teuchos::RCP<const row_matrix_type>& A)
  : impl_ (A),
    IsInitialized_ (false),
    IsComputed_ (false),
    NumInitialize_ (0),
    NumCompute_ (0),
    TimerForApply_(true),
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
  if (List.isType<bool>("timer for apply"))
    TimerForApply_ = List.get<bool>("timer for apply");
}


template<class MatrixType>
void
Chebyshev<MatrixType>::setZeroStartingSolution (bool zeroStartingSolution)
{
  impl_.setZeroStartingSolution(zeroStartingSolution);
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
Teuchos::RCP<const Tpetra::CrsMatrix<typename MatrixType::scalar_type,
                                     typename MatrixType::local_ordinal_type,
                                     typename MatrixType::global_ordinal_type,
                                     typename MatrixType::node_type> >
Chebyshev<MatrixType>::
getCrsMatrix() const {
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> crs_matrix_type;
  return Teuchos::rcp_dynamic_cast<const crs_matrix_type> (impl_.getMatrix ());
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
size_t Chebyshev<MatrixType>::getNodeSmootherComplexity() const {
  Teuchos::RCP<const row_matrix_type> A = impl_.getMatrix();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null (), std::runtime_error, "Ifpack2::Chevyshev::getNodeSmootherComplexity: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix, then call compute(), before calling this method.");
  // Chevyshev costs roughly one apply + one diagonal inverse per iteration
  return A->getLocalNumRows() + A->getLocalNumEntries();
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
  Teuchos::RCP<Teuchos::Time> timer;
  const std::string timerName ("Ifpack2::Chebyshev::apply");
  if (TimerForApply_) {
    timer = Teuchos::TimeMonitor::lookupCounter (timerName);
    if (timer.is_null ()) {
      timer = Teuchos::TimeMonitor::getNewCounter (timerName);
    }
  }

  Teuchos::Time time = Teuchos::Time(timerName);
  double startTime = time.wallTime();

  // Start timing here.
  {
    Teuchos::RCP<Teuchos::TimeMonitor> timeMon;
    if (TimerForApply_)
      timeMon = Teuchos::rcp(new Teuchos::TimeMonitor(*timer));

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
  ApplyTime_ += (time.wallTime() - startTime);
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

  double startTime = timer->wallTime();

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

  ComputeTime_ += (timer->wallTime() - startTime);
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
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::TypeNameTraits;
  using std::endl;

  // Default verbosity level is VERB_LOW
  const Teuchos::EVerbosityLevel vl =
    (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;

  if (vl == Teuchos::VERB_NONE) {
    return; // print NOTHING, not even the class name
  }

  // By convention, describe() starts with a tab.
  //
  // This does affect all processes on which it's valid to print to
  // 'out'.  However, it does not actually print spaces to 'out'
  // unless operator<< gets called, so it's safe to use on all
  // processes.
  Teuchos::OSTab tab0 (out);
  const int myRank = this->getComm ()->getRank ();
  if (myRank == 0) {
    // Output is a valid YAML dictionary.
    // In particular, we quote keys with colons in them.
    out << "\"Ifpack2::Chebyshev\":" << endl;
  }

  Teuchos::OSTab tab1 (out);
  if (vl >= Teuchos::VERB_LOW && myRank == 0) {
    out << "Template parameters:" << endl;
    {
      Teuchos::OSTab tab2 (out);
      out << "Scalar: " << TypeNameTraits<scalar_type>::name () << endl
          << "LocalOrdinal: " << TypeNameTraits<local_ordinal_type>::name () << endl
          << "GlobalOrdinal: " << TypeNameTraits<global_ordinal_type>::name () << endl
          << "Device: " << TypeNameTraits<device_type>::name () << endl;
    }
    out << "Initialized: " << (isInitialized () ? "true" : "false") << endl
        << "Computed: " << (isComputed () ? "true" : "false") << endl;
    impl_.describe (out, vl);

    if (impl_.getMatrix ().is_null ()) {
      out << "Matrix: null" << endl;
    }
    else {
      out << "Global matrix dimensions: ["
          << impl_.getMatrix ()->getGlobalNumRows () << ", "
          << impl_.getMatrix ()->getGlobalNumCols () << "]" << endl
          << "Global nnz: " << impl_.getMatrix ()->getGlobalNumEntries() << endl;
    }
  }
}

template<class MatrixType>
void
Chebyshev<MatrixType>::
applyImpl (const MV& X,
           MV& Y,
           Teuchos::ETransp /* mode */,
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
  if (X.aliases(Y)) {
    X_copy = rcp (new MV (X, Teuchos::Copy));
    copiedInput = true;
  } else {
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
  template class Ifpack2::Chebyshev< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_CHEBYSHEV_DEF_HPP
