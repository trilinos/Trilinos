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

#ifndef IFPACK2_HIPTMAIR_DEF_HPP
#define IFPACK2_HIPTMAIR_DEF_HPP

#include "Ifpack2_Hiptmair_decl.hpp"
#include "Ifpack2_Details_OneLevelFactory_decl.hpp"
#include "Ifpack2_Details_OneLevelFactory_def.hpp"

namespace Ifpack2 {

template <class MatrixType>
Hiptmair<MatrixType>::
Hiptmair (const Teuchos::RCP<const row_matrix_type>& A,
          const Teuchos::RCP<const row_matrix_type>& PtAP,
          const Teuchos::RCP<const row_matrix_type>& P) :
  A_ (A),
  PtAP_ (PtAP),
  P_ (P),
  // Default values
  precType1_ ("CHEBYSHEV"),
  precType2_ ("CHEBYSHEV"),
  preOrPost_ ("both"),
  ZeroStartingSolution_ (true),
  // General
  Condest_ (-STM::one ()),
  IsInitialized_ (false),
  IsComputed_ (false),
  NumInitialize_ (0),
  NumCompute_ (0),
  NumApply_ (0),
  InitializeTime_ (0.0),
  ComputeTime_ (0.0),
  ApplyTime_ (0.0)
{}


template <class MatrixType>
Hiptmair<MatrixType>::~Hiptmair() {}

template <class MatrixType>
void Hiptmair<MatrixType>::setParameters (const Teuchos::ParameterList& plist)
{
  using Teuchos::as;
  using Teuchos::ParameterList;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;

  ParameterList params = plist;

  // Get the current parameters' values.  We don't assign to the
  // instance data directly until we've gotten all the parameters.
  // This ensures "transactional" semantics, so that if attempting to
  // get some parameter throws an exception, the class' state doesn't
  // change.
  std::string precType1 = precType1_;
  std::string precType2 = precType2_;
  std::string preOrPost = preOrPost_;
  Teuchos::ParameterList precList1 = precList1_;
  Teuchos::ParameterList precList2 = precList2_;
  bool zeroStartingSolution = ZeroStartingSolution_;

  precType1 = params.get("hiptmair: smoother type 1", precType1);
  precType2 = params.get("hiptmair: smoother type 2", precType2);
  precList1 = params.get("hiptmair: smoother list 1", precList1);
  precList2 = params.get("hiptmair: smoother list 2", precList2);
  preOrPost = params.get("hiptmair: pre or post",     preOrPost);
  zeroStartingSolution = params.get("hiptmair: zero starting solution",
                                    zeroStartingSolution);

  // "Commit" the new values to the instance data.
  precType1_ = precType1;
  precType2_ = precType2;
  precList1_ = precList1;
  precList2_ = precList2;
  preOrPost_ = preOrPost;
  ZeroStartingSolution_ = zeroStartingSolution;
}


template <class MatrixType>
Teuchos::RCP<const Teuchos::Comm<int> >
Hiptmair<MatrixType>::getComm () const {
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Hiptmair::getComm: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");
  return A_->getComm ();
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Hiptmair<MatrixType>::getMatrix () const {
  return A_;
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Hiptmair<MatrixType>::getDomainMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Hiptmair::getDomainMap: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");
  return A_->getDomainMap ();
}


template <class MatrixType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Hiptmair<MatrixType>::getRangeMap () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Hiptmair::getRangeMap: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");
  return A_->getRangeMap ();
}


template <class MatrixType>
bool Hiptmair<MatrixType>::hasTransposeApply () const {
  // FIXME (mfh 17 Jan 2014) apply() does not currently work with mode
  // != NO_TRANS, so it's correct to return false here.
  return false;
}


template <class MatrixType>
int Hiptmair<MatrixType>::getNumInitialize () const {
  return NumInitialize_;
}


template <class MatrixType>
int Hiptmair<MatrixType>::getNumCompute () const {
  return NumCompute_;
}


template <class MatrixType>
int Hiptmair<MatrixType>::getNumApply () const {
  return NumApply_;
}


template <class MatrixType>
double Hiptmair<MatrixType>::getInitializeTime () const {
  return InitializeTime_;
}


template <class MatrixType>
double Hiptmair<MatrixType>::getComputeTime () const {
  return ComputeTime_;
}


template <class MatrixType>
double Hiptmair<MatrixType>::getApplyTime () const {
  return ApplyTime_;
}


template <class MatrixType>
typename Hiptmair<MatrixType>::magnitude_type
Hiptmair<MatrixType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const row_matrix_type>& matrix)
{
  if (! isComputed ()) { // cannot compute right now
    return -STM::one ();
  }
  // NOTE: this is computing the *local* condest
  if (Condest_ == -STM::one ()) {
    Condest_ = Ifpack2::Condest (*this, CT, MaxIters, Tol, matrix);
  }
  return Condest_;
}


template <class MatrixType>
void Hiptmair<MatrixType>::initialize ()
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Hiptmair::initialize: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");

  // clear any previous allocation
  IsInitialized_ = false;
  IsComputed_ = false;

  Teuchos::Time timer ("initialize");
  { // The body of code to time
    Teuchos::TimeMonitor timeMon (timer);

    Details::OneLevelFactory<MatrixType> factory;

    ifpack2_prec1_=factory.create(precType1_,A_);
    ifpack2_prec1_->initialize();
    ifpack2_prec1_->setParameters(precList1_);

    ifpack2_prec2_=factory.create(precType2_,PtAP_);
    ifpack2_prec2_->initialize();
    ifpack2_prec2_->setParameters(precList2_);

  }
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += timer.totalElapsedTime ();
}


template <class MatrixType>
void Hiptmair<MatrixType>::compute ()
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Hiptmair::compute: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");

  // Don't time the initialize(); that gets timed separately.
  if (! isInitialized ()) {
    initialize ();
  }

  Teuchos::Time timer ("compute");
  { // The body of code to time
    Teuchos::TimeMonitor timeMon (timer);
    ifpack2_prec1_->compute();
    ifpack2_prec2_->compute();
  }
  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += timer.totalElapsedTime ();
}


template <class MatrixType>
void Hiptmair<MatrixType>::
apply (const Tpetra::MultiVector<typename MatrixType::scalar_type,
       typename MatrixType::local_ordinal_type,
       typename MatrixType::global_ordinal_type,
       typename MatrixType::node_type>& X,
       Tpetra::MultiVector<typename MatrixType::scalar_type,
                           typename MatrixType::local_ordinal_type,
                           typename MatrixType::global_ordinal_type,
                           typename MatrixType::node_type>& Y,
       Teuchos::ETransp mode,
       typename MatrixType::scalar_type alpha,
       typename MatrixType::scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
                              global_ordinal_type, node_type> MV;
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! isComputed (), std::runtime_error,
    "Ifpack2::Hiptmair::apply: You must call compute() before you may call apply().");
  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors () != Y.getNumVectors (), std::invalid_argument,
    "Ifpack2::Hiptmair::apply: The MultiVector inputs X and Y do not have the "
    "same number of columns.  X.getNumVectors() = " << X.getNumVectors ()
    << " != Y.getNumVectors() = " << Y.getNumVectors () << ".");

  // Catch unimplemented cases: alpha != 1, beta != 0, mode != NO_TRANS.
  TEUCHOS_TEST_FOR_EXCEPTION(
    alpha != STS::one (), std::logic_error,
    "Ifpack2::Hiptmair::apply: alpha != 1 has not been implemented.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    beta != STS::zero (), std::logic_error,
    "Ifpack2::Hiptmair::apply: zero != 0 has not been implemented.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    mode != Teuchos::NO_TRANS, std::logic_error,
    "Ifpack2::Hiptmair::apply: mode != Teuchos::NO_TRANS has not been implemented.");

  Teuchos::Time timer ("apply");
  { // The body of code to time
    Teuchos::TimeMonitor timeMon (timer);

    // If X and Y are pointing to the same memory location,
    // we need to create an auxiliary vector, Xcopy
    RCP<const MV> Xcopy;
    if (X.getLocalMV ().getValues () == Y.getLocalMV ().getValues ()) {
      Xcopy = rcp (new MV (X, Teuchos::Copy));
    } else {
      Xcopy = rcpFromRef (X);
    }

    RCP<MV> Ycopy = rcpFromRef (Y);
    if (ZeroStartingSolution_) {
      Ycopy->putScalar (STS::zero ());
    }

    // apply Hiptmair Smoothing
    applyHiptmairSmoother (*Xcopy, *Ycopy);

  }
  ++NumApply_;
  ApplyTime_ += timer.totalElapsedTime ();
}

template <class MatrixType>
void Hiptmair<MatrixType>::
applyHiptmairSmoother(const Tpetra::MultiVector<typename MatrixType::scalar_type,
                      typename MatrixType::local_ordinal_type,
                      typename MatrixType::global_ordinal_type,
                      typename MatrixType::node_type>& X,
                      Tpetra::MultiVector<typename MatrixType::scalar_type,
                      typename MatrixType::local_ordinal_type,
                      typename MatrixType::global_ordinal_type,
                      typename MatrixType::node_type>& Y) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> MV;
  const scalar_type ZERO = STS::zero ();
  const scalar_type ONE = STS::one ();

  RCP<MV> res1 = rcp (new MV (A_->getRowMap (), X.getNumVectors ()));
  RCP<MV> vec1 = rcp (new MV (A_->getRowMap (), X.getNumVectors ()));
  RCP<MV> res2 = rcp (new MV (PtAP_->getRowMap (), X.getNumVectors ()));
  RCP<MV> vec2 = rcp (new MV (PtAP_->getRowMap (), X.getNumVectors ()));

  if (preOrPost_ == "pre" || preOrPost_ == "both") {
    // apply initial relaxation to primary space
    A_->apply (Y, *res1);
    res1->update (ONE, X, -ONE);
    vec1->putScalar (ZERO);
    ifpack2_prec1_->apply (*res1, *vec1);
    Y.update (ONE, *vec1, ONE);
  }

  // project to auxiliary space and smooth
  A_->apply (Y, *res1);
  res1->update (ONE, X, -ONE);
  P_->apply (*res1, *res2, Teuchos::TRANS);
  vec2->putScalar (ZERO);
  ifpack2_prec2_->apply (*res2, *vec2);
  P_->apply (*vec2, *vec1, Teuchos::NO_TRANS);
  Y.update (ONE,*vec1,ONE);

  if (preOrPost_ == "post" || preOrPost_ == "both") {
    // smooth again on primary space
    A_->apply (Y, *res1);
    res1->update (ONE, X, -ONE);
    vec1->putScalar (ZERO);
    ifpack2_prec1_->apply (*res1, *vec1);
    Y.update (ONE, *vec1, ONE);
  }
}

template <class MatrixType>
std::string Hiptmair<MatrixType>::description () const
{
  std::ostringstream os;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  os << "\"Ifpack2::Hiptmair\": {";
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
       << A_->getGlobalNumRows () << ", " << A_->getGlobalNumCols () << "]";
  }

  os << "}";
  return os.str ();
}


template <class MatrixType>
void Hiptmair<MatrixType>::
describe (Teuchos::FancyOStream &out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

  if (vl != VERB_NONE) {
    // describe() always starts with a tab by convention.
    Teuchos::OSTab tab0 (out);
    out << "\"Ifpack2::Hiptmair\":";

    Teuchos::OSTab tab1 (out);
    if (this->getObjectLabel () != "") {
      out << "Label: " << this->getObjectLabel () << endl;
    }
    out << "Initialized: " << (isInitialized () ? "true" : "false") << endl
        << "Computed: " << (isComputed () ? "true" : "false") << endl
        << "Global number of rows: " << A_->getGlobalNumRows () << endl
        << "Global number of columns: " << A_->getGlobalNumCols () << endl
        << "Matrix:";
    if (A_.is_null ()) {
      out << " null" << endl;
    } else {
      A_->describe (out, vl);
    }
  }
}

} // namespace Ifpack2

#define IFPACK2_HIPTMAIR_INSTANT(S,LO,GO,N) \
  template class Ifpack2::Hiptmair< Tpetra::CrsMatrix<S, LO, GO, N> >;

#endif /* IFPACK2_HIPTMAIR_DEF_HPP */
