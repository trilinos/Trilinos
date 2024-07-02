// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_HIPTMAIR_DEF_HPP
#define IFPACK2_HIPTMAIR_DEF_HPP

#include "Ifpack2_Details_OneLevelFactory.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Details_residual.hpp"
#include <Tpetra_RowMatrixTransposer.hpp>
#include <cmath>
#include <iostream>
#include <sstream>

namespace Ifpack2 {

template <class MatrixType>
Hiptmair<MatrixType>::
Hiptmair (const Teuchos::RCP<const row_matrix_type>& A,
          const Teuchos::RCP<const row_matrix_type>& PtAP,
          const Teuchos::RCP<const row_matrix_type>& P,
          const Teuchos::RCP<const row_matrix_type>& Pt) :
  A_ (A),
  PtAP_ (PtAP),
  P_ (P),
  Pt_ (Pt),
  // Default values
  precType1_ ("CHEBYSHEV"),
  precType2_ ("CHEBYSHEV"),
  preOrPost_ ("both"),
  ZeroStartingSolution_ (true),
  ImplicitTranspose_ (Pt.is_null()),
  // General
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
Hiptmair<MatrixType>::
Hiptmair (const Teuchos::RCP<const row_matrix_type>& A):
  A_ (A),
  PtAP_ (),
  P_ (),
  Pt_ (),
  // Default values
  precType1_ ("CHEBYSHEV"),
  precType2_ ("CHEBYSHEV"),
  preOrPost_ ("both"),
  ZeroStartingSolution_ (true),
  ImplicitTranspose_ (true),
  // General
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
  using Teuchos::RCP;
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
  bool implicitTranspose = ImplicitTranspose_;

  precType1 = params.get("hiptmair: smoother type 1", precType1);
  precType2 = params.get("hiptmair: smoother type 2", precType2);
  precList1 = params.get("hiptmair: smoother list 1", precList1);
  precList2 = params.get("hiptmair: smoother list 2", precList2);
  preOrPost = params.get("hiptmair: pre or post",     preOrPost);
  zeroStartingSolution = params.get("hiptmair: zero starting solution",
                                    zeroStartingSolution);
  implicitTranspose = params.get("hiptmair: implicit transpose", implicitTranspose);

  // Grab the matrices off of the parameter list if we need them
  // This will intentionally throw if they're not there and we need them
  if(PtAP_.is_null()) 
    PtAP_ = params.get<RCP<row_matrix_type> >("PtAP");
  if(P_.is_null()) 
    P_ = params.get<RCP<row_matrix_type> >("P");
  if (params.isType<RCP<row_matrix_type> >("Pt"))
    Pt_ = params.get<RCP<row_matrix_type> >("Pt");


  // "Commit" the new values to the instance data.
  precType1_ = precType1;
  precType2_ = precType2;
  precList1_ = precList1;
  precList2_ = precList2;
  preOrPost_ = preOrPost;
  ZeroStartingSolution_ = zeroStartingSolution;
  ImplicitTranspose_ = implicitTranspose;
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
Teuchos::RCP<Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> > Hiptmair<MatrixType>::getPrec1() {
  return ifpack2_prec1_;
}


template <class MatrixType>
Teuchos::RCP<Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> > Hiptmair<MatrixType>::getPrec2() {
  return ifpack2_prec2_;
}


template <class MatrixType>
void Hiptmair<MatrixType>::initialize ()
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const char methodName[] = "Ifpack2::Hiptmair::initialize";

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Hiptmair::initialize: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");

  // clear any previous allocation
  IsInitialized_ = false;
  IsComputed_ = false;

  Teuchos::RCP<Teuchos::Time> timer =
    Teuchos::TimeMonitor::getNewCounter (methodName);

  double startTime = timer->wallTime();

  { // The body of code to time
    Teuchos::TimeMonitor timeMon (*timer);

    Details::OneLevelFactory<row_matrix_type> factory;

    ifpack2_prec1_=factory.create(precType1_,A_);
    ifpack2_prec1_->initialize();
    ifpack2_prec1_->setParameters(precList1_);

    ifpack2_prec2_=factory.create(precType2_,PtAP_);
    ifpack2_prec2_->initialize();
    ifpack2_prec2_->setParameters(precList2_);

  }
  IsInitialized_ = true;
  ++NumInitialize_;
  InitializeTime_ += (timer->wallTime() - startTime);
}


template <class MatrixType>
void Hiptmair<MatrixType>::compute ()
{
  const char methodName[] = "Ifpack2::Hiptmair::initialize";

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Hiptmair::compute: "
    "The input matrix A is null.  Please call setMatrix() with a nonnull "
    "input matrix before calling this method.");

  // Don't time the initialize(); that gets timed separately.
  if (! isInitialized ()) {
    initialize ();
  }

  Teuchos::RCP<Teuchos::Time> timer =
    Teuchos::TimeMonitor::getNewCounter (methodName);

  double startTime = timer->wallTime();
  { // The body of code to time
    Teuchos::TimeMonitor timeMon (*timer);
    ifpack2_prec1_->compute();
    ifpack2_prec2_->compute();

    if (!ImplicitTranspose_ && Pt_.is_null()) {
      using crs_type = Tpetra::CrsMatrix<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type>;
      Teuchos::RCP<const crs_type> crsP = Teuchos::rcp_dynamic_cast<const crs_type>(P_);
      if (!crsP.is_null()) {
        using transposer_type = Tpetra::RowMatrixTransposer<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type>;
        Pt_ = transposer_type(crsP).createTranspose();
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Ifpack2::Hiptmair::compute: "
                                   "ImplicitTranspose == false, but no Pt was provided and transposing P was not possible.");
      }
    }
  }
  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += (timer->wallTime() - startTime);
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

  const std::string timerName ("Ifpack2::Hiptmair::apply");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }
  double startTime = timer->wallTime();
  { // The body of code to time
    Teuchos::TimeMonitor timeMon (*timer);

    // If X and Y are pointing to the same memory location,
    // we need to create an auxiliary vector, Xcopy
    RCP<const MV> Xcopy;
    {
      if (X.aliases(Y)) {
        Xcopy = rcp (new MV (X, Teuchos::Copy));
      } else {
        Xcopy = rcpFromRef (X);
      }
    }

    RCP<MV> Ycopy = rcpFromRef (Y);
    if (ZeroStartingSolution_) {
      Ycopy->putScalar (STS::zero ());
    }

    // apply Hiptmair Smoothing
    applyHiptmairSmoother (*Xcopy, *Ycopy);

  }
  ++NumApply_;
  ApplyTime_ += (timer->wallTime() - startTime);
}


template<class MatrixType>
void
Hiptmair<MatrixType>::
updateCachedMultiVectors (const Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>>& map1,
                          const Teuchos::RCP<const Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type>>& map2,
                          size_t numVecs) const
{
  // Allocate a multivector if the cached one isn't perfect.  Checking
  // for map pointer equality is much cheaper than Map::isSameAs.
  using MV = Tpetra::MultiVector<scalar_type, local_ordinal_type,
                                 global_ordinal_type, node_type>;
  if (cachedResidual1_.is_null () ||
      map1.get () != cachedResidual1_->getMap ().get () ||
      cachedResidual1_->getNumVectors () != numVecs) {

    cachedResidual1_ = Teuchos::rcp (new MV (map1, numVecs, false));
    cachedSolution1_ = Teuchos::rcp (new MV (map1, numVecs, false));
  }
  if (cachedResidual2_.is_null () ||
      map2.get () != cachedResidual2_->getMap ().get () ||
      cachedResidual2_->getNumVectors () != numVecs) {

    cachedResidual2_ = Teuchos::rcp (new MV (map2, numVecs, false));
    cachedSolution2_ = Teuchos::rcp (new MV (map2, numVecs, false));
  }
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
  const scalar_type ZERO = STS::zero ();
  const scalar_type ONE = STS::one ();

  const std::string timerName1 ("Ifpack2::Hiptmair::apply 1");
  const std::string timerName2 ("Ifpack2::Hiptmair::apply 2");

  Teuchos::RCP<Teuchos::Time> timer1 = Teuchos::TimeMonitor::lookupCounter (timerName1);
  if (timer1.is_null ()) {
    timer1 = Teuchos::TimeMonitor::getNewCounter (timerName1);
  }
  Teuchos::RCP<Teuchos::Time> timer2 = Teuchos::TimeMonitor::lookupCounter (timerName2);
  if (timer2.is_null ()) {
    timer2 = Teuchos::TimeMonitor::getNewCounter (timerName2);
  }

  //#define IFPACK2_DEBUG_SMOOTHER
#ifdef IFPACK2_DEBUG_SMOOTHER
  int mypid = X.getMap()->getComm()->getRank();
  Teuchos::Array<double> ttt(1);
   printf("\n--------------------------------\n");
   printf("Coming into matrix Hiptmair\n");
   Y.norm2(ttt());
   if (!mypid) printf("\t||x|| = %15.10e\n", ttt[0]);
   X.norm2(ttt());
   if (!mypid) printf("\t||rhs|| = %15.10e\n", ttt[0]);
   {
     double normA = A_->getFrobeniusNorm();
     if (!mypid) printf("\t||A|| = %15.10e\n", normA);     
     Tpetra::Vector<typename MatrixType::scalar_type,
                    typename MatrixType::local_ordinal_type,
                    typename MatrixType::global_ordinal_type,
                    typename MatrixType::node_type> d(A_->getRowMap());
     A_->getLocalDiagCopy(d);
     d.norm2(ttt);
     if (!mypid) printf("\t||diag(A)|| = %15.10e\n", ttt[0]);
   }
   fflush(stdout);
#endif


  updateCachedMultiVectors (A_->getRowMap (),
                            PtAP_->getRowMap (),
                            X.getNumVectors ());

  if (preOrPost_ == "pre" || preOrPost_ == "both") {
    // apply initial relaxation to primary space
    Teuchos::TimeMonitor timeMon (*timer1);
    Tpetra::Details::residual(*A_,Y,X,*cachedResidual1_);
    cachedSolution1_->putScalar (ZERO);
    ifpack2_prec1_->apply (*cachedResidual1_, *cachedSolution1_);
    Y.update (ONE, *cachedSolution1_, ONE);
  }



  {
    // project to auxiliary space and smooth
    Teuchos::TimeMonitor timeMon (*timer2);
    Tpetra::Details::residual(*A_,Y,X,*cachedResidual1_);
#ifdef IFPACK2_DEBUG_SMOOTHER
      if (!mypid) printf("  After smoothing on edges\n");
      Y.norm2(ttt());
      if (!mypid) printf("\t||x|| = %15.10e\n", ttt[0]);
      cachedResidual1_->norm2(ttt());
      if (!mypid) printf("\t||res|| = %15.10e\n", ttt[0]);
#endif



    if (!Pt_.is_null())
      Pt_->apply (*cachedResidual1_, *cachedResidual2_, Teuchos::NO_TRANS);
    else
      P_->apply (*cachedResidual1_, *cachedResidual2_, Teuchos::TRANS);
    cachedSolution2_->putScalar (ZERO);

#ifdef IFPACK2_DEBUG_SMOOTHER
      if (!mypid)printf("  Before smoothing on nodes\n");
      cachedSolution2_->norm2(ttt());
      if (!mypid)printf("\t||x_nodal|| = %15.10e\n",ttt[0]);
      cachedResidual2_->norm2(ttt());
      if (!mypid)printf("\t||rhs_nodal|| = %15.10e\n", ttt[0]);
      {
        auto An = ifpack2_prec2_->getMatrix();
        double normA = An->getFrobeniusNorm();
        if (!mypid) printf("\t||An|| = %15.10e\n", normA);     
     Tpetra::Vector<typename MatrixType::scalar_type,
                    typename MatrixType::local_ordinal_type,
                    typename MatrixType::global_ordinal_type,
                    typename MatrixType::node_type> d(An->getRowMap());
     An->getLocalDiagCopy(d);
     d.norm2(ttt);
     if (!mypid) printf("\t||diag(An)|| = %15.10e\n", ttt[0]);
   }

#endif

    ifpack2_prec2_->apply (*cachedResidual2_, *cachedSolution2_);

#ifdef IFPACK2_DEBUG_SMOOTHER
      if (!mypid)printf("  After smoothing on nodes\n");
      cachedSolution2_->norm2(ttt());
      if (!mypid)printf("\t||x_nodal|| = %15.10e\n",ttt[0]);
      cachedResidual2_->norm2(ttt());
      if (!mypid)printf("\t||rhs_nodal|| = %15.10e\n", ttt[0]);
#endif


    P_->apply (*cachedSolution2_, Y, Teuchos::NO_TRANS, ONE, ONE);
  }

  if (preOrPost_ == "post" || preOrPost_ == "both") {
    // smooth again on primary space
    Teuchos::TimeMonitor timeMon (*timer1);
    Tpetra::Details::residual(*A_,Y,X,*cachedResidual1_);
    cachedSolution1_->putScalar (ZERO);
    ifpack2_prec1_->apply (*cachedResidual1_, *cachedSolution1_);
    Y.update (ONE, *cachedSolution1_, ONE);
  }

#ifdef IFPACK2_DEBUG_SMOOTHER
  if (!mypid)printf("  After updating edge solution\n");
  Y.norm2(ttt());
  if (!mypid)printf("\t||x|| = %15.10e\n",ttt[0]);
  if (!mypid)printf("--------------------------------\n");
#endif


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
    os << "Matrix: null, ";
  }
  else {
    os << "Matrix: not null"
       << ", Global matrix dimensions: ["
       << A_->getGlobalNumRows () << ", " << A_->getGlobalNumCols () << "], ";
  }

  os << "Smoother 1: ";
  os << ifpack2_prec1_->description() << ", ";
  os << "Smoother 2: ";
  os << ifpack2_prec2_->description();

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
    out << "Smoother 1: ";
    ifpack2_prec1_->describe(out, vl);
    out << "Smoother 2: ";
    ifpack2_prec2_->describe(out, vl);
  }
}

} // namespace Ifpack2

#define IFPACK2_HIPTMAIR_INSTANT(S,LO,GO,N) \
  template class Ifpack2::Hiptmair< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif /* IFPACK2_HIPTMAIR_DEF_HPP */
