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

//-----------------------------------------------------
// Ifpack2::KRYLOV is based on the Krylov iterative
// solvers in Belos.
// written by Paul Tsuji.
//-----------------------------------------------------

#ifndef IFPACK2_KRYLOV_DEF_HPP
#define IFPACK2_KRYLOV_DEF_HPP

namespace Ifpack2 {

//Definitions for the Krylov methods:
//==============================================================================
template <class MatrixType, class PrecType>
Krylov<MatrixType,PrecType>::Krylov(const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> >& A) :
  A_(A),
  Comm_(A->getRowMap()->getComm()),
  // Default values
  IterationType_(1),
  Iterations_(5),
  ResidualTolerance_(0.001),
  BlockSize_(1),
  ZeroStartingSolution_(true),
  PreconditionerType_(1),
  // General
  Condest_ (- Teuchos::ScalarTraits<magnitude_type>::one()),
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0),
  Time_("Ifpack2::Krylov"),
  NumMyRows_(-1),
  NumGlobalNonzeros_(0)
{
  TEUCHOS_TEST_FOR_EXCEPTION(A_ == Teuchos::null, std::runtime_error,
      Teuchos::typeName(*this) << "::Krylov(): input matrix reference was null.");
}

//==========================================================================
template <class MatrixType, class PrecType>
Krylov<MatrixType,PrecType>::~Krylov() {
}

//==========================================================================
template <class MatrixType, class PrecType>
void Krylov<MatrixType,PrecType>::setParameters(const Teuchos::ParameterList& params) {
  using Teuchos::as;
  using Teuchos::Exceptions::InvalidParameterName;
  using Teuchos::Exceptions::InvalidParameterType;
  // Read in parameters
  Ifpack2::getParameter(params, "krylov: iteration type",IterationType_);
  Ifpack2::getParameter(params, "krylov: number of iterations",Iterations_);
  Ifpack2::getParameter(params, "krylov: residual tolerance",ResidualTolerance_);
  Ifpack2::getParameter(params, "krylov: block size",BlockSize_);
  Ifpack2::getParameter(params, "krylov: zero starting solution",ZeroStartingSolution_);
  Ifpack2::getParameter(params, "krylov: preconditioner type",PreconditionerType_);
  params_=params;
  // Separate preconditioner parameters into another list
  if(PreconditionerType_==1) {
    precParams_.set("relaxation: sweeps",                 params_.get("relaxation: sweeps",1));
    precParams_.set("relaxation: damping factor",         params_.get("relaxation: damping factor",(scalar_type)1.0));
    precParams_.set("relaxation: min diagonal value",     params_.get("relaxation: min diagonal value",(scalar_type)1.0));
    precParams_.set("relaxation: zero starting solution", params_.get("relaxation: zero starting solution",true));
    precParams_.set("relaxation: backward mode",          params_.get("relaxation: backward mode",false));
  }
  if(PreconditionerType_==2 || PreconditionerType_==3) {
    precParams_.set("fact: ilut level-of-fill", params_.get("fact: ilut level-of-fill",(double)1.0));
    precParams_.set("fact: absolute threshold", params_.get("fact: absolute threshold",(double)0.0));
    precParams_.set("fact: relative threshold", params_.get("fact: relative threshold",(double)1.0));
    precParams_.set("fact: relax value",        params_.get("fact: relax value",(double)0.0));
  }
  if(PreconditionerType_==3) {
    precParams_.set("schwarz: compute condest",   params_.get("schwarz: compute condest",true));
    precParams_.set("schwarz: combine mode",      params_.get("schwarz: combine mode","Zero")); // use string mode for this
    precParams_.set("schwarz: use reordering",    params_.get("schwarz: use reordering",true));
    precParams_.set("schwarz: filter singletons", params_.get("schwarz: filter singletons",false));
    precParams_.set("schwarz: overlap level",     params_.get("schwarz: overlap level",(int)0));
  }
}

//==========================================================================
template <class MatrixType, class PrecType>
Teuchos::RCP<const Teuchos::Comm<int> >
Krylov<MatrixType,PrecType>::getComm () const {
  return Comm_;
}

//==========================================================================
template <class MatrixType, class PrecType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Krylov<MatrixType,PrecType>::getMatrix () const {
  return A_;
}

//==========================================================================
template <class MatrixType, class PrecType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Krylov<MatrixType,PrecType>::getDomainMap () const
{
  return A_->getDomainMap ();
}

//==========================================================================
template <class MatrixType, class PrecType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Krylov<MatrixType,PrecType>::getRangeMap () const
{
  return A_->getRangeMap ();
}

//==============================================================================
template <class MatrixType, class PrecType>
bool Krylov<MatrixType,PrecType>::hasTransposeApply () const {
  return true;
}

//==========================================================================
template <class MatrixType, class PrecType>
int Krylov<MatrixType,PrecType>::getNumInitialize () const {
  return NumInitialize_;
}

//==========================================================================
template <class MatrixType, class PrecType>
int Krylov<MatrixType,PrecType>::getNumCompute() const {
  return(NumCompute_);
}

//==========================================================================
template <class MatrixType, class PrecType>
int Krylov<MatrixType,PrecType>::getNumApply() const {
  return(NumApply_);
}

//==========================================================================
template <class MatrixType, class PrecType>
double Krylov<MatrixType,PrecType>::getInitializeTime() const {
  return(InitializeTime_);
}

//==========================================================================
template <class MatrixType, class PrecType>
double Krylov<MatrixType,PrecType>::getComputeTime() const {
  return(ComputeTime_);
}

//==========================================================================
template <class MatrixType, class PrecType>
double Krylov<MatrixType,PrecType>::getApplyTime() const {
  return(ApplyTime_);
}

//=============================================================================
template <class MatrixType, class PrecType>
typename Krylov<MatrixType,PrecType>::magnitude_type
Krylov<MatrixType,PrecType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &matrix) {
  if (! isComputed ()) { // cannot compute right now
    return -Teuchos::ScalarTraits<magnitude_type>::one ();
  }
  // NOTE: this is computing the *local* condest
  if (Condest_ == -Teuchos::ScalarTraits<magnitude_type>::one()) {
    Condest_ = Ifpack2::Condest(*this, CT, MaxIters, Tol, matrix);
  }
  return Condest_;
}

//==========================================================================
template <class MatrixType, class PrecType>
void Krylov<MatrixType,PrecType>::initialize() {
  // clear any previous allocation
  IsInitialized_ = false;
  IsComputed_ = false;
  Time_.start(true);
  // check only in serial
  TEUCHOS_TEST_FOR_EXCEPTION(Comm_->getSize() == 1 && A_->getNodeNumRows() != A_->getNodeNumCols(), std::runtime_error, "Ifpack2::Krylov::Initialize ERROR, matrix must be square");
  NumMyRows_ = A_->getNodeNumRows();
  // Belos parameter list
  belosList_ = Teuchos::rcp( new Teuchos::ParameterList("GMRES") );
  belosList_->set("Maximum Iterations",Iterations_);
  belosList_->set("Convergence Tolerance",ResidualTolerance_);
  if(PreconditionerType_==0) {
    // no preconditioner
  }
  else if(PreconditionerType_==1) {
    ifpack2_prec_=Teuchos::rcp( new Relaxation<MatrixType> (A_) );
  }
  else if(PreconditionerType_==2) {
    ifpack2_prec_=Teuchos::rcp( new ILUT<MatrixType>(A_) );
  }
  else if(PreconditionerType_==3) {
    ifpack2_prec_ = Teuchos::rcp (new AdditiveSchwarz<MatrixType, ILUT<MatrixType> > (A_));
  }
  else if(PreconditionerType_==4) {
    ifpack2_prec_=Teuchos::rcp( new Chebyshev<MatrixType>(A_) );
  }
  if(PreconditionerType_>0) {
    ifpack2_prec_->initialize();
    ifpack2_prec_->setParameters(precParams_);
  }
  belosProblem_ = Teuchos::rcp( new Belos::LinearProblem<scalar_type,TMV,TOP> );
  belosProblem_ -> setOperator(A_);
  if(IterationType_==1) {
    belosSolver_ =
      Teuchos::rcp( new Belos::BlockGmresSolMgr<scalar_type,TMV,TOP> (belosProblem_, belosList_) );
  }
  else {
    belosSolver_ =
      Teuchos::rcp( new Belos::BlockCGSolMgr<scalar_type,TMV,TOP> (belosProblem_, belosList_) );
  }
  IsInitialized_ = true;
  ++NumInitialize_;
  Time_.stop();
  InitializeTime_ += Time_.totalElapsedTime();
}

//==========================================================================
template <class MatrixType, class PrecType>
void Krylov<MatrixType,PrecType>::compute() {
  using Teuchos::as;
  //--------------------------------------------------------------------------
  // Ifpack2::Krylov
  // Computes necessary preconditioner info
  //--------------------------------------------------------------------------
  if (!isInitialized()) {
    initialize();
  }
  Time_.start(true);
  if(PreconditionerType_>0) {
    ifpack2_prec_->compute();
    belosProblem_->setLeftPrec(ifpack2_prec_);
  }
  IsComputed_ = true;
  ++NumCompute_;
  Time_.stop();
  ComputeTime_ += Time_.totalElapsedTime();
}

//==========================================================================
template <class MatrixType, class PrecType>
void Krylov<MatrixType,PrecType>::apply(const Tpetra::MultiVector<typename MatrixType::scalar_type,
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
  TEUCHOS_TEST_FOR_EXCEPTION(!isComputed(), std::runtime_error,
    "Ifpack2::Krylov::apply() ERROR, compute() hasn't been called yet.");

  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
    "Ifpack2::Krylov::apply() ERROR, X.getNumVectors() != Y.getNumVectors().");

  Time_.start(true);

  // If X and Y are pointing to the same memory location,
  // we need to create an auxiliary vector, Xcopy
  Teuchos::RCP< const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Xcopy;
  if (X.getLocalMV().getValues() == Y.getLocalMV().getValues()) {
    Xcopy = Teuchos::rcp( new Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>(X) );
  }
  else {
    Xcopy = Teuchos::rcp( &X, false );
  }

  Teuchos::RCP< Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > Ycopy = Teuchos::rcpFromRef(Y);
  if(ZeroStartingSolution_==true) {
    Ycopy->putScalar((scalar_type) 0.0);
  }
  
  // Set left and right hand sides for Belos
  belosProblem_->setProblem(Ycopy,Xcopy);
  // iterative solve
  belosSolver_->solve();

  ++NumApply_;
  Time_.stop();
  ApplyTime_ += Time_.totalElapsedTime();
}

//=============================================================================
template <class MatrixType, class PrecType>
std::string Krylov<MatrixType,PrecType>::description() const {
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed";
    }
    else {
      oss << "{status = initialized, not computed";
    }
  }
  else {
    oss << "{status = not initialized, not computed";
  }
  oss << ", global rows = " << A_->getGlobalNumRows()
      << ", global cols = " << A_->getGlobalNumCols()
      << "}";
  return oss.str();
}

//=============================================================================
template <class MatrixType, class PrecType>
void Krylov<MatrixType,PrecType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  using std::setw;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  Teuchos::EVerbosityLevel vl = verbLevel;
  if (vl == VERB_DEFAULT) vl = VERB_LOW;
  const int myImageID = Comm_->getRank();
  Teuchos::OSTab tab(out);
  //    none: print nothing
  //     low: print O(1) info from node 0
  //  medium:
  //    high:
  // extreme:
  if (vl != VERB_NONE && myImageID == 0) {
  }
}

}//namespace Ifpack2

#endif /* IFPACK2_KRYLOV_DEF_HPP */
