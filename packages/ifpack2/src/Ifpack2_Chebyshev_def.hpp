/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_CHEBYSHEV_DEF_HPP
#define IFPACK2_CHEBYSHEV_DEF_HPP


namespace Ifpack2 {

template<class MatrixType>
Chebyshev<MatrixType>::
Chebyshev (const Teuchos::RCP<const row_matrix_type>& A)
  : impl_ (A),
    Time_ (Teuchos::rcp (new Teuchos::Time ("Ifpack2::Chebyshev"))),
    Condest_ (-1.0),
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
{this->setObjectLabel("Ifpack2::Chebyshev");}

//==========================================================================
template<class MatrixType>
Chebyshev<MatrixType>::~Chebyshev() {
}

//==========================================================================
template<class MatrixType>
void 
Chebyshev<MatrixType>::setParameters (const Teuchos::ParameterList& List) 
{
  // FIXME (mfh 25 Jan 2013) Casting away const is bad here.
  impl_.setParameters (const_cast<Teuchos::ParameterList&> (List));
}


//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & 
Chebyshev<MatrixType>::getComm() const {
  return impl_.getMatrix ()->getRowMap ()->getComm ();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const typename Chebyshev<MatrixType>::row_matrix_type>
Chebyshev<MatrixType>::
getMatrix() const {
  return impl_.getMatrix ();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const typename Chebyshev<MatrixType>::map_type>&
Chebyshev<MatrixType>::
getDomainMap() const {
  return impl_.getMatrix ()->getDomainMap ();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const typename Chebyshev<MatrixType>::map_type>&
Chebyshev<MatrixType>::
getRangeMap() const {
  return impl_.getMatrix ()->getRangeMap ();
}

//==========================================================================
template<class MatrixType>
bool Chebyshev<MatrixType>::hasTransposeApply() const {
  return impl_.hasTransposeApply ();
}

//==========================================================================
template<class MatrixType>
int Chebyshev<MatrixType>::getNumInitialize() const {
  return(NumInitialize_);
}

//==========================================================================
template<class MatrixType>
int Chebyshev<MatrixType>::getNumCompute() const {
  return(NumCompute_);
}

//==========================================================================
template<class MatrixType>
int Chebyshev<MatrixType>::getNumApply() const {
  return(NumApply_);
}

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getInitializeTime() const {
  return(InitializeTime_);
}

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getComputeTime() const {
  return(ComputeTime_);
}

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getApplyTime() const {
  return(ApplyTime_);
}

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getComputeFlops() const {
  return(ComputeFlops_);
}

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getApplyFlops() const {
  return(ApplyFlops_);
}

//==========================================================================
template<class MatrixType>
typename Chebyshev<MatrixType>::magnitude_type
Chebyshev<MatrixType>::
getCondEst() const {
  return(Condest_);
}

//==========================================================================
template<class MatrixType>
typename Chebyshev<MatrixType>::magnitude_type
Chebyshev<MatrixType>::
computeCondEst (CondestType CT,
		local_ordinal_type MaxIters, 
		magnitude_type Tol,
		const Teuchos::Ptr<const row_matrix_type>& matrix) 
{
  if (! isComputed ()) {
    return -Teuchos::ScalarTraits<magnitude_type>::one ();
  }
  else {
    // Always compute it. Call Condest() with no parameters to get
    // the previous estimate.
    Condest_ = Ifpack2::Condest (*this, CT, MaxIters, Tol, matrix);
    return Condest_;
  }
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
  {
    Teuchos::TimeMonitor timeMon (*Time_);

    // compute() calls initialize() if it hasn't already been called.
    // Thus, we only need to check isComputed().
    TEUCHOS_TEST_FOR_EXCEPTION(! isComputed(), std::runtime_error, 
      "Ifpack2::Chebyshev::apply(): You must call the compute() method before "
      "you may call apply().");
    TEUCHOS_TEST_FOR_EXCEPTION(
       X.getNumVectors() != Y.getNumVectors(), 
       std::runtime_error,
       "Ifpack2::Chebyshev::apply(): X and Y must have the same number of "
       "columns.  X.getNumVectors() = " << X.getNumVectors() << " != "
       << "Y.getNumVectors() = " << Y.getNumVectors() << ".");
#ifdef HAVE_TEUCHOS_DEBUG
    {
      // The relation 'isSameAs' is transitive.  It's also a collective,
      // so we don't have to do a "shared" test for exception (i.e., a
      // global reduction on the test value).
      TEUCHOS_TEST_FOR_EXCEPTION(
         ! X.getMap ()->isSameAs (*getDomainMap ()),
         std::runtime_error,
         "Ifpack2::Chebyshev: The domain Map of the matrix must be the same as "
	 "the Map of the input vector(s) X.");
      TEUCHOS_TEST_FOR_EXCEPTION(
         ! Y.getMap ()->isSameAs (*getRangeMap ()),
         std::runtime_error,
         "Ifpack2::Chebyshev: The range Map of the matrix must be the same as "
	 "the Map of the output vector(s) Y.");
    }
#endif // HAVE_TEUCHOS_DEBUG
    applyImpl (X, Y, mode, alpha, beta);
  }
  ++NumApply_;
  ApplyTime_ += Time_->totalElapsedTime ();
}


template<class MatrixType>
void 
Chebyshev<MatrixType>::
applyMat (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
	  Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
	  Teuchos::ETransp mode) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
   "Ifpack2::Chebyshev::applyMat(): X.getNumVectors() != Y.getNumVectors().");
  impl_.getMatrix ()->apply (X, Y, mode);
}


template<class MatrixType>
void Chebyshev<MatrixType>::initialize() {
  // This method doesn't do anything anymore, so there's no need to time it.
  IsInitialized_ = true;
  ++NumInitialize_;
  // Note to developers: Defer fetching any data that relate to the
  // structure of the matrix until compute().  That way, it will
  // always be correct to omit calling initialize(), even if the
  // number of entries in the sparse matrix has changed since the last
  // call to compute().
}


template<class MatrixType>
void Chebyshev<MatrixType>::compute()
{
  {
    Teuchos::TimeMonitor timeMon (*Time_);
    if (! isInitialized ()) {
      initialize ();
    }
    IsComputed_ = false;
    Condest_ = -1.0;  
    impl_.compute ();
  }
  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += Time_->totalElapsedTime();
}


template<class MatrixType>
void 
Chebyshev<MatrixType>::
PowerMethod (const Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Operator, 
	     const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& InvPointDiagonal, 
	     const int MaximumIterations, 
	     scalar_type& lambda_max)
{
  const scalar_type one = STS::one();
  const scalar_type zero = STS::zero();

  lambda_max = zero;
  Teuchos::Array<scalar_type> RQ_top(1), RQ_bottom(1);
  vector_type x (Operator.getDomainMap ());
  vector_type y (Operator.getRangeMap ());
  x.randomize ();
  Teuchos::Array<magnitude_type> norms (x.getNumVectors ());
  x.norm2 (norms ());
  x.scale (one / norms[0]);

  for (int iter = 0; iter < MaximumIterations; ++iter) {
    Operator.apply (x, y);
    y.elementWiseMultiply (one, InvPointDiagonal, y, zero);
    y.dot (x, RQ_top ());
    x.dot (x, RQ_bottom ());
    lambda_max = RQ_top[0] / RQ_bottom[0];
    y.norm2 (norms ());
    TEUCHOS_TEST_FOR_EXCEPTION(
      norms[0] == zero, 
      std::runtime_error, 
      "Ifpack2::Chebyshev::PowerMethod: norm == 0 at iteration " << (iter+1) 
      << " of " << MaximumIterations);
    x.update (one / norms[0], y, zero);
  }
}

//==========================================================================
template<class MatrixType>
void Chebyshev<MatrixType>::
CG(const Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Operator, 
            const Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& InvPointDiagonal, 
   const int MaximumIterations, 
   scalar_type& lambda_min, scalar_type& lambda_max)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::logic_error,
    "Ifpack2::Chebyshev::CG: Not implemented.  "
    "Please use Belos' implementation of CG with Tpetra objects.");
}

//==========================================================================
template <class MatrixType>
std::string Chebyshev<MatrixType>::description() const {
  std::ostringstream oss;
  oss << Teuchos::LabeledObject::getObjectLabel();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed, ";
    }
    else {
      oss << "{status = initialized, not computed, ";
    }
  }
  else {
    oss << "{status = not initialized, not computed, ";
  }

  oss << impl_.description();

  oss << ", global rows = " << impl_.getMatrix ()->getGlobalNumRows()
      << ", global cols = " << impl_.getMatrix ()->getGlobalNumCols()
      << "}";
  return oss.str();
}

//==========================================================================
template <class MatrixType>
void Chebyshev<MatrixType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
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
    if   (Condest_ == -1.0) { out << "Condition number estimate       = N/A" << endl; }
    else                    { out << "Condition number estimate       = " << Condest_ << endl; }
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

  // If beta != 0, then we need to keep a copy of the initial value of
  // Y, so that we can add beta*it to the Chebyshev result at the end.
  // Usually this method is called with beta == 0, so we don't have to 
  // worry about caching Y_org.
  RCP<MV> Y_orig;
  if (beta != zero) {
    Y_orig = rcp (new MV (Y));
  }

  // If X and Y point to the same memory location, we need to use a
  // copy of X (X_copy) as the input MV.  Otherwise, just let X_copy
  // point to X.
  //
  // This is hopefully an uncommon use case, so we don't bother to
  // optimize for it by caching X_copy.
  RCP<const MV> X_copy;
  bool copiedInput = false;
  if (X.getLocalMV().getValues() == Y.getLocalMV().getValues()) {
    X_copy = rcp (new MV (X));
    copiedInput = true;
  }
  else {
    X_copy = rcpFromRef (X);
  }
  
  // If alpha != 1, fold alpha into (a copy of) X.
  //
  // This is an uncommon use case, so we don't bother to optimize for
  // it by caching X_copy.  However, we do check whether we've already
  // copied X above, to avoid a second copy.
  if (alpha != one) {
    RCP<MV> X_copy_nonConst = rcp_const_cast<MV> (X_copy);
    if (! copiedInput) {
      X_copy_nonConst = rcp (new MV (X));
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

//==========================================================================
template<class MatrixType>
double Chebyshev<MatrixType>::getLambdaMaxForApply() const {
  return impl_.getLambdaMaxForApply();
}


}//namespace Ifpack2

#endif // IFPACK2_CHEBYSHEV_DEF_HPP

