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
: A_ (A),
  PolyDegree_ (1),
  EigRatio_ (30.0),
  LambdaMin_ (0.0),
  LambdaMax_ (100.0),
  MinDiagonalValue_ (0.0),
  ZeroStartingSolution_ (true),
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
{
  TEUCHOS_TEST_FOR_EXCEPTION(A_.is_null (), std::invalid_argument,
    "Ifpack2::Chebyshev: Input matrix to constructor is null.");

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_->getGlobalNumRows() != A_->getGlobalNumCols(), 
    std::invalid_argument,
    "Ifpack2::Chebyshev: The input matrix A must be square.  "
    "A has " << A_->getGlobalNumRows() << " rows and " 
    << A_->getGlobalNumCols() << " columns.");

#ifdef HAVE_TEUCHOS_DEBUG
  using Teuchos::RCP;

  RCP<const map_type> domainMap = A_->getDomainMap ();
  RCP<const map_type> rangeMap = A_->getRangeMap ();

  // The relation 'isSameAs' is transitive.  It's also a collective,
  // so we don't have to do a "shared" test for exception (i.e., a
  // global reduction on the test value).
  TEUCHOS_TEST_FOR_EXCEPTION(
     ! domainMap->isSameAs (*rangeMap),
     std::runtime_error,
     "Ifpack2::Chebyshev: The domain Map and range Map of the matrix must be "
     "the same (in the sense of isSameAs()).  We will only check for this if "
     "Trilinos was built with the CMake configuration option Teuchos_ENABLE_"
     "DEBUG set to ON.");
#endif // HAVE_TEUCHOS_DEBUG
}

//==========================================================================
template<class MatrixType>
Chebyshev<MatrixType>::~Chebyshev() {
}

//==========================================================================
template<class MatrixType>
void 
Chebyshev<MatrixType>::setParameters(const Teuchos::ParameterList& List) 
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::Export<local_ordinal_type,global_ordinal_type,node_type> export_type;

  Ifpack2::getParameter(List, "chebyshev: ratio eigenvalue", EigRatio_);
  Ifpack2::getParameter(List, "chebyshev: min eigenvalue", LambdaMin_);
  Ifpack2::getParameter(List, "chebyshev: max eigenvalue", LambdaMax_);
  Ifpack2::getParameter(List, "chebyshev: degree", PolyDegree_);
  TEUCHOS_TEST_FOR_EXCEPTION(
    PolyDegree_ <= 0, 
    std::invalid_argument,
    "Ifpack2::Chebyshev::setParameters(): The \"chebyshev: degree\" parameter "
    "must be an integer >= 1, but you specified a value of " << PolyDegree_ 
    << ".");
  Ifpack2::getParameter(List, "chebyshev: min diagonal value", MinDiagonalValue_);
  Ifpack2::getParameter(List, "chebyshev: zero starting solution", ZeroStartingSolution_);

  vector_type* globalDiag = NULL;
  Ifpack2::getParameter(List, "chebyshev: operator inv diagonal", globalDiag);

  if (globalDiag != NULL) {
    RCP<const map_type> sourceMap = globalDiag->getMap ();
    RCP<const map_type> rangeMap = A_->getRangeMap ();
    RCP<const map_type> rowMap = A_->getRowMap ();
    if (rangeMap.getRawPtr () == sourceMap.getRawPtr ()) {
      // The given vector's Map is identical to the matrix's range
      // Map.  That means we don't need to Export.  (We don't call
      // isSameAs() here, because that takes at least one global
      // reduction.  This is the fast path.  We may call isSameAs() in
      // the slow path below, to avoid an even slower path.
      userSuppliedInvDiag_ = rcp (new vector_type (*globalDiag));
    }
    else { // We need to Export.
      RCP<const export_type> exporter;
      // Making an Export object from scratch is expensive enough that
      // it's worth the O(1) global reductions to call isSameAs().
      if (sourceMap.getRawPtr () == rowMap.getRawPtr () || 
	  sourceMap->isSameAs (*rowMap)) {
	// We can reuse the matrix's Export object, if there is one.
	exporter = A_->getGraph ()->getExporter ();
      }
      else { // We have to make a new Export object.
	exporter = rcp (new export_type (sourceMap, rangeMap));
      }
      
      userSuppliedInvDiag_ = rcp (new vector_type (rangeMap));
      if (exporter.is_null ()) { 
	// Row Map and range Map are the same; no need to Export.
	*userSuppliedInvDiag_ = *globalDiag;
      }
      else {
	userSuppliedInvDiag_->doExport (*globalDiag, *exporter, Tpetra::ADD);
      }
    } // if we don't need to Export, or if we do
  } // the user gave us a vector of diagonal entries
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const Teuchos::Comm<int> > & 
Chebyshev<MatrixType>::getComm() const{
  return A_->getRowMap ()->getComm ();
}

//==========================================================================
template<class MatrixType>
Teuchos::RCP<const typename Chebyshev<MatrixType>::row_matrix_type>
Chebyshev<MatrixType>::
getMatrix() const {
  return A_;
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const typename Chebyshev<MatrixType>::map_type>&
Chebyshev<MatrixType>::
getDomainMap() const {
  return A_->getDomainMap();
}

//==========================================================================
template<class MatrixType>
const Teuchos::RCP<const typename Chebyshev<MatrixType>::map_type>&
Chebyshev<MatrixType>::
getRangeMap() const {
  return A_->getRangeMap();
}

//==========================================================================
template<class MatrixType>
bool Chebyshev<MatrixType>::hasTransposeApply() const {
  return true;
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

//==========================================================================
template<class MatrixType>
void 
Chebyshev<MatrixType>::
apply (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
       Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const 
{
  using Teuchos::ArrayRCP;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

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

  TEUCHOS_TEST_FOR_EXCEPTION(alpha != STS::one(), std::logic_error,
    "Ifpack2::Chebyshev::apply(): Not yet implemented for alpha != 1.");

  TEUCHOS_TEST_FOR_EXCEPTION(beta != STS::zero(), std::logic_error,
    "Ifpack2::Chebyshev::apply(): Not yet implemented for beta != 0.");

#ifdef HAVE_TEUCHOS_DEBUG
  {
    RCP<const map_type> domainMap = A_->getDomainMap ();
    RCP<const map_type> rangeMap = A_->getRangeMap ();

    // The relation 'isSameAs' is transitive.  It's also a collective,
    // so we don't have to do a "shared" test for exception (i.e., a
    // global reduction on the test value).
    TEUCHOS_TEST_FOR_EXCEPTION(
       ! X.getMap ()->isSameAs (*domainMap),
       std::runtime_error,
       "Ifpack2::Chebyshev: The domain Map of the matrix must be the same as the Map of the input vector(s) X.");
    TEUCHOS_TEST_FOR_EXCEPTION(
       ! Y.getMap ()->isSameAs (*rangeMap),
       std::runtime_error,
       "Ifpack2::Chebyshev: The range Map of the matrix must be the same as the Map of the output vector(s) Y.");
  }
#endif // HAVE_TEUCHOS_DEBUG

  Time_->start();

  // If X and Y are pointing to the same memory location,
  // we need to create an auxiliary vector, Xcopy
  RCP<const MV> Xcopy;
  if (X.getLocalMV().getValues() == Y.getLocalMV().getValues()) {
    Xcopy = rcp (new MV (X));
  }
  else {
    Xcopy = rcpFromRef (X);
  }

  ArrayRCP<ArrayRCP<const scalar_type> > xView = Xcopy->get2dView();
  ArrayRCP<ArrayRCP<scalar_type> > yView = Y.get2dViewNonConst();
  ArrayRCP<const scalar_type> invdiag = InvDiagonal_->get1dView();

  size_t nVecs = Y.getNumVectors();

  //--- Do a quick solve when the matrix is identity
  if ((LambdaMin_ == 1.0) && (LambdaMax_ == LambdaMin_)) {
    // Y_j = X_j .* invdiag, for all columns j.
    Y.elementWiseMultiply (STS::one(), *InvDiagonal_, X, STS::zero());
    return;
  }

  //--- initialize coefficients

  const scalar_type one = STS::one();
  const scalar_type two = one + one;
  const scalar_type one_half = one / two;
  const scalar_type lambdaMaxIncr = as<scalar_type> (1.1);

  // Note that delta stores the inverse of ML_Cheby::delta
  scalar_type alpha_cheby = LambdaMax_ / EigRatio_;
  // ML source comment: "try and bracket high".  Presumably this
  // explains that the coefficient is increased by a small factor
  // (1.1) in order for the smoother to be more effective.
  scalar_type beta_cheby = lambdaMaxIncr * LambdaMax_;
  scalar_type delta = two / (beta_cheby - alpha_cheby);
  scalar_type theta = one_half * (beta_cheby + alpha_cheby);
  scalar_type s1 = theta * delta;

  //--- Define vectors
  // In ML_Cheby, V corresponds to pAux and W to dk
  MV V (X);
  MV W (X);

  ArrayRCP<ArrayRCP<const scalar_type> > vView = V.get2dView();
  ArrayRCP<ArrayRCP<scalar_type> > wView = W.get2dViewNonConst();

  scalar_type oneOverTheta = one/theta;
  const size_t numMyRows = A_->getNodeNumRows ();

  // Do the smoothing when block scaling is turned OFF
  // --- Treat the initial guess
  if (ZeroStartingSolution_ == false) {
    // Compute W = (1/theta) D^{-1} (X - A*Y).

    applyMat (Y, V, mode); // V = A*Y
    if (nVecs == 1) {
      ArrayRCP<const scalar_type> x = xView[0];
      ArrayRCP<scalar_type> w = wView[0];
      ArrayRCP<const scalar_type> v = vView[0];
      for (size_t i = 0; i < numMyRows; ++i) {
        w[i] = invdiag[i] * (x[i] - v[i]) * oneOverTheta;
      }
    }
    else {
      for (size_t k = 0; k < nVecs; ++k) {
        ArrayRCP<scalar_type> wk = wView[k];
        ArrayRCP<const scalar_type> vk = vView[k];
        for (size_t i = 0; i < numMyRows; ++i) {
          scalar_type coeff = invdiag[i]*oneOverTheta;
          wk[i] = (xView[k][i] - (vk[i])) * coeff;
        }
      }
    } // if (nVec == 1)

    Y.update(one, W, one); // Y = Y + W
  }
  else {
    // Compute W = (1/theta) D^{-1} X and set Y = W.

    if (nVecs == 1) {
      ArrayRCP<const scalar_type> x= xView[0];
      ArrayRCP<scalar_type> w = wView[0];
      ArrayRCP<scalar_type> y = yView[0];
      for (size_t i = 0; i < numMyRows; ++i) {
        w[i] = invdiag[i] * x[i] * oneOverTheta;
        y[i] = w[i];
      }
    }
    else {
      for (size_t k = 0; k < nVecs; ++k) {
        for (size_t i = 0; i < numMyRows; ++i) {
          scalar_type coeff = invdiag[i]*oneOverTheta;
          wView[k][i] = xView[k][i] * coeff;
          yView[k][i] = wView[k][i];
        }
      }
    } // if (nVec == 1)
  } // if (ZeroStartingSolution_ == false)

  //--- Apply the polynomial
  scalar_type rhok = one/s1, rhokp1;
  scalar_type dtemp1, dtemp2;
  int degreeMinusOne = PolyDegree_ - 1;
  for (int deg = 0; deg < degreeMinusOne; ++deg) {
    applyMat (Y, V, mode); // V = A*Y
    rhokp1 = one / (two *s1 - rhok);
    dtemp1 = rhokp1 * rhok;
    dtemp2 = two * rhokp1 * delta;
    rhok = rhokp1;

    // Compute W = dtemp1*W + dtemp2 * D^{-1} * (X - V).

    W.scale (dtemp1); // W = dtemp1 * W
    for (size_t k = 0; k < nVecs; ++k) { // W = W + dtemp2 * D^{-1} * (X - V)
      for (size_t i = 0; i < numMyRows; ++i) {
        scalar_type coeff = invdiag[i]*dtemp2;
        wView[k][i] += (xView[k][i] - (vView[k][i])) * coeff;
      }
    }
    Y.update (one, W, one); // Y = Y + W
  } // for (deg = 0; deg < degreeMinusOne; ++deg)

  ++NumApply_;
  Time_->stop();
  ApplyTime_ += Time_->totalElapsedTime();
}

//==========================================================================
template<class MatrixType>
void 
Chebyshev<MatrixType>::
applyMat (const Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& X,
	  Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Y,
	  Teuchos::ETransp mode) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(! isComputed(), std::runtime_error,
   "Ifpack2::Chebyshev::applyMat(): isComputed() must be true prior to calling applyMat().");
  TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != Y.getNumVectors(), std::runtime_error,
   "Ifpack2::Chebyshev::applyMat(): X.getNumVectors() != Y.getNumVectors().");
  A_->apply (X, Y, mode);
}

//==========================================================================
template<class MatrixType>
void Chebyshev<MatrixType>::initialize() {
  IsInitialized_ = false;
  Time_->start (true);

  // mfh 15 Jan 2013: Defer fetching the global number of nonzeros in
  // the matrix until compute().  That way, it will always be correct
  // to omit calling initialize(), even if the number of entries in
  // the sparse matrix has changed since the last call to compute().

  ++NumInitialize_;
  Time_->stop ();
  InitializeTime_ += Time_->totalElapsedTime();
  IsInitialized_ = true;
}

//==========================================================================
template<class MatrixType>
void Chebyshev<MatrixType>::compute()
{
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::Export<local_ordinal_type,global_ordinal_type,node_type> export_type;

  Time_->start (true);

  if (! isInitialized ()) {
    initialize ();
  }

  // Reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  // mfh 15 Jan 2013: Defer fetching the global number of nonzeros
  // until here.  That way, it will always be correct to omit calling
  // initialize(), even if the number of entries in the sparse matrix
  // has changed since the last call to compute().
  const Tpetra::global_size_t numGlobalNonzeros = A_->getGlobalNumEntries ();
  (void) numGlobalNonzeros; // Forestall compiler warning.

  TEUCHOS_TEST_FOR_EXCEPTION(PolyDegree_ <= 0, std::runtime_error,
    "Ifpack2::Chebyshev::compute(): PolyDegree_ must be at least one");

  // If the user has not given us the inverse diagonal, compute it
  // here.  We do so on every call to compute(), in case the values in
  // the matrix have changed since the last call.
  if (! userSuppliedInvDiag_.is_null ()) {
    InvDiagonal_ = userSuppliedInvDiag_;
  }
  else {
    RCP<const map_type> rowMap = A_->getRowMap ();
    RCP<vector_type> globalDiagRowMap = rcp (new vector_type (rowMap));
    A_->getLocalDiagCopy (*globalDiagRowMap);

    // globalDiag should be distributed according to the row Map of A
    // at this point.  If necessary, Export the vector of diagonal
    // entries from the row Map to the range Map.
    RCP<vector_type> globalDiagRangeMap;
    RCP<const export_type> exporter = A_->getGraph ()->getExporter ();
    if (exporter.is_null ()) {
      globalDiagRangeMap = globalDiagRowMap; // Row & range Maps are the same.
    }
    else {
      RCP<const map_type> rangeMap = A_->getRangeMap ();
      globalDiagRangeMap = rcp (new vector_type (rangeMap));
      globalDiagRangeMap->doExport (*globalDiagRowMap, *exporter, Tpetra::ADD);
    }
    // Get rid of the row Map version of the diagonal.
    globalDiagRowMap = Teuchos::null;

    // Invert all diagonal elements no smaller in magnitude than the
    // minimum allowed diagonal value d_min.  Those smaller than d_min
    // in magnitude are replaced with d_min.  (That's why d_min itself
    // is a scalar and not a magnitude.)  Tpetra::MultiVector's
    // reciprocal() method doesn't have the threshold feature that we
    // want, but Kokkos (mfh 15 Jan 2013: as of today) does.
    typedef Kokkos::MultiVector<scalar_type, node_type> KMV;
    KMV& localDiag = globalDiagRangeMap->getLocalMVNonConst ();
    typedef Kokkos::DefaultArithmetic<KMV> KMVT;
    KMVT::ReciprocalThreshold (localDiag, MinDiagonalValue_);

    InvDiagonal_ = globalDiagRangeMap; // "Commit" the result.
    ComputeFlops_ += A_->getNodeNumRows ();
  }

  ++NumCompute_;
  Time_->stop();
  ComputeTime_ += Time_->totalElapsedTime();
  IsComputed_ = true;
}

//==========================================================================
template<class MatrixType>
void 
Chebyshev<MatrixType>::
PowerMethod (const Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& Operator, 
	     const Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type>& InvPointDiagonal, 
	     const int MaximumIterations, 
	     scalar_type& lambda_max)
{
  using Teuchos::Array;

  const scalar_type one = STS::one();
  const scalar_type zero = STS::zero();

  lambda_max = zero;
  Array<scalar_type> RQ_top(1), RQ_bottom(1);
  vector_type x (Operator.getDomainMap ());
  vector_type y (Operator.getRangeMap ());
  x.randomize ();
  Array<magnitude_type> norms (x.getNumVectors ());
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
  //
  oss << ", global rows = " << A_->getGlobalNumRows()
      << ", global cols = " << A_->getGlobalNumCols()
      << "}";
  return oss.str();
}

//==========================================================================
template <class MatrixType>
void Chebyshev<MatrixType>::describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const {
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
}

}//namespace Ifpack2

#endif // IFPACK2_CHEBYSHEV_DEF_HPP

