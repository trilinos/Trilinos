// ***********************************************************************
// 
//      Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2004) Sandia Corporation
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

#ifndef IFPACK2_DETAILS_CHEBYSHEV_DEF_HPP
#define IFPACK2_DETAILS_CHEBYSHEV_DEF_HPP

/// \file Ifpack2_Details_Chebyshev_def.hpp
/// \brief Definition of Chebyshev implementation
/// \author Mark Hoemmen
///
/// This file is meant for Ifpack2 developers only, not for users.
/// It defines a new implementation of Chebyshev iteration.

#include "Ifpack2_Details_Chebyshev_decl.hpp"
#include <cmath>

namespace Ifpack2 {
namespace Details {

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
checkConstructorInput () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(A_.is_null (), std::invalid_argument,
    "Ifpack2::Details::Chebyshev: Input matrix to constructor is null.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_->getGlobalNumRows() != A_->getGlobalNumCols(), 
    std::invalid_argument,
    "Ifpack2::Details::Chebyshev: The input matrix A must be square.  "
    "A has " << A_->getGlobalNumRows() << " rows and " 
    << A_->getGlobalNumCols() << " columns.");
#ifdef HAVE_TEUCHOS_DEBUG
  Teuchos::RCP<const map_type> domainMap = A_->getDomainMap ();
  Teuchos::RCP<const map_type> rangeMap = A_->getRangeMap ();

  // The relation 'isSameAs' is transitive.  It's also a collective,
  // so we don't have to do a "shared" test for exception (i.e., a
  // global reduction on the test value).
  TEUCHOS_TEST_FOR_EXCEPTION(
     ! domainMap->isSameAs (*rangeMap),
     std::invalid_argument,
     "Ifpack2::Details::Chebyshev: The domain Map and range Map of the matrix "
     "must be the same (in the sense of isSameAs()).  We only check for this "
     "if Trilinos was built with the CMake configuration option Teuchos_ENABLE_"
     "DEBUG set to ON.");
#endif // HAVE_TEUCHOS_DEBUG
}

template<class ScalarType, class MV, class MAT>
Chebyshev<ScalarType, MV, MAT>::
Chebyshev (Teuchos::RCP<const MAT> A) : 
  A_ (A),
  computedLambdaMax_ (STS::nan ()),
  computedLambdaMin_ (STS::nan ()),
  lambdaMax_ (STS::nan ()),
  lambdaMin_ (STS::nan ()),
  eigRatio_ (Teuchos::as<ST> (30)),
  minDiagVal_ (STS::eps ()),
  numIters_ (1),
  eigMaxIters_ (10),
  zeroStartingSolution_ (true),
  textbookAlgorithm_ (false),
  computeMaxResNorm_ (false)
{
  checkConstructorInput ();
}

template<class ScalarType, class MV, class MAT>
Chebyshev<ScalarType, MV, MAT>::
Chebyshev (Teuchos::RCP<const MAT> A, Teuchos::ParameterList& params) : 
  A_ (A), 
  computedLambdaMax_ (STS::nan ()),
  computedLambdaMin_ (STS::nan ()),
  lambdaMax_ (STS::nan ()),
  lambdaMin_ (STS::nan ()),
  eigRatio_ (Teuchos::as<ST> (30)),
  minDiagVal_ (STS::eps ()),
  numIters_ (1),
  eigMaxIters_ (10),
  zeroStartingSolution_ (true),
  textbookAlgorithm_ (false),
  computeMaxResNorm_ (false)
{
  checkConstructorInput ();
  setParameters (params);
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
setParameters (Teuchos::ParameterList& plist) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcpFromRef;

  // Note to developers: The logic for this method is complicated,
  // because we want to accept Ifpack and ML parameters whenever
  // possible, but we don't want to add their default values to the
  // user's ParameterList.  That's why we do all the isParameter()
  // checks, instead of using the two-argument version of get()
  // everywhere.  The min and max eigenvalue parameters are also a
  // special case, because we decide whether or not to do eigenvalue
  // analysis based on whether the user supplied the max eigenvalue.

  // Default values of all the parameters.
  const ST defaultLambdaMax = STS::nan ();
  const ST defaultLambdaMin = STS::nan ();
  // 30 is Ifpack::Chebyshev's default.  ML has a different default
  // eigRatio for smoothers and the coarse-grid solve (if using
  // Chebyshev for that).  The former uses 20; the latter uses 30.
  // We're testing smoothers here, so use 20.  (However, if you give
  // ML an Epetra matrix, it will use Ifpack for Chebyshev, in which
  // case it would defer to Ifpack's default settings.)
  const ST defaultEigRatio = Teuchos::as<ST> (30);
  const ST defaultMinDiagVal = STS::eps ();
  const int defaultNumIters = 1;
  const int defaultEigMaxIters = 10;
  const bool defaultZeroStartingSolution = true; // Ifpack::Chebyshev default
  const bool defaultTextbookAlgorithm = false;
  const bool defaultComputeMaxResNorm = false;

  // We'll set the instance data transactionally, after all reads
  // from the ParameterList.  That way, if any of the ParameterList
  // reads fail (e.g., due to the wrong parameter type), we will not
  // have left the instance data in a half-changed state.
  RCP<const V> userInvDiag;
  ST lambdaMax = defaultLambdaMax;
  ST lambdaMin = defaultLambdaMin;
  ST eigRatio = defaultEigRatio;
  ST minDiagVal = defaultMinDiagVal;
  int numIters = defaultNumIters;
  int eigMaxIters = defaultEigMaxIters;
  bool zeroStartingSolution = defaultZeroStartingSolution;
  bool textbookAlgorithm = defaultTextbookAlgorithm;
  bool computeMaxResNorm = defaultComputeMaxResNorm;

  // Fetch the parameters from the ParameterList.  Defer all
  // externally visible side effects until we have finished all
  // ParameterList interaction.  This makes the method satisfy the
  // strong exception guarantee.

  // Get the user-supplied inverse diagonal.
  //
  // Check for a raw pointer (const V* or V*), for Ifpack
  // compatibility, as well as for an RCP<const V> or RCP<V>.  We'll
  // copy the vector anyway, so it doesn't matter whether it's const
  // or nonconst.
  if (plist.isParameter ("chebyshev: operator inv diagonal")) {
    try { // Could the type be const V*?
      const V* rawUserInvDiag = plist.get<const V*> ("chebyshev: operator inv diagonal");
      // It's OK to have a nonowning reference; we'll copy the vector anyway.
      userInvDiag = rcp (rawUserInvDiag, false);
    } catch (Teuchos::Exceptions::InvalidParameterType&) {
    }
    if (userInvDiag.is_null ()) {
      try { // Could the type be V*?
	V* rawUserInvDiag = plist.get<V*> ("chebyshev: operator inv diagonal");
	// It's OK to have a nonowning reference; we'll copy the vector anyway.
	userInvDiag = rcp (const_cast<const V*> (rawUserInvDiag), false);
      } catch (Teuchos::Exceptions::InvalidParameterType&) {
      }
    }
    if (userInvDiag.is_null ()) {
      try { // Could the type be RCP<const V>?
	userInvDiag = plist.get<RCP<const V> > ("chebyshev: operator inv diagonal");
      } catch (Teuchos::Exceptions::InvalidParameterType&) {
      }
    }
    if (userInvDiag.is_null ()) {
      try { // Could the type be RCP<const V>?
	RCP<V> userInvDiagNonconst = plist.get<RCP<V> > ("chebyshev: operator inv diagonal");
	userInvDiag = rcp_const_cast<const V> (userInvDiagNonconst);
      } catch (Teuchos::Exceptions::InvalidParameterType&) {
      }
    }
    // We've tried all the possible types.  If we got something, then
    // make a range Map copy of it.
    if (! userInvDiag.is_null ()) {
      userInvDiag = makeRangeMapVector (*userInvDiag);
    }
  }

  // Don't fill in defaults for the max or min eigenvalue, because
  // this class uses the existence of those parameters to determine
  // whether it should do eigenanalysis.
  if (plist.isParameter ("chebyshev: max eigenvalue")) {
    lambdaMax = plist.get<ST> ("chebyshev: max eigenvalue");
  }
  if (plist.isParameter ("chebyshev: min eigenvalue")) {
    lambdaMin = plist.get<ST> ("chebyshev: min eigenvalue");
  }

  // Only fill in Ifpack2's name for the default parameter, not ML's.
  if (plist.isParameter ("smoother: Chebyshev alpha")) { // ML compatibility
    eigRatio = plist.get<ST> ("smoother: Chebyshev alpha");
  }
  // Ifpack2's name overrides ML's name.
  eigRatio = plist.get ("chebyshev: ratio eigenvalue", eigRatio);

  // Same name in Ifpack2 and Ifpack.
  minDiagVal = plist.get ("chebyshev: min diagonal value", minDiagVal);

  // Only fill in Ifpack2's name, not ML's or Ifpack's.
  if (plist.isParameter ("smoother: sweeps")) { // ML compatibility
    numIters = plist.get<int> ("smoother: sweeps");
  } // Ifpack's name overrides ML's name.
  if (plist.isParameter ("relaxation: sweeps")) { // Ifpack compatibility
    numIters = plist.get<int> ("relaxation: sweeps");
  } // Ifpack2's name overrides Ifpack's name.
  numIters = plist.get ("chebyshev: degree", numIters);

  // The last parameter name overrides the first.
  if (plist.isParameter ("eigen-analysis: iterations")) { // ML compatibility
    eigMaxIters = plist.get<int> ("eigen-analysis: iterations");
  } // Ifpack2's name overrides ML's name.
  eigMaxIters = plist.get ("chebyshev: eigenvalue max iterations", eigMaxIters);

  zeroStartingSolution = plist.get ("chebyshev: zero starting solution", zeroStartingSolution);

  // We don't want to fill these parameters in, because they shouldn't
  // be visible to Ifpack2::Chebyshev users.
  if (plist.isParameter ("chebyshev: textbook algorithm")) {
    textbookAlgorithm = plist.get<bool> ("chebyshev: textbook algorithm");
  }
  if (plist.isParameter ("chebyshev: compute max residual norm")) {
    computeMaxResNorm = plist.get<bool> ("chebyshev: compute max residual norm");
  }

  // Test for Ifpack parameters that we won't ever implement here.
  // Be careful to use the one-argument version of get(), since the
  // two-argment version adds the parameter if it's not there.
  TEUCHOS_TEST_FOR_EXCEPTION(
    plist.isParameter ("chebyshev: use block mode") && 
    ! plist.get<bool> ("chebyshev: use block mode"),
    std::invalid_argument,
    "Ifpack2::Chebyshev requires that if you supply the Ifpack parameter "
    "\"chebyshev: use block mode\", it must be set to false.  Ifpack2's "
    "Chebyshev does not implement Ifpack's block mode.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    plist.isParameter ("chebyshev: solve normal equations") && 
    ! plist.get<bool> ("chebyshev: solve normal equations"),
    std::invalid_argument,
    "Ifpack2::Chebyshev does not and will never implement the Ifpack "
    "parameter \"chebyshev: solve normal equations\".  If you want to solve "
    "the normal equations, construct a Tpetra::Operator that implements "
    "A^* A, and use Chebyshev to solve A^* A x = A^* b.");

  // Test for Ifpack parameters that we haven't implemented yet.
  //
  // For now, we only check that this ML parameter, if provided, has
  // the one value that we support.  We consider other values "invalid
  // arguments" rather than "logic errors," because Ifpack does not
  // implement eigenanalyses other than the power method.
  std::string eigenAnalysisType ("power-method");
  if (plist.isParameter ("eigen-analysis: type")) {
    eigenAnalysisType = plist.get<std::string> ("eigen-analysis: type");
    TEUCHOS_TEST_FOR_EXCEPTION(
      eigenAnalysisType != "power-method" && 
      eigenAnalysisType != "power method",
      std::invalid_argument,
      "Ifpack2::Chebyshev: This class supports the ML parameter \"eigen-"
      "analysis: type\" for backwards compatibility.  However, Ifpack2 "
      "currently only supports the \"power-method\" option for this "
      "parameter.  This imitates Ifpack, which only implements the power "
      "method for eigenanalysis.");
  }

  userInvDiag_ = userInvDiag;
  lambdaMax_ = lambdaMax;
  lambdaMin_ = lambdaMin;
  eigRatio_ = eigRatio;
  minDiagVal_ = minDiagVal;
  numIters_ = numIters;
  eigMaxIters_ = eigMaxIters;
  zeroStartingSolution_ = zeroStartingSolution;
  textbookAlgorithm_ = textbookAlgorithm;
  computeMaxResNorm_ = computeMaxResNorm;
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
compute () {
  if (userInvDiag_.is_null ()) {
    // The matrix may have changed, and the parameters also affect
    // this process, so recompute the inverse diagonal.
    D_ = makeInverseDiagonal (*A_);
  } else {
    D_ = userInvDiag_;
  }

  // The matrix may have changed, so we always compute the
  // eigenvalue estimates, even if we computed them before.
  // However, if the user gave us lambdaMax already (if it's not
  // NaN), then don't (re)compute it.
  if (STS::isnaninf (lambdaMax_)) {
    const ST computedLambdaMax = powerMethod (*A_, *D_, eigMaxIters_);
    const ST computedLambdaMin = computedLambdaMax / eigRatio_;

    // Defer "committing" results until all computations succeeded.
    computedLambdaMax_ = computedLambdaMax;
    computedLambdaMin_ = computedLambdaMin;
  }
}

template<class ScalarType, class MV, class MAT>
typename Chebyshev<ScalarType, MV, MAT>::MT
Chebyshev<ScalarType, MV, MAT>::
apply (const MV& B, MV& X) {
  ST eigRatio = eigRatio_;
  ST lambdaMax = lambdaMax_;
  ST lambdaMin = lambdaMax / eigRatio_; // lambdaMin_; (???)

  // FIXME (mfh 22 Jan 2013) We really only want to check if it's
  // NaN.  Inf means something different.  However,
  // Teuchos::ScalarTraits doesn't distinguish the two cases.
  if (STS::isnaninf (lambdaMax)) {
    lambdaMax = computedLambdaMax_;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isnaninf (lambdaMax), 
    std::runtime_error, 
    "Chebyshev::apply: Both lambdaMax_ and computedLambdaMax_ are NaN or Inf.  "
    "This means one of the following: (a) you didn't call compute() before "
    "calling apply(), "
    "(b) you didn't call compute() after calling setParameters(), or "
    "(c) there is a bug in this class (compute() has not done the "
    "eigenanalysis that it should have done).");
  if (STS::isnaninf (lambdaMin)) {
    lambdaMin = computedLambdaMin_;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isnaninf (lambdaMin), 
    std::runtime_error, 
    "Chebyshev::apply: Both lambdaMin_ and computedLambdaMin_ are NaN or Inf.  "
    "This means one of the following: (a) you didn't call compute() before "
    "calling apply(), "
    "(b) you didn't call compute() after calling setParameters(), or "
    "(c) there is a bug in this class (compute() has not done the "
    "eigenanalysis that it should have done).");
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isnaninf (eigRatio), 
    std::logic_error, 
    "Chebyshev::apply: eigRatio is NaN or Inf.");

  if (! textbookAlgorithm_) {
    const ST one = Teuchos::as<ST> (1);
    const ST diff = lambdaMax - one;
    const ST tol = Teuchos::as<ST> (1.0e-6);
    if (STS::magnitude(diff) < STS::magnitude(tol)) {
      lambdaMin = one;
      lambdaMax = lambdaMin;
      eigRatio = one; // Ifpack doesn't include this line.
    }
    ifpackApplyImpl (*A_, B, X, numIters_, lambdaMax, lambdaMin, eigRatio, *D_);
  } else {
    textbookApplyImpl (*A_, B, X, numIters_, lambdaMax, lambdaMin, eigRatio, *D_);
  }

  if (computeMaxResNorm_) {
    MV R (B.getMap (), B.getNumVectors ());
    computeResidual (R, B, *A_, X);
    Teuchos::Array<MT> norms (B.getNumVectors ());
    R.norm2 (norms ());
    return *std::max_element (norms.begin (), norms.end ());
  } else {
    return Teuchos::as<MT> (0);
  }
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
print (std::ostream& out) {
  using std::endl;
  out << "Chebyshev:" << endl
      << "- computedLambdaMax_ = " << computedLambdaMax_ << endl
      << "- computedLambdaMin_ = " << computedLambdaMin_ << endl
      << "- lambdaMax_ = " << lambdaMax_ << endl
      << "- lambdaMin_ = " << lambdaMin_ << endl
      << "- eigRatio_ = " << eigRatio_ << endl
      << "- numIters_ = " << numIters_ << endl
      << "- eigMaxIters_ = " << eigMaxIters_ << endl
      << "- textbookAlgorithm_ = " << textbookAlgorithm_ << endl;
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
computeResidual (MV& R, const MV& B, const MAT& A, const MV& X, 
		 const Teuchos::ETransp mode) 
{
  R = B;
  A.apply (X, R, mode, -STS::one(), STS::one());
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
solve (MV& Z, const V& D_inv, const MV& R) {
  Z.elementWiseMultiply (STS::one(), D_inv, R, STS::zero());
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
solve (MV& Z, const ST alpha, const V& D_inv, const MV& R) {
  Z.elementWiseMultiply (alpha, D_inv, R, STS::zero());
}

template<class ScalarType, class MV, class MAT>
Teuchos::RCP<typename Chebyshev<ScalarType, MV, MAT>::V>
Chebyshev<ScalarType, MV, MAT>::
makeInverseDiagonal (const MAT& A) const {
  using Teuchos::RCP;

  RCP<V> D_rowMap (new V (A.getGraph ()->getRowMap ()));
  A.getLocalDiagCopy (*D_rowMap);
  RCP<V> D_rangeMap = makeRangeMapVector (*D_rowMap);

  // Invert the diagonal entries, replacing entries less than
  // machine precision with machine precision.
  typedef Kokkos::MultiVector<ST, typename MAT::node_type> KMV;
  KMV& localDiag = D_rangeMap->getLocalMVNonConst ();
  typedef Kokkos::DefaultArithmetic<KMV> KMVT;
  KMVT::ReciprocalThreshold (localDiag, minDiagVal_);
  return D_rangeMap;
}

template<class ScalarType, class MV, class MAT>
Teuchos::RCP<typename Chebyshev<ScalarType, MV, MAT>::V>
Chebyshev<ScalarType, MV, MAT>::
makeRangeMapVector (const V& D) const {
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::Export<typename MV::local_ordinal_type,
			 typename MV::global_ordinal_type,
			 typename MV::node_type> export_type;
  RCP<const map_type> sourceMap = D.getMap ();
  RCP<const map_type> rangeMap = A_->getRangeMap ();
  RCP<const map_type> rowMap = A_->getRowMap ();
  RCP<V> D_out;

  if (rangeMap.getRawPtr () == sourceMap.getRawPtr ()) {
    // The given vector's Map is identical to the matrix's range Map.
    // That means we don't need to Export.  (We don't call isSameAs()
    // here, because that takes at least one global reduction.  This
    // is the fast path.  We may call isSameAs() in the slow path
    // below, to avoid an even slower path.
    D_out = rcp (new V (D));
  }
  else { // We need to Export.
    RCP<const export_type> exporter;
    // Making an Export object from scratch is expensive enough that
    // it's worth the O(1) global reductions to call isSameAs(), to
    // see if we can avoid that cost.
    if (sourceMap.getRawPtr () == rowMap.getRawPtr () || 
	sourceMap->isSameAs (*rowMap)) {
      // We can reuse the matrix's Export object, if there is one.
      exporter = A_->getGraph ()->getExporter ();
    }
    else { // We have to make a new Export object.
      exporter = rcp (new export_type (sourceMap, rangeMap));
    }
     
    D_out = rcp (new V (D));
    if (exporter.is_null ()) { 
      // Row Map and range Map are the same; no need to Export.
      *D_out = D;
    }
    else {
      D_out->doExport (D, *exporter, Tpetra::ADD);
    }
  } // if we don't need to Export, or if we do

  return D_out;
}


template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
textbookApplyImpl (const MAT& A,
		   const MV& B,
		   MV& X,
		   const int numIters,
		   const ST lambdaMax,
		   const ST lambdaMin,
		   const ST eigRatio,
		   const V& D_inv) const
{
  (void) lambdaMin; // Forestall compiler warning.
  const ST myLambdaMin = lambdaMax / eigRatio;

  const ST zero = Teuchos::as<ST> (0);
  const ST one = Teuchos::as<ST> (1);
  const ST two = Teuchos::as<ST> (2);
  const ST d = (lambdaMax + myLambdaMin) / two; // Ifpack2 calls this theta
  const ST c = (lambdaMax - myLambdaMin) / two; // Ifpack2 calls this 1/delta

  if (zeroStartingSolution_ && numIters > 0) {
    // If zero iterations, then input X is output X.
    X.putScalar (zero);
  }
  MV R (B.getMap (), B.getNumVectors (), false);
  MV P (B.getMap (), B.getNumVectors (), false);
  MV Z (B.getMap (), B.getNumVectors (), false);
  ST alpha, beta;
  for (int i = 0; i < numIters; ++i) {
    computeResidual (R, B, A, X); // R = B - A*X
    solve (Z, D_inv, R); // z = D_inv * R, that is, D \ R.
    if (i == 0) {
      P = Z;
      alpha = two / d;
    } else {
      //beta = (c * alpha / two)^2;
      //const ST sqrtBeta = c * alpha / two;
      //beta = sqrtBeta * sqrtBeta;
      beta = alpha * (c/two) * (c/two);
      alpha = one / (d - beta);
      P.update (one, Z, beta); // P = Z + beta*P
    }
    X.update (alpha, P, one); // X = X + alpha*P
    // If we compute the residual here, we could either do R = B -
    // A*X, or R = R - alpha*A*P.  Since we choose the former, we
    // can move the computeResidual call to the top of the loop.
  }
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
ifpackApplyImpl (const MAT& A,
		 const MV& B,
		 MV& X,
		 const int numIters,
		 const ST lambdaMax,
		 const ST lambdaMin,
		 const ST eigRatio,
		 const V& D_inv) 
{
  if (numIters <= 0) {
    return;
  }
  const ST one = Teuchos::as<ST> (1);
  const ST two = Teuchos::as<ST> (2);

  // Quick solve when the matrix A is the identity.
  if (lambdaMin == one && lambdaMax == lambdaMin) {
    solve (X, D_inv, B);
    return;
  }

  // Initialize coefficients
  const ST alpha = lambdaMax / eigRatio;
  const ST beta = Teuchos::as<ST> (1.1) * lambdaMax;
  const ST delta = two / (beta - alpha);
  const ST theta = (beta + alpha) / two;
  const ST s1 = theta * delta;

  // Fetch cached temporary vectors.
  Teuchos::RCP<MV> V_ptr, W_ptr;
  makeTempMultiVectors (V_ptr, W_ptr, B);
  //MV V (B.getMap (), B.getNumVectors (), false);
  //MV W (B.getMap (), B.getNumVectors (), false);
  MV& V = *V_ptr;
  MV& W = *W_ptr;
  const ST oneOverTheta = one / theta;

  // Special case for the first iteration.
  if (! zeroStartingSolution_) {
    computeResidual (V, B, A, X); // V = B - A*X
    solve (W, one/theta, D_inv, V); // W = (1/theta)*D_inv*(B-A*X)
    X.update (one, W, one); // X = X + W
  } else {
    solve (W, one/theta, D_inv, B); // W = (1/theta)*D_inv*B
    X = W; // X = 0 + W
  }

  // The rest of the iterations.
  ST rhok = one / s1;
  ST rhokp1, dtemp1, dtemp2;
  for (int deg = 1; deg < numIters; ++deg) {
    computeResidual (V, B, A, X); // V = B - A*X

    rhokp1 = one / (two * s1 - rhok);
    dtemp1 = rhokp1 * rhok;
    dtemp2 = two * rhokp1 * delta;
    rhok = rhokp1;

    W.scale (dtemp1);
    W.elementWiseMultiply (dtemp2, D_inv, V, one);
    X.update (one, W, one);
  }
}

template<class ScalarType, class MV, class MAT>
typename Chebyshev<ScalarType, MV, MAT>::ST
Chebyshev<ScalarType, MV, MAT>::
powerMethod (const MAT& A, const V& D_inv, const int numIters) {
  const ST zero = Teuchos::as<ST> (0);
  const ST one = Teuchos::as<ST> (1);
  ST lambdaMax = zero;
  ST RQ_top, RQ_bottom, norm;

  V x (A.getDomainMap ());
  V y (A.getRangeMap ());
  x.randomize ();
  norm = x.norm2 ();
  TEUCHOS_TEST_FOR_EXCEPTION(norm == zero, std::runtime_error, 
    "Chebyshev::powerMethod: Tpetra::Vector::randomize filled the vector "
    "with zeros.  This is not impossible, but is unlikely.");

  x.scale (one / norm);
  for (int iter = 0; iter < numIters; ++iter) {
    A.apply (x, y);
    solve (y, D_inv, y);
    RQ_top = y.dot (x);
    RQ_bottom = x.dot (x);
    lambdaMax = RQ_top / RQ_bottom;
    norm = y.norm2 ();
    if (norm == zero) { // Return something reasonable.
      return zero;
    }
    x.update (one / norm, y, zero);
  }
  return lambdaMax;
}

template<class ScalarType, class MV, class MAT>
Teuchos::RCP<const MAT>
Chebyshev<ScalarType, MV, MAT>::
getMatrix() const {
  return A_;
}

template<class ScalarType, class MV, class MAT>
bool
Chebyshev<ScalarType, MV, MAT>::
hasTransposeApply() const {
  // Technically, this is true, because the matrix must be symmetric.
  return true;
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
makeTempMultiVectors (Teuchos::RCP<MV>& V,
		      Teuchos::RCP<MV>& W,
		      const MV& X)
{
  // Don't fill the vectors with zeros.
  if (V_.is_null ()) {
    V_ = Teuchos::rcp (new MV (X.getMap (), X.getNumVectors (), false));
  }
  if (W_.is_null ()) {
    W_ = Teuchos::rcp (new MV (X.getMap (), X.getNumVectors (), false));
  }
  V = V_;
  W = W_;
}

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_CHEBYSHEV_DEF_HPP
