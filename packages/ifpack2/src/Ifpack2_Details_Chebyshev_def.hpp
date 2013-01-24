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
Chebyshev<ScalarType, MV, MAT>::
Chebyshev (Teuchos::RCP<const MAT> A) : 
  A_ (A),
  computedLambdaMax_ (STS::nan ()),
  computedLambdaMin_ (STS::nan ()),
  lambdaMax_ (STS::nan ()),
  lambdaMin_ (STS::nan ()),
  eigRatio_ (Teuchos::as<ST> (30)),
  numIters_ (1),
  eigMaxIters_ (10),
  zeroStartingSolution_ (true),
  textbookAlgorithm_ (false)
{}

template<class ScalarType, class MV, class MAT>
Chebyshev<ScalarType, MV, MAT>::
Chebyshev (Teuchos::RCP<const MAT> A, Teuchos::ParameterList& params) : 
  A_ (A), 
  computedLambdaMax_ (STS::nan ()),
  computedLambdaMin_ (STS::nan ()),
  lambdaMax_ (STS::nan ()),
  lambdaMin_ (STS::nan ()),
  eigRatio_ (Teuchos::as<ST> (30)),
  numIters_ (1),
  eigMaxIters_ (10),
  zeroStartingSolution_ (true),
  textbookAlgorithm_ (false)
{
  setParameters (params);
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
setParameters (Teuchos::ParameterList& plist) {
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
  const int defaultNumIters = 1;
  const int defaultEigMaxIters = 10;
  const bool defaultZeroStartingSolution = true; // Ifpack::Chebyshev default
  const bool defaultTextbookAlgorithm = false;

  // We'll set the instance data transactionally, after all reads
  // from the ParameterList.  That way, if any of the ParameterList
  // reads fail (e.g., due to the wrong parameter type), we will not
  // have left the instance data in a half-changed state.
  ST lambdaMax = defaultLambdaMax;
  ST lambdaMin = defaultLambdaMin;
  ST eigRatio = defaultEigRatio;
  int numIters = defaultNumIters;
  int eigMaxIters = defaultEigMaxIters;
  bool zeroStartingSolution = defaultZeroStartingSolution;
  bool textbookAlgorithm = defaultTextbookAlgorithm;

  // Defer all externally visible side effects until we have
  // finished all ParameterList interaction.  This makes the method
  // satisfy the strong exception guarantee.
  lambdaMax = plist.get ("chebyshev: max eigenvalue", lambdaMax);
  lambdaMin = plist.get ("chebyshev: min eigenvalue", lambdaMin);

  // The last parameter name overrides the first.
  eigRatio = plist.get ("smoother: Chebyshev alpha", eigRatio); // ML compatibility
  eigRatio = plist.get ("chebyshev: ratio eigenvalue", eigRatio);

  // Each parameter name overrides the one above it.
  numIters = plist.get ("smoother: sweeps", numIters); // ML compatibility
  numIters = plist.get ("relaxation: sweeps", numIters); // Ifpack compatibility
  numIters = plist.get ("chebyshev: degree", numIters);

  // The last parameter name overrides the first.
  eigMaxIters = plist.get ("eigen-analysis: iterations", eigMaxIters); // ML compatibility
  eigMaxIters = plist.get ("chebyshev: eigenvalue max iterations", eigMaxIters);

  zeroStartingSolution = plist.get ("chebyshev: zero starting solution", zeroStartingSolution),
    textbookAlgorithm = plist.get ("chebyshev: textbook algorithm", textbookAlgorithm);

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
  TEUCHOS_TEST_FOR_EXCEPTION(
    plist.isParameter ("chebyshev: min diagonal value"),
    std::logic_error,
    "Ifpack2::Chebyshev: The Ifpack parameter \"chebyshev: min diagonal "
    "value\" is not currently supported.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    plist.isParameter ("chebyshev: operator inv diagonal"),
    std::logic_error,
    "Ifpack2::Chebyshev: The Ifpack parameter \"chebyshev: operator inv "
    "diagonal\" is not currently supported.");
  // For now, we only check that this ML parameter, if provided, has
  // the one value that we support.
  std::string eigenAnalysisType ("power-method");
  if (plist.isParameter ("eigen-analysis: type")) {
    eigenAnalysisType = plist.get<std::string> ("eigen-analysis: type");
    TEUCHOS_TEST_FOR_EXCEPTION(
      eigenAnalysisType != "power-method" && eigenAnalysisType != "power method",
      std::invalid_argument,
      "Ifpack2::Chebyshev: This class does support the ML parameter \"eigen-"
      "analysis: type\" for backwards compatibility.  However, Ifpack2 "
      "currently only supports the \"power-method\" option for this "
      "parameter.  This imitates Ifpack, which only implements the power "
      "method for eigenanalysis.");
  }

  lambdaMax_ = lambdaMax;
  lambdaMin_ = lambdaMin;
  eigRatio_ = eigRatio;
  numIters_ = numIters;
  eigMaxIters_ = eigMaxIters;
  zeroStartingSolution_ = zeroStartingSolution;
  textbookAlgorithm_ = textbookAlgorithm;
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
compute () {
  // The matrix may have changed, and the parameters also affect
  // this process, so recompute the inverse diagonal.
  D_ = getDiagonal (*A_);

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
    if (STS::magnitude (lambdaMax - one) < Teuchos::as<ST> (1.0e-6)) {
      lambdaMin = one;
      lambdaMax = lambdaMin;
      eigRatio = one; // Ifpack doesn't include this line.
    }
    ifpackApplyImpl (*A_, B, X, numIters_, lambdaMax, lambdaMin, eigRatio, *D_);
  } else {
    applyImpl (*A_, B, X, numIters_, lambdaMax, lambdaMin, eigRatio, *D_);
  }

  MV R (B.getMap (), B.getNumVectors ());
  computeResidual (R, B, *A_, X);
  Teuchos::Array<MT> norms (B.getNumVectors ());
  R.norm2 (norms ());
  return *std::max_element (norms.begin (), norms.end ());
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
getDiagonal (const MAT& A) {
  Teuchos::RCP<V> D (new V (A.getGraph ()->getRowMap ()));
  A.getLocalDiagCopy (*D);
  // Invert the diagonal entries, replacing entries less than
  // machine precision with machine precision.
  typedef Kokkos::MultiVector<ST, typename MAT::node_type> KMV;
  KMV& localDiag = D->getLocalMVNonConst ();
  typedef Kokkos::DefaultArithmetic<KMV> KMVT;
  KMVT::ReciprocalThreshold (localDiag, STS::eps ());
  return D;
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
applyImpl (const MAT& A,
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
		 const V& D_inv) const
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

  // Define vectors
  MV V (B.getMap (), B.getNumVectors (), false);
  MV W (B.getMap (), B.getNumVectors (), false);
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

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_CHEBYSHEV_DEF_HPP
