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

// Uncommit the #define line below if you want Chebyshev to do extra
// debug checking and generate a lot of debugging output to stderr (on
// all processes in the communicator).  Even if you uncomment this
// line, the debugging code will only be enabled if the CMake option
// Teuchos_ENABLE_DEBUG was set to ON when configuring Trilinos.
//#define IFPACK_DETAILS_CHEBYSHEV_DEBUG 1

namespace Ifpack2 {
namespace Details {

namespace {
  // We use this text a lot in error messages.
  const char computeBeforeApplyReminder[] =
    "This means one of the following:\n"
    "  - you have not yet called compute() on this instance, or \n"
    "  - you didn't call compute() after calling setParameters().\n\n"
    "After creating an Ifpack2::Chebyshev instance,\n"
    "you must _always_ call compute() at least once before calling apply().\n"
    "After calling compute() once, you do not need to call it again,\n"
    "unless the matrix has changed or you have changed parameters\n"
    "(by calling setParameters()).";
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
checkConstructorInput () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(STS::isComplex, std::logic_error,
    "Ifpack2::Details::Chebyshev: This class' implementation of Chebyshev "
    "iteration only works for real-valued, symmetric positive definite "
    "matrices.  However, you instantiated this class for ScalarType="
    << Teuchos::TypeNameTraits<ScalarType>::name () << ", which is a complex-"
    "valued type.  While this may be algorithmically correct if all of the "
    "complex numbers in the matrix have zero imaginary part, we forbid using "
    "complex ScalarType altogether in order to remind you of the limitations "
    "of our implementation (and of the algorithm itself).");
  TEUCHOS_TEST_FOR_EXCEPTION(A_.is_null (), std::invalid_argument,
    "Ifpack2::Chebyshev: Input matrix to constructor is null.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_->getGlobalNumRows() != A_->getGlobalNumCols(),
    std::invalid_argument,
    "Ifpack2::Chebyshev: The input matrix A must be square.  "
    "A has " << A_->getGlobalNumRows() << " rows and "
    << A_->getGlobalNumCols() << " columns.");

  // In a debug build, test that the domain and range Maps of the
  // matrix are the same.
#ifdef HAVE_TEUCHOS_DEBUG
  Teuchos::RCP<const map_type> domainMap = A_->getDomainMap ();
  Teuchos::RCP<const map_type> rangeMap = A_->getRangeMap ();

  // isSameAs is a collective, but if the two pointers are the same,
  // isSameAs will assume that they are the same on all processes, and
  // return true without an all-reduce.
  TEUCHOS_TEST_FOR_EXCEPTION(
     ! domainMap->isSameAs (*rangeMap), std::invalid_argument,
     "Ifpack2::Chebyshev: The domain Map and range Map of the matrix must be "
     "the same (in the sense of isSameAs())." << std::endl << "We only check "
     "for this if Trilinos was built with the CMake configuration option "
     "Teuchos_ENABLE_DEBUG set to ON.");
#endif // HAVE_TEUCHOS_DEBUG
}

template<class ScalarType, class MV, class MAT>
Chebyshev<ScalarType, MV, MAT>::
Chebyshev (Teuchos::RCP<const MAT> A) :
  A_ (A),
  savedDiagOffsets_ (false),
  computedLambdaMax_ (STS::nan ()),
  computedLambdaMin_ (STS::nan ()),
  lambdaMaxForApply_ (STS::nan ()),
  lambdaMinForApply_ (STS::nan ()),
  eigRatioForApply_ (STS::nan ()),
  userLambdaMax_ (STS::nan ()),
  userLambdaMin_ (STS::nan ()),
  userEigRatio_ (Teuchos::as<ST> (30)),
  minDiagVal_ (STS::eps ()),
  numIters_ (1),
  eigMaxIters_ (10),
  zeroStartingSolution_ (true),
  assumeMatrixUnchanged_ (false),
  textbookAlgorithm_ (false),
  computeMaxResNorm_ (false)
{
  checkConstructorInput ();
}

template<class ScalarType, class MV, class MAT>
Chebyshev<ScalarType, MV, MAT>::
Chebyshev (Teuchos::RCP<const MAT> A, Teuchos::ParameterList& params) :
  A_ (A),
  savedDiagOffsets_ (false),
  computedLambdaMax_ (STS::nan ()),
  computedLambdaMin_ (STS::nan ()),
  lambdaMaxForApply_ (STS::nan ()),
  lambdaMinForApply_ (STS::nan ()),
  eigRatioForApply_ (STS::nan ()),
  userLambdaMax_ (STS::nan ()),
  userLambdaMin_ (STS::nan ()),
  userEigRatio_ (Teuchos::as<ST> (30)),
  minDiagVal_ (STS::eps ()),
  numIters_ (1),
  eigMaxIters_ (10),
  zeroStartingSolution_ (true),
  assumeMatrixUnchanged_ (false),
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
  const bool defaultAssumeMatrixUnchanged = false;
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
  bool assumeMatrixUnchanged = defaultAssumeMatrixUnchanged;
  bool textbookAlgorithm = defaultTextbookAlgorithm;
  bool computeMaxResNorm = defaultComputeMaxResNorm;

  // Fetch the parameters from the ParameterList.  Defer all
  // externally visible side effects until we have finished all
  // ParameterList interaction.  This makes the method satisfy the
  // strong exception guarantee.

  // Get the user-supplied inverse diagonal.
  //
  // Check for a raw pointer (const V* or V*), for Ifpack
  // compatibility, as well as for RCP<const V>, RCP<V>, const V, or
  // V.  We'll copy the vector anyway, so it doesn't matter whether
  // it's const or nonconst.
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
      try { // Could the type be RCP<V>?
        RCP<V> userInvDiagNonconst = plist.get<RCP<V> > ("chebyshev: operator inv diagonal");
        userInvDiag = rcp_const_cast<const V> (userInvDiagNonconst);
      } catch (Teuchos::Exceptions::InvalidParameterType&) {
      }
    }
    if (userInvDiag.is_null ()) {
      try { // Could the type be const V?
        // The line below does a deep copy (V::operator=).
        userInvDiag = rcp (new V (plist.get<const V> ("chebyshev: operator inv diagonal")));
      } catch (Teuchos::Exceptions::InvalidParameterType&) {
      }
    }
    if (userInvDiag.is_null ()) {
      try { // Could the type be V?
        V userInvDiagNonconst = plist.get<V> ("chebyshev: operator inv diagonal");
        userInvDiag = rcp (new V (userInvDiagNonconst)); // deep copy
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
    TEUCHOS_TEST_FOR_EXCEPTION(
      STS::isnaninf (lambdaMax), std::invalid_argument,
      "Ifpack2::Chebyshev::setParameters: \"chebyshev: max eigenvalue\" "
      "parameter is NaN or Inf.  This parameter is optional, but if you "
      "choose to supply it, it must have a finite value.");
  }
  if (plist.isParameter ("chebyshev: min eigenvalue")) {
    lambdaMin = plist.get<ST> ("chebyshev: min eigenvalue");
    TEUCHOS_TEST_FOR_EXCEPTION(
      STS::isnaninf (lambdaMin), std::invalid_argument,
      "Ifpack2::Chebyshev::setParameters: \"chebyshev: min eigenvalue\" "
      "parameter is NaN or Inf.  This parameter is optional, but if you "
      "choose to supply it, it must have a finite value.");
  }

  // Only fill in Ifpack2's name for the default parameter, not ML's.
  if (plist.isParameter ("smoother: Chebyshev alpha")) { // ML compatibility
    eigRatio = plist.get<ST> ("smoother: Chebyshev alpha");
  }
  // Ifpack2's name overrides ML's name.
  eigRatio = plist.get ("chebyshev: ratio eigenvalue", eigRatio);
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isnaninf (eigRatio), std::invalid_argument,
    "Ifpack2::Chebyshev::setParameters: \"chebyshev: ratio eigenvalue\" "
    "parameter (also called \"smoother: Chebyshev alpha\") is NaN or Inf.  "
    "This parameter is optional, but if you choose to supply it, it must have "
    "a finite value.");
  // mfh 11 Feb 2013: This class is currently only correct for real
  // Scalar types, but we still want it to build for complex Scalar
  // type so that users of Ifpack2::Factory can build their
  // executables for real or complex Scalar type.  Thus, we take the
  // real parts here, which are always less-than comparable.
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::real (eigRatio) < STS::real (STS::one ()),
    std::invalid_argument,
    "Ifpack2::Chebyshev::setParameters: \"chebyshev: ratio eigenvalue\""
    "parameter (also called \"smoother: Chebyshev alpha\") must be >= 1, "
    "but you supplied the value " << eigRatio << ".");

  // Same name in Ifpack2 and Ifpack.
  minDiagVal = plist.get ("chebyshev: min diagonal value", minDiagVal);
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isnaninf (minDiagVal), std::invalid_argument,
    "Ifpack2::Chebyshev::setParameters: \"chebyshev: min diagonal value\" "
    "parameter is NaN or Inf.  This parameter is optional, but if you choose "
    "to supply it, it must have a finite value.");

  // Only fill in Ifpack2's name, not ML's or Ifpack's.
  if (plist.isParameter ("smoother: sweeps")) { // ML compatibility
    numIters = plist.get<int> ("smoother: sweeps");
  } // Ifpack's name overrides ML's name.
  if (plist.isParameter ("relaxation: sweeps")) { // Ifpack compatibility
    numIters = plist.get<int> ("relaxation: sweeps");
  } // Ifpack2's name overrides Ifpack's name.
  numIters = plist.get ("chebyshev: degree", numIters);
  TEUCHOS_TEST_FOR_EXCEPTION(
    numIters < 0, std::invalid_argument,
    "Ifpack2::Chebyshev::setParameters: \"chebyshev: degree\" parameter (also "
    "called \"smoother: sweeps\" or \"relaxation: sweeps\") must be a "
    "nonnegative integer.  You gave a value of " << numIters << ".");

  // The last parameter name overrides the first.
  if (plist.isParameter ("eigen-analysis: iterations")) { // ML compatibility
    eigMaxIters = plist.get<int> ("eigen-analysis: iterations");
  } // Ifpack2's name overrides ML's name.
  eigMaxIters = plist.get ("chebyshev: eigenvalue max iterations", eigMaxIters);
  TEUCHOS_TEST_FOR_EXCEPTION(
    eigMaxIters < 0, std::invalid_argument,
    "Ifpack2::Chebyshev::setParameters: \"chebyshev: eigenvalue max iterations"
    "\" parameter (also called \"eigen-analysis: iterations\") must be a "
    "nonnegative integer.  You gave a value of " << eigMaxIters << ".");

  zeroStartingSolution = plist.get ("chebyshev: zero starting solution",
                                    zeroStartingSolution);
  assumeMatrixUnchanged = plist.get ("chebyshev: assume matrix does not change",
                                     assumeMatrixUnchanged);

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

  // We've validated all the parameters, so it's safe now to "commit" them.
  userInvDiag_ = userInvDiag;
  userLambdaMax_ = lambdaMax;
  userLambdaMin_ = lambdaMin;
  userEigRatio_ = eigRatio;
  minDiagVal_ = minDiagVal;
  numIters_ = numIters;
  eigMaxIters_ = eigMaxIters;
  zeroStartingSolution_ = zeroStartingSolution;
  assumeMatrixUnchanged_ = assumeMatrixUnchanged;
  textbookAlgorithm_ = textbookAlgorithm;
  computeMaxResNorm_ = computeMaxResNorm;
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
compute () {
  using std::endl;

  // Some of the optimizations below only work if A_ is a
  // Tpetra::CrsMatrix.  We'll make our best guess about its type
  // here, since we have no way to get back the original fifth
  // template parameter.
  typedef Tpetra::CrsMatrix<typename MV::scalar_type,
    typename MV::local_ordinal_type,
    typename MV::global_ordinal_type,
    typename MV::node_type> crs_matrix_type;

  // If A_ is a CrsMatrix and its graph is constant, we presume that
  // the user plans to reuse the structure of A_, but possibly change
  // A_'s values before each compute() call.  This is the intended use
  // case for caching the offsets of the diagonal entries of A_, to
  // speed up extraction of diagonal entries on subsequent compute()
  // calls.

  // FIXME (mfh 22 Jan 2013, 10 Feb 2013) In all cases when we use
  // isnaninf() in this method, we really only want to check if the
  // number is NaN.  Inf means something different.  However,
  // Teuchos::ScalarTraits doesn't distinguish the two cases.

  if (userInvDiag_.is_null ()) {
    Teuchos::RCP<const crs_matrix_type> A_crsMat =
      Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_);

    if (D_.is_null ()) { // We haven't computed D_ before
      if (! A_crsMat.is_null () && A_crsMat->isStaticGraph ()) {
        // It's a CrsMatrix with a const graph; cache diagonal offsets.
        A_crsMat->getLocalDiagOffsets (diagOffsets_);
        savedDiagOffsets_ = true;
        D_ = makeInverseDiagonal (*A_, true);
      }
      else { // either A_ is not a CrsMatrix, or its graph is nonconst
        D_ = makeInverseDiagonal (*A_);
      }
    }
    else if (! assumeMatrixUnchanged_) { // D_ exists but A_ may have changed
      if (! A_crsMat.is_null () && A_crsMat->isStaticGraph ()) {
        // It's a CrsMatrix with a const graph; cache diagonal offsets
        // if we haven't already.
        if (! savedDiagOffsets_) {
          A_crsMat->getLocalDiagOffsets (diagOffsets_);
          savedDiagOffsets_ = true;
        }
        // Now we're guaranteed to have cached diagonal offsets.
        D_ = makeInverseDiagonal (*A_, true);
      }
      else { // either A_ is not a CrsMatrix, or its graph is nonconst
        D_ = makeInverseDiagonal (*A_);
      }
    }
  }
  else { // the user provided an inverse diagonal
    D_ = userInvDiag_;
  }

  // Have we estimated eigenvalues before?
  const bool computedEigenvalueEstimates =
    STS::isnaninf (computedLambdaMax_) || STS::isnaninf (computedLambdaMin_);

  // Only recompute the eigenvalue estimates if
  // - we are supposed to assume that the matrix may have changed, or
  // - they haven't been computed before, and the user hasn't given
  //   us at least an estimate of the max eigenvalue.
  //
  // We at least need an estimate of the max eigenvalue.  This is the
  // most important one if using Chebyshev as a smoother.
  if (! assumeMatrixUnchanged_ ||
      (! computedEigenvalueEstimates && STS::isnaninf (userLambdaMax_))) {
    const ST computedLambdaMax = powerMethod (*A_, *D_, eigMaxIters_);
    TEUCHOS_TEST_FOR_EXCEPTION(
      STS::isnaninf (computedLambdaMax),
      std::runtime_error,
      "Ifpack2::Chebyshev::compute: Estimation of the max eigenvalue "
      "of D^{-1} A failed, by producing Inf or NaN.  This probably means that "
      "the matrix contains Inf or NaN values, or that it is badly scaled.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      STS::isnaninf (userEigRatio_),
      std::logic_error,
      "Ifpack2::Chebyshev::compute: userEigRatio_ is Inf or NaN."
      << endl << "This should be impossible." << endl <<
      "Please report this bug to the Ifpack2 developers.");

    // The power method doesn't estimate the min eigenvalue, so we
    // do our best to provide an estimate.  userEigRatio_ has a
    // reasonable default value, and if the user provided it, we
    // have already checked that its value is finite and >= 1.
    const ST computedLambdaMin = computedLambdaMax / userEigRatio_;

    // Defer "committing" results until all computations succeeded.
    computedLambdaMax_ = computedLambdaMax;
    computedLambdaMin_ = computedLambdaMin;
  } else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      STS::isnaninf (userLambdaMax_) && STS::isnaninf (computedLambdaMax_),
      std::logic_error,
      "Ifpack2::Chebyshev::compute: " << endl <<
      "Both userLambdaMax_ and computedLambdaMax_ are Inf or NaN."
      << endl << "This should be impossible." << endl <<
      "Please report this bug to the Ifpack2 developers.");
  }

  ////////////////////////////////////////////////////////////////////
  // Figure out the eigenvalue estimates that apply() will use.
  ////////////////////////////////////////////////////////////////////

  // Always favor the user's max eigenvalue estimate, if provided.
  if (STS::isnaninf (userLambdaMax_)) {
    lambdaMaxForApply_ = computedLambdaMax_;
  } else {
    lambdaMaxForApply_ = userLambdaMax_;
  }
  // mfh 11 Feb 2013: For now, we imitate Ifpack by ignoring the
  // user's min eigenvalue estimate, and using the given eigenvalue
  // ratio to estimate the min eigenvalue.  We could instead do this:
  // favor the user's eigenvalue ratio estimate, but if it's not
  // provided, use lambdaMax / lambdaMin.  However, we want Chebyshev
  // to have sensible smoother behavior if the user did not provide
  // eigenvalue estimates.  Ifpack's behavior attempts to push down
  // the error terms associated with the largest eigenvalues, while
  // expecting that users will only want a small number of iterations,
  // so that error terms associated with the smallest eigenvalues
  // won't grow too much.  This is sensible behavior for a smoother.
  lambdaMinForApply_ = lambdaMaxForApply_ / userEigRatio_;
  eigRatioForApply_ = userEigRatio_;

  if (! textbookAlgorithm_) {
    // Ifpack has a special-case modification of the eigenvalue bounds
    // for the case where the max eigenvalue estimate is close to one.
    const ST one = Teuchos::as<ST> (1);
    if (STS::magnitude (lambdaMaxForApply_ - one) < Teuchos::as<MT> (1.0e-6)) {
      lambdaMinForApply_ = one;
      lambdaMaxForApply_ = lambdaMinForApply_;
      eigRatioForApply_ = one; // Ifpack doesn't include this line.
    }
  }
} //compute()

template<class ScalarType, class MV, class MAT>
ScalarType
Chebyshev<ScalarType, MV, MAT>::
getLambdaMaxForApply() const {
  return lambdaMaxForApply_;
}

template<class ScalarType, class MV, class MAT>
typename Chebyshev<ScalarType, MV, MAT>::MT
Chebyshev<ScalarType, MV, MAT>::
apply (const MV& B, MV& X) {
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isnaninf (lambdaMaxForApply_), std::runtime_error,
    "Ifpack2::Chebyshev::apply: There is no estimate for the max eigenvalue."
    << std::endl << computeBeforeApplyReminder);
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isnaninf (lambdaMinForApply_), std::runtime_error,
    "Ifpack2::Chebyshev::apply: There is no estimate for the min eigenvalue."
    << std::endl << computeBeforeApplyReminder);
  TEUCHOS_TEST_FOR_EXCEPTION(
    STS::isnaninf (eigRatioForApply_), std::runtime_error,
    "Ifpack2::Chebyshev::apply: There is no estimate for the ratio of the max "
    "eigenvalue to the min eigenvalue."
    << std::endl << computeBeforeApplyReminder);
  TEUCHOS_TEST_FOR_EXCEPTION(
    D_.is_null (), std::runtime_error,
    "Ifpack2::Chebyshev::apply: "
    "The vector of inverse diagonal entries of the matrix has not yet been "
    "computed." << std::endl << computeBeforeApplyReminder);

  if (textbookAlgorithm_) {
    textbookApplyImpl (*A_, B, X, numIters_, lambdaMaxForApply_,
                       lambdaMinForApply_, eigRatioForApply_, *D_);
  } else {
    ifpackApplyImpl (*A_, B, X, numIters_, lambdaMaxForApply_,
                     lambdaMinForApply_, eigRatioForApply_, *D_);
  }

  if (computeMaxResNorm_ && B.getNumVectors () > 0) {
    MV R (B.getMap (), B.getNumVectors ());
    computeResidual (R, B, *A_, X);
    Teuchos::Array<MT> norms (B.getNumVectors ());
    R.norm2 (norms ());
    return *std::max_element (norms.begin (), norms.end ());
  } else {
    return Teuchos::ScalarTraits<MT>::zero ();
  }
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
print (std::ostream& out) {
  this->describe (* (Teuchos::getFancyOStream (Teuchos::rcpFromRef (out))),
                  Teuchos::VERB_MEDIUM);
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
makeInverseDiagonal (const MAT& A, const bool useDiagOffsets) const {
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using Teuchos::rcp_dynamic_cast;

  RCP<V> D_rowMap (new V (A.getGraph ()->getRowMap ()));
  if (useDiagOffsets) {
    // The optimizations below only work if A_ is a Tpetra::CrsMatrix.
    // We'll make our best guess about its type here, since we have no
    // way to get back the original fifth template parameter.
    typedef Tpetra::CrsMatrix<typename MV::scalar_type,
      typename MV::local_ordinal_type,
      typename MV::global_ordinal_type,
      typename MV::node_type> crs_matrix_type;
    RCP<const crs_matrix_type> A_crsMat =
      rcp_dynamic_cast<const crs_matrix_type> (rcpFromRef (A));
    if (! A_crsMat.is_null ()) {
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! savedDiagOffsets_, std::logic_error,
        "Ifpack2::Details::Chebyshev::makeInverseDiagonal: "
        "It is not allowed to call this method with useDiagOffsets=true, "
        "if you have not previously saved offsets of diagonal entries.  "
        "This situation should never arise if this class is used properly.  "
        "Please report this bug to the Ifpack2 developers.");
      A_crsMat->getLocalDiagCopy (*D_rowMap, diagOffsets_ ());
    }
  }
  else {
    A.getLocalDiagCopy (*D_rowMap);
  }
  RCP<V> D_rangeMap = makeRangeMapVector (*D_rowMap);

  // Invert the diagonal entries, replacing entries less (in
  // magnitude) than the user-specified value with that value.
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
typename Chebyshev<ScalarType, MV, MAT>::MT
Chebyshev<ScalarType, MV, MAT>::maxNormInf (const MV& X) {
  std::vector<MT> norms (X.getNumVectors ());
  X.normInf (norms);
  return *std::max_element (norms.begin (), norms.end ());
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
#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
  using std::cerr;
  using std::endl;
  cerr << "\\|B\\|_{\\infty} = " << maxNormInf (B) << endl;
  cerr << "\\|X\\|_{\\infty} = " << maxNormInf (X) << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

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

#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
  cerr << "alpha = " << alpha << endl
       << "beta = " << beta << endl
       << "delta = " << delta << endl
       << "theta = " << theta << endl
       << "s1 = " << s1 << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

  // Fetch cached temporary vectors.
  Teuchos::RCP<MV> V_ptr, W_ptr;
  makeTempMultiVectors (V_ptr, W_ptr, B);

  // mfh 28 Jan 2013: We write V1 instead of V, so as not to confuse
  // the multivector V with the typedef V (for Tpetra::Vector).
  //MV V1 (B.getMap (), B.getNumVectors (), false);
  //MV W (B.getMap (), B.getNumVectors (), false);
  MV& V1 = *V_ptr;
  MV& W = *W_ptr;

#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
  cerr << "Iteration " << 1 << ":" << endl
       << "- \\|D\\|_{\\infty} = " << D_->normInf () << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

  // Special case for the first iteration.
  if (! zeroStartingSolution_) {
    computeResidual (V1, B, A, X); // V1 = B - A*X

#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
    cerr << "- \\|B - A*X\\|_{\\infty} = " << maxNormInf (V1) << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

    solve (W, one/theta, D_inv, V1); // W = (1/theta)*D_inv*(B-A*X)

#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
    cerr << "- \\|W\\|_{\\infty} = " << maxNormInf (W) << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

    X.update (one, W, one); // X = X + W
  } else {
    solve (W, one/theta, D_inv, B); // W = (1/theta)*D_inv*B

#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
    cerr << "- \\|W\\|_{\\infty} = " << maxNormInf (W) << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

    X = W; // X = 0 + W
  }
#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
  cerr << "- \\|X\\|_{\\infty} = " << maxNormInf (X) << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

  // The rest of the iterations.
  ST rhok = one / s1;
  ST rhokp1, dtemp1, dtemp2;
  for (int deg = 1; deg < numIters; ++deg) {

#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
    cerr << "Iteration " << deg+1 << ":" << endl;
    cerr << "- \\|D\\|_{\\infty} = " << D_->normInf () << endl;
    cerr << "- \\|B\\|_{\\infty} = " << maxNormInf (B) << endl;
    cerr << "- \\|A\\|_{\\text{frob}} = " << A_->getFrobeniusNorm () << endl;
    cerr << "- rhok = " << rhok << endl;
    V1.putScalar (STS::zero ());
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

    computeResidual (V1, B, A, X); // V1 = B - A*X

#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
    cerr << "- \\|B - A*X\\|_{\\infty} = " << maxNormInf (V1) << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

    rhokp1 = one / (two * s1 - rhok);
    dtemp1 = rhokp1 * rhok;
    dtemp2 = two * rhokp1 * delta;
    rhok = rhokp1;

#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
    cerr << "- dtemp1 = " << dtemp1 << endl
         << "- dtemp2 = " << dtemp2 << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG

    W.scale (dtemp1);
    W.elementWiseMultiply (dtemp2, D_inv, V1, one);
    X.update (one, W, one);

#ifdef HAVE_TEUCHOS_DEBUG
#ifdef IFPACK_DETAILS_CHEBYSHEV_DEBUG
    cerr << "- \\|W\\|_{\\infty} = " << maxNormInf (W) << endl;
    cerr << "- \\|X\\|_{\\infty} = " << maxNormInf (X) << endl;
#endif // IFPACK_DETAILS_CHEBYSHEV_DEBUG
#endif // HAVE_TEUCHOS_DEBUG
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
    "Ifpack2::Chebyshev::powerMethod: "
    "Tpetra::Vector's randomize() method filled the vector "
    "with zeros.  This is not impossible, but is unlikely.  "
    "It's far more likely that there is a bug in Tpetra.");

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
makeTempMultiVectors (Teuchos::RCP<MV>& V1,
                      Teuchos::RCP<MV>& W,
                      const MV& X)
{
  if (V_.is_null ()) {
    V_ = Teuchos::rcp (new MV (X.getMap (), X.getNumVectors (), false));
  }
  //W must be initialized to zero when it is used as a multigrid smoother.
  if (W_.is_null ()) {
    W_ = Teuchos::rcp (new MV (X.getMap (), X.getNumVectors (), true));
  }
  V1 = V_;
  W = W_;
}

template<class ScalarType, class MV, class MAT>
std::string
Chebyshev<ScalarType, MV, MAT>::
description() const {
  std::ostringstream oss;
  oss << "Ifpack2::Details::Chebyshev : "
      << "degree = " << numIters_
      << ", lambdaMax = " << lambdaMaxForApply_
      << ", alpha = " << eigRatioForApply_
      << ", lambdaMin = " << lambdaMinForApply_;

  return oss.str();
}

template<class ScalarType, class MV, class MAT>
void
Chebyshev<ScalarType, MV, MAT>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::TypeNameTraits;
  using std::endl;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;
  if (vl > Teuchos::VERB_NONE) {
    if (vl == Teuchos::VERB_LOW) {
      out << "Ifpack2::Details::Chebyshev" << endl;
    } else { // vl > Teuchos::VERB_LOW
      // YAML flow mapping style (with curly braces) lets the key (the
      // above "description") be arbitrarily long.  Simple indenting
      // limits the key length to 1024 characters.  It's not an issue
      // right now, but it's better to be more general.
      out << "Ifpack2::Details::Chebyshev {" << endl;
      Teuchos::OSTab tab1 (Teuchos::rcpFromRef (out));
      out << "Template parameters: {" << endl;
      {
        Teuchos::OSTab tab2 (Teuchos::rcpFromRef (out));
        out << "ScalarType: \"" << TypeNameTraits<ScalarType>::name () << "\"" << endl
            << "MV: \"" << TypeNameTraits<MV>::name () << "\"" << endl
            << "MAT: \"" << TypeNameTraits<MAT>::name () << "\"" << endl;
      }
      // "Computed parameters" literally means "parameters whose
      // values were computed by compute()."
      out << "}" << endl << "Computed parameters: {" << endl;
      {
        Teuchos::OSTab tab2 (Teuchos::rcpFromRef (out));
        // Users might want to see the values in the computed inverse
        // diagonal, so we print them out at the highest verbosity.
        out << "D_: ";
        if (D_.is_null ()) {
          out << "unset" << endl;
        } else if (vl <= Teuchos::VERB_HIGH) {
          out << "set" << endl;
        } else { // D_ not null and vl > Teuchos::VERB_HIGH
          out << "{" << endl;
          {
            Teuchos::OSTab tab3 (Teuchos::rcpFromRef (out));
            D_->describe (out, vl);
          }
          out << "}" << endl;
        }
        // V_ and W_ are scratch space; their values are irrelevant.
        // All that matters is whether or not they have been set.
        out << "V_: " << (V_.is_null () ? "unset" : "set") << endl
            << "W_: " << (W_.is_null () ? "unset" : "set") << endl
            << "computedLambdaMax_: " << computedLambdaMax_ << endl
            << "computedLambdaMin_: " << computedLambdaMin_ << endl
            << "lambdaMaxForApply_: " << lambdaMaxForApply_ << endl
            << "lambdaMinForApply_: " << lambdaMinForApply_ << endl
            << "eigRatioForApply_: " << eigRatioForApply_ << endl;
      }
      out << "}" << endl << "User parameters: {" << endl;
      {
        Teuchos::OSTab tab2 (Teuchos::rcpFromRef (out));
        out << "userInvDiag_: ";
        if (userInvDiag_.is_null ()) {
          out << "unset" << endl;
        } else if (vl <= Teuchos::VERB_HIGH) {
          out << "set" << endl;
        } else { // userInvDiag_ not null and vl > Teuchos::VERB_HIGH
          out << "{" << endl;
          {
            Teuchos::OSTab tab3 (Teuchos::rcpFromRef (out));
            userInvDiag_->describe (out, vl);
          }
          out << "}" << endl;
        }
        out << "userLambdaMax_: " << userLambdaMax_ << endl
            << "userLambdaMin_: " << userLambdaMin_ << endl
            << "userEigRatio_: " << userEigRatio_ << endl
            << "numIters_: " << numIters_ << endl
            << "eigMaxIters_: " << eigMaxIters_ << endl
            << "zeroStartingSolution_: " << zeroStartingSolution_ << endl
            << "assumeMatrixUnchanged_: " << assumeMatrixUnchanged_ << endl
            << "textbookAlgorithm_: " << textbookAlgorithm_ << endl
            << "computeMaxResNorm_: " << computeMaxResNorm_ << endl;
      }
      out << "}" << endl;
    }
  }
}

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_CHEBYSHEV_DEF_HPP
