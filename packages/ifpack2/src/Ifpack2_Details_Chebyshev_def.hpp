/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_DETAILS_CHEBYSHEV_DEF_HPP
#define IFPACK2_DETAILS_CHEBYSHEV_DEF_HPP

/// \file Ifpack2_Details_Chebyshev_def.hpp
/// \brief Definition of Chebyshev implementation
/// \author Mark Hoemmen
///
/// This file is meant for Ifpack2 developers only, not for users.
/// It defines a new implementation of Chebyshev iteration.

#include "Ifpack2_Details_Chebyshev_decl.hpp"
#include "Kokkos_ArithTraits.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include <cmath>
#include <iostream>

namespace Ifpack2 {
namespace Details {

namespace { // (anonymous)

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

} // namespace (anonymous)

// ReciprocalThreshold stuff below needs to be in a namspace visible outside
// of this file
template<class XV, class SizeType = typename XV::size_type>
struct V_ReciprocalThresholdSelfFunctor
{
  typedef typename XV::execution_space execution_space;
  typedef typename XV::non_const_value_type value_type;
  typedef SizeType size_type;
  typedef Kokkos::Details::ArithTraits<value_type> KAT;
  typedef typename KAT::mag_type mag_type;

  XV X_;
  const value_type minVal_;
  const mag_type minValMag_;

  V_ReciprocalThresholdSelfFunctor (const XV& X,
                                    const value_type& min_val) :
    X_ (X),
    minVal_ (min_val),
    minValMag_ (KAT::abs (min_val))
  {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const size_type& i) const
  {
    const mag_type X_i_abs = KAT::abs (X_[i]);

    if (X_i_abs < minValMag_) {
      X_[i] = minVal_;
    }
    else {
      X_[i] = KAT::one () / X_[i];
    }
  }
};

template<class XV, class SizeType = typename XV::size_type>
struct LocalReciprocalThreshold {
  static void
  compute (const XV& X,
           const typename XV::non_const_value_type& minVal)
  {
    typedef typename XV::execution_space execution_space;
    Kokkos::RangePolicy<execution_space, SizeType> policy (0, X.extent (0));
    V_ReciprocalThresholdSelfFunctor<XV, SizeType> op (X, minVal);
    Kokkos::parallel_for (policy, op);
  }
};

template <class TpetraVectorType,
          const bool classic = TpetraVectorType::node_type::classic>
struct GlobalReciprocalThreshold {};

template <class TpetraVectorType>
struct GlobalReciprocalThreshold<TpetraVectorType, true> {
  static void
  compute (TpetraVectorType& V,
           const typename TpetraVectorType::scalar_type& min_val)
  {
    typedef typename TpetraVectorType::scalar_type scalar_type;
    typedef typename TpetraVectorType::mag_type mag_type;
    typedef Kokkos::Details::ArithTraits<scalar_type> STS;

    const scalar_type ONE = STS::one ();
    const mag_type min_val_abs = STS::abs (min_val);

    Teuchos::ArrayRCP<scalar_type> V_0 = V.getDataNonConst (0);
    const size_t lclNumRows = V.getLocalLength ();

    for (size_t i = 0; i < lclNumRows; ++i) {
      const scalar_type V_0i = V_0[i];
      if (STS::abs (V_0i) < min_val_abs) {
        V_0[i] = min_val;
      } else {
        V_0[i] = ONE / V_0i;
      }
    }
  }
};

template <class TpetraVectorType>
struct GlobalReciprocalThreshold<TpetraVectorType, false> {
  static void
  compute (TpetraVectorType& X,
           const typename TpetraVectorType::scalar_type& minVal)
  {
    typedef typename TpetraVectorType::impl_scalar_type value_type;
    typedef typename TpetraVectorType::device_type::memory_space memory_space;

    X.template sync<memory_space> ();
    X.template modify<memory_space> ();

    const value_type minValS = static_cast<value_type> (minVal);
    auto X_0 = Kokkos::subview (X.template getLocalView<memory_space> (),
                                Kokkos::ALL (), 0);
    LocalReciprocalThreshold<decltype (X_0) >::compute (X_0, minValS);
  }
};

// Utility function for inverting diagonal with threshold.
template <typename S, typename L, typename G, typename N>
void
reciprocal_threshold (Tpetra::Vector<S,L,G,N>& V, const S& minVal)
{
  GlobalReciprocalThreshold<Tpetra::Vector<S,L,G,N> >::compute (V, minVal);
}

namespace { // (anonymous)

// Functor for making sure the real parts of all entries of a vector
// are nonnegative.  We use this in computeInitialGuessForPowerMethod
// below.
template<class OneDViewType,
         class LocalOrdinal = typename OneDViewType::size_type>
class PositivizeVector {
  static_assert (Kokkos::Impl::is_view<OneDViewType>::value,
                 "OneDViewType must be a 1-D Kokkos::View.");
  static_assert (static_cast<int> (OneDViewType::rank) == 1,
                 "This functor only works with a 1-D View.");
  static_assert (std::is_integral<LocalOrdinal>::value,
                 "The second template parameter, LocalOrdinal, "
                 "must be an integer type.");
public:
  PositivizeVector (const OneDViewType& x) : x_ (x) {}

  KOKKOS_INLINE_FUNCTION void
  operator () (const LocalOrdinal& i) const
  {
    typedef typename OneDViewType::non_const_value_type IST;
    typedef Kokkos::Details::ArithTraits<IST> STS;
    typedef Kokkos::Details::ArithTraits<typename STS::mag_type> STM;

    if (STS::real (x_(i)) < STM::zero ()) {
      x_(i) = -x_(i);
    }
  }

private:
  OneDViewType x_;
};

} // namespace (anonymous)

template<class ScalarType, class MV>
void Chebyshev<ScalarType, MV>::checkInputMatrix () const
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! A_.is_null () && A_->getGlobalNumRows () != A_->getGlobalNumCols (),
    std::invalid_argument,
    "Ifpack2::Chebyshev: The input matrix A must be square.  "
    "A has " << A_->getGlobalNumRows () << " rows and "
    << A_->getGlobalNumCols () << " columns.");

  // In debug mode, test that the domain and range Maps of the matrix
  // are the same.
  if (debug_ && ! A_.is_null ()) {
    Teuchos::RCP<const map_type> domainMap = A_->getDomainMap ();
    Teuchos::RCP<const map_type> rangeMap = A_->getRangeMap ();

    // isSameAs is a collective, but if the two pointers are the same,
    // isSameAs will assume that they are the same on all processes, and
    // return true without an all-reduce.
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! domainMap->isSameAs (*rangeMap), std::invalid_argument,
      "Ifpack2::Chebyshev: The domain Map and range Map of the matrix must be "
      "the same (in the sense of isSameAs())." << std::endl << "We only check "
      "for this in debug mode.");
  }
}


template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
checkConstructorInput () const
{
  // mfh 12 Aug 2016: The if statement avoids an "unreachable
  // statement" warning for the checkInputMatrix() call, when
  // STS::isComplex is false.
  if (STS::isComplex) {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::logic_error, "Ifpack2::Chebyshev: This class' implementation "
       "of Chebyshev iteration only works for real-valued, symmetric positive "
       "definite matrices.  However, you instantiated this class for ScalarType"
       " = " << Teuchos::TypeNameTraits<ScalarType>::name () << ", which is a "
       "complex-valued type.  While this may be algorithmically correct if all "
       "of the complex numbers in the matrix have zero imaginary part, we "
       "forbid using complex ScalarType altogether in order to remind you of "
       "the limitations of our implementation (and of the algorithm itself).");
  }
  else {
    checkInputMatrix ();
  }
}

template<class ScalarType, class MV>
Chebyshev<ScalarType, MV>::
Chebyshev (Teuchos::RCP<const row_matrix_type> A) :
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
  computeMaxResNorm_ (false),
  debug_ (false)
{
  checkConstructorInput ();
}

template<class ScalarType, class MV>
Chebyshev<ScalarType, MV>::
Chebyshev (Teuchos::RCP<const row_matrix_type> A, Teuchos::ParameterList& params) :
  A_ (A),
  savedDiagOffsets_ (false),
  computedLambdaMax_ (STS::nan ()),
  computedLambdaMin_ (STS::nan ()),
  lambdaMaxForApply_ (STS::nan ()),
  boostFactor_ (static_cast<MT> (1.1)),
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
  computeMaxResNorm_ (false),
  debug_ (false)
{
  checkConstructorInput ();
  setParameters (params);
}

template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
setParameters (Teuchos::ParameterList& plist)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;

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
  const MT defaultBoostFactor = static_cast<MT> (1.1);
  const ST defaultMinDiagVal = STS::eps ();
  const int defaultNumIters = 1;
  const int defaultEigMaxIters = 10;
  const bool defaultZeroStartingSolution = true; // Ifpack::Chebyshev default
  const bool defaultAssumeMatrixUnchanged = false;
  const bool defaultTextbookAlgorithm = false;
  const bool defaultComputeMaxResNorm = false;
  const bool defaultDebug = false;

  // We'll set the instance data transactionally, after all reads
  // from the ParameterList.  That way, if any of the ParameterList
  // reads fail (e.g., due to the wrong parameter type), we will not
  // have left the instance data in a half-changed state.
  RCP<const V> userInvDiagCopy; // if nonnull: deep copy of user's Vector
  ST lambdaMax = defaultLambdaMax;
  ST lambdaMin = defaultLambdaMin;
  ST eigRatio = defaultEigRatio;
  MT boostFactor = defaultBoostFactor;
  ST minDiagVal = defaultMinDiagVal;
  int numIters = defaultNumIters;
  int eigMaxIters = defaultEigMaxIters;
  bool zeroStartingSolution = defaultZeroStartingSolution;
  bool assumeMatrixUnchanged = defaultAssumeMatrixUnchanged;
  bool textbookAlgorithm = defaultTextbookAlgorithm;
  bool computeMaxResNorm = defaultComputeMaxResNorm;
  bool debug = defaultDebug;

  // Fetch the parameters from the ParameterList.  Defer all
  // externally visible side effects until we have finished all
  // ParameterList interaction.  This makes the method satisfy the
  // strong exception guarantee.

  if (plist.isParameter ("debug")) {
    try { // a bool
      debug = plist.get<bool> ("debug");
    }
    catch (Teuchos::Exceptions::InvalidParameterType&) {}

    try { // an int instead of a bool
      int debugInt = plist.get<bool> ("debug");
      debug = debugInt != 0;
    }
    catch (Teuchos::Exceptions::InvalidParameterType&) {}
  }

  // Get the user-supplied inverse diagonal.
  //
  // Check for a raw pointer (const V* or V*), for Ifpack
  // compatibility, as well as for RCP<const V>, RCP<V>, const V, or
  // V.  We'll make a deep copy of the vector at the end of this
  // method anyway, so its const-ness doesn't matter.  We handle the
  // latter two cases ("const V" or "V") specially (copy them into
  // userInvDiagCopy first, which is otherwise null at the end of the
  // long if-then chain) to avoid an extra copy.
  if (plist.isParameter ("chebyshev: operator inv diagonal")) {
    // Pointer to the user's Vector, if provided.
    RCP<const V> userInvDiag;

    try { // Could the type be const V*?
      const V* rawUserInvDiag =
        plist.get<const V*> ("chebyshev: operator inv diagonal");
      // Nonowning reference (we'll make a deep copy below)
      userInvDiag = rcp (rawUserInvDiag, false);
    } catch (Teuchos::Exceptions::InvalidParameterType&) {
    }
    if (userInvDiag.is_null ()) {
      try { // Could the type be V*?
        V* rawUserInvDiag = plist.get<V*> ("chebyshev: operator inv diagonal");
        // Nonowning reference (we'll make a deep copy below)
        userInvDiag = rcp (const_cast<const V*> (rawUserInvDiag), false);
      } catch (Teuchos::Exceptions::InvalidParameterType&) {
      }
    }
    if (userInvDiag.is_null ()) {
      try { // Could the type be RCP<const V>?
        userInvDiag =
          plist.get<RCP<const V> > ("chebyshev: operator inv diagonal");
      } catch (Teuchos::Exceptions::InvalidParameterType&) {
      }
    }
    if (userInvDiag.is_null ()) {
      try { // Could the type be RCP<V>?
        RCP<V> userInvDiagNonConst =
          plist.get<RCP<V> > ("chebyshev: operator inv diagonal");
        userInvDiag = rcp_const_cast<const V> (userInvDiagNonConst);
      } catch (Teuchos::Exceptions::InvalidParameterType&) {
      }
    }
    if (userInvDiag.is_null ()) {
#ifndef _MSC_VER
      try { // Could the type be const V?
        // ParameterList::get() returns by reference.  Thus, we don't
        // have to invoke Vector's copy constructor here.  It's good
        // practice not to make an RCP to this reference, even though
        // it should be valid as long as the ParameterList that holds
        // it is valid.  Thus, we make our deep copy here, rather than
        // waiting to do it below.
        const V& userInvDiagRef =
          plist.get<const V> ("chebyshev: operator inv diagonal");
        userInvDiagCopy = rcp (new V (userInvDiagRef, Teuchos::Copy));
        // Tell the if-chain below not to keep trying.
        userInvDiag = userInvDiagCopy;
      } catch (Teuchos::Exceptions::InvalidParameterType&) {
      }
#else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error,
      "Ifpack2::Chebyshev::setParameters: \"chebyshev: operator inv diagonal\" "
      "plist.get<const V> does not compile around return held == other_held "
      "in Teuchos::any in Visual Studio.  Can't fix it now, so throwing "
      "in case someone builds there.");
#endif
    }
    if (userInvDiag.is_null ()) {
#ifndef _MSC_VER
      try { // Could the type be V?
        // ParameterList::get() returns by reference.  Thus, we don't
        // have to invoke Vector's copy constructor here.  It's good
        // practice not to make an RCP to this reference, even though
        // it should be valid as long as the ParameterList that holds
        // it is valid.  Thus, we make our deep copy here, rather than
        // waiting to do it below.
        V& userInvDiagNonConstRef =
          plist.get<V> ("chebyshev: operator inv diagonal");
        const V& userInvDiagRef = const_cast<const V&> (userInvDiagNonConstRef);
        userInvDiagCopy = rcp (new V (userInvDiagRef, Teuchos::Copy));
        // Tell the if-chain below not to keep trying.
        userInvDiag = userInvDiagCopy;
      } catch (Teuchos::Exceptions::InvalidParameterType&) {
      }
#else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error,
      "Ifpack2::Chebyshev::setParameters: \"chebyshev: operator inv diagonal\" "
      "plist.get<V> does not compile around return held == other_held "
      "in Teuchos::any in Visual Studio.  Can't fix it now, so throwing "
      "in case someone builds there.");
#endif
    }

    // NOTE: If the user's parameter has some strange type that we
    // didn't test above, userInvDiag might still be null.  You may
    // want to add an error test for this condition.  Currently, we
    // just say in this case that the user didn't give us a Vector.

    // If we have userInvDiag but don't have a deep copy yet, make a
    // deep copy now.
    if (! userInvDiag.is_null () && userInvDiagCopy.is_null ()) {
      userInvDiagCopy = rcp (new V (*userInvDiag, Teuchos::Copy));
    }

    // NOTE: userInvDiag, if provided, is a row Map version of the
    // Vector.  We don't necessarily have a range Map yet.  compute()
    // would be the proper place to compute the range Map version of
    // userInvDiag.
  }

  // Don't fill in defaults for the max or min eigenvalue, because
  // this class uses the existence of those parameters to determine
  // whether it should do eigenanalysis.
  if (plist.isParameter ("chebyshev: max eigenvalue")) {
    if (plist.isType<double>("chebyshev: max eigenvalue"))
      lambdaMax = plist.get<double> ("chebyshev: max eigenvalue");
    else
      lambdaMax = plist.get<ST> ("chebyshev: max eigenvalue");
    TEUCHOS_TEST_FOR_EXCEPTION(
      STS::isnaninf (lambdaMax), std::invalid_argument,
      "Ifpack2::Chebyshev::setParameters: \"chebyshev: max eigenvalue\" "
      "parameter is NaN or Inf.  This parameter is optional, but if you "
      "choose to supply it, it must have a finite value.");
  }
  if (plist.isParameter ("chebyshev: min eigenvalue")) {
    if (plist.isType<double>("chebyshev: min eigenvalue"))
      lambdaMin = plist.get<double> ("chebyshev: min eigenvalue");
    else
      lambdaMin = plist.get<ST> ("chebyshev: min eigenvalue");
    TEUCHOS_TEST_FOR_EXCEPTION(
      STS::isnaninf (lambdaMin), std::invalid_argument,
      "Ifpack2::Chebyshev::setParameters: \"chebyshev: min eigenvalue\" "
      "parameter is NaN or Inf.  This parameter is optional, but if you "
      "choose to supply it, it must have a finite value.");
  }

  // Only fill in Ifpack2's name for the default parameter, not ML's.
  if (plist.isParameter ("smoother: Chebyshev alpha")) { // ML compatibility
    if (plist.isType<double>("smoother: Chebyshev alpha"))
      eigRatio = plist.get<double> ("smoother: Chebyshev alpha");
    else
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

  // See Github Issue #234.  This parameter may be either MT
  // (preferred) or double.  We check both.
  {
    const char paramName[] = "chebyshev: boost factor";

    if (plist.isParameter (paramName)) {
      if (plist.isType<MT> (paramName)) { // MT preferred
        boostFactor = plist.get<MT> (paramName);
      }
      else if (! std::is_same<double, MT>::value &&
               plist.isType<double> (paramName)) {
        const double dblBF = plist.get<double> (paramName);
        boostFactor = static_cast<MT> (dblBF);
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument,
           "Ifpack2::Chebyshev::setParameters: \"chebyshev: boost factor\""
           "parameter must have type magnitude_type (MT) or double.");
      }
    }
    else { // parameter not in the list
      // mfh 12 Aug 2016: To preserve current behavior (that fills in
      // any parameters not in the input ParameterList with their
      // default values), we call set() here.  I don't actually like
      // this behavior; I prefer the Belos model, where the input
      // ParameterList is a delta from current behavior.  However, I
      // don't want to break things.
      plist.set (paramName, defaultBoostFactor);
    }
    TEUCHOS_TEST_FOR_EXCEPTION
      (boostFactor < Teuchos::ScalarTraits<MT>::one (), std::invalid_argument,
       "Ifpack2::Chebyshev::setParameters: \"" << paramName << "\" parameter "
       "must be >= 1, but you supplied the value " << boostFactor << ".");
  }

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
  userInvDiag_ = userInvDiagCopy;
  userLambdaMax_ = lambdaMax;
  userLambdaMin_ = lambdaMin;
  userEigRatio_ = eigRatio;
  boostFactor_ = static_cast<MT> (boostFactor);
  minDiagVal_ = minDiagVal;
  numIters_ = numIters;
  eigMaxIters_ = eigMaxIters;
  zeroStartingSolution_ = zeroStartingSolution;
  assumeMatrixUnchanged_ = assumeMatrixUnchanged;
  textbookAlgorithm_ = textbookAlgorithm;
  computeMaxResNorm_ = computeMaxResNorm;
  debug_ = debug;

  if (debug_) {
    // Only print if myRank == 0.
    int myRank = -1;
    if (A_.is_null () || A_->getComm ().is_null ()) {
      // We don't have a communicator (yet), so we assume that
      // everybody can print.  Revise this expectation in setMatrix().
      myRank = 0;
    }
    else {
      myRank = A_->getComm ()->getRank ();
    }

    if (myRank == 0) {
      out_ = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
    }
    else {
      Teuchos::RCP<Teuchos::oblackholestream> blackHole (new Teuchos::oblackholestream ());
      out_ = Teuchos::getFancyOStream (blackHole); // print nothing on other processes
    }
  }
  else { // NOT debug
    // free the "old" output stream, if there was one
    out_ = Teuchos::null;
  }
}


template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::reset ()
{
  D_ = Teuchos::null;
  diagOffsets_ = offsets_type ();
  savedDiagOffsets_ = false;
  V_ = Teuchos::null;
  W_ = Teuchos::null;
  computedLambdaMax_ = STS::nan ();
  computedLambdaMin_ = STS::nan ();
}


template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A.getRawPtr () != A_.getRawPtr ()) {
    if (! assumeMatrixUnchanged_) {
      reset ();
    }
    A_ = A;

    // The communicator may have changed, or we may not have had a
    // communicator before.  Thus, we may have to reset the debug
    // output stream.
    if (debug_) {
      // Only print if myRank == 0.
      int myRank = -1;
      if (A.is_null () || A->getComm ().is_null ()) {
        // We don't have a communicator (yet), so we assume that
        // everybody can print.  Revise this expectation in setMatrix().
        myRank = 0;
      }
      else {
        myRank = A->getComm ()->getRank ();
      }

      if (myRank == 0) {
        out_ = Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr));
      }
      else {
        Teuchos::RCP<Teuchos::oblackholestream> blackHole (new Teuchos::oblackholestream ());
        out_ = Teuchos::getFancyOStream (blackHole); // print nothing on other processes
      }
    }
    else { // NOT debug
      // free the "old" output stream, if there was one
      out_ = Teuchos::null;
    }
  }
}


template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::compute ()
{
  using std::endl;
  // Some of the optimizations below only work if A_ is a
  // Tpetra::CrsMatrix.  We'll make our best guess about its type
  // here, since we have no way to get back the original fifth
  // template parameter.
  typedef Tpetra::CrsMatrix<typename MV::scalar_type,
    typename MV::local_ordinal_type,
    typename MV::global_ordinal_type,
    typename MV::node_type> crs_matrix_type;

  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::runtime_error, "Ifpack2::Chebyshev::compute: The input "
    "matrix A is null.  Please call setMatrix() with a nonnull input matrix "
    "before calling this method.");

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

  // makeInverseDiagonal() returns a range Map Vector.
  if (userInvDiag_.is_null ()) {
    Teuchos::RCP<const crs_matrix_type> A_crsMat =
      Teuchos::rcp_dynamic_cast<const crs_matrix_type> (A_);

    if (D_.is_null ()) { // We haven't computed D_ before
      if (! A_crsMat.is_null () && A_crsMat->isStaticGraph ()) {
        // It's a CrsMatrix with a const graph; cache diagonal offsets.
        const size_t lclNumRows = A_crsMat->getNodeNumRows ();
        if (diagOffsets_.extent (0) < lclNumRows) {
          diagOffsets_ = offsets_type (); // clear first to save memory
          diagOffsets_ = offsets_type ("offsets", lclNumRows);
        }
        A_crsMat->getCrsGraph ()->getLocalDiagOffsets (diagOffsets_);
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
          const size_t lclNumRows = A_crsMat->getNodeNumRows ();
          if (diagOffsets_.extent (0) < lclNumRows) {
            diagOffsets_ = offsets_type (); // clear first to save memory
            diagOffsets_ = offsets_type ("offsets", lclNumRows);
          }
          A_crsMat->getCrsGraph ()->getLocalDiagOffsets (diagOffsets_);
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
    D_ = makeRangeMapVectorConst (userInvDiag_);
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
    // FIXME (mfh 20 Nov 2013) Should scale this 1.0e-6 term
    // appropriately for MT's machine precision.
    if (STS::magnitude (lambdaMaxForApply_ - one) < Teuchos::as<MT> (1.0e-6)) {
      lambdaMinForApply_ = one;
      lambdaMaxForApply_ = lambdaMinForApply_;
      eigRatioForApply_ = one; // Ifpack doesn't include this line.
    }
  }
}


template<class ScalarType, class MV>
ScalarType
Chebyshev<ScalarType, MV>::
getLambdaMaxForApply() const {
  return lambdaMaxForApply_;
}


template<class ScalarType, class MV>
typename Chebyshev<ScalarType, MV>::MT
Chebyshev<ScalarType, MV>::apply (const MV& B, MV& X)
{
  const char prefix[] = "Ifpack2::Chebyshev::apply: ";

  if (debug_) {
    *out_ << "apply: " << std::endl;
  }
  TEUCHOS_TEST_FOR_EXCEPTION
    (A_.is_null (), std::runtime_error, prefix << "The input matrix A is null. "
     " Please call setMatrix() with a nonnull input matrix, and then call "
     "compute(), before calling this method.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (STS::isnaninf (lambdaMaxForApply_), std::runtime_error,
     prefix << "There is no estimate for the max eigenvalue."
     << std::endl << computeBeforeApplyReminder);
  TEUCHOS_TEST_FOR_EXCEPTION
    (STS::isnaninf (lambdaMinForApply_), std::runtime_error,
     prefix << "There is no estimate for the min eigenvalue."
     << std::endl << computeBeforeApplyReminder);
  TEUCHOS_TEST_FOR_EXCEPTION
    (STS::isnaninf (eigRatioForApply_), std::runtime_error,
     prefix <<"There is no estimate for the ratio of the max "
     "eigenvalue to the min eigenvalue."
     << std::endl << computeBeforeApplyReminder);
  TEUCHOS_TEST_FOR_EXCEPTION
    (D_.is_null (), std::runtime_error, prefix << "The vector of inverse "
     "diagonal entries of the matrix has not yet been computed."
     << std::endl << computeBeforeApplyReminder);

  if (textbookAlgorithm_) {
    textbookApplyImpl (*A_, B, X, numIters_, lambdaMaxForApply_,
                       lambdaMinForApply_, eigRatioForApply_, *D_);
  }
  else {
    ifpackApplyImpl (*A_, B, X, numIters_, lambdaMaxForApply_,
                     lambdaMinForApply_, eigRatioForApply_, *D_);
  }

  if (computeMaxResNorm_ && B.getNumVectors () > 0) {
    MV R (B.getMap (), B.getNumVectors ());
    computeResidual (R, B, *A_, X);
    Teuchos::Array<MT> norms (B.getNumVectors ());
    R.norm2 (norms ());
    return *std::max_element (norms.begin (), norms.end ());
  }
  else {
    return Teuchos::ScalarTraits<MT>::zero ();
  }
}

template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
print (std::ostream& out) {
  this->describe (* (Teuchos::getFancyOStream (Teuchos::rcpFromRef (out))),
                  Teuchos::VERB_MEDIUM);
}

template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
computeResidual (MV& R, const MV& B, const op_type& A, const MV& X,
                 const Teuchos::ETransp mode)
{
  Tpetra::deep_copy(R, B);
  A.apply (X, R, mode, -STS::one(), STS::one());
}

template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
solve (MV& Z, const V& D_inv, const MV& R) {
  Z.elementWiseMultiply (STS::one(), D_inv, R, STS::zero());
}

template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
solve (MV& Z, const ST alpha, const V& D_inv, const MV& R) {
  Z.elementWiseMultiply (alpha, D_inv, R, STS::zero());
}

template<class ScalarType, class MV>
Teuchos::RCP<const typename Chebyshev<ScalarType, MV>::V>
Chebyshev<ScalarType, MV>::
makeInverseDiagonal (const row_matrix_type& A, const bool useDiagOffsets) const
{
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
      A_crsMat->getLocalDiagCopy (*D_rowMap, diagOffsets_);
    }
  }
  else {
    // This always works for a Tpetra::RowMatrix, even if it is not a
    // Tpetra::CrsMatrix.  We just don't have offsets in this case.
    A.getLocalDiagCopy (*D_rowMap);
  }
  RCP<V> D_rangeMap = makeRangeMapVector (D_rowMap);

  if (debug_) {
    // In debug mode, make sure that all diagonal entries are
    // positive, on all processes.  Note that *out_ only prints on
    // Process 0 of the matrix's communicator.
    D_rangeMap->template sync<Kokkos::HostSpace> ();
    auto D_lcl = D_rangeMap->template getLocalView<Kokkos::HostSpace> ();
    auto D_lcl_1d = Kokkos::subview (D_lcl, Kokkos::ALL (), 0);

    typedef typename MV::impl_scalar_type IST;
    typedef typename MV::local_ordinal_type LO;
    typedef Kokkos::Details::ArithTraits<IST> STS;
    typedef Kokkos::Details::ArithTraits<typename STS::mag_type> STM;

    const LO lclNumRows = static_cast<LO> (D_rangeMap->getLocalLength ());
    bool foundNonpositiveValue = false;
    for (LO i = 0; i < lclNumRows; ++i) {
      if (STS::real (D_lcl_1d(i)) <= STM::zero ()) {
        foundNonpositiveValue = true;
        break;
      }
    }

    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

    const int lclSuccess = foundNonpositiveValue ? 0 : 1;
    int gblSuccess = lclSuccess; // to be overwritten
    if (! D_rangeMap->getMap ().is_null () && D_rangeMap->getMap ()->getComm ().is_null ()) {
      const Teuchos::Comm<int>& comm = * (D_rangeMap->getMap ()->getComm ());
      reduceAll<int, int> (comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    }
    if (gblSuccess == 1) {
      *out_ << "makeInverseDiagonal: The matrix's diagonal entries all have "
        "positive real part (this is good for Chebyshev)." << std::endl;
    }
    else {
      *out_ << "makeInverseDiagonal: The matrix's diagonal has at least one "
        "entry with nonpositive real part, on at least one process in the "
        "matrix's communicator.  This is bad for Chebyshev." << std::endl;
    }
  }

  // Invert the diagonal entries, replacing entries less (in
  // magnitude) than the user-specified value with that value.
  reciprocal_threshold (*D_rangeMap, minDiagVal_);
  return Teuchos::rcp_const_cast<const V> (D_rangeMap);
}


template<class ScalarType, class MV>
Teuchos::RCP<const typename Chebyshev<ScalarType, MV>::V>
Chebyshev<ScalarType, MV>::
makeRangeMapVectorConst (const Teuchos::RCP<const V>& D) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef Tpetra::Export<typename MV::local_ordinal_type,
                         typename MV::global_ordinal_type,
                         typename MV::node_type> export_type;
  // This throws logic_error instead of runtime_error, because the
  // methods that call makeRangeMapVector should all have checked
  // whether A_ is null before calling this method.
  TEUCHOS_TEST_FOR_EXCEPTION(
    A_.is_null (), std::logic_error, "Ifpack2::Details::Chebyshev::"
    "makeRangeMapVector: The input matrix A is null.  Please call setMatrix() "
    "with a nonnull input matrix before calling this method.  This is probably "
    "a bug in Ifpack2; please report this bug to the Ifpack2 developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    D.is_null (), std::logic_error, "Ifpack2::Details::Chebyshev::"
    "makeRangeMapVector: The input Vector D is null.  This is probably "
    "a bug in Ifpack2; please report this bug to the Ifpack2 developers.");

  RCP<const map_type> sourceMap = D->getMap ();
  RCP<const map_type> rangeMap = A_->getRangeMap ();
  RCP<const map_type> rowMap = A_->getRowMap ();

  if (rangeMap->isSameAs (*sourceMap)) {
    // The given vector's Map is the same as the matrix's range Map.
    // That means we don't need to Export.  This is the fast path.
    return D;
  }
  else { // We need to Export.
    RCP<const export_type> exporter;
    // Making an Export object from scratch is expensive enough that
    // it's worth the O(1) global reductions to call isSameAs(), to
    // see if we can avoid that cost.
    if (sourceMap->isSameAs (*rowMap)) {
      // We can reuse the matrix's Export object, if there is one.
      exporter = A_->getGraph ()->getExporter ();
    }
    else { // We have to make a new Export object.
      exporter = rcp (new export_type (sourceMap, rangeMap));
    }

    if (exporter.is_null ()) {
      return D; // Row Map and range Map are the same; no need to Export.
    }
    else { // Row Map and range Map are _not_ the same; must Export.
      RCP<V> D_out = rcp (new V (*D, Teuchos::Copy));
      D_out->doExport (*D, *exporter, Tpetra::ADD);
      return Teuchos::rcp_const_cast<const V> (D_out);
    }
  }
}


template<class ScalarType, class MV>
Teuchos::RCP<typename Chebyshev<ScalarType, MV>::V>
Chebyshev<ScalarType, MV>::
makeRangeMapVector (const Teuchos::RCP<V>& D) const
{
  using Teuchos::rcp_const_cast;
  return rcp_const_cast<V> (makeRangeMapVectorConst (rcp_const_cast<V> (D)));
}


template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
textbookApplyImpl (const op_type& A,
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

template<class ScalarType, class MV>
typename Chebyshev<ScalarType, MV>::MT
Chebyshev<ScalarType, MV>::maxNormInf (const MV& X) {
  Teuchos::Array<MT> norms (X.getNumVectors ());
  X.normInf (norms());
  return *std::max_element (norms.begin (), norms.end ());
}

template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
ifpackApplyImpl (const op_type& A,
                 const MV& B,
                 MV& X,
                 const int numIters,
                 const ST lambdaMax,
                 const ST lambdaMin,
                 const ST eigRatio,
                 const V& D_inv)
{
  using std::endl;
#ifdef HAVE_IFPACK2_DEBUG
  const bool debug = debug_;
#else
  const bool debug = false;
#endif

  if (debug) {
    *out_ << " \\|B\\|_{\\infty} = " << maxNormInf (B) << endl;
    *out_ << " \\|X\\|_{\\infty} = " << maxNormInf (X) << endl;
  }

  if (numIters <= 0) {
    return;
  }
  const ST one = static_cast<ST> (1.0);
  const ST two = static_cast<ST> (2.0);

  // Quick solve when the matrix A is the identity.
  if (lambdaMin == one && lambdaMax == lambdaMin) {
    solve (X, D_inv, B);
    return;
  }

  // Initialize coefficients
  const ST alpha = lambdaMax / eigRatio;
  const ST beta = boostFactor_ * lambdaMax;
  const ST delta = two / (beta - alpha);
  const ST theta = (beta + alpha) / two;
  const ST s1 = theta * delta;

  if (debug) {
    *out_ << " alpha = " << alpha << endl
          << " beta = " << beta << endl
          << " delta = " << delta << endl
          << " theta = " << theta << endl
          << " s1 = " << s1 << endl;
  }

  // Fetch cached temporary vectors.
  Teuchos::RCP<MV> V_ptr, W_ptr;
  makeTempMultiVectors (V_ptr, W_ptr, B);

  // mfh 28 Jan 2013: We write V1 instead of V, so as not to confuse
  // the multivector V with the typedef V (for Tpetra::Vector).
  //MV V1 (B.getMap (), B.getNumVectors (), false);
  //MV W (B.getMap (), B.getNumVectors (), false);
  MV& V1 = *V_ptr;
  MV& W = *W_ptr;

  if (debug) {
    *out_ << " Iteration " << 1 << ":" << endl
          << " - \\|D\\|_{\\infty} = " << D_->normInf () << endl;
  }

  // Special case for the first iteration.
  if (! zeroStartingSolution_) {
    computeResidual (V1, B, A, X); // V1 = B - A*X
    if (debug) {
      *out_ << " - \\|B - A*X\\|_{\\infty} = " << maxNormInf (V1) << endl;
    }

    solve (W, one/theta, D_inv, V1); // W = (1/theta)*D_inv*(B-A*X)
    if (debug) {
      *out_ << " - \\|W\\|_{\\infty} = " << maxNormInf (W) << endl;
    }
    X.update (one, W, one); // X = X + W
  }
  else {
    solve (W, one/theta, D_inv, B); // W = (1/theta)*D_inv*B
    if (debug) {
      *out_ << " - \\|W\\|_{\\infty} = " << maxNormInf (W) << endl;
    }
    Tpetra::deep_copy (X, W); // X = 0 + W
  }
  if (debug) {
    *out_ << " - \\|X\\|_{\\infty} = " << maxNormInf (X) << endl;
  }

  // The rest of the iterations.
  ST rhok = one / s1;
  ST rhokp1, dtemp1, dtemp2;
  for (int deg = 1; deg < numIters; ++deg) {
    if (debug) {
      *out_ << " Iteration " << deg+1 << ":" << endl
            << " - \\|D\\|_{\\infty} = " << D_->normInf () << endl
            << " - \\|B\\|_{\\infty} = " << maxNormInf (B) << endl
            << " - \\|A\\|_{\\text{frob}} = " << A_->getFrobeniusNorm () << endl
            << " - rhok = " << rhok << endl;
      V1.putScalar (STS::zero ()); // ???????
    }

    computeResidual (V1, B, A, X); // V1 = B - A*X
    if (debug) {
      *out_ << " - \\|B - A*X\\|_{\\infty} = " << maxNormInf (V1) << endl;
    }

    rhokp1 = one / (two * s1 - rhok);
    dtemp1 = rhokp1 * rhok;
    dtemp2 = two * rhokp1 * delta;
    rhok = rhokp1;

    if (debug) {
      *out_ << " - dtemp1 = " << dtemp1 << endl
            << " - dtemp2 = " << dtemp2 << endl;
    }


    W.elementWiseMultiply (dtemp2, D_inv, V1, dtemp1);
    X.update (one, W, one);

    if (debug) {
      *out_ << " - \\|W\\|_{\\infty} = " << maxNormInf (W) << endl;
      *out_ << " - \\|X\\|_{\\infty} = " << maxNormInf (X) << endl;
    }
  }
}

template<class ScalarType, class MV>
typename Chebyshev<ScalarType, MV>::ST
Chebyshev<ScalarType, MV>::
powerMethodWithInitGuess (const op_type& A,
                          const V& D_inv,
                          const int numIters,
                          V& x)
{
  using std::endl;
  if (debug_) {
    *out_ << " powerMethodWithInitGuess:" << endl;
  }

  const ST zero = static_cast<ST> (0.0);
  const ST one = static_cast<ST> (1.0);
  ST lambdaMax = zero;
  ST RQ_top, RQ_bottom, norm;

  V y (A.getRangeMap ());
  norm = x.norm2 ();
  TEUCHOS_TEST_FOR_EXCEPTION
    (norm == zero, std::runtime_error,
    "Ifpack2::Chebyshev::powerMethodWithInitGuess: "
    "The initial guess has zero norm.  This could be either because Tpetra::"
    "Vector::randomize() filled the vector with zeros (if that was used to "
    "compute the initial guess), or because the norm2 method has a bug.  The "
    "first is not impossible, but unlikely.");

  if (debug_) {
    *out_ << "  Original norm1(x): " << x.norm1 () << ", norm2(x): " << norm << endl;
  }

  x.scale (one / norm);

  if (debug_) {
    *out_ << "  norm1(x.scale(one/norm)): " << x.norm1 () << endl;
  }

  for (int iter = 0; iter < numIters; ++iter) {
    if (debug_) {
      *out_ << "  Iteration " << iter << endl;
    }
    A.apply (x, y);
    solve (y, D_inv, y);
    RQ_top = y.dot (x);
    RQ_bottom = x.dot (x);
    if (debug_) {
      *out_ << "   RQ_top: " << RQ_top << ", RQ_bottom: " << RQ_bottom << endl;
    }
    lambdaMax = RQ_top / RQ_bottom;
    norm = y.norm2 ();
    if (norm == zero) { // Return something reasonable.
      if (debug_) {
        *out_ << "   norm is zero; returning zero" << endl;
      }
      return zero;
    }
    x.update (one / norm, y, zero);
  }
  if (debug_) {
    *out_ << "  lambdaMax: " << lambdaMax << endl;
  }
  return lambdaMax;
}

template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
computeInitialGuessForPowerMethod (V& x, const bool nonnegativeRealParts) const
{
  typedef typename MV::device_type::execution_space dev_execution_space;
  typedef typename MV::device_type::memory_space dev_memory_space;
  typedef typename MV::local_ordinal_type LO;

  x.randomize ();

  if (nonnegativeRealParts) {
    x.template modify<dev_memory_space> ();
    auto x_lcl = x.template getLocalView<dev_memory_space> ();
    auto x_lcl_1d = Kokkos::subview (x_lcl, Kokkos::ALL (), 0);

    const LO lclNumRows = static_cast<LO> (x.getLocalLength ());
    Kokkos::RangePolicy<dev_execution_space, LO> range (0, lclNumRows);
    PositivizeVector<decltype (x_lcl_1d), LO> functor (x_lcl_1d);
    Kokkos::parallel_for (range, functor);
  }
}

template<class ScalarType, class MV>
typename Chebyshev<ScalarType, MV>::ST
Chebyshev<ScalarType, MV>::
powerMethod (const op_type& A, const V& D_inv, const int numIters)
{
  using std::endl;
  if (debug_) {
    *out_ << "powerMethod:" << endl;
  }

  const ST zero = static_cast<ST> (0.0);
  V x (A.getDomainMap ());
  // For the first pass, just let the pseudorandom number generator
  // fill x with whatever values it wants; don't try to make its
  // entries nonnegative.
  computeInitialGuessForPowerMethod (x, false);

  ST lambdaMax = powerMethodWithInitGuess (A, D_inv, numIters, x);

  // mfh 07 Jan 2015: Taking the real part here is only a concession
  // to the compiler, so that this class can build with ScalarType =
  // std::complex<T>.  Our Chebyshev implementation only works with
  // real, symmetric positive definite matrices.  The right thing to
  // do would be what Belos does, which is provide a partial
  // specialization for ScalarType = std::complex<T> with a stub
  // implementation (that builds, but whose constructor throws).
  if (STS::real (lambdaMax) < STS::real (zero)) {
    if (debug_) {
      *out_ << "real(lambdaMax) = " << STS::real (lambdaMax) << " < 0; try "
        "again with a different random initial guess" << endl;
    }
    // Max eigenvalue estimate was negative.  Perhaps we got unlucky
    // with the random initial guess.  Try again with a different (but
    // still random) initial guess.  Only try again once, so that the
    // run time is bounded.

    // For the second pass, make all the entries of the initial guess
    // vector have nonnegative real parts.
    computeInitialGuessForPowerMethod (x, true);
    lambdaMax = powerMethodWithInitGuess (A, D_inv, numIters, x);
  }
  return lambdaMax;
}

template<class ScalarType, class MV>
Teuchos::RCP<const typename Chebyshev<ScalarType, MV>::row_matrix_type>
Chebyshev<ScalarType, MV>::getMatrix () const {
  return A_;
}

template<class ScalarType, class MV>
bool
Chebyshev<ScalarType, MV>::
hasTransposeApply () const {
  // Technically, this is true, because the matrix must be symmetric.
  return true;
}

template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
makeTempMultiVectors (Teuchos::RCP<MV>& V1,
                      Teuchos::RCP<MV>& W,
                      const MV& X)
{
  // ETP 02/08/17:  We must check not only if the temporary vectors are
  // null, but also if the number of columns match, since some multi-RHS
  // solvers (e.g., Belos) may call apply() with different numbers of columns.
  if (V_.is_null () || V_->getNumVectors () != X.getNumVectors ()) {
    V_ = Teuchos::rcp (new MV (X.getMap (), X.getNumVectors (), false));
  }
  //W must be initialized to zero when it is used as a multigrid smoother.
  if (W_.is_null () || W_->getNumVectors () != X.getNumVectors ()) {
    W_ = Teuchos::rcp (new MV (X.getMap (), X.getNumVectors (), true));
  }
  V1 = V_;
  W = W_;
}

template<class ScalarType, class MV>
std::string
Chebyshev<ScalarType, MV>::
description () const {
  std::ostringstream oss;
  // YAML requires quoting the key in this case, to distinguish
  // key's colons from the colon that separates key from value.
  oss << "\"Ifpack2::Details::Chebyshev\":"
      << "{"
      << "degree: " << numIters_
      << ", lambdaMax: " << lambdaMaxForApply_
      << ", alpha: " << eigRatioForApply_
      << ", lambdaMin: " << lambdaMinForApply_
      << ", boost factor: " << boostFactor_
      << "}";
  return oss.str();
}

template<class ScalarType, class MV>
void
Chebyshev<ScalarType, MV>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::TypeNameTraits;
  using std::endl;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;
  if (vl == Teuchos::VERB_NONE) {
    return; // print NOTHING
  }

  // By convention, describe() starts with a tab.
  //
  // This does affect all processes on which it's valid to print to
  // 'out'.  However, it does not actually print spaces to 'out'
  // unless operator<< gets called, so it's safe to use on all
  // processes.
  Teuchos::OSTab tab0 (out);

  // We only print on Process 0 of the matrix's communicator.  If
  // the matrix isn't set, we don't have a communicator, so we have
  // to assume that every process can print.
  int myRank = -1;
  if (A_.is_null () || A_->getComm ().is_null ()) {
    myRank = 0;
  }
  else {
    myRank = A_->getComm ()->getRank ();
  }
  if (myRank == 0) {
    // YAML requires quoting the key in this case, to distinguish
    // key's colons from the colon that separates key from value.
    out << "\"Ifpack2::Details::Chebyshev\":" << endl;
  }
  Teuchos::OSTab tab1 (out);

  if (vl == Teuchos::VERB_LOW) {
    if (myRank == 0) {
      out << "degree: " << numIters_ << endl
          << "lambdaMax: " << lambdaMaxForApply_ << endl
          << "alpha: " << eigRatioForApply_ << endl
          << "lambdaMin: " << lambdaMinForApply_ << endl
          << "boost factor: " << boostFactor_ << endl;
    }
    return;
  }

  // vl > Teuchos::VERB_LOW

  if (myRank == 0) {
    out << "Template parameters:" << endl;
    {
      Teuchos::OSTab tab2 (out);
      out << "ScalarType: " << TypeNameTraits<ScalarType>::name () << endl
          << "MV: " << TypeNameTraits<MV>::name () << endl;
    }

    // "Computed parameters" literally means "parameters whose
    // values were computed by compute()."
    if (myRank == 0) {
      out << "Computed parameters:" << endl;
    }
  }

  // Print computed parameters
  {
    Teuchos::OSTab tab2 (out);
    // Users might want to see the values in the computed inverse
    // diagonal, so we print them out at the highest verbosity.
    if (myRank == 0) {
      out << "D_: ";
    }
    if (D_.is_null ()) {
      if (myRank == 0) {
        out << "unset" << endl;
      }
    }
    else if (vl <= Teuchos::VERB_HIGH) {
      if (myRank == 0) {
        out << "set" << endl;
      }
    }
    else { // D_ not null and vl > Teuchos::VERB_HIGH
      if (myRank == 0) {
        out << endl;
      }
      // By convention, describe() first indents, then prints.
      // We can rely on other describe() implementations to do that.
      D_->describe (out, vl);
    }
    if (myRank == 0) {
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
  } // print computed parameters

  if (myRank == 0) {
    out << "User parameters:" << endl;
  }

  // Print user parameters
  {
    Teuchos::OSTab tab2 (out);
    out << "userInvDiag_: ";
    if (userInvDiag_.is_null ()) {
      out << "unset" << endl;
    }
    else if (vl <= Teuchos::VERB_HIGH) {
      out << "set" << endl;
    }
    else { // userInvDiag_ not null and vl > Teuchos::VERB_HIGH
      if (myRank == 0) {
        out << endl;
      }
      userInvDiag_->describe (out, vl);
    }
    if (myRank == 0) {
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
  } // print user parameters
}

} // namespace Details
} // namespace Ifpack2

#define IFPACK2_DETAILS_CHEBYSHEV_INSTANT(S,LO,GO,N) \
  template class Ifpack2::Details::Chebyshev< S, Tpetra::MultiVector<S, LO, GO, N> >;

#endif // IFPACK2_DETAILS_CHEBYSHEV_DEF_HPP
