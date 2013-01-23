// ***********************************************************************
// 
//      Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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


/*! \file Ifpack2_UnitTestChebyshev2.cpp
\brief A convergence test for Ifpack2::Chebyshev.
\author Mark Hoemmen

This test compares Ifpack2's implementation of Chebyshev iteration
against the following implementations:
1. A textbook version of the algorithm
2. A direct imitation of Ifpack's implementation
3. A direct imitation of ML's implementation
4. A textbook implementation of CG

All three do the same left diagonal scaling that Ifpack2 does.
"Textbook" implementations (#1 and #4) come from "Templates for the
Solution of Linear Systems," 2nd edition.  #2 imitates
Ifpack::Chebyshev, both in how it sets parameters and in the actual
iteration (ApplyInverse()).  #3 imitates ML_Cheby
(packages/ml/src/Smoother/ml_smoother.c).  We include CG just to give
us a rough measure of how fast the methods "should" converge.

The test exercises all four algorithms with a 1-D Poisson equation.
We know the eigenvalues of the matrix exactly as a function of its
dimensions, so we can give perfect eigenvalue bounds.  We also
experiment with changing the max eigenvalue bound so that it's either
an overestimate or an underestimate.

This test has the following command-line arguments:
- numIters: The number of iterations of Chebyshev or CG.
- localNumRows: The number of rows of the matrix on each process.

Currently, Ifpack2's Chebyshev and the Ifpack imitation <i>almost</i>
agree, but not quite.  The difference is almost always less than an
order of magnitude.  The textbook implementation of Chebyshev
converges much faster if the eigenvalue bounds are good, but both
Ifpack2's Chebyshev and its imitation are much less sensitive to an
incorrect upper bound on the eigenvalues.  This gives me confidence
that Ifpack2's version is correct.

I haven't yet done the analysis to understand why the different
implementations produce different results.  However, I feel better
now about Ifpack2's version.  I would still prefer to use the textbook
version of the algorithm, though, just to remove the puzzle of whether 
what you see on the screen is actually what it claims to be.
*/

#include <Ifpack2_ConfigDefs.hpp>
#include <Ifpack2_Chebyshev.hpp>
#include <Ifpack2_UnitTestHelpers.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_UnitTestRepository.hpp>
#include <cmath>

namespace {
/// \class Chebyshev
/// \brief Left-scaled Chebyshev iteration.
/// \tparam ScalarType The type of entries in the matrix and vectors.
/// \tparam MV Specialization of Tpetra::MultiVector.
/// \tparam MAT Corresponding specialization of Tpetra::CrsMatrix.
///
/// This class implements the version of Chebyshev iteration in
/// "Templates for the Solution of Linear Systems," 2nd edition.  It
/// uses the diagonal of the matrix to precondition the linear system
/// on the left.  Diagonal entries less than machine precision are
/// replaced with machine precision.
///
/// We require that the matrix A be real valued and symmetric positive
/// definite.  If users could provide the ellipse parameters ("d" and
/// "c" in the literature, where d is the real-valued center of the
/// ellipse, and d-c and d+c the two foci), the iteration itself would
/// work fine with nonsymmetric A, as long as the eigenvalues of A can
/// be bounded in an ellipse that is entirely to the right of the
/// origin.
///
/// This class implements several variants of Chebyshev iteration:
/// 1. A textbook version of the algorithm
/// 2. A direct imitation of Ifpack's implementation
/// 3. A direct imitation of ML's implementation
///
/// The "textbook" is "Templates for the Solution of Linear Systems,"
/// 2nd edition.  #2 imitates Ifpack::Chebyshev, both in how it sets
/// parameters and in the actual iteration (ApplyInverse()).  #3
/// imitates ML_Cheby (packages/ml/src/Smoother/ml_smoother.c).
template<class ScalarType, class MV, class MAT>
class Chebyshev {
public:
  typedef ScalarType ST;
  typedef Teuchos::ScalarTraits<ScalarType> STS;
  typedef typename STS::magnitudeType MT;
  typedef Tpetra::Operator<typename MV::scalar_type,
			   typename MV::local_ordinal_type,
			   typename MV::global_ordinal_type,
			   typename MV::node_type> OP;
  typedef Tpetra::Vector<typename MV::scalar_type,
			 typename MV::local_ordinal_type,
			 typename MV::global_ordinal_type,
			 typename MV::node_type> V;
  /// Constructor that takes a MAT and sets default parameters.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  ///   A must be real-valued and symmetric positive definite.
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
    imitateML_ (false),
    imitateIfpack_ (false)
  {}

  /// Constructor that sets the user's parameters.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  ///   A must be real-valued and symmetric positive definite.
  /// \param params [in/out] On input: the parameters.  On output:
  ///   filled with the current parameter settings.
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
    imitateML_ (false),
    imitateIfpack_ (false)
  {
    setParameters (params);
  }

  /// \brief Set (or reset) parameters.
  ///
  /// This method fills in the input ParameterList with missing
  /// parameters set to their default values.  You may call this
  /// method as many times as you want.  On each call, the input
  /// ParameterList is treated as a complete list of the desired
  /// parameters, not as a "delta" or change list from the current set
  /// of parameters.  (That is, if you remove parameters from the list
  /// that were there in the last call to setParameters() and call
  /// setParameters() again with the revised list, this method will
  /// use default values for the removed parameters, rather than
  /// letting the current settings remain.)  However, since the method
  /// fills in missing parameters, you may keep calling it with the
  /// ParameterList used in the previous call in order to get the same
  /// behavior as before.
  ///
  /// Parameters that govern spectral bounds of the matrix:
  /// - "chebyshev: max eigenvalue" (\c ScalarType): lambdaMax, an
  ///   upper bound of the bounding ellipse of the eigenvalues of the
  ///   matrix A.  If you do not set this parameter, we will compute
  ///   an approximation.  See "Parameters that govern eigenvalue
  ///   analysis" to control this approximation process.
  /// - "chebyshev: ratio eigenvalue" (\c ScalarType): eigRatio, the
  ///   ratio of lambdaMax to the lower bound of the bounding ellipse
  ///   of the eigenvalues of A.  We use lambdaMax and eigRatio to
  ///   determine the Chebyshev iteration coefficients.  This
  ///   parameter is optional and defaults to 30.
  /// - "chebyshev: min eigenvalue" (\c ScalarType): lambdaMin, a
  ///   lower bound of real part of bounding ellipse of eigenvalues of
  ///   the matrix A.  This parameter is optional and only used for
  ///   sanity checks.
  ///
  /// Parameters that govern the number of Chebyshev iterations:
  /// - "chebyshev: degree" (\c int): numIters, the number of
  ///   iterations.  This overrides "relaxation: sweeps" and
  ///   "smoother: sweeps" (see below).
  /// - "relaxation: sweeps" (\c int): numIters, the number of
  ///   iterations.  We include this for compatibility with Ifpack.
  ///   This overrides "smoother: sweeps" (see below).
  /// - "smoother: sweeps" (\c int): numIters, as above.
  ///   We include this for compatibility with ML.
  ///
  /// Parameters that govern eigenvalue analysis:
  /// - "chebyshev: eigenvalue max iterations" (\c int): eigMaxIters,
  ///   the number of power method iterations used to compute the
  ///   maximum eigenvalue.  This overrides "eigen-analysis:
  ///   iterations" (see below).
  /// - "eigen-analysis: iterations" (\c int): eigMaxIters, as above.
  ///   We include this parameter for compatibility with ML.
  /// - "eigen-analysis: type" (<tt>std::string</tt>): The algorithm
  ///   to use for estimating the max eigenvalue.  This parameter is
  ///   optional.  Currently, we only support "power-method" (or
  ///   "power method"), which is what Ifpack::Chebyshev uses for
  ///   eigenanalysis.  We include this parameter for compatibility
  ///   with ML.
  ///
  /// Parameters that govern backwards compatibility:
  /// - "chebyshev: imitate ML" (\c bool): If true, imitate
  ///   ML's implementation of Chebyshev iteration (ML_Cheby).
  ///   This overwrites the "chebyshev: imitate Ifpack" parameter.
  /// - "chebyshev: imitate Ifpack" (\c bool): If true, imitate
  ///   Ifpack's implementation of Chebyshev iteration.
  ///
  /// \pre lambdaMin, lambdaMax, and eigRatio are real
  /// \pre 0 < lambdaMin <= lambdaMax
  /// \pre numIters >= 0
  /// \pre eigMaxIters >= 0
  ///
  /// Default settings for parameters relating to spectral bounds come
  /// from Ifpack.
  void setParameters (Teuchos::ParameterList& plist) {
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
    const bool defaultImitateML = false;
    const bool defaultImitateIfpack = false;

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
    bool imitateML = defaultImitateML;
    bool imitateIfpack = defaultImitateIfpack;

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
    imitateML = plist.get ("chebyshev: imitate ML", imitateML);
    imitateIfpack = plist.get ("chebyshev: imitate Ifpack", imitateIfpack);

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
    imitateML_ = imitateML;
    imitateIfpack_ = imitateIfpack;
  }

  /// \brief (Re)compute the left scaling, and (if applicable)
  ///   estimate max and min eigenvalues of D_inv * A.
  ///
  /// You must call this method before calling apply(),
  /// - if you have not yet called this method,
  /// - if the matrix (either its values or its structure) has changed, or
  /// - any time after you call setParameters().
  ///
  /// Advanced users may omit calling compute() after calling
  /// setParameters(), as long as none of the changed parameters
  /// affect either computation of the inverse diagonal, or estimation
  /// of the max or min eigenvalues.
  ///
  /// If estimation of the eigenvalues is required, this method may
  /// take as long as several Chebyshev iterations.
  void compute () {
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

  /// Solve Ax=b for x with Chebyshev iteration with left diagonal scaling.
  ///
  /// \param B [in] Right-hand side(s) in the linear system to solve.
  /// \param X [in] Initial guess(es) for the linear system to solve.
  ///
  /// \return Max (over all columns) absolute residual 2-norm after iterating.
  MT apply (const MV& B, MV& X) {
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
      "This means one of the following: (a) you didn't call compute() before calling apply(), "
      "(b) you didn't call compute() after calling setParameters(), or "
      "(c) there is a bug in this class (compute() has not done the eigenanalysis that it should have done).");
    if (STS::isnaninf (lambdaMin)) {
      lambdaMin = computedLambdaMin_;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      STS::isnaninf (lambdaMin), 
      std::runtime_error, 
      "Chebyshev::apply: Both lambdaMin_ and computedLambdaMin_ are NaN or Inf.  "
      "This means one of the following: (a) you didn't call compute() before calling apply(), "
      "(b) you didn't call compute() after calling setParameters(), or "
      "(c) there is a bug in this class (compute() has not done the eigenanalysis that it should have done).");
    TEUCHOS_TEST_FOR_EXCEPTION(
      STS::isnaninf (eigRatio), 
      std::logic_error, 
      "Chebyshev::apply: eigRatio is NaN or Inf.");

    if (imitateML_) {
      mlApplyImpl (*A_, B, X, numIters_, lambdaMax, lambdaMin, eigRatio, *D_);
    } else if (imitateIfpack_) {
      const ST one = Teuchos::as<ST> (1);
      if (STS::magnitude (lambdaMax - one) < Teuchos::as<ST> (1.0e-6)) {
	lambdaMin = one;
	lambdaMax = lambdaMin;
	eigRatio = one; // Ifpack doesn't include this line.
      }
      ifpackApplyImpl (*A_, B, X, numIters_, lambdaMax, lambdaMin, eigRatio, *D_);
    } else {
      // mfh 17,18 Jan 2013: For convergence in 50 iterations with the
      // 1-D Poisson equation, when you use the exact min and max
      // eigenvalues, the best thing is to use the actual max
      // eigenvalue (not to multiply it by 1.1, as in ML), but to use
      // lambdaMax/30 for the "min eigenvalue," instead of the actual
      // min eigenvalue.
      applyImpl (*A_, B, X, numIters_, lambdaMax, lambdaMin, eigRatio, *D_);
    }

    MV R (B.getMap (), B.getNumVectors ());
    computeResidual (R, B, *A_, X);
    Teuchos::Array<MT> norms (B.getNumVectors ());
    R.norm2 (norms ());
    return *std::max_element (norms.begin (), norms.end ());
  }

  //! Print instance data to the given output stream.
  void print (std::ostream& out) {
    using std::endl;
    out << "Chebyshev:" << endl
	<< "- computedLambdaMax_ = " << computedLambdaMax_ << endl
	<< "- computedLambdaMin_ = " << computedLambdaMin_ << endl
	<< "- lambdaMax_ = " << lambdaMax_ << endl
	<< "- lambdaMin_ = " << lambdaMin_ << endl
	<< "- eigRatio_ = " << eigRatio_ << endl
	<< "- numIters_ = " << numIters_ << endl
	<< "- eigMaxIters_ = " << eigMaxIters_ << endl
	<< "- imitateML_ = " << imitateML_ << endl
	<< "- imitateIfpack_ = " << imitateIfpack_ << endl;
  }

private:
  Teuchos::RCP<const MAT> A_; //!< The sparse matrix A.
  Teuchos::RCP<const V> D_; //!< The inverse of the diagonal entries of A.

  // Computed values

  /// Estimate that we compute for maximum eigenvalue of A.
  /// compute() will always recompute this. 
  /// This is set to NaN if it hasn't been computed yet.
  ST computedLambdaMax_;
  /// Estimate that we compute for minimum eigenvalue of A.
  /// compute() will always recompute this.
  /// This is set to NaN if it hasn't been computed yet.
  ST computedLambdaMin_;

  // Parameters

  /// User-provided estimate for maximum eigenvalue of A.
  /// This is set to NaN if the user did not provide this.
  ST lambdaMax_; 
  /// User-provided estimate for minimum eigenvalue of A.
  /// This is set to NaN if the user did not provide this.
  ST lambdaMin_; 
  /// Estimate for ratio of max to max eigenvalue of A.
  /// Not necessarily equal to lambdaMax_/lambdaMin_.
  ST eigRatio_;  
  //! Number of Chebyshev iterations to run on each call to apply().
  int numIters_;
  //! Number of iterations of the power method for estimating lambdaMax_.
  int eigMaxIters_;
  //! Whether to assume that the X input to apply() is always zero.
  bool zeroStartingSolution_;
  /// Whether to imitate ML's Chebyshev implementation.
  /// This overrides imitateIfpack_, and enables the power method.
  bool imitateML_;
  //! Whether to imitate Ifpack's Chebyshev implementation.
  bool imitateIfpack_;

  //! R = B - Op(A) * X, where Op(A) is either A, \f$A^T\f$, or \f$A^H\f$.
  static void 
  computeResidual (MV& R, const MV& B, const MAT& A, const MV& X, 
		   const Teuchos::ETransp mode = Teuchos::NO_TRANS) 
  {
    R = B;
    A.apply (X, R, mode, -STS::one(), STS::one());
  }

  /// \brief Z = D_inv * R, = D \ R.
  ///
  /// \param D_inv [in] A vector representing a diagonal matrix.
  /// \param R [in] Input multivector.
  /// \param Z [out] Result of multiplying the diagonal matrix D_inv with R.
  static void solve (MV& Z, const V& D_inv, const MV& R) {
    Z.elementWiseMultiply (STS::one(), D_inv, R, STS::zero());
  }

  /// \brief Z = alpha * D_inv * R, = alpha * (D \ R).
  ///
  /// \param D_inv [in] A vector representing a diagonal matrix.
  /// \param R [in] Input multivector.
  /// \param Z [out] Result of multiplying the diagonal matrix D_inv with R.
  static void solve (MV& Z, const ST alpha, const V& D_inv, const MV& R) {
    Z.elementWiseMultiply (alpha, D_inv, R, STS::zero());
  }

  //! Get a copy of the diagonal of the matrix, as a row Map vector.
  static Teuchos::RCP<V> getDiagonal (const MAT& A) {
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

  /// Solve AX=B for X with Chebyshev iteration with left diagonal scaling.
  ///
  /// \pre A must be real-valued and symmetric positive definite.
  /// \pre numIters >= 0.
  /// \pre 0 < lambdaMin <= lambdaMax
  /// \pre All entries of D_inv are positive.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  /// \param B [in] Right-hand side(s) in the linear system to solve.
  /// \param X [in] Initial guess(es) for the linear system to solve.
  /// \param numIters [in] Number of Chebyshev iterations.
  /// \param lambdaMax [in] Estimate of max eigenvalue of A.
  /// \param lambdaMin [in] Estimate of min eigenvalue of A.
  /// \param eigRatio [in] Estimate of ratio of max to min eigenvalue of A.
  ///   This need not be the same as lambdaMax / lambdaMin.
  /// \param D_inv [in] Vector of diagonal entries of A.  It must have
  ///   the same distribution as B.
  void
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

  /// \brief Solve AX=B for X with Chebyshev iteration with left
  ///   diagonal scaling, imitating Ifpack's implementation.
  ///
  /// \pre A must be real-valued and symmetric positive definite.
  /// \pre numIters >= 0
  /// \pre eigRatio >= 1
  /// \pre 0 < lambdaMax
  /// \pre All entries of D_inv are positive.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  /// \param B [in] Right-hand side(s) in the linear system to solve.
  /// \param X [in] Initial guess(es) for the linear system to solve.
  /// \param numIters [in] Number of Chebyshev iterations.
  /// \param lambdaMax [in] Estimate of max eigenvalue of D_inv*A.
  /// \param lambdaMin [in] Estimate of min eigenvalue of D_inv*A.  We
  ///   only use this to determine if A is the identity matrix.
  /// \param eigRatio [in] Estimate of max / min eigenvalue ratio of
  ///   D_inv*A.  We use this along with lambdaMax to compute the
  ///   Chebyshev coefficients.  This need not be the same as
  ///   lambdaMax/lambdaMin.
  /// \param D_inv [in] Vector of diagonal entries of A.  It must have
  ///   the same distribution as b.
  void
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

  /// \brief Use numIters power method iterations to estimate the
  ///   maximum eigenvalue of A*D_inv.
  static ST
  powerMethod (const MAT& A, const V& D_inv, const int numIters)
  {
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

  /// \brief Solve AX=B for X with Chebyshev iteration with left
  ///   diagonal scaling, imitating ML's implementation.
  ///
  /// \pre A must be real-valued and symmetric positive definite.
  /// \pre numIters >= 0
  /// \pre eigRatio >= 1
  /// \pre 0 < lambdaMax
  /// \pre All entries of D_inv are positive.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  /// \param B [in] Right-hand side(s) in the linear system to solve.
  /// \param X [in] Initial guess(es) for the linear system to solve.
  /// \param numIters [in] Number of Chebyshev iterations.
  /// \param lambdaMax [in] Estimate of max eigenvalue of D_inv*A.
  /// \param lambdaMin [in] Estimate of min eigenvalue of D_inv*A.  We
  ///   only use this to determine if A is the identity matrix.
  /// \param eigRatio [in] Estimate of max / min eigenvalue ratio of
  ///   D_inv*A.  We use this along with lambdaMax to compute the
  ///   Chebyshev coefficients.  This need not be the same as
  ///   lambdaMax/lambdaMin.
  /// \param D_inv [in] Vector of diagonal entries of A.  It must have
  ///   the same distribution as b.
  void
  mlApplyImpl (const MAT& A,
	       const MV& B,
	       MV& X,
	       const int numIters,
	       const ST lambdaMax,
	       const ST lambdaMin,
	       const ST eigRatio,
	       const V& D_inv) const
  {
    const ST zero = Teuchos::as<ST> (0);
    const ST one = Teuchos::as<ST> (1);
    const ST two = Teuchos::as<ST> (2);

    MV pAux (B.getMap (), B.getNumVectors ()); // Result of A*X
    MV dk (B.getMap (), B.getNumVectors ()); // Solution update
    MV R (B.getMap (), B.getNumVectors ()); // Not in original ML; need for B - pAux

    ST beta = Teuchos::as<ST> (1.1) * lambdaMax;
    ST alpha = lambdaMax / eigRatio;

    ST delta = (beta - alpha) / two;
    ST theta = (beta + alpha) / two;
    ST s1 = theta / delta;
    ST rhok = one / s1;

    // Diagonal: ML replaces entries containing 0 with 1.  We
    // shouldn't have any entries like that in typical test problems,
    // so it's OK not to do that here.

    // The (scaled) matrix is the identity: set X = D_inv * B.  (If A
    // is the identity, then certainly D_inv is too.  D_inv comes from
    // A, so if D_inv * A is the identity, then we still need to apply
    // the "preconditioner" D_inv to B as well, to get X.)
    if (lambdaMin == one && lambdaMin == lambdaMax) {
      solve (X, D_inv, B);
      return;
    }

    // The next bit of code is a direct translation of code from ML's
    // ML_Cheby function, in the "normal point scaling" section, which
    // is in lines 7365-7392 of ml_smoother.c.

    if (! zeroStartingSolution_) {
      // dk = (1/theta) * D_inv * (B - (A*X))
      A.apply (X, pAux); // pAux = A * X
      R = B;
      R.update (-one, pAux, one); // R = B - pAux
      dk.elementWiseMultiply (one/theta, D_inv, R, zero); // dk = (1/theta)*D_inv*R
      X.update (one, dk, one); // X = X + dk
    } else {
      dk.elementWiseMultiply (one/theta, D_inv, B, zero); // dk = (1/theta)*D_inv*B
      X = dk;
    }

    ST rhokp1, dtemp1, dtemp2;
    for (int k = 0; k < numIters-1; ++k) {
      A.apply (X, pAux);
      rhokp1 = one / (two*s1 - rhok);
      dtemp1 = rhokp1*rhok;
      dtemp2 = two*rhokp1/delta;
      rhok = rhokp1;

      R = B;
      R.update (-one, pAux, one); // R = B - pAux
      // dk = dtemp1 * dk + dtemp2 * D_inv * (B - pAux)
      dk.elementWiseMultiply (dtemp2, D_inv, B, dtemp1);
      X.update (one, dk, one); // X = X + dk
    }
  }
};



/// \class CG
/// \brief Method of conjugate gradients, with left-scaling preconditioning.
/// \tparam ScalarType The type of entries in the matrix and vectors.
/// \tparam MV Specialization of Tpetra::MultiVector.
/// \tparam MAT Corresponding specialization of Tpetra::CrsMatrix.
///
/// This class requires that the matrix A be symmetric
/// (resp. Hermitian) positive definite.
template<class ScalarType, class MV, class MAT>
class CG {
public:
  typedef ScalarType ST;
  typedef Teuchos::ScalarTraits<ScalarType> STS;
  typedef typename STS::magnitudeType MT;
  typedef Tpetra::Vector<typename MV::scalar_type,
			 typename MV::local_ordinal_type,
			 typename MV::global_ordinal_type,
			 typename MV::node_type> V;
  /// Constructor.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  ///   A must be symmetric (resp. Hermitian) positive definite.
  CG (Teuchos::RCP<const MAT> A) : 
    A_ (A), 
    D_ (getDiagonal (*A)),
    numIters_ (1)
  {}

  /// Constructor with parameters.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  ///   A must be symmetric (resp. Hermitian) positive definite.
  /// \param params [in/out] On input: the parameters.  On output:
  ///   filled with the current parameter settings.
  CG (Teuchos::RCP<const MAT> A, Teuchos::ParameterList& params) : 
    A_ (A), 
    D_ (getDiagonal (*A)),
    numIters_ (1)
  {
    setParameters (params);
  }

  /// \brief Set (or reset) parameters.
  ///
  /// This method accepts the following parameters:
  /// - "relaxation: sweeps" (\c int): numIters, the number of iterations.
  ///
  /// \pre numIters >= 0
  void setParameters (Teuchos::ParameterList& plist) {
    int numIters = numIters_;
    if (plist.isParameter ("chebyshev: degree")) {
      numIters = plist.get<int> ("chebyshev: degree");
    } else if (plist.isParameter ("CG: iterations")) {
      numIters = plist.get<int> ("CG: iterations");
    } else {
      numIters = plist.get ("relaxation: sweeps", numIters);
    }
    numIters_ = numIters;
  }

  /// Solve Ax=b for x with Chebyshev iteration, using diagonal left preconditioning.
  ///
  /// \param b [in] Right-hand side(s) in the linear system to solve.
  /// \param x [in] Initial guess(es) for the linear system to solve.
  ///
  /// \return Max (over all columns) absolute residual 2-norm after iterating.
  MT apply (const MV& b, MV& x) {
    return leftScaledCG (*A_, b, x, numIters_, *D_);
  }

private:
  Teuchos::RCP<const MAT> A_;
  Teuchos::RCP<const V> D_;
  int numIters_;

  //! r = b - A * x
  static void 
  computeResidual (MV& r, const MV& b, const MAT& A, const MV& x, 
		   const Teuchos::ETransp mode = Teuchos::NO_TRANS) 
  {
    r = b;
    A.apply (x, r, mode, -STS::one(), STS::one());
  }

  //! z = D_inv * r, = D \ r.
  static void solve (MV& z, const V& D_inv, const MV& r) {
    z.elementWiseMultiply (STS::one(), D_inv, r, STS::zero());
  }

  //! Get a copy of the diagonal of the matrix, as a row Map vector.
  static Teuchos::RCP<V> getDiagonal (const MAT& A) {
    Teuchos::RCP<V> D (new V (A.getGraph ()->getRowMap ()));
    A.getLocalDiagCopy (*D);

    typedef Kokkos::MultiVector<ST, typename MAT::node_type> KMV;
    KMV& localDiag = D->getLocalMVNonConst ();
    typedef Kokkos::DefaultArithmetic<KMV> KMVT;
    KMVT::ReciprocalThreshold (localDiag, STS::eps ());

    return D;
  }

  /// Solve Ax=b for x with CG, using diagonal left preconditioning.
  ///
  /// \pre A must be real-valued and symmetric positive definite.
  /// \pre iterNum >= 0.
  /// \pre All entries of D_inv are positive.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  /// \param b [in] Right-hand side(s) in the linear system to solve.
  /// \param x [in] Initial guess(es) for the linear system to solve.
  /// \param numIters [in] Number of iterations.
  /// \param D_inv [in] Vector of diagonal entries of A.  It must have
  ///   the same distribution as b.
  ///
  /// \return Max (over all columns) absolute residual 2-norm after iterating.
  static MT
  leftScaledCG (const MAT& A,
		const MV& B,
		MV& X,
		const int numIters,
		const V& D_inv)
  {
    Teuchos::Array<MT> norms (B.getNumVectors ());
    for (size_t j = 0; j < B.getNumVectors (); ++j) {
      Teuchos::RCP<const V> b_j = B.getVector (j);
      Teuchos::RCP<V> x_j = X.getVectorNonConst (j);
      norms[j] = oneVecLeftScaledCG (A, *b_j, *x_j, numIters, D_inv);
    }
    return *std::max_element (norms.begin (), norms.end ());
  }

  /// Solve Ax=b for x with CG, using diagonal left preconditioning.
  ///
  /// \pre A must be real-valued and symmetric positive definite.
  /// \pre iterNum >= 0.
  /// \pre All entries of D_inv are positive.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  /// \param b [in] Right-hand side(s) in the linear system to solve.
  /// \param x [in] Initial guess(es) for the linear system to solve.
  /// \param numIters [in] Number of iterations.
  /// \param D_inv [in] Vector of diagonal entries of A.  It must have
  ///   the same distribution as b.
  ///
  /// \return Max (over all columns) absolute residual 2-norm after iterating.
  static MT
  oneVecLeftScaledCG (const MAT& A,
		      const V& b,
		      V& x,
		      const int numIters,
		      const V& D_inv)
  {
    const ST one = STS::one ();
    V r (b.getMap ());
    V p (b.getMap ());
    V q (b.getMap ());
    V z (b.getMap ());

    ST alpha, beta, rho, rho_prev;
    computeResidual (r, b, A, x); // r = b - A*x
    for (int i = 0; i < numIters; ++i) {
      solve (z, D_inv, r); // z = D_inv * r, that is, D \ r.
      rho = r.dot (z); // rho = r^T z; not sure if the order is right for complex arithmetic.
      if (i == 0) {
	p = z;
      } else {
	beta = rho / rho_prev;
	p.update (one, z, beta); // p = z + beta*p
      }
      A.apply (p, q);
      const ST p_dot_q = p.dot (q); // p_dot_q = p^T q; not sure if the order is right for complex arithmetic.
      alpha = rho / p_dot_q;
      x.update (+alpha, p, one); // x = x + alpha*p
      r.update (-alpha, q, one); // r = r - alpha*q
      rho_prev = rho;
    }

    const bool computeResidualNorm = true;
    if (computeResidualNorm) {
      computeResidual (r, b, A, x);
      return r.norm2 ();
    }
  }
};

//////////////////////////////////////////////////////////////////////
// Command-line arguments
//////////////////////////////////////////////////////////////////////

// They have long names so I don't confuse them with the shorter-named
// actual options in the body of the test.
int numberOfIterations = 50;
int numberOfEigenanalysisIterations = 15;
int localNumberOfRows = 10000;
bool perturbMaxEigenvalue = false;

} // namespace (anonymous) 


TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP ();
  clp.setOption ("numIters", &numberOfIterations, 
		 "Number of Chebyshev iterations");
  clp.setOption ("numEigIters", &numberOfEigenanalysisIterations, 
		 "Number of iterations of eigenvalue analysis (e.g., power method)");
  clp.setOption ("localNumRows", &localNumberOfRows, 
		 "Number of rows per process in the sparse matrix.");
  clp.setOption ("perturbMaxEigenvalue", "dontPerturbMaxEigenvalue", 
		 &perturbMaxEigenvalue, "If true, test sensitivity of the "
		 "Chebyshev implementations to changes in the provided max "
		 "eigenvalue estimate.");
}


TEUCHOS_UNIT_TEST(Ifpack2Chebyshev, Convergence)
{
  // We are now in a class method declared by the above macro, and
  // that method has these input arguments:
  // Teuchos::FancyOStream& out, bool& success

  using Tpetra::global_size_t;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using std::cout;
  using std::endl;

  // Typedefs for basic Tpetra template parameters.
  typedef double ST;
  typedef int LO;
  //typedef long GO;
  typedef int GO;
  //typedef Kokkos::SerialNode NT;
  typedef Kokkos::DefaultNode::DefaultNodeType NT;

  // Convenience typedefs.
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;
  typedef Tpetra::MultiVector<ST, LO, GO, NT> MV;
  typedef Tpetra::Vector<ST, LO, GO, NT> V;
  typedef Ifpack2::Chebyshev<crs_matrix_type> prec_type;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef STS::magnitudeType MT;

  const ST zero = STS::zero ();
  const ST one = STS::one ();
  const ST two = one + one;

  // Any real-valued, symmetric positive definite matrix is spectrally
  // equivalent to a diagonal matrix.  Thus, it suffices to test
  // Chebyshev with a diagonal matrix.

  // Prepare arguments for creating the Map.
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm ();
  RCP<NT> node;
  {
    ParameterList junk;
    node = rcp (new NT (junk));
  }
  const size_t localNumRows = as<size_t> (localNumberOfRows);
  const global_size_t globalNumRows = localNumRows * comm->getSize ();
  const GO indexBase = 0;
  const Tpetra::LocalGlobal lg = Tpetra::GloballyDistributed;

  // Create the row Map of the matrix, and the matrix's other Maps.
  RCP<const map_type> rowMap (new map_type (globalNumRows, indexBase, comm, lg, node));
  RCP<const map_type> rangeMap = rowMap;
  RCP<const map_type> domainMap = rowMap;

  // Create the matrix, with static profile.
  RCP<crs_matrix_type> A (new crs_matrix_type (rowMap, 3, Tpetra::StaticProfile));

  // Fill the matrix.
  Array<GO> cols (3);
  Array<ST> vals (3);
  for (GO globalRow = rowMap->getMinGlobalIndex (); 
       globalRow <= rowMap->getMaxGlobalIndex (); ++globalRow) {
    size_t numEntries = 3;
    if (globalRow == rowMap->getMinAllGlobalIndex ()) {
      numEntries = 2;
      cols[0] = globalRow;
      cols[1] = globalRow+1;
      vals[0] = two;
      vals[1] = -one;
    }
    else if (globalRow == rowMap->getMaxAllGlobalIndex ()) {
      numEntries = 2;
      cols[0] = globalRow-1;
      cols[1] = globalRow;
      vals[0] = -one;
      vals[1] = two;
    }
    else {
      numEntries = 3;
      cols[0] = globalRow-1;
      cols[1] = globalRow;
      cols[2] = globalRow+1;
      vals[0] = -one;
      vals[1] = two;
      vals[2] = -one;
    }
    ArrayView<const GO> colsView = cols.view (0, numEntries);
    ArrayView<const ST> valsView = vals.view (0, numEntries);
    A->insertGlobalValues (globalRow, colsView, valsView);
  }
  A->fillComplete (domainMap, rangeMap);

  // See James Demmel, "Applied Numerical Linear Algebra," SIAM,
  // pp. 267-8.  The eigenvalue approximations apply to the N x N
  // matrix with typical row (-1, 2, -1), representing the
  // discretization of the 1-D Poisson equation with Dirichlet
  // boundary conditions.
  const ST pi = acos (-1.0);
  const ST N = as<ST> (globalNumRows);
  const ST lambdaMax = two * (one - cos ((pi*N) / (N+one)));
  const ST lambdaMin = (pi / N) * (pi / N);
  ST eigRatio = lambdaMax / lambdaMin;
  const int numIters = numberOfIterations;
  int numEigIters = numberOfEigenanalysisIterations;

  // Set up the linear system to solve.
  V x_exact (domainMap), x (domainMap), b (rangeMap);
  x_exact.randomize ();
  A->apply (x_exact, b);
  x.putScalar (zero);

  V r (rangeMap); // for storing residual vector(s).
  Array<MT> norms (b.getNumVectors ());

  // Compute max initial absolute residual 2-norm.
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  const MT maxInitResNorm = *std::max_element (norms.begin (), norms.end ());

  cout << std::scientific;
  cout << endl
       << "numIters: " << numIters << endl
       << "localNumRows: " << localNumRows << endl
       << "globalNumRows: " << globalNumRows << endl
       << "lambdaMin: " << lambdaMin << endl
       << "lambdaMax: " << lambdaMax << endl
       << "eigRatio: " << eigRatio << endl
       << "Initial residual norm: " << maxInitResNorm << endl 
       << endl;

  Teuchos::ParameterList params;
  // Set parameters for the various Chebyshev implementations.  The
  // above Chebyshev class understands many of the same parameters as
  // Ifpack2, Ifpack, and ML.  For this first pass, we only set the
  // max eigenvalue.  Below, we'll experiment with also setting the
  // min eigenvalue and the min / max eigenvalue ratio.
  params.set ("chebyshev: eigenvalue max iterations", numEigIters);
  params.set ("chebyshev: degree", numIters);
  params.set ("chebyshev: max eigenvalue", lambdaMax);
  params.set ("chebyshev: imitate Ifpack", false);  
  params.set ("chebyshev: imitate ML", false);  

  // Create the operators: Ifpack2, custom Chebyshev, custom CG.
  prec_type ifpack2Cheby (A);
  Chebyshev<ST, MV, crs_matrix_type> myCheby (A);
  CG<ST, MV, crs_matrix_type> cg (A);

  // Residual 2-norms for comparison.
  MT maxResNormM1, maxResNormM2, maxResNormM3, maxResNormM4, maxResNormM5;

  ////////////////////////////////////////////////////////////////////
  // Test 1: set lambdaMax exactly, use default values of eigRatio and
  // lambdaMin.  Run each version of Chebyshev and compare results.
  ////////////////////////////////////////////////////////////////////

  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormM1 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormM2 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate Ifpack.
  params.set ("chebyshev: imitate Ifpack", true);
  params.set ("chebyshev: imitate ML", false);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM3 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate ML.
  params.set ("chebyshev: imitate Ifpack", false);
  params.set ("chebyshev: imitate ML", true);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM4 = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormM5 = cg.apply (b, x);

  cout << "Results with lambdaMax = " << lambdaMax 
       << ", default lambdaMin and eigRatio:" << endl
       << "- Ifpack2::Chebyshev:         " << maxResNormM1 / maxInitResNorm << endl
       << "- Chebyshev imitating Ifpack: " << maxResNormM3 / maxInitResNorm << endl
       << "- Textbook Chebyshev:         " << maxResNormM2 / maxInitResNorm << endl
       << "- Chebyshev imitating ML:     " << maxResNormM4 / maxInitResNorm << endl
       << "- CG:                         " << maxResNormM5 / maxInitResNorm << endl;

  // At least the first four or so digits of Ifpack2::Chebyshev and
  // "Chebyshev imitating Ifpack" should be the same.  (We're not so
  // dreadfully picky, but generally many more digits than this should
  // match.)
  {
    const MT tol = Teuchos::as<MT> (1.0e-4);
    const MT relDiff = STS::magnitude (maxResNormM1 - maxResNormM3) / maxResNormM3;
    TEUCHOS_TEST_FOR_EXCEPTION(
      relDiff > tol, 
      std::runtime_error, 
      "After " << numIters << " iterations of Chebyshev, with lambdaMax = " 
      << lambdaMax << " and default lambdaMin and eigRatio, Ifpack2::Chebyshev "
      "and an imitation of Ifpack's implementation have a relative difference "
      "of " << relDiff << " > tol = " << tol << ".");
  }

  // Reset the parameters to their original values.
  params.set ("chebyshev: imitate Ifpack", false);
  params.set ("chebyshev: imitate ML", false);

  ////////////////////////////////////////////////////////////////////
  // Test 2: set lambdaMax and lambdaMin exactly, and set eigRatio =
  // lambdaMax / lambdaMin.
  ////////////////////////////////////////////////////////////////////

  params.set ("chebyshev: min eigenvalue", lambdaMin);
  params.set ("chebyshev: ratio eigenvalue", eigRatio);

  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormM1 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormM2 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate Ifpack.
  params.set ("chebyshev: imitate Ifpack", true);
  params.set ("chebyshev: imitate ML", false);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM3 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate ML.
  params.set ("chebyshev: imitate Ifpack", false);
  params.set ("chebyshev: imitate ML", true);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM4 = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormM5 = cg.apply (b, x);

  cout << "Results with lambdaMax = " << lambdaMax 
       << ", lambdaMin = " << lambdaMin << ", eigRatio = " << eigRatio << endl
       << "- Ifpack2::Chebyshev:         " << maxResNormM1 / maxInitResNorm << endl
       << "- Chebyshev imitating Ifpack: " << maxResNormM3 / maxInitResNorm << endl
       << "- Textbook Chebyshev:         " << maxResNormM2 / maxInitResNorm << endl
       << "- Chebyshev imitating ML:     " << maxResNormM4 / maxInitResNorm << endl
       << "- CG:                         " << maxResNormM5 / maxInitResNorm << endl;

  // Reset parameters to their first values.
  params.set ("chebyshev: imitate Ifpack", false);
  params.set ("chebyshev: imitate ML", false);

  // ParameterList is NOT a delta.  That is, if we remove these
  // parameters from the list, setParameters() will use default
  // values, rather than letting the current settings remain.
  params.remove ("chebyshev: min eigenvalue", false);
  params.remove ("chebyshev: ratio eigenvalue", false);

  ////////////////////////////////////////////////////////////////////
  // Test 3: set lambdaMax exactly, and set eigRatio = 20 (ML's
  // default).
  ////////////////////////////////////////////////////////////////////

  // Now try again, but set the max / min eigenvalue ratio to 20 (ML's
  // default value for smoothers).
  eigRatio = Teuchos::as<ST> (20);
  params.set ("chebyshev: ratio eigenvalue", eigRatio);

  //
  // Run each version of Chebyshev and compare their results.
  //
  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormM1 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormM2 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate Ifpack.
  params.set ("chebyshev: imitate Ifpack", true);
  params.set ("chebyshev: imitate ML", false);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM3 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate ML.
  params.set ("chebyshev: imitate Ifpack", false);
  params.set ("chebyshev: imitate ML", true);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM4 = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormM5 = cg.apply (b, x);

  cout << "Results with lambdaMax = " << lambdaMax 
       << ", default lambdaMin, eigRatio = " << eigRatio << endl
       << "- Ifpack2::Chebyshev:         " << maxResNormM1 / maxInitResNorm << endl
       << "- Chebyshev imitating Ifpack: " << maxResNormM3 / maxInitResNorm << endl
       << "- Textbook Chebyshev:         " << maxResNormM2 / maxInitResNorm << endl
       << "- Chebyshev imitating ML:     " << maxResNormM4 / maxInitResNorm << endl
       << "- CG:                         " << maxResNormM5 / maxInitResNorm << endl;

  params.set ("chebyshev: imitate Ifpack", false);
  params.set ("chebyshev: imitate ML", false);

  ////////////////////////////////////////////////////////////////////
  // Test 4: set lambdaMax exactly, and set eigRatio = 30 (Ifpack's
  // default).
  ////////////////////////////////////////////////////////////////////
  eigRatio = Teuchos::as<ST> (30);
  params.set ("chebyshev: ratio eigenvalue", eigRatio);

  //
  // Run each version of Chebyshev and compare their results.
  //
  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormM1 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormM2 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate Ifpack.
  params.set ("chebyshev: imitate Ifpack", true);
  params.set ("chebyshev: imitate ML", false);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM3 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate ML.
  params.set ("chebyshev: imitate Ifpack", false);
  params.set ("chebyshev: imitate ML", true);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM4 = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormM5 = cg.apply (b, x);

  cout << "Results with lambdaMax = " << lambdaMax 
       << ", default lambdaMin, eigRatio = " << eigRatio << endl
       << "- Ifpack2::Chebyshev:         " << maxResNormM1 / maxInitResNorm << endl
       << "- Chebyshev imitating Ifpack: " << maxResNormM3 / maxInitResNorm << endl
       << "- Textbook Chebyshev:         " << maxResNormM2 / maxInitResNorm << endl
       << "- Chebyshev imitating ML:     " << maxResNormM4 / maxInitResNorm << endl
       << "- CG:                         " << maxResNormM5 / maxInitResNorm << endl;

  // Reset parameters to their original values.
  params.set ("chebyshev: imitate Ifpack", false);
  params.set ("chebyshev: imitate ML", false);
  params.remove ("chebyshev: ratio eigenvalue", false);
  eigRatio = lambdaMax / lambdaMin;

  ////////////////////////////////////////////////////////////////////
  // Test 5: Clear lambdaMax, lambdaMin, and eigRatio.  Let the
  // smoother do eigenanalysis to estimate lambdaMax.
  ////////////////////////////////////////////////////////////////////

  params.remove ("chebyshev: min eigenvalue", false);
  params.remove ("chebyshev: max eigenvalue", false);
  params.remove ("chebyshev: ratio eigenvalue", false);

  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormM1 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormM2 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate Ifpack.
  params.set ("chebyshev: imitate Ifpack", true);
  params.set ("chebyshev: imitate ML", false);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM3 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate ML.
  params.set ("chebyshev: imitate Ifpack", false);
  params.set ("chebyshev: imitate ML", true);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM4 = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormM5 = cg.apply (b, x);

  cout << "Results with default lambdaMax, lambdaMin, and eigRatio, "
    "with numEigIters = " << numEigIters << ":" << endl
       << "- Ifpack2::Chebyshev:         " << maxResNormM1 / maxInitResNorm << endl
       << "- Chebyshev imitating Ifpack: " << maxResNormM3 / maxInitResNorm << endl
       << "- Textbook Chebyshev:         " << maxResNormM2 / maxInitResNorm << endl
       << "- Chebyshev imitating ML:     " << maxResNormM4 / maxInitResNorm << endl
       << "- CG:                         " << maxResNormM5 / maxInitResNorm << endl;

  // Print the computed max and min eigenvalues, and other details.
  cout << endl;
  myCheby.print (cout);

  ////////////////////////////////////////////////////////////////////
  // Test 6: Clear lambdaMax, lambdaMin, and eigRatio.  Let the
  // smoother do eigenanalysis to estimate lambdaMax, with more
  // iterations.
  ////////////////////////////////////////////////////////////////////

  params.remove ("chebyshev: min eigenvalue", false);
  params.remove ("chebyshev: max eigenvalue", false);
  params.remove ("chebyshev: ratio eigenvalue", false);
  numEigIters = 2 * numEigIters;
  params.set ("chebyshev: eigenvalue max iterations", numEigIters);

  // Run Ifpack2's version of Chebyshev.
  ifpack2Cheby.setParameters (params);
  ifpack2Cheby.initialize ();
  ifpack2Cheby.compute ();
  ifpack2Cheby.apply (b, x);
  r = b;
  A->apply (x, r, Teuchos::NO_TRANS, -one, one);
  r.norm2 (norms ());
  maxResNormM1 = *std::max_element (norms.begin (), norms.end ());

  // Run our custom version of Chebyshev.
  x.putScalar (zero); // Reset the initial guess(es).
  myCheby.setParameters (params);
  myCheby.compute ();
  maxResNormM2 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate Ifpack.
  params.set ("chebyshev: imitate Ifpack", true);
  params.set ("chebyshev: imitate ML", false);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM3 = myCheby.apply (b, x);

  // Run our custom version of Chebyshev, but imitate ML.
  params.set ("chebyshev: imitate Ifpack", false);
  params.set ("chebyshev: imitate ML", true);
  myCheby.setParameters (params);
  x.putScalar (zero); // Reset the initial guess(es).
  maxResNormM4 = myCheby.apply (b, x);

  // Run CG, just to compare.
  x.putScalar (zero); // Reset the initial guess(es).
  cg.setParameters (params);
  maxResNormM5 = cg.apply (b, x);

  cout << "Results with default lambdaMax, lambdaMin, and eigRatio, "
    "with numEigIters = " << numEigIters << ":" << endl
       << "- Ifpack2::Chebyshev:         " << maxResNormM1 / maxInitResNorm << endl
       << "- Chebyshev imitating Ifpack: " << maxResNormM3 / maxInitResNorm << endl
       << "- Textbook Chebyshev:         " << maxResNormM2 / maxInitResNorm << endl
       << "- Chebyshev imitating ML:     " << maxResNormM4 / maxInitResNorm << endl
       << "- CG:                         " << maxResNormM5 / maxInitResNorm << endl;

  // Print the computed max and min eigenvalues, and other details.
  cout << endl;
  myCheby.print (cout);
}


