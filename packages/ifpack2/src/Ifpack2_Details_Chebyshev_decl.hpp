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

#ifndef IFPACK2_DETAILS_CHEBYSHEV_DECL_HPP
#define IFPACK2_DETAILS_CHEBYSHEV_DECL_HPP

/// \file Ifpack2_Details_Chebyshev_decl.hpp
/// \brief Declaration of Chebyshev implementation
/// \author Mark Hoemmen
///
/// This file is meant for Ifpack2 developers only, not for users.
/// It declares a new implementation of Chebyshev iteration.

#include <Ifpack2_ConfigDefs.hpp>
#include <Teuchos_VerbosityLevel.hpp>
#include <Teuchos_Describable.hpp>
#include <Tpetra_CrsMatrix.hpp>

namespace Ifpack2 {
/// \namespace Details
/// \brief Ifpack2 implementation details
///
/// This namespace contains implementation details of Ifpack2.
/// It is <i>not</i> meant for users.  Users should not rely on 
/// anything in this namespace.
namespace Details {
/// \class Chebyshev
/// \brief Left-scaled Chebyshev iteration.
/// \tparam ScalarType The type of entries in the matrix and vectors.
/// \tparam MV Specialization of Tpetra::MultiVector.
/// \tparam MAT Corresponding specialization of Tpetra::RowMatrix.
///
/// This class implements two variants of Chebyshev iteration:
/// 1. A direct imitation of Ifpack's implementation
/// 2. A textbook version of the algorithm
///
/// All implemented variants use the diagonal of the matrix to
/// precondition the linear system on the left.  Diagonal entries less
/// than machine precision are replaced with machine precision.
///
/// The first version imitates Ifpack::Chebyshev, both in how it sets
/// parameters and in the actual iteration (ApplyInverse()).  The
/// "textbook" in variant #2 above is "Templates for the Solution of
/// Linear Systems," 2nd edition.  Experiments show that the Ifpack
/// imitation is much less sensitive to the eigenvalue bounds than the
/// textbook version, so users should prefer it.  (In fact, it is the
/// default.)
///
/// We require that the matrix A be real valued and symmetric positive
/// definite.  If users could provide the ellipse parameters ("d" and
/// "c" in the literature, where d is the real-valued center of the
/// ellipse, and d-c and d+c the two foci), the iteration itself would
/// work fine with nonsymmetric real-valued A, as long as the
/// eigenvalues of A can be bounded in an ellipse that is entirely to
/// the right of the origin.
///
/// There is also dead code for imitating ML's Chebyshev
/// implementation (ML_Cheby(), in
/// packages/ml/src/Smoother/ml_smoother.c).  I couldn't get it to
/// converge in time to be useful for testing, so it is disabled.
template<class ScalarType, class MV, class MAT>
class Chebyshev : public Teuchos::Describable {
public:
  //! \name Typedefs
  //@{ 
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
  typedef Tpetra::Map<typename MV::local_ordinal_type,
		      typename MV::global_ordinal_type,
		      typename MV::node_type> map_type;
  //@}

  /// Constructor that takes a sparse matrix and sets default parameters.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  ///   A must be real-valued and symmetric positive definite.
  Chebyshev (Teuchos::RCP<const MAT> A);

  /// Constructor that takes a sparse matrix and sets the user's parameters.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  ///   A must be real-valued and symmetric positive definite.
  /// \param params [in/out] On input: the parameters.  On output:
  ///   filled with the current parameter settings.
  Chebyshev (Teuchos::RCP<const MAT> A, Teuchos::ParameterList& params);

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
  ///   the matrix A.  This parameter is optional and only used for a
  ///   quick check if the matrix is the identity matrix (if lambdaMax
  ///   == lambdaMin == 1).
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
  /// Parameters that govern other algorithmic details:
  /// - "chebyshev: operator inv diagonal" (<tt>RCP<const V></tt> or
  ///   <tt>const V*</tt>): If nonnull, we will use a deep copy of
  ///   this vector for left scaling as the inverse diagonal of the
  ///   matrix A, instead of computing the inverse diagonal ourselves.
  ///   We will make a copy every time you call setParameters().  If
  ///   you ever call setParameters() without this parameter, we will
  ///   clear our copy and compute the inverse diagonal ourselves
  ///   again.  You are responsible for updating this if the matrix
  ///   has changed.
  /// - "chebyshev: min diagonal value" (\c ST): minDiagVal.  If any
  ///   entry of the diagonal of the matrix is less than this in
  ///   magnitude, it will be replaced with this value in the inverse
  ///   diagonal used for left scaling.
  /// - "chebyshev: zero starting solution" (\c bool): If true, then
  ///   always use the zero vector(s) as the initial guess(es).  If
  ///   false, then apply() will use X on input as the initial
  ///   guess(es).
  ///
  /// Parameters that govern backwards compatibility:
  /// - "chebyshev: textbook algorithm" (\c bool): If true, use the
  ///   textbook version of Chebyshev iteration.  We recommend against
  ///   this, since the default algorithm is less sensitive to the
  ///   quality of the eigenvalue bounds.
  /// - "chebyshev: compute max residual norm" (\c bool): If true,
  ///   apply() will compute and return the max (absolute) residual
  ///   norm.  Otherwise, apply() returns 0.  This defaults to false.
  ///
  /// \pre lambdaMin, lambdaMax, and eigRatio are real
  /// \pre 0 < lambdaMin <= lambdaMax
  /// \pre numIters >= 0
  /// \pre eigMaxIters >= 0
  ///
  /// Default settings for parameters relating to spectral bounds come
  /// from Ifpack.
  void setParameters (Teuchos::ParameterList& plist);

  /// \brief (Re)compute the left scaling D_inv, and estimate min and
  ///   max eigenvalues of D_inv * A.
  ///
  /// You must call this method before calling apply(),
  /// - if you have not yet called this method,
  /// - if the matrix (either its values or its structure) has changed, or
  /// - any time after you call setParameters(). 
  ///
  /// Users have the option to supply the left scaling vector D_inv
  /// and estimates of the min and max eigenvalues of D_inv * A as
  /// parameters to setParameters().  If users did <i>not</i> supply a
  /// left scaling, then this method will compute it by default (if
  /// assumeMatrixUnchanged is false).  Likewise, if users did
  /// <i>not</i> supply at least an estimate of the max eigenvalue,
  /// this method will estimate it by default.  If estimation of the
  /// eigenvalues is required, this method may take as long as several
  /// Chebyshev iterations.
  ///
  /// Advanced users may avoid recomputing the left scaling and max
  /// eigenvalue estimate by setting assumeMatrixUnchanged to \c true,
  /// as long as none of the changed parameters affect either
  /// computation of the inverse diagonal, or estimation of the max or
  /// min eigenvalues.
  ///
  /// \param assumeMatrixUnchanged [in] If false, always compute the
  ///   left scaling D_inv and estimate the max eigenvalue of D_inv *
  ///   A.  If true, do not recompute these quantities.  (This method
  ///   will always compute them if they have not yet been computed
  ///   and if the user has not provided them to setParameters().)
  void compute (const bool assumeMatrixUnchanged=false);

  /// \brief Solve Ax=b for x with Chebyshev iteration with left diagonal scaling.
  ///
  /// \param B [in] Right-hand side(s) in the linear system to solve.
  /// \param X [in] Initial guess(es) for the linear system to solve.
  ///
  /// If the "chebyshev: compute max residual norm" parameter is true
  /// (not the default), then this method returns the maximum (over
  /// all columns) absolute residual 2-norm after iterating.
  /// Otherwise, it returns zero.
  ///
  /// \warning If you did not set the "chebyshev: zero starting
  ///   solution" parameter to true, then this method will use X as
  ///   the starting guess for Chebyshev iteration.  If you did not
  ///   initialize X before calling this method, then the resulting
  ///   solution will be undefined, since it will be computed using
  ///   uninitialized data.
  MT apply (const MV& B, MV& X);

  //! Get the matrix given to the constructor.
  Teuchos::RCP<const MAT> getMatrix () const;

  //! Whether it's possible to apply the transpose of this operator.
  bool hasTransposeApply () const;

  //! Print instance data to the given output stream.
  void print (std::ostream& out);

  //@}
  //! \name Implementation of Teuchos::Describable
  //@{ 

  //! A single-line description of the Chebyshev solver.
  std::string description() const;

  //! Print a description of the Chebyshev solver to \c out.
  void 
  describe (Teuchos::FancyOStream& out, 
	    const Teuchos::EVerbosityLevel verbLevel = 
	    Teuchos::Describable::verbLevel_default) const;
  //@}
private:
  //! \name The sparse matrix, and other related data.
  //@{

  Teuchos::RCP<const MAT> A_; //!< The sparse matrix A.

  /// The inverse of the diagonal entries of A.
  ///
  /// This is distributed using the range Map of the matrix.  If the
  /// user has not supplied the inverse diagonal (in setParameters(),
  /// as userInvDiag_), we compute this each time compute() is called.
  /// This ensures that compute() will respect changes to the values
  /// of the matrix.
  /// 
  /// If the user <i>has</i> supplied the inverse diagonal elements,
  /// compute() sets this to point to userInvDiag_.
  Teuchos::RCP<const V> D_;

  //@}
  //! \name Cached computed data
  //@{

  /// In ifpackApplyImpl(): the result of A*Y.
  /// We cache this multivector here to avoid creating it on each call.
  Teuchos::RCP<MV> V_;
  /// In ifpackApplyImpl(): Iteration update multivector.
  /// We cache this multivector here to avoid creating it on each call.
  Teuchos::RCP<MV> W_;

  /// Estimate that we compute for maximum eigenvalue of A.
  /// compute() will always recompute this. 
  /// This is set to NaN if it hasn't been computed yet.
  ST computedLambdaMax_;
  /// Estimate that we compute for minimum eigenvalue of A.
  /// compute() will always recompute this.
  /// This is set to NaN if it hasn't been computed yet.
  ST computedLambdaMin_;

  //@}
  //! \name Eigenvalue estimates to be used by apply().
  //@{

  /// Estimate for maximum eigenvalue of A.
  /// This is the value actually used by ifpackApplyImpl().
  ST lambdaMaxForApply_; 
  /// Estimate for minimum eigenvalue of A.
  /// This is the value actually used by ifpackApplyImpl().
  ST lambdaMinForApply_; 
  /// Estimate for ratio of max to min eigenvalue of A.
  /// This is the ratio actually used by ifpackApplyImpl().
  ST eigRatioForApply_;

  //@}
  //! \name Parameters given by the user to setParameters().
  //@{

  /// Range Map version of user-supplied inverse diagonal of the matrix A.
  /// This is null if the user did not provide it.
  Teuchos::RCP<const V> userInvDiag_;
  /// User-provided estimate for maximum eigenvalue of A.
  /// This is NaN if the user did not provide this.
  ST userLambdaMax_; 
  /// User-provided estimate for minimum eigenvalue of A.
  /// This is NaN if the user did not provide this.
  ST userLambdaMin_; 
  /// User-provided estimate for ratio of max to min eigenvalue of A.
  /// Not necessarily equal to userLambdaMax_ / userLambdaMin_.
  ST userEigRatio_;  
  /// Minimum allowed value on the diagonal of the matrix.  
  /// When computing the inverse diagonal, values less than this in
  /// magnitude are replaced with 1.
  ST minDiagVal_;
  //! Number of Chebyshev iterations to run on each call to apply().
  int numIters_;
  //! Number of power method iterations for estimating the max eigenvalue.
  int eigMaxIters_;
  //! Whether to assume that the X input to apply() is always zero.
  bool zeroStartingSolution_;
  //! Whether to use the textbook version of the algorithm.
  bool textbookAlgorithm_;
  //! Whether apply() will compute and return the max residual norm.
  bool computeMaxResNorm_;

  //@}
  //! \name Computational helper methods
  //@{

  //! Called by constructors to verify their input.
  void checkConstructorInput () const;

  /// \brief Set V and W to temporary multivectors with the same Map as X.
  ///
  /// \param V [out] 
  /// \param W [out]
  /// \param X [in] Multivector, whose Map to use when making V and W.
  ///
  /// This is an optimization for apply().  This method caches the
  /// created multivectors in the class instance as V_ resp. W_.
  /// Caching optimizes the common case of calling apply() many times.
  ///
  /// We call the first argument V1 instead of V, so as not to confuse
  /// the multivector V with the typedef V (for Tpetra::Vector).
  void
  makeTempMultiVectors (Teuchos::RCP<MV>& V1,
			Teuchos::RCP<MV>& W,
			const MV& X);

  //! R = B - Op(A) * X, where Op(A) is either A, \f$A^T\f$, or \f$A^H\f$.
  static void 
  computeResidual (MV& R, const MV& B, const MAT& A, const MV& X, 
		   const Teuchos::ETransp mode = Teuchos::NO_TRANS);

  /// \brief Z = D_inv * R, = D \ R.
  ///
  /// \param D_inv [in] A vector representing a diagonal matrix.
  /// \param R [in] Input multivector.
  /// \param Z [out] Result of multiplying the diagonal matrix D_inv with R.
  static void solve (MV& Z, const V& D_inv, const MV& R);

  /// \brief Z = alpha * D_inv * R, = alpha * (D \ R).
  ///
  /// \param D_inv [in] A vector representing a diagonal matrix.
  /// \param R [in] Input multivector.
  /// \param Z [out] Result of multiplying the diagonal matrix D_inv with R.
  static void solve (MV& Z, const ST alpha, const V& D_inv, const MV& R);

  //! Compute the inverse diagonal of the matrix, as a range Map vector.
  Teuchos::RCP<V> makeInverseDiagonal (const MAT& A) const;

  /// Return a range Map copy of the vector D.
  ///
  /// If *D is a range Map vector, return a copy of D.  Otherwise,
  /// Export D to a range Map vector and return the result.
  Teuchos::RCP<V> makeRangeMapVector (const V& D) const;

  /// \brief Solve AX=B for X with Chebyshev iteration with left
  ///   diagonal scaling, using the textbook version of the algorithm.
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
  textbookApplyImpl (const MAT& A,
		     const MV& B,
		     MV& X,
		     const int numIters,
		     const ST lambdaMax,
		     const ST lambdaMin,
		     const ST eigRatio,
		     const V& D_inv) const;

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
		   const V& D_inv);

  /// \brief Use numIters power method iterations to estimate the
  ///   maximum eigenvalue of A*D_inv.
  static ST
  powerMethod (const MAT& A, const V& D_inv, const int numIters);

  //! The maximum infinity norm of all the columns of X.
  static MT maxNormInf (const MV& X);

  // mfh 22 Jan 2013: This code builds correctly, but does not
  // converge.  I'm leaving it in place in case someone else wants to
  // work on it.
#if 0
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
	       const V& D_inv) 
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
#endif // 0
  //@}
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_CHEBYSHEV_DECL_HPP
