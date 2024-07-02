// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_CHEBYSHEV_DECL_HPP
#define IFPACK2_DETAILS_CHEBYSHEV_DECL_HPP

/// \file Ifpack2_Details_Chebyshev_decl.hpp
/// \brief Declaration of Chebyshev implementation
/// \author Mark Hoemmen
///
/// This file is meant for Ifpack2 developers only, not for users.
/// It declares a new implementation of Chebyshev iteration.

#include "Ifpack2_ConfigDefs.hpp"
#include "Teuchos_VerbosityLevel.hpp"
#include "Teuchos_Describable.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace Ifpack2 {
namespace Details {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<class TpetraOperatorType>
class ChebyshevKernel; // forward declaration
#endif // DOXYGEN_SHOULD_SKIP_THIS

/// \class Chebyshev
/// \brief Left-scaled Chebyshev iteration.
/// \tparam ScalarType The type of entries in the matrix and vectors.
/// \tparam MV Specialization of Tpetra::MultiVector.
///
/// \warning This class is NOT for users; it is an implementation
///   detail of Ifpack2.  Users should use Ifpack2::Chebyshev instead.
///
/// This class implements two variants of Chebyshev iteration:
/// <ol>
/// <li> A direct imitation of Ifpack's implementation </li>
/// <li> A textbook version of the algorithm </li>
/// <li> Chebyshev polynomials of the 4th kind, using optimal coefficients </li>
/// </ol>
///
/// All implemented variants use the diagonal of the matrix to
/// precondition the linear system on the left.  %Diagonal entries less
/// than machine precision are replaced with machine precision.
///
/// The first version imitates Ifpack_Chebyshev, both in how it sets
/// parameters and in the actual iteration (ApplyInverse()).  The
/// "textbook" in variant #2 above is "Templates for the Solution of
/// Linear Systems," 2nd edition.  Experiments show that the Ifpack
/// imitation is much less sensitive to the eigenvalue bounds than the
/// textbook version, so users should prefer it.  (In fact, it is the
/// default.) Variant #3 is an experimental implementation of Chebyshev
/// polynomials of the 4th kind with optimal coefficients,
/// from https://arxiv.org/pdf/2202.08830.pdf.
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
template<class ScalarType, class MV>
class Chebyshev : public Teuchos::Describable {
public:
  //! \name Typedefs
  //@{

  //! The type of entries in the matrix and vectors.
  typedef ScalarType ST;
  //! Traits class for ST.
  typedef Teuchos::ScalarTraits<ScalarType> STS;
  //! The type of the absolute value of a ScalarType.
  typedef typename STS::magnitudeType MT;
  //! Specialization of Tpetra::Operator.
  typedef Tpetra::Operator<typename MV::scalar_type,
                           typename MV::local_ordinal_type,
                           typename MV::global_ordinal_type,
                           typename MV::node_type> op_type;
  //! Specialization of Tpetra::RowMatrix.
  typedef Tpetra::RowMatrix<typename MV::scalar_type,
                           typename MV::local_ordinal_type,
                           typename MV::global_ordinal_type,
                           typename MV::node_type> row_matrix_type;
  //! Specialization of Tpetra::Vector.
  typedef Tpetra::Vector<typename MV::scalar_type,
                         typename MV::local_ordinal_type,
                         typename MV::global_ordinal_type,
                         typename MV::node_type> V;
  //! Specialization of Tpetra::Map.
  typedef Tpetra::Map<typename MV::local_ordinal_type,
                      typename MV::global_ordinal_type,
                      typename MV::node_type> map_type;
  //@}

  /// Constructor that takes a sparse matrix and sets default parameters.
  ///
  /// \param A [in] The matrix A in the linear system to solve.  If A
  ///   is nonnull, it must be real-valued and symmetric positive
  ///   definite.  The input A may be null.  In that case, you must
  ///   call setMatrix() with a nonnull input before you may call
  ///   compute() or apply().
  Chebyshev (Teuchos::RCP<const row_matrix_type> A);

  /// Constructor that takes a sparse matrix and sets the user's parameters.
  ///
  /// \param A [in] The matrix A in the linear system to solve.  If A
  ///   is nonnull, it must be real-valued and symmetric positive
  ///   definite.  The input A may be null.  In that case, you must
  ///   call setMatrix() with a nonnull input before you may call
  ///   compute() or apply().
  /// \param params [in/out] On input: the parameters.  On output:
  ///   filled with the current parameter settings.
  Chebyshev (Teuchos::RCP<const row_matrix_type> A, Teuchos::ParameterList& params);

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
  /// \section Ifpack2_Details_Chebyshev_setParameters_List List of parameters
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
  /// - "chebyshev: assume matrix does not change": Whether compute()
  ///   should always assume that the matrix has not changed since the
  ///   last call to compute().  The default is false.  If true,
  ///   compute() will not recompute the inverse diagonal or the
  ///   estimates of the max and min eigenvalues.  compute() will
  ///   always compute any quantity which the user did not provide and
  ///   which we have not yet computed before.
  /// - "chebyshev: operator inv diagonal" (<tt>RCP<const V></tt> or
  ///   <tt>const V*</tt>): If nonnull, we will use a deep copy of
  ///   this vector for left scaling as the inverse diagonal of the
  ///   matrix A, instead of computing the inverse diagonal ourselves.
  ///   We will make a copy every time you call setParameters().  If
  ///   you ever call setParameters() without this parameter, we will
  ///   clear our copy and compute the inverse diagonal ourselves
  ///   again.  If you choose to provide this parameter, you are
  ///   responsible for updating this if the matrix has changed.
  /// - "chebyshev: min diagonal value" (\c ST): minDiagVal.  If any
  ///   entry of the diagonal of the matrix is less than this in
  ///   magnitude, it will be replaced with this value in the inverse
  ///   diagonal used for left scaling.
  /// - "chebyshev: zero starting solution" (\c bool): If true, then
  ///   always use the zero vector(s) as the initial guess(es).  If
  ///   false, then apply() will use X on input as the initial
  ///   guess(es).
  /// - "chebyshev: compute spectral radius" (\c bool): If true, the
  ///   power method will compute the spectral radius of the operator.
  ///   If false, it will compute the dominant eigenvalue.
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
  /// The above compatibility parameters are not exposed in the
  /// documentation of Ifpack2::Chebyshev, because they are more
  /// useful to Ifpack2 developers than to users.
  ///
  /// \pre lambdaMin, lambdaMax, and eigRatio are real
  /// \pre 0 < lambdaMin <= lambdaMax
  /// \pre numIters >= 0
  /// \pre eigMaxIters >= 0
  ///
  /// Default settings for parameters relating to spectral bounds come
  /// from Ifpack.
  void setParameters (Teuchos::ParameterList& plist);

  void setZeroStartingSolution (bool zeroStartingSolution) { zeroStartingSolution_ = zeroStartingSolution; }

  /// \brief (Re)compute the left scaling D_inv, and estimate min and
  ///   max eigenvalues of D_inv * A.
  ///
  /// You must call this method before calling apply(),
  /// - if you have not yet called this method,
  /// - if the matrix (either its values or its structure) has changed, or
  /// - any time after you call setParameters().
  ///
  /// The input matrix must be nonnull before you may call this
  /// method.  If the input matrix is null, you must first call
  /// setMatrix() with a nonnull input matrix before you may call this
  /// method.
  ///
  /// Users have the option to supply the left scaling vector \c D_inv
  /// and estimates of the min and max eigenvalues of <tt>D_inv * A</tt>
  /// as parameters to setParameters().  If users did <i>not</i>
  /// supply a left scaling, then this method will compute it by
  /// default (if assumeMatrixUnchanged is false).  Likewise, if users
  /// did <i>not</i> supply at least an estimate of the max
  /// eigenvalue, this method will estimate it by default.  If
  /// estimation of the eigenvalues is required, this method may take
  /// as long as several Chebyshev iterations.
  ///
  /// Advanced users may avoid recomputing the left scaling vector and
  /// eigenvalue estimates by setting the "chebyshev: assume matrix
  /// does not change" parameter of setParameters() to \c true.  The
  /// left scaling vector and eigenvalue estimates will always be
  /// computed if the user did not provide them and we have not yet
  /// computed them.  Any changes to parameters that affect
  /// computation of the inverse diagonal or estimation of the
  /// eigenvalue bounds will not affect subsequent apply() operations,
  /// until the "chebyshev: assume matrix does not change" parameter
  /// is set back to \c false (its default value).
  void compute ();

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

  ST getLambdaMaxForApply() const;

  //! Get the matrix given to the constructor.
  Teuchos::RCP<const row_matrix_type> getMatrix () const;

  /// \brief Set the matrix.
  ///
  /// It's legal to call this method with a null input.  In that case,
  /// one must then call this method with a nonnull input before one
  /// may call compute() or apply().
  void setMatrix (const Teuchos::RCP<const row_matrix_type>& A);

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
  //! \name The sparse matrix, and data related to its diagonal.
  //@{

  /// \brief The input sparse matrix A.
  ///
  /// This comes either from the constructor, or from the last call to
  /// setMatrix().  It may be null, in which case it is not legal to
  /// call compute() or apply() until setMatrix() is called with a
  /// nonnull input.
  Teuchos::RCP<const row_matrix_type> A_;

  //! "Operator" implementing W := alpha*D_inv*(B-A*X) + beta*W and X := X+W.
  Teuchos::RCP<ChebyshevKernel<op_type> > ck_;

  /// \brief The inverse of the diagonal entries of A.
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

  //! The type of diagOffsets_ (see below).
  typedef Kokkos::View<size_t*, typename MV::node_type::device_type> offsets_type;

  /// \brief Precomputed offsets of local diagonal entries of the matrix.
  ///
  /// These are only used if the matrix has a const ("static") graph.
  /// In that case, the offsets of the diagonal entries will never
  /// change, even if the values of the diagonal entries change.
  offsets_type diagOffsets_;

  /// \brief Whether we have precomputed offsets of diagonal entries.
  ///
  /// We need this flag because it is not enough just to test if
  /// diagOffsets_ has size zero.  It is perfectly legitimate for the
  /// matrix to have zero rows on the calling process, in which case
  /// diagOffsets_ would have length zero (and would "equal
  /// Teuchos::null").
  bool savedDiagOffsets_;

  //@}
  //! \name Cached computed data
  //@{

  /// \brief In ifpackApplyImpl(): Iteration update MultiVector.
  ///
  /// We cache these multivectors here to avoid creating them on each call.
  Teuchos::RCP<MV> W_;
  Teuchos::RCP<MV> W2_;

  /// \brief Estimate that we compute for maximum eigenvalue of A.
  ///
  /// compute() will always recompute this, unless
  /// assumeMatrixUnchanged_ is true.  This is set to NaN if it hasn't
  /// been computed yet.
  ST computedLambdaMax_;

  /// Estimate that we compute for minimum eigenvalue of A.
  ///
  /// compute() will always recompute this, unless
  /// assumeMatrixUnchanged_ is true.  This is set to NaN if it hasn't
  /// been computed yet.
  ST computedLambdaMin_;

  //@}
  //! \name Eigenvalue estimates to be used by apply().
  //@{

  /// Estimate for maximum eigenvalue of A.
  /// This is the value actually used by ifpackApplyImpl().
  ST lambdaMaxForApply_;

  /// \brief Factor used to increase estimate of A's maximum eigenvalue.
  ///
  /// ifpackApplyImpl() multiplies lambdaMaxForApply_ (which see) by
  /// this factor. The idea is to ensure that A's maximum eigenvalue
  /// is less than the result. Otherwise the smoother could actually
  /// magnify high-energy error modes.  The default value is 1.1.
  MT boostFactor_;
  /// Estimate for minimum eigenvalue of A.
  /// This is the value actually used by ifpackApplyImpl().
  ST lambdaMinForApply_;
  /// Estimate for ratio of max to min eigenvalue of A.
  /// This is the ratio actually used by ifpackApplyImpl().
  ST eigRatioForApply_;

  //@}
  //! \name Parameters given by the user to setParameters().
  //@{

  /// \brief User-supplied inverse diagonal of the matrix A.
  ///
  /// This is null if the user did not provide it as a parameter to
  /// setParameters().
  Teuchos::RCP<const V> userInvDiag_;

  /// \brief User-provided estimate for maximum eigenvalue of A.
  ///
  /// This is NaN if the user did not provide this.
  ST userLambdaMax_;

  /// \brief User-provided estimate for minimum eigenvalue of A.
  ///
  /// This is NaN if the user did not provide this.
  ST userLambdaMin_;

  /// \brief User-provided estimate for ratio of max to min eigenvalue of A.
  ///
  /// Not necessarily equal to userLambdaMax_ / userLambdaMin_.
  ST userEigRatio_;

  /// \brief Minimum allowed value on the diagonal of the matrix.
  ///
  /// When computing the inverse diagonal, values less than this in
  /// magnitude are replaced with 1.
  ST minDiagVal_;

  //! Number of Chebyshev iterations to run on each call to apply().
  int numIters_;

  //! Number of power method iterations for estimating the max eigenvalue.
  int eigMaxIters_;

  //! Relative tolerance for power method iterations for estimating the max eigenvalue.
  MT eigRelTolerance_;

  //! Whether the iteration vectors of the power method should be saved.
  bool eigKeepVectors_;

  //! Type of eigen-analysis, i.e. power method or cg.
  std::string eigenAnalysisType_;

  //! Iteration vectors of the power method. Can be saved for the purpose of multiple setups.
  Teuchos::RCP<V> eigVector_;
  Teuchos::RCP<V> eigVector2_;

  //! Frequency of normalization in the power method.
  int eigNormalizationFreq_;

  //! Whether to assume that the X input to apply() is always zero.
  bool zeroStartingSolution_;

  /// Whether compute() should assume that the matrix has not changed.
  ///
  /// If true, compute() will not recompute the inverse diagonal or
  /// the estimates of the max and min eigenvalues.  compute() will
  /// always compute any quantity which the user did not provide and
  /// which we have not yet computed before.
  bool assumeMatrixUnchanged_;

  //! Chebyshev type
  std::string chebyshevAlgorithm_;

  //! Whether apply() will compute and return the max residual norm.
  bool computeMaxResNorm_;

  //! Whether the power method will compute the spectral radius or the dominant eigenvalue.
  bool computeSpectralRadius_;

  /// If true, the ChebyshevKernel operator will not to use a fused kernel
  /// and insead use native blas/SpMV operators
  bool ckUseNativeSpMV_;

  /// \brief Output stream for debug output ONLY.
  ///
  /// This is ONLY valid if debug_ is true.
  Teuchos::RCP<Teuchos::FancyOStream> out_;

  //! Whether to print copious debug output.
  bool debug_;

  //@}
  //! \name Computational helper methods
  //@{

  //! Called by constructors to verify their input.
  void checkConstructorInput () const;

  //! Called by constructors and setMatrix() to verify the input matrix.
  void checkInputMatrix () const;

  /// \brief Reset internal state dependent on the matrix.
  ///
  /// Calling this method forces recomputation of diagonal entries (if
  /// not provided by the user) and offsets, cached internal Vectors,
  /// and estimates of the min and eigenvalues.  It must be called in
  /// setMatrix(), unless assumeMatrixChanged_ is false or the input
  /// matrix to setMatrix() is the same object as A_.
  void reset ();

  /// \brief Set W to temporary MultiVector with the same Map as B.
  ///
  /// This is an optimization for apply().  This method caches the
  /// created MultiVector as W_.  Caching optimizes the common case of
  /// calling apply() many times.
  Teuchos::RCP<MV> makeTempMultiVector (const MV& B);

  /// \brief Set W2 to temporary MultiVector with the same Map as B.
  ///
  /// This is an optimization for apply(). This method caches the
  /// created MultiVector as W2_. Caching optimizes the common case of
  /// calling apply() many times. This is used by fourth kind
  /// Chebyshev as two temporary multivectors are needed.
  Teuchos::RCP<MV> makeSecondTempMultiVector (const MV& B);

  //! W = alpha*D_inv*B and X = 0 + W.
  void
  firstIterationWithZeroStartingSolution
  (MV& W,
   const ScalarType& alpha,
   const V& D_inv,
   const MV& B,
   MV& X);

  //! R = B - Op(A) * X, where Op(A) is either A, \f$A^T\f$, or \f$A^H\f$.
  static void
  computeResidual (MV& R, const MV& B, const op_type& A, const MV& X,
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

  /// \brief Compute the inverse diagonal of the matrix, as a range Map vector.
  ///
  /// \param A [in] The sparse matrix for which to compute the inverse
  ///   diagonal.
  /// \param useDiagOffsets [in] If true, use previously computed
  ///   offsets of diagonal entries (diagOffsets_) to speed up
  ///   extracting the diagonal entries of the sparse matrix A.
  Teuchos::RCP<const V>
  makeInverseDiagonal (const row_matrix_type& A, const bool useDiagOffsets=false) const;

  /// Return a range Map copy of the vector D.
  ///
  /// If *D is a range Map Vector, return a shallow copy of D.
  /// Otherwise, Export D to a new range Map Vector and return the
  /// result.  (The latter case might happen if D is a row Map Vector,
  /// for example.)  This method takes D as an RCP so that it can
  /// return a shallow copy of D if appropriate.
  ///
  /// \param D [in] Nonnull Vector.
  Teuchos::RCP<V> makeRangeMapVector (const Teuchos::RCP<V>& D) const;

  //! Version of makeRangeMapVector() that takes const input and returns const output.
  Teuchos::RCP<const V>
  makeRangeMapVectorConst (const Teuchos::RCP<const V>& D) const;

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
  textbookApplyImpl (const op_type& A,
                     const MV& B,
                     MV& X,
                     const int numIters,
                     const ST lambdaMax,
                     const ST lambdaMin,
                     const ST eigRatio,
                     const V& D_inv) const;

  /// \brief Solve AX=B for X with Chebyshev iteration with left
  ///   diagonal scaling, using the fourth kind Chebyshev polynomials
  ///   (optionally) with optimal weights, see:
  ///    https://arxiv.org/pdf/2202.08830.pdf
  ///   for details.
  ///
  /// \pre A must be real-valued and symmetric positive definite.
  /// \pre numIters <= 16 if using the opt. 4th-kind Chebyshev smoothers
  ///      -- this is an arbitrary distinction,
  ///      but the weights are currently hard-coded from the
  ///      MATLAB scripts from https://arxiv.org/pdf/2202.08830.pdf.
  /// \pre 0 < lambdaMax
  /// \pre All entries of D_inv are positive.
  ///
  /// \param A [in] The matrix A in the linear system to solve.
  /// \param B [in] Right-hand side(s) in the linear system to solve.
  /// \param X [in] Initial guess(es) for the linear system to solve.
  /// \param numIters [in] Number of Chebyshev iterations.
  /// \param lambdaMax [in] Estimate of max eigenvalue of A.
  /// \param D_inv [in] Vector of diagonal entries of A.  It must have
  ///   the same distribution as B.
  void
  fourthKindApplyImpl (const op_type& A,
                       const MV& B,
                       MV& X,
                       const int numIters,
                       const ST lambdaMax,
                       const V& D_inv);

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
  ifpackApplyImpl (const op_type& A,
                   const MV& B,
                   MV& X,
                   const int numIters,
                   const ST lambdaMax,
                   const ST lambdaMin,
                   const ST eigRatio,
                   const V& D_inv);

  /// \brief Use the cg method to estimate the maximum eigenvalue
  ///   of A*D_inv, given an initial guess vector x.
  ///
  /// \param A [in] The Operator to use.
  /// \param D_inv [in] Vector to use as implicit right scaling of A.
  /// \param numIters [in] Maximum number of iterations of the power
  ///   method.
  /// \param x [in/out] On input: Initial guess Vector for the cg
  ///   method.  Its Map must be the same as that of the domain Map of
  ///   A.  This method may use this Vector as scratch space.
  ///
  /// \return Estimate of the maximum eigenvalue of A*D_inv.
  ST
  cgMethodWithInitGuess (const op_type& A, const V& D_inv, const int numIters, V& x);

  /// \brief Use the cg method to estimate the maximum eigenvalue
  ///   of A*D_inv.
  ///
  /// \param A [in] The Operator to use.
  /// \param D_inv [in] Vector to use as implicit right scaling of A.
  /// \param numIters [in] Maximum number of iterations of the power
  ///   method.
  ///
  /// \return Estimate of the maximum eigenvalue of A*D_inv.
  ST
  cgMethod (const op_type& A, const V& D_inv, const int numIters);

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
  mlApplyImpl (const op_type& A,
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
