/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
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

#ifndef IFPACK2_CHEBYSHEV_DECL_HPP
#define IFPACK2_CHEBYSHEV_DECL_HPP

/// \file Ifpack2_Chebyshev_decl.hpp
/// \brief Chebyshev iteration

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Ifpack2_Details_Chebyshev_decl.hpp"
#include "Ifpack2_Details_Chebyshev_def.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_Time.hpp>

#include <iostream>
#include <string>
#include <sstream>


namespace Ifpack2 {

/// \class Chebyshev 
/// \brief Diagonally scaled Chebyshev iteration for Tpetra sparse matrices.
/// \tparam MatrixType A specialization of Tpetra::RowMatrix or Tpetra::CrsMatrix.
///
/// \section Ifpack_Chebyshev_Summary Summary
///
/// This class implements a Chebyshev polynomial preconditioner or
/// smoother for a Tpetra sparse matrix.  Given a matrix A, it applies
/// Chebyshev iteration to the left-scaled matrix \f$D^{-1} A\f$,
/// where D = diag(A) is the matrix of the diagonal entries of A.
/// This class' constructor accepts either a Tpetra::RowMatrix or a
/// Tpetra::CrsMatrix.  Its template parameter may be a specialization
/// of either class.
///
/// Chebyshev is derived from Preconditioner, which itself is derived
/// from Tpetra::Operator.  Therefore, a Chebyshev instance may be
/// used as an operator in any code that invokes the operator as
/// apply().
///
/// \warning Our implementation currently <i>only</i> works with a
///   real symmetric positive definite (SPD) matrix.  Results for
///   matrices that are not SPD, or for complex-valued Scalar types,
///   are not defined.
///
/// \section Ifpack_Chebyshev_Algorithm Algorithm
///
/// Given a matrix A, a right-hand side X, and an initial guess Y,
/// this class computes an approximate solution to \f$AY=X\f$ via
/// Chebyshev iteration using the left-scaled matrix \f$D^{-1} A\f$,
/// where D is the matrix of the diagonal elements of A.  (You may
/// control left scaling yourself if you wish, by providing an
/// optional vector of the entries of \f$D^{-1}\f$.)  While Chebyshev
/// iteration works for any matrix, we have chosen only to allow
/// real-valued, symmetric positive definite matrices.
///
/// Chebyshev iteration was originally intended as an iterative solver
/// for linear systems.  See the following publication (the spelling
/// of "Chebyshev" in Latin characters differs in some publications):
///
/// T. Manteuffel, "The Tchebychev iteration for nonsymmetric linear
/// systems," Numer. Math., 28 (1977), pp. 307-327.
///
/// It also works as a smoother for algebraic multigrid, which is the
/// target use case of this implementation.
///
/// \section Ifpack_Chebyshev_Eig Eigenvalue bounds
///
/// We require that the input matrix A be real-valued and symmetric
/// positive definite.  Thus, all of its eigenvalues must lie in a
/// positive interval on the real line.  Furthermore, if D is the
/// matrix of the diagonal elements of A, then the same is true of
/// \f$D^{-1} A\f$.  
///
/// Suppose \f$[\lambda_{min}, \lambda_{max}]\f$ is the interval of
/// the eigenvalues of \f$D^{-1} A\f$.  Users may either give us an
/// estimate of the maximum eigenvalue \f$\lambda_{max}\f$, or let us
/// compute it (which we do with a few power iterations).  They may
/// optionally also give us the (estimated) ratio \f$\eta =
/// \lambda_{max} / \lambda_{min}\f$, or (an estimate of) the minimum
/// eigenvalue \f$\lambda_{min}\f$.  The \f$\eta\f$ parameter
/// corresponds to the "smoother: Chebyshev alpha" parameter of ML.
/// (We use "eta" instead of "alpha" to avoid confusion with the
/// "alpha" argument of the apply() method of Tpetra::Operator.)
///
/// When using Chebyshev iteration by itself to solve linear systems,
/// it is important to have good estimates of both the minimum and
/// maximum eigenvalues.  However, when using a small number of
/// Chebyshev iterations as a smoother in multigrid, the maximum
/// eigenvalue estimate is more important.  (The point of a smoother
/// is to smooth out the high-frequency components of the error, that
/// is, those that correspond to the largest eigenvalues.  The coarser
/// grids below the current grid will take care of the lower-frequency
/// components of the error.)  This is why we use a ratio \f$\eta =
/// \lambda_{max} / \lambda_{min}$, rather than requiring a guess for
/// $\lambda_{min}$.  In fact, we only use \f$\lambda_{min}\f$ for
/// error checking, not when determining the Chebyshev coefficients.
/// Often, if users give us \f$\lambda_{max}\f$, our default value of
/// \f$\eta\f$ suffices.
///
/// Underestimating \f$\lambda_{min}\f$ may make Chebyshev fail to
/// converge, or fail to reduce the highest-frequency components of
/// the error, if used as a smoother.  Thus, we always multiply the
/// given \f$\lambda_{min}\f$ by a small factor (1.1).  This heuristic
/// accounts for the fact that typical methods for estimating extremal
/// eigenvalues (like Lanczos or CG) underestimate them.
///
/// If you do not give us an estimate for the maximum eigenvalue, we
/// estimate it using a few iterations of the power method in the
/// compute() method.  We do not attempt to refine the eigenvalue
/// bounds over Chebyshev iterations, as the typical smoother case
/// does not use very many iterations.  For an example of a Chebyshev
/// implementation that updates eigenvalue bound estimates, see Steve
/// Ashby's CHEBYCODE:
///
/// S. ASHBY, "CHEBYCODE: A Fortran implementation of Manteuffel's
/// adaptive Chebyshev algorithm," Tech. Rep. UIUCDCS-R-85-1203,
/// University of Illinois, 1985.
///
/// \section Ifpack_Chebyshev_Params Setting parameters
///
/// Call the setParameters() method to give this instance your
/// eigenvalue bound estimates (if you have them), as well as to set
/// other options controlling the behavior of Chebyshev iteration.
/// The documentation of setParameters() lists all the parameters that
/// this class accepts.  Where possible, we list comparable parameters
/// in the Ifpack package and the ML multigrid package.
///
/// \section Ifpack_Chebyshev_Performance Performance
///
/// Chebyshev should spend most of its time in Tpetra's native sparse
/// matrix-vector multiply kernel.  This should give good performance,
/// since we have spent a lot of effort tuning that kernel.  Depending
/// on the Kokkos Node type of your Tpetra matrix, the kernel may also
/// exploit threads for additional parallelism within each MPI process
/// ("hybrid parallelism" a.k.a. "MPI + X").  If your application
/// depends on hybrid parallelism for performance, you should favor
/// Chebyshev smoothers whenever possible over "serial within a
/// process" smoothers like Gauss-Seidel or SOR (Symmetric
/// Over-Relaxation).
///
/// \section Ifpack_Chebyshev_History History
///
/// The original implementation of this class was an adaptation of
/// ML's ML_Cheby routine.  The original author was Ulrich Hetmaniuk,
/// a Sandia employee in what was then (2006) Org 1416.  Ifpack2 has
/// seen significant development since then.
template<class MatrixType>
class Chebyshev : 
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
					   typename MatrixType::local_ordinal_type,
					   typename MatrixType::global_ordinal_type,
					   typename MatrixType::node_type> 
{
public:
  //! \name Typedefs
  //@{ 

  //! The template parameter of this class.
  typedef MatrixType matrix_type;

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! Preserved only for backwards compatibility.  Please use "scalar_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::scalar_type Scalar;


  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "local_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::local_ordinal_type LocalOrdinal;


  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "global_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::global_ordinal_type GlobalOrdinal;


  //! The type of the Kokkos Node used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! Preserved only for backwards compatibility.  Please use "node_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::node_type Node;


  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Preserved only for backwards compatibility.  Please use "magnitude_type".
  TEUCHOS_DEPRECATED typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitudeType;


  /// \brief The Tpetra::RowMatrix specialization matching MatrixType.
  ///
  /// MatrixType may be either a Tpetra::CrsMatrix specialization or a
  /// Tpetra::RowMatrix specialization.  This typedef will always be a
  /// Tpetra::RowMatrix specialization which is either the same as
  /// MatrixType, or the parent class of MatrixType.
  typedef Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> row_matrix_type;

  //! The Tpetra::Map specialization matching MatrixType.
  typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> map_type;

  /// \brief The Tpetra::Vector specialization matching MatrixType.
  ///
  /// If you wish to supply setParameters() a precomputed vector of
  /// diagonal entries of the matrix, use a pointer to an object of
  /// this type.
  typedef Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> vector_type;

  //@}
  // \name Constructors and destructors
  //@{

  /// \brief Constructor.
  ///
  /// \param A [in] The sparse matrix to which to apply Chebyshev
  ///   iteration.  The matrix A must be square, and its domain Map
  ///   and range Map must be the same.  The latter means that the
  ///   vectors x and y in the sparse matrix-vector product y = A*x
  ///   must both have the same distribution over process(es).
  ///
  /// We do <i>not</i> require that the row Map and the range Map of A
  /// be the same.  However, set-up will take less time if they are
  /// identical (in terms of pointer equality).  This is because we
  /// have to extract the diagonal entries of A as a row Map vector:
  /// if the row and range Maps are not identical, we have to
  /// redistribute the vector from the row Map to the range Map.
  ///
  /// The constructor will only check the requirements on the various
  /// Maps of A if the CMake configuration option
  /// <tt>Teuchos_ENABLE_DEBUG</tt> was set to <tt>ON</tt> before
  /// building Trilinos.  The checks require \f$O(1)\f$ global
  /// reductions over all processes in A's communicator, so we prefer
  /// to avoid them if we can.
  explicit Chebyshev (const Teuchos::RCP<const row_matrix_type>& A);

  //! Destructor.
  virtual ~Chebyshev ();

  //@}
  //! \name Preconditioner computation methods
  //@{ 

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
  /// \section Ifpack2_Chebyshev_setParameters_List List of parameters
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
  ///
  /// \pre lambdaMin, lambdaMax, and eigRatio are real
  /// \pre 0 < lambdaMin <= lambdaMax
  /// \pre numIters >= 0
  /// \pre eigMaxIters >= 0
  ///
  /// \section Ifpack2_Chebyshev_setParameters_compat Note on compatibility with Ifpack and ML
  ///
  /// Both the Ifpack and ML packages implement a Chebyshev smoother.
  /// We accept Ifpack and ML names for parameters whenever Ifpack2
  /// has an equivalent parameter.  Default settings for parameters
  /// relating to spectral bounds come from Ifpack.
  ///
  /// The following list maps from an ML parameter to its
  /// corresponding Ifpack2 parameter.
  /// - "smoother: Chebyshev alpha": "chebyshev: ratio eigenvalue"
  /// - "smoother: sweeps": "chebyshev: degree"
  ///
  /// ML does not have a parameter corresponding to "chebyshev: max
  /// eigenvalue", because ML estimates the spectral radius
  /// automatically.  Ifpack and Ifpack2 both can estimate this
  /// automatically, but also let the user provide an estimate.
  /// Similarly, ML does not have a parameter corresponding to
  /// "chebyshev: min eigenvalue".
  ///
  /// The following list maps from an Ifpack parameter to its
  /// corresponding Ifpack2 parameter.  Many of the parameters have
  /// the same names, in which case we simply write <i>same</i>.
  /// - "chebyshev: max eigenvalue": same
  /// - "chebyshev: ratio eigenvalue": same
  /// - "chebyshev: min eigenvalue": same
  /// - "chebyshev: degree": same
  /// - "relaxation: sweeps": "chebyshev: degree"
  /// - "chebyshev: min diagonal value": same
  /// - "relaxation: min diagonal value": "chebyshev: min diagonal value"
  /// - "chebyshev: zero starting solution": same
  /// - "relaxation: zero starting solution": "chebyshev: zero starting solution"
  /// - "chebyshev: operator inv diagonal": same
  ///
  /// \section Ifpack2_Chebyshev_setParameters_details Details on parameters
  ///
  /// The optional user-provided vector of diagonal entries of the
  /// matrix may have any distribution for which an Export to the
  /// range Map of the matrix is legal.  However, if the vector is
  /// already distributed according to the range Map, that saves us
  /// the communication cost of an Export.  We also avoid the Export
  /// in case the row Map and the range Map of the matrix are the
  /// same.  If they are not the same, and if the vector is
  /// distributed according to the row Map, we will reuse the Export
  /// from the matrix.  Otherwise, we have to make a fresh Export
  /// object, which is more expensive.  To avoid this cost, you should
  /// always provide a row Map or range Map vector for this parameter.
  void setParameters (const Teuchos::ParameterList& params);

  /// \brief Initialize the preconditioner.
  ///
  /// The compute() method will call initialize() automatically if it
  /// has not yet been called, so you do not normally need to call
  /// this.  However, it is correct to call initialize() yourself, and
  /// compute() will not call it again if it already has been called.
  void initialize();

  /// Whether the preconditioner has been successfully initialized
  /// (by calling initialize()).
  inline bool isInitialized() const {
    return IsInitialized_;
  }

  /// \brief (Re)compute the left scaling, and (if applicable)
  ///   estimate max and min eigenvalues of D_inv * A.
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
  ///
  /// This method will call initialize() if it has not already been
  /// called.  However, you may call initialize() before calling this
  /// method if you wish.
  void compute ();

  /// Whether compute() has been called at least once.
  ///
  /// Note that you must <i>always</i> call compute() if the matrix
  /// has changed, if you have called setParameters(), or if you have
  /// not yet called compute().  This method only tells you if
  /// compute() has been called at least once, not if you need to call
  /// compute().
  inline bool isComputed() const {
    return IsComputed_;
  }

  //@}
  //! @name Implementation of Tpetra::Operator
  //@{ 

  /// \brief Apply the preconditioner to X, returning the result in Y.
  ///
  /// This method actually computes Y = beta*Y + alpha*(M*X), where
  /// M*X represents the result of Chebyshev iteration on X, using the
  /// matrix Op(A).  Op(A) is either A itself, its transpose
  /// \f$A^T\f$, or its Hermitian transpose \f$A^H\f$, depending on
  /// the <tt>mode</tt> argument.  Since this class currently requires
  /// A to be real and symmetric positive definite, it should always
  /// be the case that \f$A = A^T = A^H\f$, but we will still respect
  /// the <tt>mode</tt> argument.
  ///
  /// \warning If you did not set the "chebyshev: zero starting
  ///   solution" parameter to true, then this method will use X as
  ///   the starting guess for Chebyshev iteration.  If you did not
  ///   initialize X before calling this method, then the resulting
  ///   solution will be undefined, since it will be computed using
  ///   uninitialized data.
  ///
  /// \param X [in] A (multi)vector to which to apply the preconditioner.
  /// \param Y [in/out] A (multi)vector containing the result of
  ///   applying the preconditioner to X.
  /// \param mode [in] If <tt>Teuchos::NO_TRANS</tt>, apply the matrix
  ///   A.  If <tt>mode</tt> is <tt>Teuchos::NO_TRANS</tt>, apply its
  ///   transpose \f$A^T\f$.  If <tt>Teuchos::CONJ_TRANS</tt>, apply
  ///   its Hermitian transpose \f$A^H\f$.
  /// \param alpha [in] Scaling factor for the result of Chebyshev
  ///   iteration.  The default is 1.
  /// \param beta [in] Scaling factor for Y.  The default is 0.
  void
  apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
	 Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
	 Teuchos::ETransp mode = Teuchos::NO_TRANS,
	 scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
	 scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! The Tpetra::Map representing the domain of this operator.
  const Teuchos::RCP<const map_type>& getDomainMap() const;

  //! The Tpetra::Map representing the range of this operator.
  const Teuchos::RCP<const map_type>& getRangeMap() const;

  //! Whether it's possible to apply the transpose of this operator.
  bool hasTransposeApply() const;

  /// \brief Compute Y = Op(A)*X, where Op(A) is either A, \f$A^T\f$, or \f$A^H\f$.
  ///
  /// \param X [in] Input (multi)vector of sparse matrix-vector
  ///   multiply.  If mode == Teuchos::NO_TRANS, X must be in the
  ///   domain Map of the matrix A.  Otherwise, X must be in the range
  ///   Map of A.
  /// \param Y [out] Output (multi)vector of sparse matrix-vector
  ///   multiply.  If mode == Teuchos::NO_TRANS, Y must be in the
  ///   range Map of the matrix A.  Otherwise, Y must be in the domain
  ///   Map of A.
  /// \param mode [in] Whether to apply the matrix A, its transpose
  ///   \f$A^T\f$, or its conjugate transpose \f$A^H\f$.  This method
  ///   applies A if <tt>mode</tt> is <tt>Teuchos::NO_TRANS</tt>,
  ///   \f$A^T\f$ if <tt>mode</tt> is <tt>Teuchos::TRANS</tt>, and
  ///   \f$A^H\f$ (the Hermitian transpose) if <tt>mode</tt> is
  ///   <tt>Teuchos::CONJ_TRANS</tt>.
  ///
  /// Since this class currently requires A to be real and symmetric
  /// positive definite, setting <tt>mode</tt> should not affect the
  /// result.
  void 
  applyMat (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
	    Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
	    Teuchos::ETransp mode = Teuchos::NO_TRANS) const;

  //@}
  //! \name Mathematical functions
  //@{

  //! Compute and return the estimated condition number.
  magnitude_type
  computeCondEst (CondestType CT = Cheap, 
		  local_ordinal_type MaxIters = 1550,
		  magnitude_type Tol = 1e-9,
		  const Teuchos::Ptr<const row_matrix_type>& matrix = Teuchos::null);

  //@}
  //! \name Attribute accessor methods
  //@{ 

  //! The estimated condition number, or -1.0 if it has not yet been computed.
  magnitude_type getCondEst() const;

  //! The communicator over which the matrix is distributed.
  const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

  //! The matrix for which this is a preconditioner.
  Teuchos::RCP<const row_matrix_type> getMatrix() const;

  //! The total number of floating-point operations taken by all calls to compute().
  double getComputeFlops() const;

  //! The total number of floating-point operations taken by all calls to apply().
  double getApplyFlops() const;

  //! The total number of successful calls to initialize().
  int getNumInitialize() const;

  //! The total number of successful calls to compute().
  int getNumCompute() const;

  //! The total number of successful calls to apply().
  int getNumApply() const;

  //! The total time spent in all calls to initialize().
  double getInitializeTime() const;

  //! The total time spent in all calls to compute().
  double getComputeTime() const;

  //! The total time spent in all calls to apply().
  double getApplyTime() const;

  //! The estimate of the maximum eigenvalue used in the apply().
  typename MatrixType::scalar_type getLambdaMaxForApply() const;

  //@}
  //! @name Implementation of Teuchos::Describable 
  //@{

  //! A simple one-line description of this object.
  std::string description() const;

  //! Print the object with some verbosity level to a Teuchos::FancyOStream.
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}
  //! \name Utility methods
  //@{

  /// Simple power method to compute lambda_max.
  ///
  /// This method is DEPRECATED.  Please don't call it any more.  You
  /// don't normally need to, because this class now automatically
  /// uses the power method (a different implementation) to estimate
  /// the max eigenvalue, if you don't give it an estimate yourself.
  static void TEUCHOS_DEPRECATED 
  PowerMethod (const Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Operator,
	       const vector_type& InvPointDiagonal,
	       const int MaximumIterations, 
	       scalar_type& LambdaMax);

  //! Not currently implemented: Use CG to estimate lambda_min and lambda_max.
  static void TEUCHOS_DEPRECATED 
  CG (const Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Operator, 
      const vector_type& InvPointDiagonal, 
      const int MaximumIterations, 
      scalar_type& lambda_min, scalar_type& lambda_max);

  //@}

private:
  //! Abbreviation for the Teuchos::ScalarTraits specialization for scalar_type.
  typedef Teuchos::ScalarTraits<typename MatrixType::scalar_type> STS;

  //! Abbreviation for the Tpetra::MultiVector specialization used in methods like apply().
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MV;
  
  //! Copy constructor (use is syntactically forbidden)
  Chebyshev (const Chebyshev<MatrixType>&);

  //! Assignment operator (use is syntactically forbidded)
  Chebyshev<MatrixType>& operator= (const Chebyshev<MatrixType>&);

  /// \brief Y := beta*Y + alpha*M*X.
  ///
  /// M*X represents the result of Chebyshev iteration with right-hand
  /// side(s) X and initial guess(es) Y, using the matrix Op(A).  Op(A)
  /// is A if mode is <tt>Teuchos::NO_TRANS</tt>, \f$A^T\f$ if mode is
  /// <tt>Teuchos::TRANS</tt>, and \f$A^H\f$ if mode is
  /// <tt>Teuchos::CONJ_TRANS</tt>.
  void 
  applyImpl (const MV& X,
	     MV& Y,
	     Teuchos::ETransp mode,
	     scalar_type alpha,
	     scalar_type beta) const;

  //! \name Internal state
  //@{

  /// The actual implementation of Chebyshev iteration.
  ///
  /// This is declared "mutable" because the implementation caches
  /// things in its apply() method, which makes it "morally" nonconst.
  /// I prefer that morals and syntax go together, so I didn't declare
  /// this class' apply() method const.  Hence, we have to declare the
  /// whole thing mutable here.
  mutable Details::Chebyshev<scalar_type, MV, row_matrix_type> impl_;

  //! Time object to track timing.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! The estimated condition number.
  magnitude_type Condest_;
  //! If \c true, initialize() has completed successfully.
  bool IsInitialized_;
  //! If \c true, compute() has completed successfully.
  bool IsComputed_;
  //! The total number of successful calls to initialize().
  int NumInitialize_;
  //! The total number of successful calls to compute().
  int NumCompute_;
  /// \brief The total number of successful calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable int NumApply_;
  //! The total time in seconds over all calls to initialize().
  double InitializeTime_;
  //! The total time in seconds over all calls to compute().
  double ComputeTime_;
  /// \brief The total time in seconds over all calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable double ApplyTime_;
  //! The total number of floating-point operations over all calls to compute().
  double ComputeFlops_;
  /// \brief The total number of floating-point operations over all calls to apply().
  ///
  /// This is "mutable" because apply() is a const method; apply() is
  /// const because it is declared this way in Tpetra::Operator.
  mutable double ApplyFlops_;

  //@}
}; // class Chebyshev

}//namespace Ifpack2

#endif // IFPACK2_CHEBYSHEV_DECL_HPP

