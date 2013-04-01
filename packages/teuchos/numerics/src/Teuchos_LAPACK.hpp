// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#ifndef _TEUCHOS_LAPACK_HPP_
#define _TEUCHOS_LAPACK_HPP_

/*! \file Teuchos_LAPACK.hpp
    \brief Templated interface class to LAPACK routines.
*/
/** \example LAPACK/cxx_main.cpp
    This is an example of how to use the Teuchos::LAPACK class.
*/

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_ScalarTraits.hpp"

/*! \class Teuchos::LAPACK
    \brief The Templated LAPACK Wrapper Class.

    The Teuchos::LAPACK class is a wrapper that encapsulates LAPACK
    (Linear Algebra Package).  LAPACK provides portable, high-
    performance implementations of linear, eigen, SVD, etc solvers.

    The standard LAPACK interface is Fortran-specific.  Unfortunately, the
    interface between C++ and Fortran is not standard across all computer
    platforms.  The Teuchos::LAPACK class provides C++ wrappers for the LAPACK
    kernels in order to insulate the rest of Teuchos from the details of C++ to Fortran
    translation.  A Teuchos::LAPACK object is essentially nothing, but allows access to
    the LAPACK wrapper functions.

    Teuchos::LAPACK is a serial interface only.  This is appropriate since the standard
    LAPACK are only specified for serial execution (or shared memory parallel).

    \note
	<ol>
		<li>These templates are specialized to use the Fortran LAPACK routines for
		scalar types \c float and \c double.

		<li>If Teuchos is configured with \c -DTeuchos_ENABLE_COMPLEX:BOOL=ON then these templates
		are specialized for scalar types \c std::complex<float> and \c std::complex<double> also.

		<li>A short description is given for each method.  For more detailed documentation, see the
		LAPACK website (\c http://www.netlib.org/lapack/ ).
	</ol>
*/

namespace Teuchos
{

  template<class T>
  struct UndefinedLAPACKRoutine
  {
    // This function should not compile if there is an attempt to instantiate!
    static inline T notDefined() { return T::LAPACK_routine_not_defined_for_this_type(); }
  };

  template<typename OrdinalType, typename ScalarType>
  class LAPACK
  {
  public:

    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

    //! @name Constructors/Destructors.
    //@{

    //! Default Constructor.
    inline LAPACK(void) {}

    //! Copy Constructor.
    inline LAPACK(const LAPACK<OrdinalType, ScalarType>& lapack) {}

    //! Destructor.
    inline virtual ~LAPACK(void) {}
    //@}

    //! @name Symmetric Positive Definite Linear System Routines.
    //@{

    //! Computes the \c L*D*L' factorization of a Hermitian/symmetric positive definite tridiagonal matrix \c A.
    void PTTRF(const OrdinalType n, ScalarType* d, ScalarType* e, OrdinalType* info) const;

    //! Solves a tridiagonal system \c A*X=B using the \L*D*L' factorization of \c A computed by PTTRF.
    void PTTRS(const OrdinalType n, const OrdinalType nrhs, const ScalarType* d, const ScalarType* e, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes Cholesky factorization of a real symmetric positive definite matrix \c A.
    void POTRF(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const;

    //! Solves a system of linear equations \c A*X=B, where \c A is a symmetric positive definite matrix factored by POTRF and the \c nrhs solutions are returned in \c B.
    void POTRS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes the inverse of a real symmetric positive definite matrix \c A using the Cholesky factorization \c A from POTRF.
    void POTRI(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const;

    //! Estimates the reciprocal of the condition number (1-norm) of a real symmetric positive definite matrix \c A using the Cholesky factorization from POTRF.

    void POCON(const char UPLO, const OrdinalType n, const ScalarType* A, const OrdinalType lda, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Computes the solution to a real system of linear equations \c A*X=B, where \c A is a symmetric positive definite matrix and the \c nrhs solutions are returned in \c B.
    void POSV(const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes row and column scalings intended to equilibrate a symmetric positive definite matrix \c A and reduce its condition number (w.r.t. 2-norm).
    void POEQU(const OrdinalType n, const ScalarType* A, const OrdinalType lda, MagnitudeType* S, MagnitudeType* scond, MagnitudeType* amax, OrdinalType* info) const;

    //! Improves the computed solution to a system of linear equations when the coefficient matrix is symmetric positive definite, and provides error bounds and backward error estimates for the solution.
    void PORFS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Uses the Cholesky factorization to compute the solution to a real system of linear equations \c A*X=B, where \c A is symmetric positive definite.  System can be equilibrated by POEQU and iteratively refined by PORFS, if requested.
    void POSVX(const char FACT, const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* AF, const OrdinalType ldaf, char EQUED, ScalarType* S, ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* rcond, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;
    //@}

    //! @name General Linear System Routines.
    //@{

    //! Solves an over/underdetermined real \c m by \c n linear system \c A using QR or LQ factorization of A.
    void GELS(const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /// \brief Use the SVD to solve a possibly rank-deficient linear least-squares problem.
    ///
    /// GELSS uses the singular value decomposition (SVD) to compute
    /// the minimum-norm solution to a possibly rank-deficient linear
    /// least-squares problem.  The problem may be under- or
    /// overdetermined.
    ///
    /// LAPACK's _GELSS routines take different arguments, depending
    /// on whether they are for real or complex arithmetic.  This is
    /// because _GELSS imitates the interface of LAPACK's SVD routine.
    /// LAPACK's SVD routine takes an additional RWORK workspace array
    /// argument for COMPLEX*8 (CGELSS) and COMPLEX*16 (ZGELSS).
    /// LAPACK's real SVD routines (SGELSS and DGELSS) do not take the
    /// RWORK argument.
    ///
    /// This class had already exposed GELSS for ScalarType = float
    /// and double that does <i>not</i> include an RWORK argument.
    /// Backwards compatibility requirements prevent us from simply
    /// changing that interface.  We could provide a different
    /// interface for LAPACK specializations with ScalarType =
    /// std::complex<T>, but that would make the GELSS interface not
    /// generic at compile time.  This would make using GELSS in
    /// generic code harder (for example, you would need to specialize
    /// code that <i>uses</i> GELSS on a Boolean, which specifies
    /// whether ScalarType is complex).
    ///
    /// We fix this problem by providing an overloaded generic GELSS
    /// interface that does take an RWORK argument.  This does not
    /// change the existing interface, but provides the additional
    /// capability to solve complex-valued least-squares problems.
    /// The RWORK argument is ignored when ScalarType is real, and may
    /// therefore be set to NULL in that case.
    ///
    void GELSS(const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, MagnitudeType* S, const MagnitudeType rcond, OrdinalType* rank, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const;

    //! Legacy GELSS interface for real-valued ScalarType.
    void GELSS(const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* S, const ScalarType rcond, OrdinalType* rank, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    //! Solves the linear equality-constrained least squares (LSE) problem where \c A is an \c m by \c n matrix,\c B is a \c p by \c n matrix \c C is a given \c m-vector, and D is a given \c p-vector.
    void GGLSE(const OrdinalType m, const OrdinalType n, const OrdinalType p, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* C, ScalarType* D, ScalarType* X, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    //! Computes a QR factorization of a general \c m by \c n matrix \c A.
    void GEQRF( const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    //! Computes an LU factorization of a general \c m by \c n matrix \c A using partial pivoting with row interchanges.
    void GETRF(const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const;

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B with a general \c n by \c n matrix \c A using the LU factorization computed by GETRF.
    void GETRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Multiplies the \c m by \c n matrix \c A by the real scalar \c cto/cfrom.
    void LASCL(const char TYPE, const OrdinalType kl, const OrdinalType ku, const MagnitudeType cfrom, const MagnitudeType cto, const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const;

    //! Apply a series of row interchanges to the matrix A.
    void
    LASWP (const OrdinalType N,
	   ScalarType A[],
	   const OrdinalType LDA,
	   const OrdinalType K1,
	   const OrdinalType K2,
	   const OrdinalType IPIV[],
	   const OrdinalType INCX) const;

    //! Computes an LU factorization of a general banded \c m by \c n matrix \c A using partial pivoting with row interchanges.
    void GBTRF(const OrdinalType m, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const;

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B with a general banded \c n by \c n matrix \c A using the LU factorization computed by GBTRF.
    void GBTRS(const char TRANS, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes an LU factorization of a \c n by \c n tridiagonal matrix \c A using partial pivoting with row interchanges.
    void GTTRF(const OrdinalType n, ScalarType* dl, ScalarType* d, ScalarType* du, ScalarType* du2, OrdinalType* IPIV, OrdinalType* info) const;

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B or \c A^H*X=B with a tridiagonal matrix \c A using the LU factorization computed by GTTRF.
    void GTTRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* dl, const ScalarType* d, const ScalarType* du, const ScalarType* du2, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes the inverse of a matrix \c A using the LU factorization computed by GETRF.
    void GETRI(const OrdinalType n, ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /// \brief Robustly solve a possibly singular triangular linear system.
    ///
    /// \note This routine is slower than the BLAS' TRSM, but can
    ///   detect possible singularity of A.
    void
    LATRS (const char UPLO,
	   const char TRANS,
	   const char DIAG,
	   const char NORMIN,
	   const OrdinalType N,
	   ScalarType* A,
	   const OrdinalType LDA,
	   ScalarType* X,
	   MagnitudeType* SCALE,
	   MagnitudeType* CNORM,
	   OrdinalType* INFO) const;

    //! Estimates the reciprocal of the condition number of a general real matrix \c A, in either the 1-norm or the infinity-norm, using the LU factorization computed by GETRF.
    void GECON(const char NORM, const OrdinalType n, const ScalarType* A, const OrdinalType lda, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Estimates the reciprocal of the condition number of a general banded real matrix \c A, in either the 1-norm or the infinity-norm, using the LU factorization computed by GETRF.
    void GBCON(const char NORM, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Returns the value of the one norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of an \c n by \c n band matrix \c A, with \c kl sub-diagonals and \c ku super-diagonals.
    typename ScalarTraits<ScalarType>::magnitudeType LANGB(const char NORM, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const ScalarType* A, const OrdinalType lda, MagnitudeType* WORK) const;

    //! Computes the solution to a real system of linear equations \c A*X=B, where \c A is factored through GETRF and the \c nrhs solutions are computed through GETRS.
    void GESV(const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes row and column scalings intended to equilibrate an \c m by \c n matrix \c A and reduce its condition number.
    void GEEQU(const OrdinalType m, const OrdinalType n, const ScalarType* A, const OrdinalType lda, ScalarType* R, ScalarType* C, ScalarType* rowcond, ScalarType* colcond, ScalarType* amax, OrdinalType* info) const;

    //! Improves the computed solution to a system of linear equations and provides error bounds and backward error estimates for the solution.  Use after GETRF/GETRS.
    void GERFS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Computes row and column scalings intended to equilibrate an \c m by \c n banded matrix \c A and reduce its condition number.
    void GBEQU(const OrdinalType m, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const ScalarType* A, const OrdinalType lda, MagnitudeType* R, MagnitudeType* C, MagnitudeType* rowcond, MagnitudeType* colcond, MagnitudeType* amax, OrdinalType* info) const;

    //! Improves the computed solution to a banded system of linear equations and provides error bounds and backward error estimates for the solution.  Use after GBTRF/GBTRS.
    void GBRFS(const char TRANS, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    //! Uses the LU factorization to compute the solution to a real system of linear equations \c A*X=B, returning error bounds on the solution and a condition estimate.
    void GESVX(const char FACT, const char TRANS, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* AF, const OrdinalType ldaf, OrdinalType* IPIV, char EQUED, ScalarType* R, ScalarType* C, ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* rcond, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const;

    /*! \brief Reduces a real symmetric matrix \c A to tridiagonal form by orthogonal similarity transformations.
	\note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    void SYTRD(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* D, ScalarType* E, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    //! Reduces a real general matrix \c A to upper Hessenberg form by orthogonal similarity transformations.
    void GEHRD(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* A, const OrdinalType lda, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    //! Solves a triangular linear system of the form \c A*X=B or \c A**T*X=B, where \c A is a triangular matrix.
    void TRTRS(const char UPLO, const char TRANS, const char DIAG, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const;

    //! Computes the inverse of an upper or lower triangular matrix \c A.
    void TRTRI(const char UPLO, const char DIAG, const OrdinalType n, const ScalarType* A, const OrdinalType lda, OrdinalType* info) const;
    //@}

    //! @name Symmetric Eigenproblem Routines
    //@{
    /*! \brief Computes the eigenvalues and, optionally, eigenvectors of a symmetric \c n by \c n matrix \c A in packed storage.
	\note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    void SPEV(const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* AP, ScalarType* W, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, OrdinalType* info) const;

    /*! \brief Computes all the eigenvalues and, optionally, eigenvectors of a symmetric \c n by \c n matrix A.
	\note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    void SYEV(const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* W, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /*! \brief Computes all the eigenvalues and, optionally, eigenvectors of a symmetric \c n by \c n matrix pencil \c {A,B}, where \c A is symmetric and \c B is symmetric positive-definite.
	\note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    void SYGV(const OrdinalType itype, const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* W, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /*! \brief Computes all the eigenvalues and, optionally, eigenvectors of a Hermitian \c n by \c n matrix A.
	\note This method will call SYEV when ScalarType is \c float or \c double.
    */
    void HEEV(const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, MagnitudeType* W, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const;

    /*! \brief Computes all the eigenvalues and, optionally, eigenvectors of a generalized Hermitian-definite \c n by \c n matrix pencil \c {A,B}, where \c A is Hermitian and \c B is Hermitian positive-definite.
	\note This method will call SYGV when ScalarType is \c float or \c double.
    */
    void HEGV(const OrdinalType itype, const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, MagnitudeType* W, ScalarType* WORK, const OrdinalType lwork, MagnitudeType *RWORK, OrdinalType* info) const;

    //! Computes the eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal \c n by \c n matrix \c A using implicit QL/QR.  The eigenvectors can only be computed if \c A was reduced to tridiagonal form by SYTRD.
    void STEQR(const char COMPZ, const OrdinalType n, ScalarType* D, ScalarType* E, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, OrdinalType* info) const;
    //@}

    //! @name Non-Hermitian Eigenproblem Routines
    //@{
    //! Computes the eigenvalues of a real upper Hessenberg matrix \c H and, optionally, the matrices \c T and \c Z from the Schur decomposition, where T is an upper quasi-triangular matrix and Z contains the Schur vectors.
    void HSEQR(const char JOB, const char COMPZ, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* H, const OrdinalType ldh, ScalarType* WR, ScalarType* WI, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /*! Computes for an \c n by \c n nonsymmetric matrix \c A, the eigenvalues, the Schur form \c T, and, optionally, the matrix of Schur vectors \c Z. When \c ScalarType is \c float or \c double, the real Schur form is computed.
       \note (This is the version used for \c float and \c double, where \c select requires two arguments to represent a complex eigenvalue.)
    */
    void GEES(const char JOBVS, const char SORT, OrdinalType (*ptr2func)(ScalarType*, ScalarType*), const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, ScalarType* WR, ScalarType* WI, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, OrdinalType* BWORK, OrdinalType* info) const;

    /*! Computes for an \c n by \c n nonsymmetric matrix \c A, the eigenvalues, the Schur form \c T, and, optionally, the matrix of Schur vectors \c Z. When \c ScalarType is \c float or \c double, the real Schur form is computed.
       \note (This is the version used for \c std::complex<float> and \c std::complex<double>, where \c select requires one arguments to represent a complex eigenvalue.)
    */
    void GEES(const char JOBVS, const char SORT, OrdinalType (*ptr2func)(ScalarType*), const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, ScalarType* W, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* BWORK, OrdinalType* info) const;

    /*! Computes for an \c n by \c n nonsymmetric matrix \c A, the eigenvalues, the Schur form \c T, and, optionally, the matrix of Schur vectors \c Z. When \c ScalarType is \c float or \c double, the real Schur form is computed.
       \note (This is the version used for any \c ScalarType, when the user doesn't want to enable the sorting functionality.)
    */
    void GEES(const char JOBVS, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, MagnitudeType* WR, MagnitudeType* WI, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* BWORK, OrdinalType* info) const;

    //! Computes for an \c n by \c n real nonsymmetric matrix \c A, the eigenvalues and, optionally, the left and/or right eigenvectors.
    void GEEV(const char JOBVL, const char JOBVR, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* WR, ScalarType* WI, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /*! Computes for an \c n by \c n real nonsymmetric matrix \c A, the eigenvalues and, optionally, the left and/or right eigenvectors.
	Optionally, it can compute a balancing transformation to improve the conditioning of the eigenvalues and eigenvectors.
	\note (This is the function is only defined for \c ScalarType = \c float or \c double.)
    */
    void GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* WR, ScalarType* WI, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, OrdinalType* ilo, OrdinalType* ihi, MagnitudeType* SCALE, MagnitudeType* abnrm, MagnitudeType* RCONDE, MagnitudeType* RCONDV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* IWORK, OrdinalType* info) const;

    /*! Computes for a pair of \c n by \c n nonsymmetric matrices (\c A,\c B) the generalized eigenvalues, and optionally, the left and/or right generalized eigenvectors.
	Optionally, it can compute a balancing transformation to improve the conditioning of the eigenvalues and eigenvectors.
	\note (This is the function is only defined for \c ScalarType = \c float or \c double.)
    */
    void GGEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, MagnitudeType* ALPHAR, MagnitudeType* ALPHAI, ScalarType* BETA, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, OrdinalType* ilo, OrdinalType* ihi, MagnitudeType* LSCALE, MagnitudeType* RSCALE, MagnitudeType* abnrm, MagnitudeType* bbnrm, MagnitudeType* RCONDE, MagnitudeType* RCONDV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* IWORK, OrdinalType* BWORK, OrdinalType* info) const;

    /*! Computes for a pair of \c n by \c n nonsymmetric matrices (\c A,\c B) the generalized eigenvalues, and optionally, the left and/or right generalized eigenvectors.
       \note (This is the function is only defined for \c ScalarType = \c float or \c double.)
    */
    void GGEV(const char JOBVL, const char JOBVR, const OrdinalType n, ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb, MagnitudeType *ALPHAR, MagnitudeType *ALPHAI, ScalarType *BETA, ScalarType *VL, const OrdinalType ldvl, ScalarType *VR, const OrdinalType ldvr, ScalarType *WORK, const OrdinalType lwork, OrdinalType *info) const;


    /*! Reorders the real Schur factorization of a real matrix so that a selected cluster of eigenvalues appears in the leading diagonal blocks of the upper quasi-triangular matrix \c T, and the leading columns of \c Q form an orthonormal basis of the corresponding right invariant subspace.
       \note (This function is only defined for \c ScalarType = \c float or \c double.)
    */
  void TRSEN(const char JOB, const char COMPQ, const OrdinalType *SELECT, const OrdinalType n, ScalarType *T, const OrdinalType ldt, ScalarType *Q, const OrdinalType ldq, MagnitudeType *WR, MagnitudeType *WI, OrdinalType *M, ScalarType *S, MagnitudeType *SEP, ScalarType *WORK, const OrdinalType lwork, OrdinalType *IWORK, const OrdinalType liwork, OrdinalType *info ) const;


    /*! Reorders the generalized real Schur decomposition of a real matrix pair (\c A, \c B), so that a selected cluster of eigenvalues appears in the leading diagonal blocks of the upper quasi-triangular matrix \c A and the upper triangular \c B.
       \note (This function is only defined for \c ScalarType = \c float or \c double.)
    */
  void TGSEN(const OrdinalType ijob, const OrdinalType wantq, const OrdinalType wantz, const OrdinalType *SELECT, const OrdinalType n, ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb, MagnitudeType *ALPHAR, MagnitudeType *ALPHAI, MagnitudeType *BETA, ScalarType *Q, const OrdinalType ldq, ScalarType *Z, const OrdinalType ldz, OrdinalType *M, MagnitudeType *PL, MagnitudeType *PR, MagnitudeType *DIF, ScalarType *WORK, const OrdinalType lwork, OrdinalType *IWORK, const OrdinalType liwork, OrdinalType *info ) const;


    /*! Computes for a pair of \c n by \c n nonsymmetric matrices (\c A,\c B) the generalized eigenvalues, the generalized real Schur form (\c S,\c T), optionally, the left and/or right matrices of Schur vectors.
       \note (This is the function is only defined for \c ScalarType = \c float or \c double.)
    */
    void GGES(const char JOBVL, const char JOBVR, const char SORT, OrdinalType (*ptr2func)(ScalarType *, ScalarType *, ScalarType *), const OrdinalType n, ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb, OrdinalType *sdim, MagnitudeType *ALPHAR, MagnitudeType *ALPHAI, MagnitudeType *BETA, ScalarType *VL, const OrdinalType ldvl, ScalarType *VR, const OrdinalType ldvr, ScalarType *WORK, const OrdinalType lwork, OrdinalType *BWORK, OrdinalType *info ) const;

    //@}


    //! @name Singular Value Decompositon Routines
    //@{
    //! Computes the singular values (and optionally, vectors) of a real matrix \c A.
    void GESVD(const char JOBU, const char JOBVT, const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, MagnitudeType* S, ScalarType* U, const OrdinalType ldu, ScalarType* V, const OrdinalType ldv, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const;
    //@}


    //! @name Orthogonal matrix routines
    //@{

    /// Apply Householder reflectors (real case).
    ///
    /// Overwrite the general real \c m by \c n matrix \c C with the
    /// product of \c Q and \c C, whiere Q is the product of \c k
    /// elementary (Householder) reflectors as returned by \c GEQRF.
    ///
    /// \note This method is not defined when ScalarType is complex.
    /// Call \c UNMQR in that case.  ("OR" stands for "orthogonal";
    /// "UN" stands for "unitary.")
    void ORMQR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /// \brief Apply Householder reflectors (complex case).
    ///
    /// Overwrite the general complex \c m by \c n matrix \c C with
    /// the product of \c Q and \c C, where Q is the product of \c k
    /// elementary (Householder) reflectors as returned by \c GEQRF.
    ///
    /// \note This method will call \c ORMQR when ScalarType is real.
    /// (Unitary real matrices are orthogonal.)
    void UNMQR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /// \brief Compute explicit Q factor from QR factorization (GEQRF) (real case).
    ///
    /// Generate the \c m by \c n matrix Q with orthonormal columns
    /// corresponding to the first \c n columns of a product of \c k
    /// elementary reflectors of order \c m, as returned by \c GEQRF.
    ///
    /// \note This method is not defined when ScalarType is complex.
    /// Call \c UNGQR in that case.  ("OR" stands for "orthogonal";
    /// "UN" stands for "unitary.")
    void ORGQR(const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /// \brief Compute explicit QR factor from QR factorization (GEQRF) (complex case).
    ///
    /// Generate the \c m by \c n matrix Q with orthonormal columns
    /// corresponding tothe first \c n columns of a product of \c k
    /// elementary reflectors of order \c m, as returned by \c GEQRF.
    ///
    /// \note This method will call \c ORGQR when ScalarType is real.
    /// (Unitary real matrices are orthogonal.)
    void UNGQR(const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /*! \brief Generates a real orthogonal matrix \c Q which is the product of \c ihi-ilo elementary reflectors of order \c n, as returned by GEHRD.  On return \c Q is stored in \c A.
    \note This method is not defined when ScalarType is complex.
    */
    void ORGHR(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;

    /*! \brief Overwrites the general real \c m by \c n matrix \c C with the product of \c C and \c Q, which is a product of \c ihi-ilo elementary reflectors, as returned by GEHRD.
    \note This method is not defined when ScalarType is complex.
    */
    void ORMHR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const;
    //@}

    //! @name Triangular Matrix Routines
    //@{

    /*! Computes some or all of the right and/or left eigenvectors of an upper triangular matrix \c T. If ScalarType is \c float or \c double, then the matrix is quasi-triangular and arugments \c RWORK is ignored.
    */
    void TREVC(const char SIDE, const char HOWMNY, OrdinalType* select, const OrdinalType n, const ScalarType* T, const OrdinalType ldt, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, ScalarType* WORK, OrdinalType* info) const;

    /*! Computes some or all of the right and/or left eigenvectors of an upper triangular matrix \c T. If ScalarType is \c float or \c double, then the matrix is quasi-triangular and arugments \c RWORK is ignored.
       \note (This is the version used for any \c ScalarType, when the user doesn't want to enable the selecting functionality, with HOWMNY='A'.)
    */
    void TREVC(const char SIDE, const OrdinalType n, const ScalarType* T, const OrdinalType ldt, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, ScalarType* WORK, MagnitudeType* RWORK, OrdinalType* info) const;

    /*! Reorders the Schur factorization of a matrix \c T via unitary similarity transformations so that the diagonal element of \c T with row index \c ifst is moved to row \c ilst. If \c ScalarType is \c float or \c double, then \c T should be in real Schur form and the operation affects the diagonal block referenced by \c ifst.
      \note This method will ignore the WORK vector when ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    void TREXC(const char COMPQ, const OrdinalType n, ScalarType* T, const OrdinalType ldt, ScalarType* Q, const OrdinalType ldq, OrdinalType ifst, OrdinalType ilst, ScalarType* WORK, OrdinalType* info) const;

    /*! Computes some or all of the right and/or left eigenvectors of a pair of real matrices ( \c S, \c P ), where \c S is a quasi-triangular matrix and \c P is upper triangular.
       \note This method is only defined for \c ScalarType = \c float or \c double.
    */
    void TGEVC(const char SIDE, const char HOWMNY, const OrdinalType *SELECT, const OrdinalType n, ScalarType *S, const OrdinalType lds, ScalarType *P, const OrdinalType ldp, ScalarType *VL, const OrdinalType ldvl, ScalarType *VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType *M, ScalarType *WORK, OrdinalType *info) const;


    //@}

    //! @name Rotation/Reflection generators
    //@{

    //! Gnerates a plane rotation that zeros out the second component of the input vector.
    void LARTG( const ScalarType f, const ScalarType g, MagnitudeType* c, ScalarType* s, ScalarType* r ) const;

    //! Generates an elementary reflector of order \c n that zeros out the last \c n-1 components of the input vector.
    void LARFG( const OrdinalType n, ScalarType* alpha, ScalarType* x, const OrdinalType incx, ScalarType* tau ) const;

    //@}

    //! @name Matrix Balancing Routines
    //@{

    //! Balances a general matrix A, through similarity transformations to make the rows and columns as close in norm as possible.
    void GEBAL(const char JOBZ, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType ilo, OrdinalType ihi, MagnitudeType* scale, OrdinalType* info) const;

    //! Forms the left or right eigenvectors of a general matrix that has been balanced by GEBAL by backward transformation of the computed eigenvectors \c V.
    void GEBAK(const char JOBZ, const char SIDE, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const MagnitudeType* scale , const OrdinalType m, ScalarType* V, const OrdinalType ldv, OrdinalType* info) const;

    //@}

    //! @name Random number generators
    //@{
    //! Returns a random number from a uniform or normal distribution.
    ScalarType LARND( const OrdinalType idist, OrdinalType* seed ) const;

    //! Returns a vector of random numbers from a chosen distribution.
    void LARNV( const OrdinalType idist, OrdinalType* seed, const OrdinalType n, ScalarType* v ) const;
    //@}

    //! @name Machine Characteristics Routines.
    //@{
    /*! \brief Determines machine parameters for floating point characteristics.
	\note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    ScalarType LAMCH(const char CMACH) const;

    /*! \brief Chooses problem-dependent parameters for the local environment.
	\note This method should give parameters for good, but not optimal, performance on many currently
	available computers.
    */
    OrdinalType ILAENV( const OrdinalType ispec, const std::string& NAME, const std::string& OPTS, const OrdinalType N1 = -1, const OrdinalType N2 = -1, const OrdinalType N3 = -1, const OrdinalType N4 = -1 ) const;
    //@}

    //! @name Miscellaneous Utilities.
    //@{
    /*! \brief Computes x^2 + y^2 safely, to avoid overflow.
	\note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    ScalarType LAPY2(const ScalarType x, const ScalarType y) const;
    //@}
  };

  // END GENERAL TEMPLATE DECLARATION //

  // BEGIN GENERAL TEMPLATE IMPLEMENTATION //


  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::PTTRF(const OrdinalType n, ScalarType* d, ScalarType* e, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::PTTRS(const OrdinalType n, const OrdinalType nrhs, const ScalarType* d, const ScalarType* e, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POTRF(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POTRS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POTRI(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POCON(const char UPLO, const OrdinalType n, const ScalarType* A, const OrdinalType lda, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POSV(const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POEQU(const OrdinalType n, const ScalarType* A, const OrdinalType lda, MagnitudeType* S, MagnitudeType* scond, MagnitudeType* amax, OrdinalType* info) const
  {
    // Test the input parameters
    info = 0;
    if (n < 0) {
      info = -1;
    } else if (lda < TEUCHOS_MAX(1, n)) {
      info = -3;
    }
    if (info != 0) {
      return;
    }

    ScalarType sZero = ScalarTraits<ScalarType>::zero();
    ScalarType sOne  = ScalarTraits<ScalarType>::one();
    MagnitudeType mZero = ScalarTraits<ScalarType>::magnitude(sZero);
    MagnitudeType mOne = ScalarTraits<ScalarType>::magnitude(sOne);

    // Quick return
    if (n == 0) {
      *scond = mOne;
      *amax = mZero;
      return;
    }

    // Find the minimum and maximum diagonal elements
    S[0] = ScalarTraits<ScalarType>::magnitude( A[0] );
    MagnitudeType smin = S[0];
    amax = S[0];
    for (OrdinalType i=0; i<n; ++i) {
      S[i] = ScalarTraits<ScalarType>::magnitude( A[i*lda + i] );
      smin = TEUCHOS_MIN( smin, S[i] );
      amax = TEUCHOS_MAX( amax, S[i] );
    }

    if (smin < mZero) {
      // Find the first non-positve diagonal element and return an error code
      for (OrdinalType i=0; i<n; ++i) {
	if (S[i] < mZero)
	  *info = i;
      }
    } else {
      // Set the scale factors to the reciprocals of the diagonal elements
      for (OrdinalType i=0; i<n; ++i) {
	S[i] = mOne / ScalarTraits<ScalarType>::squareroot( S[i] );
      }
      // Compute scond = min(S(i)) / max(S(i))
      *scond = ScalarTraits<ScalarType>::squareroot( smin ) / ScalarTraits<ScalarType>::squareroot( amax );
    }
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::PORFS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::POSVX(const char FACT, const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* AF, const OrdinalType ldaf, char EQUED, ScalarType* S, ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* rcond, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GELS(const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GELSS(const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, MagnitudeType* S, const MagnitudeType rcond, OrdinalType* rank, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GELSS(const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* S, const ScalarType rcond, OrdinalType* rank, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GGLSE(const OrdinalType m, const OrdinalType n, const OrdinalType p, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* C, ScalarType* D, ScalarType* X, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GEQRF( const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GETRF(const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GETRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::LASCL(const char TYPE, const OrdinalType kl, const OrdinalType ku, const MagnitudeType cfrom, const MagnitudeType cto, const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const
  {
    MagnitudeType safeMin = ScalarTraits<ScalarType>::sfmin();
    ScalarType sZero = ScalarTraits<ScalarType>::zero();
    ScalarType sOne  = ScalarTraits<ScalarType>::one();
    MagnitudeType mZero = ScalarTraits<ScalarType>::magnitude(sZero);
    MagnitudeType mOne = ScalarTraits<ScalarType>::magnitude(sOne);

    MagnitudeType smlnum = ScalarTraits<ScalarType>::magnitude(safeMin);
    MagnitudeType bignum = ScalarTraits<ScalarType>::magnitude(sOne/smlnum);

    OrdinalType i, j;
    ScalarType* ptr;
    MagnitudeType mul;
    bool done = false;

    MagnitudeType cfromc = cfrom;
    MagnitudeType ctoc = cto;
    MagnitudeType cfrom1;
    MagnitudeType cto1;

    while (!done) {

      cfrom1 = cfromc*smlnum;
      if (cfrom1 == cfromc) {
	// cfromc is an inf. Multiply by a correctly signed zero for finite ctoc, or a NaN if ctoc is infinite.
	mul = ctoc / cfromc;
	done = true;
	cto1 = ctoc;
    } else {
	cto1 = ctoc / bignum;
	if (cto1 == ctoc) {
	  // ctoc is either 0 or an inf. In both cases, ctoc itself serves as the correct multiplication factor.
	  mul = ctoc;
	  done = true;
	  cfromc = mOne;
	} else if (ScalarTraits<ScalarType>::magnitude(cfrom1) > ScalarTraits<ScalarType>::magnitude(ctoc) && ctoc != mZero) {
	  mul = smlnum;
	  done = false;
	  cfromc = cfrom1;
	} else if (ScalarTraits<ScalarType>::magnitude(cto1) > ScalarTraits<ScalarType>::magnitude(cfromc)) {
	  mul = bignum;
	  done = false;
	  ctoc = cto1;
	} else {
	  mul = ctoc / cfromc;
	  done = true;
	}
      }

      for (j=0; j<n; j++) {
	ptr = A + j*lda;
	for (i=0; i<m; i++) { *ptr = mul * (*ptr); ptr++; }
      }
    }

  }

  template<typename OrdinalType, typename ScalarType>
  void
  LAPACK<OrdinalType, ScalarType>::
  LASWP (const OrdinalType N,
	 ScalarType A[],
	 const OrdinalType LDA,
	 const OrdinalType K1,
	 const OrdinalType K2,
	 const OrdinalType IPIV[],
	 const OrdinalType INCX) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GBTRF(const OrdinalType m, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GBTRS(const char TRANS, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GTTRF(const OrdinalType n, ScalarType* dl, ScalarType* d, ScalarType* du, ScalarType* du2, OrdinalType* IPIV, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GTTRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* dl, const ScalarType* d, const ScalarType* du, const ScalarType* du2, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GETRI(const OrdinalType n, ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void
  LAPACK<OrdinalType,ScalarType>::
  LATRS (const char UPLO,
	 const char TRANS,
	 const char DIAG,
	 const char NORMIN,
	 const OrdinalType N,
	 ScalarType* A,
	 const OrdinalType LDA,
	 ScalarType* X,
	 MagnitudeType* SCALE,
	 MagnitudeType* CNORM,
	 OrdinalType* INFO) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GECON(const char NORM, const OrdinalType n, const ScalarType* A, const OrdinalType lda, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GBCON(const char NORM, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  typename ScalarTraits<ScalarType>::magnitudeType LAPACK<OrdinalType,ScalarType>::LANGB(const char NORM, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const ScalarType* A, const OrdinalType lda, MagnitudeType* WORK) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GESV(const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GEEQU(const OrdinalType m, const OrdinalType n, const ScalarType* A, const OrdinalType lda, ScalarType* R, ScalarType* C, ScalarType* rowcond, ScalarType* colcond, ScalarType* amax, OrdinalType* info) const
  {

    // Test the input parameters
    info = 0;
    if (m < 0) {
      info = -1;
    } else if (n < 0) {
      info = -2;
    } else if (lda < TEUCHOS_MAX(1, m)) {
      info = -4;
    }
    if (info != 0) {
      return;
    }

    ScalarType sZero = ScalarTraits<ScalarType>::zero();
    ScalarType sOne  = ScalarTraits<ScalarType>::one();
    MagnitudeType mZero = ScalarTraits<ScalarType>::magnitude(sZero);
    MagnitudeType mOne = ScalarTraits<ScalarType>::magnitude(sOne);

    // Quick return
    if (m == 0 || n == 0) {
      *rowcond = mOne;
      *colcond = mOne;
      *amax = mZero;
      return;
    }

    MagnitudeType safeMin = ScalarTraits<ScalarType>::sfmin();
    MagnitudeType smlnum = ScalarTraits<ScalarType>::magnitude(safeMin);
    MagnitudeType bignum = ScalarTraits<ScalarType>::magnitude(sOne/smlnum);

    // Compute the row scale factors
    for (OrdinalType i=0; i<m; i++) {
      R[i] = mZero;
    }

    // Find the maximum element in each row
    for (OrdinalType j=0; j<n; j++) {
      for (OrdinalType i=0; i<m; i++) {
	R[i] = TEUCHOS_MAX( R[i], ScalarTraits<ScalarType>::magnitude( A[j*lda + i] ) );
      }
    }

    // Find the maximum and minimum scale factors
    MagnitudeType rcmin = bignum;
    MagnitudeType rcmax = mZero;
    for (OrdinalType i=0; i<m; i++) {
      rcmax = TEUCHOS_MAX( rcmax, R[i] );
      rcmin = TEUCHOS_MIN( rcmin, R[i] );
    }
    *amax = rcmax;

    if (rcmin == mZero) {
      // Find the first zero scale factor and return an error code
      for (OrdinalType i=0; i<m; i++) {
	if (R[i] == mZero)
	  *info = i;
      }
    } else {
      // Invert the scale factors
      for (OrdinalType i=0; i<m; i++) {
	R[i] = mOne / TEUCHOS_MIN( TEUCHOS_MAX( R[i], smlnum ), bignum );
      }
      // Compute rowcond = min(R(i)) / max(R(i))
      *rowcond = TEUCHOS_MAX( rcmin, smlnum ) / TEUCHOS_MIN( rcmax, bignum );
    }

    // Compute the column scale factors
    for (OrdinalType j=0; j<n; j++) {
      C[j] = mZero;
    }

    // Find the maximum element in each column, assuming the row scaling computed above
    for (OrdinalType j=0; j<n; j++) {
      for (OrdinalType i=0; i<m; i++) {
	C[j] = TEUCHOS_MAX( C[j], R[i]*ScalarTraits<ScalarType>::magnitude( A[j*lda + i] ) );
      }
    }

    // Find the maximum and minimum scale factors
    rcmin = bignum;
    rcmax = mZero;
    for (OrdinalType j=0; j<n; j++) {
      rcmax = TEUCHOS_MAX( rcmax, C[j] );
      rcmin = TEUCHOS_MIN( rcmin, C[j] );
    }

    if (rcmin == mZero) {
      // Find the first zero scale factor and return an error code
      for (OrdinalType j=0; j<n; j++) {
	if (C[j] == mZero)
	  *info = m+j;
      }
    } else {
      // Invert the scale factors
      for (OrdinalType j=0; j<n; j++) {
	C[j] = mOne / TEUCHOS_MIN( TEUCHOS_MAX( C[j], smlnum ), bignum );
      }
      // Compute colcond = min(C(j)) / max(C(j))
      *colcond = TEUCHOS_MAX( rcmin, smlnum ) / TEUCHOS_MIN( rcmax, bignum );
    }
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GERFS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GBEQU(const OrdinalType m, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const ScalarType* A, const OrdinalType lda, MagnitudeType* R, MagnitudeType* C, MagnitudeType* rowcond, MagnitudeType* colcond, MagnitudeType* amax, OrdinalType* info) const
  {

    // Test the input parameters
    info = 0;
    if (m < 0) {
      info = -1;
    } else if (n < 0) {
      info = -2;
    } else if (kl < 0) {
      info = -3;
    } else if (ku < 0) {
      info = -4;
    } else if (lda < kl+ku+1) {
      info = -6;
    }
    if (info != 0) {
      return;
    }

    ScalarType sZero = ScalarTraits<ScalarType>::zero();
    ScalarType sOne  = ScalarTraits<ScalarType>::one();
    MagnitudeType mZero = ScalarTraits<ScalarType>::magnitude(sZero);
    MagnitudeType mOne = ScalarTraits<ScalarType>::magnitude(sOne);

    // Quick return
    if (m == 0 || n == 0) {
      *rowcond = mOne;
      *colcond = mOne;
      *amax = mZero;
      return;
    }

    MagnitudeType safeMin = ScalarTraits<ScalarType>::sfmin();
    MagnitudeType smlnum = ScalarTraits<ScalarType>::magnitude(safeMin);
    MagnitudeType bignum = ScalarTraits<ScalarType>::magnitude(sOne/smlnum);

    // Compute the row scale factors
    for (OrdinalType i=0; i<m; i++) {
      R[i] = mZero;
    }

    // Find the maximum element in each row
    for (OrdinalType j=0; j<n; j++) {
      for (OrdinalType i=TEUCHOS_MAX(j-ku,0); i<TEUCHOS_MIN(j+kl,m-1); i++) {
	R[i] = TEUCHOS_MAX( R[i], ScalarTraits<ScalarType>::magnitude( A[j*lda + ku+i-j] ) );
      }
    }

    // Find the maximum and minimum scale factors
    MagnitudeType rcmin = bignum;
    MagnitudeType rcmax = mZero;
    for (OrdinalType i=0; i<m; i++) {
      rcmax = TEUCHOS_MAX( rcmax, R[i] );
      rcmin = TEUCHOS_MIN( rcmin, R[i] );
    }
    *amax = rcmax;

    if (rcmin == mZero) {
      // Find the first zero scale factor and return an error code
      for (OrdinalType i=0; i<m; i++) {
	if (R[i] == mZero)
	  *info = i;
      }
    } else {
      // Invert the scale factors
      for (OrdinalType i=0; i<m; i++) {
	R[i] = mOne / TEUCHOS_MIN( TEUCHOS_MAX( R[i], smlnum ), bignum );
      }
      // Compute rowcond = min(R(i)) / max(R(i))
      *rowcond = TEUCHOS_MAX( rcmin, smlnum ) / TEUCHOS_MIN( rcmax, bignum );
    }

    // Compute the column scale factors
    for (OrdinalType j=0; j<n; j++) {
      C[j] = mZero;
    }

    // Find the maximum element in each column, assuming the row scaling computed above
    for (OrdinalType j=0; j<n; j++) {
      for (OrdinalType i=TEUCHOS_MAX(j-ku,0); i<TEUCHOS_MIN(j+kl,m-1); i++) {
	C[j] = TEUCHOS_MAX( C[j], R[i]*ScalarTraits<ScalarType>::magnitude( A[j*lda + ku+i-j] ) );
      }
    }

    // Find the maximum and minimum scale factors
    rcmin = bignum;
    rcmax = mZero;
    for (OrdinalType j=0; j<n; j++) {
      rcmax = TEUCHOS_MAX( rcmax, C[j] );
      rcmin = TEUCHOS_MIN( rcmin, C[j] );
    }

    if (rcmin == mZero) {
      // Find the first zero scale factor and return an error code
      for (OrdinalType j=0; j<n; j++) {
	if (C[j] == mZero)
	  *info = m+j;
      }
    } else {
      // Invert the scale factors
      for (OrdinalType j=0; j<n; j++) {
	C[j] = mOne / TEUCHOS_MIN( TEUCHOS_MAX( C[j], smlnum ), bignum );
      }
      // Compute colcond = min(C(j)) / max(C(j))
      *colcond = TEUCHOS_MAX( rcmin, smlnum ) / TEUCHOS_MIN( rcmax, bignum );
    }
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GBRFS(const char TRANS, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GESVX(const char FACT, const char TRANS, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* AF, const OrdinalType ldaf, OrdinalType* IPIV, char EQUED, ScalarType* R, ScalarType* C, ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* rcond, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::SYTRD(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* D, ScalarType* E, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::GEHRD(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* A, const OrdinalType lda, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::TRTRS(const char UPLO, const char TRANS, const char DIAG, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::TRTRI(const char UPLO, const char DIAG, const OrdinalType n, const ScalarType* A, const OrdinalType lda, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::SPEV(const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* AP, ScalarType* W, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::SYEV(const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* W, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::SYGV(const OrdinalType itype, const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* W, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::HEEV(const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, MagnitudeType* W, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::HEGV(const OrdinalType itype, const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, MagnitudeType* W, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::STEQR(const char COMPZ, const OrdinalType n, ScalarType* D, ScalarType* E, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::HSEQR(const char JOB, const char COMPZ, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* H, const OrdinalType ldh, ScalarType* WR, ScalarType* WI, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEES(const char JOBVS, const char SORT, OrdinalType (*ptr2func)(ScalarType*, ScalarType*), const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, ScalarType* WR, ScalarType* WI, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, OrdinalType* BWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEES(const char JOBVS, const char SORT, OrdinalType (*ptr2func)(ScalarType*), const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, ScalarType* W, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, MagnitudeType *RWORK, OrdinalType* BWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEES(const char JOBVS, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, MagnitudeType* WR, MagnitudeType* WI, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, MagnitudeType *RWORK, OrdinalType* BWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEEV(const char JOBVL, const char JOBVR, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* WR, ScalarType* WI, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* WR, ScalarType* WI, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, OrdinalType* ilo, OrdinalType* ihi, MagnitudeType* SCALE, MagnitudeType* abnrm, MagnitudeType* RCONDE, MagnitudeType* RCONDV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* IWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GESVD(const char JOBU, const char JOBVT, const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, MagnitudeType* S, ScalarType* U, const OrdinalType ldu, ScalarType* V, const OrdinalType ldv, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GGEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, MagnitudeType* ALPHAR, MagnitudeType* ALPHAI, ScalarType* BETA, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, OrdinalType* ilo, OrdinalType* ihi, MagnitudeType* LSCALE, MagnitudeType* RSCALE, MagnitudeType* abnrm, MagnitudeType* bbnrm, MagnitudeType* RCONDE, MagnitudeType* RCONDV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* IWORK, OrdinalType* BWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GGEV(const char JOBVL, const char JOBVR, const OrdinalType n, ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb, MagnitudeType *ALPHAR, MagnitudeType *ALPHAI, ScalarType *BETA, ScalarType *VL, const OrdinalType ldvl, ScalarType *VR, const OrdinalType ldvr, ScalarType *WORK, const OrdinalType lwork, OrdinalType *info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }


  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::TRSEN(const char JOB, const char COMPQ, const OrdinalType *SELECT, const OrdinalType n, ScalarType *T, const OrdinalType ldt, ScalarType *Q, const OrdinalType ldq, MagnitudeType *WR, MagnitudeType *WI, OrdinalType *M, ScalarType *S, MagnitudeType *SEP, ScalarType *WORK, const OrdinalType lwork, OrdinalType *IWORK, const OrdinalType liwork, OrdinalType *info ) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }


  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType,ScalarType>::TGSEN(const OrdinalType ijob, const OrdinalType wantq, const OrdinalType wantz, const OrdinalType *SELECT, const OrdinalType n, ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb, MagnitudeType *ALPHAR, MagnitudeType *ALPHAI, MagnitudeType *BETA, ScalarType *Q, const OrdinalType ldq, ScalarType *Z, const OrdinalType ldz, OrdinalType *M, MagnitudeType *PL, MagnitudeType *PR, MagnitudeType *DIF, ScalarType *WORK, const OrdinalType lwork, OrdinalType *IWORK, const OrdinalType liwork, OrdinalType *info ) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }


  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GGES(const char JOBVL, const char JOBVR, const char SORT, OrdinalType (*ptr2func)(ScalarType*, ScalarType*, ScalarType*), const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* sdim, MagnitudeType* ALPHAR, MagnitudeType* ALPHAI, MagnitudeType* BETA, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, ScalarType* WORK, const OrdinalType lwork, OrdinalType *BWORK, OrdinalType* info ) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::ORMQR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::UNMQR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::ORGQR(const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::UNGQR(const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::ORGHR(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::ORMHR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::TREVC(const char SIDE, const char HOWMNY, OrdinalType* select, const OrdinalType n, const ScalarType* T, const OrdinalType ldt, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, ScalarType* WORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::TREVC(const char SIDE, const OrdinalType n, const ScalarType* T, const OrdinalType ldt, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, ScalarType* WORK, MagnitudeType* RWORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::TREXC(const char COMPQ, const OrdinalType n, ScalarType* T, const OrdinalType ldt, ScalarType* Q, const OrdinalType ldq, OrdinalType ifst, OrdinalType ilst, ScalarType* WORK, OrdinalType* info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }


  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::TGEVC(const char SIDE, const char HOWMNY, const OrdinalType *SELECT, const OrdinalType n, ScalarType *S, const OrdinalType lds, ScalarType *P, const OrdinalType ldp, ScalarType *VL, const OrdinalType ldvl, ScalarType *VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType *M, ScalarType *WORK, OrdinalType *info) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }


  template<typename OrdinalType, typename ScalarType>
  ScalarType LAPACK<OrdinalType, ScalarType>::LAMCH(const char CMACH) const
  {
    return UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  OrdinalType LAPACK<OrdinalType, ScalarType>::ILAENV( const OrdinalType ispec, const std::string& NAME, const std::string& OPTS, const OrdinalType N1, const OrdinalType N2, const OrdinalType N3, const OrdinalType N4 ) const
  {
    return UndefinedLAPACKRoutine<OrdinalType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType LAPACK<OrdinalType, ScalarType>::LAPY2(const ScalarType x, const ScalarType y) const
  {
    return UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::LARTG( const ScalarType f, const ScalarType g, MagnitudeType* c, ScalarType* s, ScalarType* r ) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::LARFG( const OrdinalType n, ScalarType* alpha, ScalarType* x, const OrdinalType incx, ScalarType* tau ) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEBAL( const char JOBZ, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType ilo, OrdinalType ihi, MagnitudeType* scale, OrdinalType* info ) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }


  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::GEBAK( const char JOBZ, const char SIDE, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const MagnitudeType* scale, const OrdinalType m, ScalarType* V, const OrdinalType ldv, OrdinalType* info ) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  ScalarType LAPACK<OrdinalType, ScalarType>::LARND( const OrdinalType idist, OrdinalType* seed ) const
  {
    return UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  template<typename OrdinalType, typename ScalarType>
  void LAPACK<OrdinalType, ScalarType>::LARNV( const OrdinalType idist, OrdinalType* seed, const OrdinalType n, ScalarType* v ) const
  {
    UndefinedLAPACKRoutine<ScalarType>::notDefined();
  }

  // END GENERAL TEMPLATE IMPLEMENTATION //

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  // BEGIN INT, FLOAT SPECIALIZATION DECLARATION //

  template<>
  class TEUCHOSNUMERICS_LIB_DLL_EXPORT LAPACK<int, float>
  {
  public:
    inline LAPACK(void) {}
    inline LAPACK(const LAPACK<int, float>& lapack) {}
    inline virtual ~LAPACK(void) {}

    // Symmetric positive definite linear system routines
    void POTRF(const char UPLO, const int n, float* A, const int lda, int * info) const;
    void POTRS(const char UPLO, const int n, const int nrhs, const float* A, const int lda, float* B, const int ldb, int* info) const;
    void PTTRF(const int n, float* d, float* e, int* info) const;
    void PTTRS(const int n, const int nrhs, const float* d, const float* e, float* B, const int ldb, int* info) const;
    void POTRI(const char UPLO, const int n, float* A, const int lda, int* info) const;
    void POCON(const char UPLO, const int n, const float* A, const int lda, const float anorm, float* rcond, float* WORK, int* IWORK, int* info) const;
    void POSV(const char UPLO, const int n, const int nrhs, float* A, const int lda, float* B, const int ldb, int* info) const;
    void POEQU(const int n, const float* A, const int lda, float* S, float* scond, float* amax, int* info) const;
    void PORFS(const char UPLO, const int n, const int nrhs, float* A, const int lda, const float* AF, const int ldaf, const float* B, const int ldb, float* X, const int ldx, float* FERR, float* BERR, float* WORK, int* IWORK, int* info) const;
    void POSVX(const char FACT, const char UPLO, const int n, const int nrhs, float* A, const int lda, float* AF, const int ldaf, char EQUED, float* S, float* B, const int ldb, float* X, const int ldx, float* rcond, float* FERR, float* BERR, float* WORK, int* IWORK, int* info) const;

    // General Linear System Routines
    void GELS(const char TRANS, const int m, const int n, const int nrhs, float* A, const int lda, float* B, const int ldb, float* WORK, const int lwork, int* info) const;
    void GELSS(const int m, const int n, const int nrhs, float* A, const int lda, float* B, const int ldb, float* S, const float rcond, int* rank, float* WORK, const int lwork, float* RWORK, int* info) const;
    void GELSS(const int m, const int n, const int nrhs, float* A, const int lda, float* B, const int ldb, float* S, const float rcond, int* rank, float* WORK, const int lwork, int* info) const;
    void GGLSE(const int m, const int n, const int p, float* A, const int lda, float* B, const int ldb, float* C, float* D, float* X, float* WORK, const int lwork, int* info) const;
    void GEQRF( const int m, const int n, float* A, const int lda, float* TAU, float* WORK, const int lwork, int* info) const;
    void GETRF(const int m, const int n, float* A, const int lda, int* IPIV, int* info) const;
    void GETRS(const char TRANS, const int n, const int nrhs, const float* A, const int lda, const int* IPIV, float* B, const int ldb, int* info) const;
    void LASCL(const char TYPE, const int kl, const int ku, const float cfrom, const float cto, const int m, const int n, float* A, const int lda, int* info) const;


    void LASWP (const int N,
		float A[],
		const int LDA,
		const int K1,
		const int K2,
		const int IPIV[],
		const int INCX) const;

    void GBTRF(const int m, const int n, const int kl, const int ku, float* A, const int lda, int* IPIV, int* info) const;
    void GBTRS(const char TRANS, const int n, const int kl, const int ku, const int nrhs, const float* A, const int lda, const int* IPIV, float* B, const int ldb, int* info) const;
    void GTTRF(const int n, float* dl, float* d, float* du, float* du2, int* IPIV, int* info) const;
    void GTTRS(const char TRANS, const int n, const int nrhs, const float* dl, const float* d, const float* du, const float* du2, const int* IPIV, float* B, const int ldb, int* info) const;


    void GETRI(const int n, float* A, const int lda, const int* IPIV, float* WORK, const int lwork, int* info) const;
    void LATRS (const char UPLO, const char TRANS, const char DIAG, const char NORMIN, const int N, float* A, const int LDA, float* X, float* SCALE, float* CNORM, int* INFO) const;
    void GECON(const char NORM, const int n, const float* A, const int lda, const float anorm, float* rcond, float* WORK, int* IWORK, int* info) const;
    void GBCON(const char NORM, const int n, const int kl, const int ku, const float* A, const int lda, int* IPIV, const float anorm, float* rcond, float* WORK, int* IWORK, int* info) const;
    float LANGB(const char NORM, const int n, const int kl, const int ku, const float* A, const int lda, float* WORK) const;
    void GESV(const int n, const int nrhs, float* A, const int lda, int* IPIV, float* B, const int ldb, int* info) const;
    void GEEQU(const int m, const int n, const float* A, const int lda, float* R, float* C, float* rowcond, float* colcond, float* amax, int* info) const;
    void GERFS(const char TRANS, const int n, const int nrhs, const float* A, const int lda, const float* AF, const int ldaf, const int* IPIV, const float* B, const int ldb, float* X, const int ldx, float* FERR, float* BERR, float* WORK, int* IWORK, int* info) const;
    void GBEQU(const int m, const int n, const int kl, const int ku, const float* A, const int lda, float* R, float* C, float* rowcond, float* colcond, float* amax, int* info) const;
    void GBRFS(const char TRANS, const int n, const int kl, const int ku, const int nrhs, const float* A, const int lda, const float* AF, const int ldaf, const int* IPIV, const float* B, const int ldb, float* X, const int ldx, float* FERR, float* BERR, float* WORK, int* IWORK, int* info) const;
    void GESVX(const char FACT, const char TRANS, const int n, const int nrhs, float* A, const int lda, float* AF, const int ldaf, int* IPIV, char EQUED, float* R, float* C, float* B, const int ldb, float* X, const int ldx, float* rcond, float* FERR, float* BERR, float* WORK, int* IWORK, int* info) const;
    void SYTRD(const char UPLO, const int n, float* A, const int lda, float* D, float* E, float* TAU, float* WORK, const int lwork, int* info) const;
    void GEHRD(const int n, const int ilo, const int ihi, float* A, const int lda, float* TAU, float* WORK, const int lwork, int* info) const;
    void TRTRS(const char UPLO, const char TRANS, const char DIAG, const int n, const int nrhs, const float* A, const int lda, float* B, const int ldb, int* info) const;
    void TRTRI(const char UPLO, const char DIAG, const int n, const float* A, const int lda, int* info) const;

    // Symmetric eigenvalue routines.
    void SPEV(const char JOBZ, const char UPLO, const int n, float* AP, float* W, float* Z, const int ldz, float* WORK, int* info) const;
    void SYEV(const char JOBZ, const char UPLO, const int n, float* A, const int lda, float* W, float* WORK, const int lwork, int* info) const;
    void SYGV(const int itype, const char JOBZ, const char UPLO, const int n, float* A, const int lda, float* B, const int ldb, float* W, float* WORK, const int lwork, int* info) const;
    void HEEV(const char JOBZ, const char UPLO, const int n, float* A, const int lda, float* W, float* WORK, const int lwork, float* RWORK, int* info) const;
    void HEGV(const int itype, const char JOBZ, const char UPLO, const int n, float* A, const int lda, float* B, const int ldb, float* W, float* WORK, const int lwork, float *RWORK, int* info) const;
    void STEQR(const char COMPZ, const int n, float* D, float* E, float* Z, const int ldz, float* WORK, int* info) const;

    // Non-Hermitian eigenvalue routines.
    void HSEQR(const char JOB, const char COMPZ, const int n, const int ilo, const int ihi, float* H, const int ldh, float* WR, float* WI, float* Z, const int ldz, float* WORK, const int lwork, int* info) const;
    void GEES(const char JOBVS, const char SORT, int (*ptr2func)(float*, float*), const int n, float* A, const int lda, int* sdim, float* WR, float* WI, float* VS, const int ldvs, float* WORK, const int lwork, int* BWORK, int* info) const;
    void GEES(const char JOBVS, const int n, float* A, const int lda, int* sdim, float* WR, float* WI, float* VS, const int ldvs, float* WORK, const int lwork, float* RWORK, int* BWORK, int* info) const;
    void GEEV(const char JOBVL, const char JOBVR, const int n, float* A, const int lda, float* WR, float* WI, float* VL, const int ldvl, float* VR, const int ldvr, float* WORK, const int lwork, int* info) const;
    void GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const int n, float* A, const int lda, float* WR, float* WI, float* VL, const int ldvl, float* VR, const int ldvr, int* ilo, int* ihi, float* SCALE, float* abnrm, float* RCONDE, float* RCONDV, float* WORK, const int lwork, int* IWORK, int* info) const;
    void GGEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const int n, float* A, const int lda, float* B, const int ldb, float* ALPHAR, float* ALPHAI, float* BETA, float* VL, const int ldvl, float* VR, const int ldvr, int* ilo, int* ihi, float* LSCALE, float* RSCALE, float* abnrm, float* bbnrm, float* RCONDE, float* RCONDV, float* WORK, const int lwork, int* IWORK, int* BWORK, int* info) const;
    void GGEV(const char JOBVL, const char JOBVR, const int n, float *A, const int lda, float *B, const int ldb, float *ALPHAR, float *ALPHAI, float *BETA, float *VL, const int ldvl, float *VR, const int ldvr, float *WORK, const int lwork, int *info) const;
    void TRSEN(const char JOB, const char COMPQ, const int *SELECT, const int n, float *T, const int ldt, float *Q, const int ldq, float *WR, float *WI, int *M, float *S, float *SEP, float *WORK, const int lwork, int *IWORK, const int liwork, int *info ) const;
    void TGSEN(const int ijob, const int wantq, const int wantz, const int *SELECT, const int n, float *A, const int lda, float *B, const int ldb, float *ALPHAR, float *ALPHAI, float *BETA, float *Q, const int ldq, float *Z, const int ldz, int *M, float *PL, float *PR, float *DIF, float *WORK, const int lwork, int *IWORK, const int liwork, int *info ) const;
    void GGES(const char JOBVL, const char JOBVR, const char SORT, int (*ptr2func)(float*, float*, float*), const int n, float* A, const int lda, float* B, const int ldb, int* sdim, float* ALPHAR, float* ALPHAI, float* BETA, float* VL, const int ldvl, float* VR, const int ldvr, float* WORK, const int lwork, int *bwork, int* info ) const;

    // SVD routine
    void GESVD(const char JOBU, const char JOBVT, const int m, const int n, float* A, const int lda, float* S, float* U, const int ldu, float* V, const int ldv, float* WORK, const int lwork, float* RWORK, int* info) const;

    // Orthogonal matrix routines.
    void ORMQR(const char SIDE, const char TRANS, const int m, const int n, const int k, float* A, const int lda, const float* TAU, float* C, const int ldc, float* WORK, const int lwork, int* info) const;
    void UNMQR(const char SIDE, const char TRANS, const int m, const int n, const int k, float* A, const int lda, const float* TAU, float* C, const int ldc, float* WORK, const int lwork, int* info) const;

    void ORGQR(const int m, const int n, const int k, float* A, const int lda, const float* TAU, float* WORK, const int lwork, int* info) const;
    void UNGQR(const int m, const int n, const int k, float* A, const int lda, const float* TAU, float* WORK, const int lwork, int* info) const;
    void ORGHR(const int n, const int ilo, const int ihi, float* A, const int lda, const float* TAU, float* WORK, const int lwork, int* info) const;
    void ORMHR(const char SIDE, const char TRANS, const int m, const int n, const int ilo, const int ihi, const float* A, const int lda, const float* TAU, float* C, const int ldc, float* WORK, const int lwork, int* info) const;

    // Triangular matrix routines.
    void TREVC(const char SIDE, const char HOWMNY, int* select, const int n, const float* T, const int ldt, float* VL, const int ldvl, float* VR, const int ldvr, const int mm, int* m, float* WORK, int* info) const;
    void TREVC(const char SIDE, const int n, const float* T, const int ldt, float* VL, const int ldvl, float* VR, const int ldvr, const int mm, int* m, float* WORK, float *RWORK, int* info) const;
    void TREXC(const char COMPQ, const int n, float* T, const int ldt, float* Q, const int ldq, int ifst, int ilst, float* WORK, int* info) const;
    void TGEVC(const char SIDE, const char HOWMNY, const int *SELECT, const int n, float *S, const int lds, float *P, const int ldp, float *VL, const int ldvl, float *VR, const int ldvr, const int mm, int *M, float *WORK, int *info) const;

    // Rotation/reflection generators
    void LARTG( const float f, const float g, float* c, float* s, float* r ) const;
    void LARFG( const int n, float* alpha, float* x, const int incx, float* tau ) const;

    // Matrix balancing routines.
    void GEBAL(const char JOBZ, const int n, float* A, const int lda, int ilo, int ihi, float* scale, int* info) const;
    void GEBAK(const char JOBZ, const char SIDE, const int n, const int ilo, const int ihi, const float* scale, const int m, float* V, const int ldv, int* info) const;

    // Random number generators
    float LARND( const int idist, int* seed ) const;
    void LARNV( const int idist, int* seed, const int n, float* v ) const;

    // Machine characteristics.
    float LAMCH(const char CMACH) const;
    int ILAENV( const int ispec, const std::string& NAME, const std::string& OPTS, const int N1 = -1, const int N2 = -1, const int N3 = -1, const int N4 = -1 ) const;

    // Miscellaneous routines.
    float LAPY2(const float x, const float y) const;

  };

  // END INT, FLOAT SPECIALIZATION DECLARATION //

  // BEGIN INT, DOUBLE SPECIALIZATION DECLARATION //

  template<>
  class TEUCHOSNUMERICS_LIB_DLL_EXPORT LAPACK<int, double>
  {
  public:
    inline LAPACK(void) {}
    inline LAPACK(const LAPACK<int, double>& lapack) {}
    inline virtual ~LAPACK(void) {}

    // Symmetric positive definite linear system routines
    void PTTRF(const int n, double* d, double* e, int* info) const;
    void PTTRS(const int n, const int nrhs, const double* d, const double* e, double* B, const int ldb, int* info) const;
    void POTRF(const char UPLO, const int n, double* A, const int lda, int* info) const;
    void POTRS(const char UPLO, const int n, const int nrhs, const double* A, const int lda, double* B, const int ldb, int* info) const;
    void POTRI(const char UPLO, const int n, double* A, const int lda, int* info) const;
    void POCON(const char UPLO, const int n, const double* A, const int lda, const double anorm, double* rcond, double* WORK, int* IWORK, int* info) const;
    void POSV(const char UPLO, const int n, const int nrhs, double* A, const int lda, double* B, const int ldb, int* info) const;
    void POEQU(const int n, const double* A, const int lda, double* S, double* scond, double* amax, int* info) const;
    void PORFS(const char UPLO, const int n, const int nrhs, double* A, const int lda, const double* AF, const int ldaf, const double* B, const int ldb, double* X, const int ldx, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const;
    void POSVX(const char FACT, const char UPLO, const int n, const int nrhs, double* A, const int lda, double* AF, const int ldaf, char EQUED, double* S, double* B, const int ldb, double* X, const int ldx, double* rcond, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const;

    // General linear system routines
    void GELS(const char TRANS, const int m, const int n, const int nrhs, double* A, const int lda, double* B, const int ldb, double* WORK, const int lwork, int* info) const;
    void GELSS(const int m, const int n, const int nrhs, double* A, const int lda, double* B, const int ldb, double* S, const double rcond, int* rank, double* WORK, const int lwork, double* RWORK, int* info) const;
    void GELSS(const int m, const int n, const int nrhs, double* A, const int lda, double* B, const int ldb, double* S, const double rcond, int* rank, double* WORK, const int lwork, int* info) const;
    void GGLSE(const int m, const int n, const int p, double* A, const int lda, double* B, const int ldb, double* C, double* D, double* X, double* WORK, const int lwork, int* info) const;
    void GEQRF( const int m, const int n, double* A, const int lda, double* TAU, double* WORK, const int lwork, int* info) const;
    void GETRF(const int m, const int n, double* A, const int lda, int* IPIV, int* info) const;
    void GETRS(const char TRANS, const int n, const int nrhs, const double* A, const int lda, const int* IPIV, double* B, const int ldb, int* info) const;
    void LASCL(const char TYPE, const int kl, const int ku, const double cfrom, const double cto, const int m, const int n, double* A, const int lda, int* info) const;

    void LASWP (const int N,
		double A[],
		const int LDA,
		const int K1,
		const int K2,
		const int IPIV[],
		const int INCX) const;

    void GBTRF(const int m, const int n, const int kl, const int ku, double* A, const int lda, int* IPIV, int* info) const;
    void GBTRS(const char TRANS, const int n, const int kl, const int ku, const int nrhs, const double* A, const int lda, const int* IPIV, double* B, const int ldb, int* info) const;
    void GTTRF(const int n, double* dl, double* d, double* du, double* du2, int* IPIV, int* info) const;
    void GTTRS(const char TRANS, const int n, const int nrhs, const double* dl, const double* d, const double* du, const double* du2, const int* IPIV, double* B, const int ldb, int* info) const;
    void GETRI(const int n, double* A, const int lda, const int* IPIV, double* WORK, const int lwork, int* info) const;
    void LATRS (const char UPLO, const char TRANS, const char DIAG, const char NORMIN, const int N, double* A, const int LDA, double* X, double* SCALE, double* CNORM, int* INFO) const;
    void GECON(const char NORM, const int n, const double* A, const int lda, const double anorm, double* rcond, double* WORK, int* IWORK, int* info) const;
    void GBCON(const char NORM, const int n, const int kl, const int ku, const double* A, const int lda, int* IPIV, const double anorm, double* rcond, double* WORK, int* IWORK, int* info) const;
    double LANGB(const char NORM, const int n, const int kl, const int ku, const double* A, const int lda, double* WORK) const;
    void GESV(const int n, const int nrhs, double* A, const int lda, int* IPIV, double* B, const int ldb, int* info) const;
    void GEEQU(const int m, const int n, const double* A, const int lda, double* R, double* C, double* rowcond, double* colcond, double* amax, int* info) const;
    void GERFS(const char TRANS, const int n, const int nrhs, const double* A, const int lda, const double* AF, const int ldaf, const int* IPIV, const double* B, const int ldb, double* X, const int ldx, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const;
    void GBEQU(const int m, const int n, const int kl, const int ku, const double* A, const int lda, double* R, double* C, double* rowcond, double* colcond, double* amax, int* info) const;
    void GBRFS(const char TRANS, const int n, const int kl, const int ku, const int nrhs, const double* A, const int lda, const double* AF, const int ldaf, const int* IPIV, const double* B, const int ldb, double* X, const int ldx, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const;
    void GESVX(const char FACT, const char TRANS, const int n, const int nrhs, double* A, const int lda, double* AF, const int ldaf, int* IPIV, char EQUED, double* R, double* C, double* B, const int ldb, double* X, const int ldx, double* rcond, double* FERR, double* BERR, double* WORK, int* IWORK, int* info) const;
    void SYTRD(const char UPLO, const int n, double* A, const int lda, double* D, double* E, double* TAU, double* WORK, const int lwork, int* info) const;
    void GEHRD(const int n, const int ilo, const int ihi, double* A, const int lda, double* TAU, double* WORK, const int lwork, int* info) const;
    void TRTRS(const char UPLO, const char TRANS, const char DIAG, const int n, const int nrhs, const double* A, const int lda, double* B, const int ldb, int* info) const;
    void TRTRI(const char UPLO, const char DIAG, const int n, const double* A, const int lda, int* info) const;

    // Symmetric eigenproblem routines.
    void SPEV(const char JOBZ, const char UPLO, const int n, double* AP, double* W, double* Z, const int ldz, double* WORK, int* info) const;
    void SYEV(const char JOBZ, const char UPLO, const int n, double* A, const int lda, double* W, double* WORK, const int lwork, int* info) const;
    void SYGV(const int itype, const char JOBZ, const char UPLO, const int n, double* A, const int lda, double* B, const int ldb, double* W, double* WORK, const int lwork, int* info) const;
    void HEEV(const char JOBZ, const char UPLO, const int n, double* A, const int lda, double* W, double* WORK, const int lwork, double* RWORK, int* info) const;
    void HEGV(const int itype, const char JOBZ, const char UPLO, const int n, double* A, const int lda, double* B, const int ldb, double* W, double* WORK, const int lwork, double *RWORK, int* info) const;
    void STEQR(const char COMPZ, const int n, double* D, double* E, double* Z, const int ldz, double* WORK, int* info) const;

    // Non-Hermitian eigenproblem routines.
    void HSEQR(const char JOB, const char COMPZ, const int n, const int ilo, const int ihi, double* H, const int ldh, double* WR, double* WI, double* Z, const int ldz, double* WORK, const int lwork, int* info) const;
    void GEES(const char JOBVS, const char SORT, int (*ptr2func)(double*, double*), const int n, double* A, const int lda, int* sdim, double* WR, double* WI, double* VS, const int ldvs, double* WORK, const int lwork, int* BWORK, int* info) const;
    void GEES(const char JOBVS, const int n, double* A, const int lda, int* sdim, double* WR, double* WI, double* VS, const int ldvs, double* WORK, const int lwork, double* RWORK, int* BWORK, int* info) const;
    void GEEV(const char JOBVL, const char JOBVR, const int n, double* A, const int lda, double* WR, double* WI, double* VL, const int ldvl, double* VR, const int ldvr, double* WORK, const int lwork, int* info) const;
    void GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const int n, double* A, const int lda, double* WR, double* WI, double* VL, const int ldvl, double* VR, const int ldvr, int* ilo, int* ihi, double* SCALE, double* abnrm, double* RCONDE, double* RCONDV, double* WORK, const int lwork, int* IWORK, int* info) const;
    void GGEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const int n, double* A, const int lda, double* B, const int ldb, double* ALPHAR, double* ALPHAI, double* BETA, double* VL, const int ldvl, double* VR, const int ldvr, int* ilo, int* ihi, double* LSCALE, double* RSCALE, double* abnrm, double* bbnrm, double* RCONDE, double* RCONDV, double* WORK, const int lwork, int* IWORK, int* BWORK, int* info) const;
    void GGEV(const char JOBVL, const char JOBVR, const int n, double *A, const int lda, double *B, const int ldb, double *ALPHAR, double *ALPHAI, double *BETA, double *VL, const int ldvl, double *VR, const int ldvr, double *WORK, const int lwork, int *info) const;
    void TRSEN(const char JOB, const char COMPQ, const int *SELECT, const int n, double *T, const int ldt, double *Q, const int ldq, double *WR, double *WI, int *M, double *S, double *SEP, double *WORK, const int lwork, int *IWORK, const int liwork, int *info ) const;
    void TGSEN(const int ijob, const int wantq, const int wantz, const int *SELECT, const int n, double *A, const int lda, double *B, const int ldb, double *ALPHAR, double *ALPHAI, double *BETA, double *Q, const int ldq, double *Z, const int ldz, int *M, double *PL, double *PR, double *DIF, double *WORK, const int lwork, int *IWORK, const int liwork, int *info ) const;
    void GGES(const char JOBVL, const char JOBVR, const char SORT, int (*ptr2func)(double*, double*, double*), const int n, double* A, const int lda, double* B, const int ldb, int* sdim, double* ALPHAR, double* ALPHAI, double* BETA, double* VL, const int ldvl, double* VR, const int ldvr, double* WORK, const int lwork, int *bwork, int* info ) const;


    // SVD routine
    void GESVD(const char JOBU, const char JOBVT, const int m, const int n, double* A, const int lda, double* S, double* U, const int ldu, double* V, const int ldv, double* WORK, const int lwork, double* RWORK, int* info) const;

    // Orthogonal matrix routines.
    void ORMQR(const char SIDE, const char TRANS, const int m, const int n, const int k, double* A, const int lda, const double* TAU, double* C, const int ldc, double* WORK, const int lwork, int* info) const;
    void UNMQR(const char SIDE, const char TRANS, const int m, const int n, const int k, double* A, const int lda, const double* TAU, double* C, const int ldc, double* WORK, const int lwork, int* info) const;
    void ORGQR(const int m, const int n, const int k, double* A, const int lda, const double* TAU, double* WORK, const int lwork, int* info) const;
    void UNGQR(const int m, const int n, const int k, double* A, const int lda, const double* TAU, double* WORK, const int lwork, int* info) const;
    void ORGHR(const int n, const int ilo, const int ihi, double* A, const int lda, const double* TAU, double* WORK, const int lwork, int* info) const;
    void ORMHR(const char SIDE, const char TRANS, const int m, const int n, const int ilo, const int ihi, const double* A, const int lda, const double* TAU, double* C, const int ldc, double* WORK, const int lwork, int* info) const;

    // Triangular matrix routines.
    void TREVC(const char SIDE, const char HOWMNY, int* select, const int n, const double* T, const int ldt, double* VL, const int ldvl, double* VR, const int ldvr, const int mm, int* m, double* WORK, int* info) const;
    void TREVC(const char SIDE, const int n, const double* T, const int ldt, double* VL, const int ldvl, double* VR, const int ldvr, const int mm, int* m, double* WORK, double* RWORK, int* info) const;
    void TREXC(const char COMPQ, const int n, double* T, const int ldt, double* Q, const int ldq, int ifst, int ilst, double* WORK, int* info) const;
    void TGEVC(const char SIDE, const char HOWMNY, const int *SELECT, const int n, double *S, const int lds, double *P, const int ldp, double *VL, const int ldvl, double *VR, const int ldvr, const int mm, int *M, double *WORK, int *info) const;

    // Rotation/reflection generators
    void LARTG( const double f, const double g, double* c, double* s, double* r ) const;
    void LARFG( const int n, double* alpha, double* x, const int incx, double* tau ) const;

    // Matrix balancing routines.
    void GEBAL(const char JOBZ, const int n, double* A, const int lda, int ilo, int ihi, double* scale, int* info) const;
    void GEBAK(const char JOBZ, const char SIDE, const int n, const int ilo, const int ihi, const double* scale, const int m, double* V, const int ldv, int* info) const;

    // Random number generators
    double LARND( const int idist, int* seed ) const;
    void LARNV( const int idist, int* seed, const int n, double* v ) const;

    // Machine characteristic routines.
    double LAMCH(const char CMACH) const;
    int ILAENV( const int ispec, const std::string& NAME, const std::string& OPTS, const int N1 = -1, const int N2 = -1, const int N3 = -1, const int N4 = -1 ) const;

    // Miscellaneous routines.
    double LAPY2(const double x, const double y) const;

  };

  // END INT, DOUBLE SPECIALIZATION DECLARATION //

#ifdef HAVE_TEUCHOS_COMPLEX

  // BEGIN INT, COMPLEX<FLOAT> SPECIALIZATION DECLARATION //

  template<>
  class TEUCHOSNUMERICS_LIB_DLL_EXPORT LAPACK<int, std::complex<float> >
  {
  public:
    inline LAPACK(void) {}
    inline LAPACK(const LAPACK<int, std::complex<float> >& lapack) {}
    inline virtual ~LAPACK(void) {}

    // Symmetric positive definite linear system routines
    void PTTRF(const int n, std::complex<float>* d, std::complex<float>* e, int* info) const;
    void PTTRS(const int n, const int nrhs, const std::complex<float>* d, const std::complex<float>* e, std::complex<float>* B, const int ldb, int* info) const;
    void POTRF(const char UPLO, const int n, std::complex<float>* A, const int lda, int* info) const;
    void POTRS(const char UPLO, const int n, const int nrhs, const std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb, int* info) const;
    void POTRI(const char UPLO, const int n, std::complex<float>* A, const int lda, int* info) const;
    void POCON(const char UPLO, const int n, const std::complex<float>* A, const int lda, const float anorm, float* rcond, std::complex<float>* WORK, float* rwork, int* info) const;
    void POSV(const char UPLO, const int n, const int nrhs, std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb, int* info) const;
    void POEQU(const int n, const std::complex<float>* A, const int lda, float* S, float* scond, float* amax, int* info) const;
    void PORFS(const char UPLO, const int n, const int nrhs, std::complex<float>* A, const int lda, const std::complex<float>* AF, const int ldaf, const std::complex<float>* B, const int ldb, std::complex<float>* X, const int ldx, float* FERR, float* BERR, std::complex<float>* WORK, float* RWORK, int* info) const;
    void POSVX(const char FACT, const char UPLO, const int n, const int nrhs, std::complex<float>* A, const int lda, std::complex<float>* AF, const int ldaf, char EQUED, float* S, std::complex<float>* B, const int ldb, std::complex<float>* X, const int ldx, float* rcond, float* FERR, float* BERR, std::complex<float>* WORK, float* RWORK, int* info) const;

    // General Linear System Routines
    void GELS(const char TRANS, const int m, const int n, const int nrhs, std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb, std::complex<float>* WORK, const int lwork, int* info) const;
    void GELSS(const int m, const int n, const int nrhs, std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb, float* S, const float rcond, int* rank, std::complex<float>* WORK, const int lwork, float* RWORK, int* info) const;
    void GEQRF( const int m, const int n, std::complex<float>* A, const int lda, std::complex<float>* TAU, std::complex<float>* WORK, const int lwork, int* info) const;
    void UNGQR(const int m, const int n, const int k, std::complex<float>* A, const int lda, const std::complex<float>* TAU, std::complex<float>* WORK, const int lwork, int* info) const;
    void UNMQR(const char SIDE, const char TRANS, const int m, const int n, const int k, std::complex<float>* A, const int lda, const std::complex<float>* TAU, std::complex<float>* C, const int ldc, std::complex<float>* WORK, const int lwork, int* info) const;
    void GETRF(const int m, const int n, std::complex<float>* A, const int lda, int* IPIV, int* info) const;
    void GETRS(const char TRANS, const int n, const int nrhs, const std::complex<float>* A, const int lda, const int* IPIV, std::complex<float>* B, const int ldb, int* info) const;
    void LASCL(const char TYPE, const int kl, const int ku, const float cfrom, const float cto, const int m, const int n, std::complex<float>* A, const int lda, int* info) const;

    void LASWP (const int N,
		std::complex<float> A[],
		const int LDA,
		const int K1,
		const int K2,
		const int IPIV[],
		const int INCX) const;

    void GBTRF(const int m, const int n, const int kl, const int ku, std::complex<float>* A, const int lda, int* IPIV, int* info) const;
    void GBTRS(const char TRANS, const int n, const int kl, const int ku, const int nrhs, const std::complex<float>* A, const int lda, const int* IPIV, std::complex<float>* B, const int ldb, int* info) const;
    void GTTRF(const int n, std::complex<float>* dl, std::complex<float>* d, std::complex<float>* du, std::complex<float>* du2, int* IPIV, int* info) const;
    void GTTRS(const char TRANS, const int n, const int nrhs, const std::complex<float>* dl, const std::complex<float>* d, const std::complex<float>* du, const std::complex<float>* du2, const int* IPIV, std::complex<float>* B, const int ldb, int* info) const;
    void GETRI(const int n, std::complex<float>* A, const int lda, const int* IPIV, std::complex<float>* WORK, const int lwork, int* info) const;
    void LATRS (const char UPLO, const char TRANS, const char DIAG, const char NORMIN, const int N, std::complex<float>* A, const int LDA, std::complex<float>* X, float* SCALE, float* CNORM, int* INFO) const;
    void GECON(const char NORM, const int n, const std::complex<float>* A, const int lda, const float anorm, float* rcond, std::complex<float>* WORK, float* RWORK, int* info) const;
    void GBCON(const char NORM, const int n, const int kl, const int ku, const std::complex<float>* A, const int lda, int* IPIV, const float anorm, float* rcond, std::complex<float>* WORK, float* RWORK, int* info) const;
    float LANGB(const char NORM, const int n, const int kl, const int ku, const std::complex<float>* A, const int lda, float* WORK) const;
    void GESV(const int n, const int nrhs, std::complex<float>* A, const int lda, int* IPIV, std::complex<float>* B, const int ldb, int* info) const;
    void GEEQU(const int m, const int n, const std::complex<float>* A, const int lda, float* R, float* C, float* rowcond, float* colcond, float* amax, int* info) const;
    void GERFS(const char TRANS, const int n, const int nrhs, const std::complex<float>* A, const int lda, const std::complex<float>* AF, const int ldaf, const int* IPIV, const std::complex<float>* B, const int ldb, std::complex<float>* X, const int ldx, float* FERR, float* BERR, std::complex<float>* WORK, float* RWORK, int* info) const;
    void GBEQU(const int m, const int n, const int kl, const int ku, const std::complex<float>* A, const int lda, float* R, float* C, float* rowcond, float* colcond, float* amax, int* info) const;
    void GBRFS(const char TRANS, const int n, const int kl, const int ku, const int nrhs, const std::complex<float>* A, const int lda, const std::complex<float>* AF, const int ldaf, const int* IPIV, const std::complex<float>* B, const int ldb, std::complex<float>* X, const int ldx, float* FERR, float* BERR, std::complex<float>* WORK, float* RWORK, int* info) const;
    void GESVX(const char FACT, const char TRANS, const int n, const int nrhs, std::complex<float>* A, const int lda, std::complex<float>* AF, const int ldaf, int* IPIV, char EQUED, float* R, float* C, std::complex<float>* B, const int ldb, std::complex<float>* X, const int ldx, float* rcond, float* FERR, float* BERR, std::complex<float>* WORK, float* RWORK, int* info) const;
    void GEHRD(const int n, const int ilo, const int ihi, std::complex<float>* A, const int lda, std::complex<float>* TAU, std::complex<float>* WORK, const int lwork, int* info) const;
    void TRTRS(const char UPLO, const char TRANS, const char DIAG, const int n, const int nrhs, const std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb, int* info) const;
    void TRTRI(const char UPLO, const char DIAG, const int n, const std::complex<float>* A, const int lda, int* info) const;

    // Symmetric eigenvalue routines.
    void STEQR(const char COMPZ, const int n, float* D, float* E, std::complex<float>* Z, const int ldz, float* WORK, int* info) const;
    void HEEV(const char JOBZ, const char UPLO, const int n, std::complex<float>* A, const int lda, float* W, std::complex<float>* WORK, const int lwork, float* RWORK, int* info) const;
    void HEGV(const int itype, const char JOBZ, const char UPLO, const int n, std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb, float* W, std::complex<float>* WORK, const int lwork, float *RWORK, int* info) const;

    // Non-Hermitian eigenvalue routines.
    void HSEQR(const char JOB, const char COMPZ, const int n, const int ilo, const int ihi, std::complex<float>* H, const int ldh, std::complex<float>* W, std::complex<float>* Z, const int ldz, std::complex<float>* WORK, const int lwork, int* info) const;
    void GEES(const char JOBVS, const char SORT, int (*ptr2func)(std::complex<float>*), const int n, std::complex<float>* A, const int lda, int* sdim, std::complex<float>* W, std::complex<float>* VS, const int ldvs, std::complex<float>* WORK, const int lwork, float* RWORK, int* BWORK, int* info) const;
    void GEES(const char JOBVS, const int n, std::complex<float>* A, const int lda, int* sdim, float* WR, float* WI, std::complex<float>* VS, const int ldvs, std::complex<float>* WORK, const int lwork, float* RWORK, int* BWORK, int* info) const;
    void GEEV(const char JOBVL, const char JOBVR, const int n, std::complex<float>* A, const int lda, std::complex<float>* W, std::complex<float>* VL, const int ldvl, std::complex<float>* VR, const int ldvr, std::complex<float>* WORK, const int lwork, float* RWORK, int* info) const;
    void GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const int n, std::complex<float>* A, const int lda, std::complex<float>* W, std::complex<float>* VL, const int ldvl, std::complex<float>* VR, const int ldvr, int* ilo, int* ihi, float* SCALE, float* abnrm, float* RCONDE, float* RCONDV, std::complex<float>* WORK, const int lwork, float* RWORK, int* info) const;
    void GGEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const int n, std::complex<float>* A, const int lda, std::complex<float>* B, const int ldb, std::complex<float>* ALPHA, std::complex<float>* BETA, std::complex<float>* VL, const int ldvl, std::complex<float>* VR, const int ldvr, int* ilo, int* ihi, float* LSCALE, float* RSCALE, float* abnrm, float* bbnrm, float* RCONDE, float* RCONDV, std::complex<float>* WORK, const int lwork, float * RWORK, int* IWORK, int* BWORK, int* info) const;
    void GGEV(const char JOBVL, const char JOBVR, const int n, std::complex<float> *A, const int lda, std::complex<float> *B, const int ldb, std::complex<float>* ALPHA, std::complex<float>* BETA, std::complex<float>* VL, const int ldvl, std::complex<float>* VR, const int ldvr, std::complex<float> *WORK, const int lwork, float* RWORK, int* info) const;

    // SVD routine
    void GESVD(const char JOBU, const char JOBVT, const int m, const int n, std::complex<float>* A, const int lda, float* S, std::complex<float>* U, const int ldu, std::complex<float>* V, const int ldv, std::complex<float>* WORK, const int lwork, float* RWORK, int* info) const;

    // Triangular matrix routines.
    void TREVC(const char SIDE, const char HOWMNY, int* select, const int n, const std::complex<float>* T, const int ldt, std::complex<float>* VL, const int ldvl, std::complex<float>* VR, const int ldvr, const int mm, int* m, std::complex<float>* WORK, float* RWORK, int* info) const;
    void TREVC(const char SIDE, const int n, const std::complex<float>* T, const int ldt, std::complex<float>* VL, const int ldvl, std::complex<float>* VR, const int ldvr, const int mm, int* m, std::complex<float>* WORK, float* RWORK, int* info) const;
    void TREXC(const char COMPQ, const int n, std::complex<float>* T, const int ldt, std::complex<float>* Q, const int ldq, int ifst, int ilst, std::complex<float>* WORK, int* info) const;

    // Rotation/reflection generators
    void LARTG( const std::complex<float> f, const std::complex<float> g, float* c, std::complex<float>* s, std::complex<float>* r ) const;
    void LARFG( const int n, std::complex<float>* alpha, std::complex<float>* x, const int incx, std::complex<float>* tau ) const;

    // Matrix balancing routines.
    void GEBAL(const char JOBZ, const int n, std::complex<float>* A, const int lda, int ilo, int ihi, float* scale, int* info) const;
    void GEBAK(const char JOBZ, const char SIDE, const int n, const int ilo, const int ihi, const float* scale, const int m, std::complex<float>* V, const int ldv, int* info) const;

    // Random number generators
    std::complex<float> LARND( const int idist, int* seed ) const;
    void LARNV( const int idist, int* seed, const int n, std::complex<float>* v ) const;

    // Machine characteristics
    int ILAENV( const int ispec, const std::string& NAME, const std::string& OPTS, const int N1 = -1, const int N2 = -1, const int N3 = -1, const int N4 = -1 ) const;

  };

  // END INT, COMPLEX<FLOAT> SPECIALIZATION DECLARATION //

  // BEGIN INT, COMPLEX<DOUBLE> SPECIALIZATION DECLARATION //

  template<>
  class TEUCHOSNUMERICS_LIB_DLL_EXPORT LAPACK<int, std::complex<double> >
  {
  public:
    inline LAPACK(void) {}
    inline LAPACK(const LAPACK<int, std::complex<double> >& lapack) {}
    inline virtual ~LAPACK(void) {}

    // Symmetric positive definite linear system routines
    void PTTRF(const int n, std::complex<double>* d, std::complex<double>* e, int* info) const;
    void PTTRS(const int n, const int nrhs, const std::complex<double>* d, const std::complex<double>* e, std::complex<double>* B, const int ldb, int* info) const;
    void POTRF(const char UPLO, const int n, std::complex<double>* A, const int lda, int* info) const;
    void POTRS(const char UPLO, const int n, const int nrhs, const std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb, int* info) const;
    void POTRI(const char UPLO, const int n, std::complex<double>* A, const int lda, int* info) const;
    void POCON(const char UPLO, const int n, const std::complex<double>* A, const int lda, const double anorm, double* rcond, std::complex<double>* WORK, double* RWORK, int* info) const;
    void POSV(const char UPLO, const int n, const int nrhs, std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb, int* info) const;
    void POEQU(const int n, const std::complex<double>* A, const int lda, double* S, double* scond, double* amax, int* info) const;
    void PORFS(const char UPLO, const int n, const int nrhs, std::complex<double>* A, const int lda, const std::complex<double>* AF, const int ldaf, const std::complex<double>* B, const int ldb, std::complex<double>* X, const int ldx, double* FERR, double* BERR, std::complex<double>* WORK, double* RWORK, int* info) const;
    void POSVX(const char FACT, const char UPLO, const int n, const int nrhs, std::complex<double>* A, const int lda, std::complex<double>* AF, const int ldaf, char EQUED, double* S, std::complex<double>* B, const int ldb, std::complex<double>* X, const int ldx, double* rcond, double* FERR, double* BERR, std::complex<double>* WORK, double* RWORK, int* info) const;

    // General Linear System Routines
    void GELS(const char TRANS, const int m, const int n, const int nrhs, std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb, std::complex<double>* WORK, const int lwork, int* info) const;
    void GELSS(const int m, const int n, const int nrhs, std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb, double* S, const double rcond, int* rank, std::complex<double>* WORK, const int lwork, double* RWORK, int* info) const;
    void GEQRF( const int m, const int n, std::complex<double>* A, const int lda, std::complex<double>* TAU, std::complex<double>* WORK, const int lwork, int* info) const;
    void UNGQR(const int m, const int n, const int k, std::complex<double>* A, const int lda, const std::complex<double>* TAU, std::complex<double>* WORK, const int lwork, int* info) const;
    void UNMQR(const char SIDE, const char TRANS, const int m, const int n, const int k, std::complex<double>* A, const int lda, const std::complex<double>* TAU, std::complex<double>* C, const int ldc, std::complex<double>* WORK, const int lwork, int* info) const;
    void GETRF(const int m, const int n, std::complex<double>* A, const int lda, int* IPIV, int* info) const;
    void GETRS(const char TRANS, const int n, const int nrhs, const std::complex<double>* A, const int lda, const int* IPIV, std::complex<double>* B, const int ldb, int* info) const;
    void LASCL(const char TYPE, const int kl, const int ku, const double cfrom, const double cto, const int m, const int n, std::complex<double>* A, const int lda, int* info) const;

    void LASWP (const int N,
		std::complex<double> A[],
		const int LDA,
		const int K1,
		const int K2,
		const int IPIV[],
		const int INCX) const;

    void GBTRF(const int m, const int n, const int kl, const int ku, std::complex<double>* A, const int lda, int* IPIV, int* info) const;
    void GBTRS(const char TRANS, const int n, const int kl, const int ku, const int nrhs, const std::complex<double>* A, const int lda, const int* IPIV, std::complex<double>* B, const int ldb, int* info) const;
    void GTTRF(const int n, std::complex<double>* dl, std::complex<double>* d, std::complex<double>* du, std::complex<double>* du2, int* IPIV, int* info) const;
    void GTTRS(const char TRANS, const int n, const int nrhs, const std::complex<double>* dl, const std::complex<double>* d, const std::complex<double>* du, const std::complex<double>* du2, const int* IPIV, std::complex<double>* B, const int ldb, int* info) const;
    void GETRI(const int n, std::complex<double>* A, const int lda, const int* IPIV, std::complex<double>* WORK, const int lwork, int* info) const;
    void LATRS (const char UPLO, const char TRANS, const char DIAG, const char NORMIN, const int N, std::complex<double>* A, const int LDA, std::complex<double>* X, double* SCALE, double* CNORM, int* INFO) const;
    void GECON(const char NORM, const int n, const std::complex<double>* A, const int lda, const double anorm, double* rcond, std::complex<double>* WORK, double* RWORK, int* info) const;
    void GBCON(const char NORM, const int n, const int kl, const int ku, const std::complex<double>* A, const int lda, int* IPIV, const double anorm, double* rcond, std::complex<double>* WORK, double* RWORK, int* info) const;
    double LANGB(const char NORM, const int n, const int kl, const int ku, const std::complex<double>* A, const int lda, double* WORK) const;
    void GESV(const int n, const int nrhs, std::complex<double>* A, const int lda, int* IPIV, std::complex<double>* B, const int ldb, int* info) const;
    void GEEQU(const int m, const int n, const std::complex<double>* A, const int lda, double* R, double* C, double* rowcond, double* colcond, double* amax, int* info) const;
    void GERFS(const char TRANS, const int n, const int nrhs, const std::complex<double>* A, const int lda, const std::complex<double>* AF, const int ldaf, const int* IPIV, const std::complex<double>* B, const int ldb, std::complex<double>* X, const int ldx, double* FERR, double* BERR, std::complex<double>* WORK, double* RWORK, int* info) const;
    void GBEQU(const int m, const int n, const int kl, const int ku, const std::complex<double>* A, const int lda, double* R, double* C, double* rowcond, double* colcond, double* amax, int* info) const;
    void GBRFS(const char TRANS, const int n, const int kl, const int ku, const int nrhs, const std::complex<double>* A, const int lda, const std::complex<double>* AF, const int ldaf, const int* IPIV, const std::complex<double>* B, const int ldb, std::complex<double>* X, const int ldx, double* FERR, double* BERR, std::complex<double>* WORK, double* RWORK, int* info) const;
    void GESVX(const char FACT, const char TRANS, const int n, const int nrhs, std::complex<double>* A, const int lda, std::complex<double>* AF, const int ldaf, int* IPIV, char EQUED, double* R, double* C, std::complex<double>* B, const int ldb, std::complex<double>* X, const int ldx, double* rcond, double* FERR, double* BERR, std::complex<double>* WORK, double* RWORK, int* info) const;
    void GEHRD(const int n, const int ilo, const int ihi, std::complex<double>* A, const int lda, std::complex<double>* TAU, std::complex<double>* WORK, const int lwork, int* info) const;
    void TRTRS(const char UPLO, const char TRANS, const char DIAG, const int n, const int nrhs, const std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb, int* info) const;
    void TRTRI(const char UPLO, const char DIAG, const int n, const std::complex<double>* A, const int lda, int* info) const;

    // Symmetric eigenvalue routines.
    void STEQR(const char COMPZ, const int n, double* D, double* E, std::complex<double>* Z, const int ldz, double* WORK, int* info) const;
    void HEEV(const char JOBZ, const char UPLO, const int n, std::complex<double>* A, const int lda, double* W, std::complex<double>* WORK, const int lwork, double* RWORK, int* info) const;
    void HEGV(const int itype, const char JOBZ, const char UPLO, const int n, std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb, double* W, std::complex<double>* WORK, const int lwork, double *RWORK, int* info) const;

    // Non-hermitian eigenvalue routines.
    void HSEQR(const char JOB, const char COMPZ, const int n, const int ilo, const int ihi, std::complex<double>* H, const int ldh, std::complex<double>* W, std::complex<double>* Z, const int ldz, std::complex<double>* WORK, const int lwork, int* info) const;
    void GEES(const char JOBVS, const char SORT, int (*ptr2func)(std::complex<double>*), const int n, std::complex<double>* A, const int lda, int* sdim, std::complex<double>* W, std::complex<double>* VS, const int ldvs, std::complex<double>* WORK, const int lwork, double* RWORK, int* BWORK, int* info) const;
    void GEES(const char JOBVS, const int n, std::complex<double>* A, const int lda, int* sdim, double* WR, double* WI, std::complex<double>* VS, const int ldvs, std::complex<double>* WORK, const int lwork, double* RWORK, int* BWORK, int* info) const;
    void GEEV(const char JOBVL, const char JOBVR, const int n, std::complex<double>* A, const int lda, std::complex<double>* W, std::complex<double>* VL, const int ldvl, std::complex<double>* VR, const int ldvr, std::complex<double>* WORK, const int lwork, double* RWORK, int* info) const;
    void GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const int n, std::complex<double>* A, const int lda, std::complex<double>* W, std::complex<double>* VL, const int ldvl, std::complex<double>* VR, const int ldvr, int* ilo, int* ihi, double* SCALE, double* abnrm, double* RCONDE, double* RCONDV, std::complex<double>* WORK, const int lwork, double* RWORK, int* info) const;
    void GGEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const int n, std::complex<double>* A, const int lda, std::complex<double>* B, const int ldb, std::complex<double>* ALPHA, std::complex<double>* BETA, std::complex<double>* VL, const int ldvl, std::complex<double>* VR, const int ldvr, int* ilo, int* ihi, double* LSCALE, double* RSCALE, double* abnrm, double* bbnrm, double* RCONDE, double* RCONDV, std::complex<double>* work, const int lwork, double* RWORK, int* IWORK, int* BWORK, int* info) const;
    void GGEV(const char JOBVL, const char JOBVR, const int n, std::complex<double> *A, const int lda, std::complex<double> *B, const int ldb, std::complex<double>* ALPHA, std::complex<double>* BETA, std::complex<double>* VL, const int ldvl, std::complex<double>*VR, const int ldvr, std::complex<double> *WORK, const int lwork, double* RWORK, int* info) const;

    // SVD routine
    void GESVD(const char JOBU, const char JOBVT, const int m, const int n, std::complex<double>* A, const int lda, double* S, std::complex<double>* U, const int ldu, std::complex<double>* V, const int ldv, std::complex<double>* WORK, const int lwork, double* RWORK, int* info) const;

    // Triangular matrix routines.
    void TREVC(const char SIDE, const char HOWMNY, int* select, const int n, const std::complex<double>* T, const int ldt, std::complex<double>* VL, const int ldvl, std::complex<double>* VR, const int ldvr, const int mm, int* m, std::complex<double>* WORK, double* RWORK, int* info) const;
    void TREVC(const char SIDE, const int n, const std::complex<double>* T, const int ldt, std::complex<double>* VL, const int ldvl, std::complex<double>* VR, const int ldvr, const int mm, int* m, std::complex<double>* WORK, double* RWORK, int* info) const;
    void TREXC(const char COMPQ, const int n, std::complex<double>* T, const int ldt, std::complex<double>* Q, const int ldq, int ifst, int ilst, std::complex<double>* WORK, int* info) const;

    // Rotation/reflection generators
    void LARTG( const std::complex<double> f, const std::complex<double> g, double* c, std::complex<double>* s, std::complex<double>* r ) const;
    void LARFG( const int n, std::complex<double>* alpha, std::complex<double>* x, const int incx, std::complex<double>* tau ) const;

    // Matrix balancing routines.
    void GEBAL(const char JOBZ, const int n, std::complex<double>* A, const int lda, int ilo, int ihi, double* scale, int* info) const;
    void GEBAK(const char JOBZ, const char SIDE, const int n, const int ilo, const int ihi, const double* scale, const int m, std::complex<double>* V, const int ldv, int* info) const;

    // Random number generators
    std::complex<double> LARND( const int idist, int* seed ) const;
    void LARNV( const int idist, int* seed, const int n, std::complex<double>* v ) const;

    // Machine characteristics
    int ILAENV( const int ispec, const std::string& NAME, const std::string& OPTS, const int N1 = -1, const int N2 = -1, const int N3 = -1, const int N4 = -1 ) const;

  };

  // END INT, COMPLEX<DOUBLE> SPECIALIZATION DECLARATION //

#endif // HAVE_TEUCHOS_COMPLEX

#endif // DOXYGEN_SHOULD_SKIP_THIS

} // namespace Teuchos

#endif // _TEUCHOS_LAPACK_HPP_
