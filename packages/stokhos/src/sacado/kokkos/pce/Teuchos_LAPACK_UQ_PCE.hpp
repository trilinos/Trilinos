// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _TEUCHOS_LAPACK_UQ_PCE_HPP_
#define _TEUCHOS_LAPACK_UQ_PCE_HPP_

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_TestForException.hpp"
#include "Sacado_UQ_PCE.hpp"

/*!
 * \class Teuchos::LAPACK
 * \brief Specialization for Sacado::UQ::PCE< Storage<...> >
 *
 * Currently a shell implementation just so things that use LAPACK
 * will compile.  However each method throws a run-time error.
*/
namespace Teuchos {

  template<typename OrdinalType, typename Storage>
  class LAPACK<OrdinalType, Sacado::UQ::PCE<Storage> >
  {
  public:

    typedef Sacado::UQ::PCE<Storage> ScalarType;
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
    void PTTRF(const OrdinalType n, ScalarType* d, ScalarType* e, OrdinalType* info) const
      { throw_error("PTTRF"); }

    //! Solves a tridiagonal system \c A*X=B using the \L*D*L' factorization of \c A computed by PTTRF.
    void PTTRS(const OrdinalType n, const OrdinalType nrhs, const ScalarType* d, const ScalarType* e, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
      { throw_error("PTTRS"); }

    //! Computes Cholesky factorization of a real symmetric positive definite matrix \c A.
    void POTRF(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const
      { throw_error("POTRF"); }

    //! Solves a system of linear equations \c A*X=B, where \c A is a symmetric positive definite matrix factored by POTRF and the \c nrhs solutions are returned in \c B.
    void POTRS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
      { throw_error("POTRS"); }

    //! Computes the inverse of a real symmetric positive definite matrix \c A using the Cholesky factorization \c A from POTRF.
    void POTRI(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const
      { throw_error("POTRI"); }

    //! Estimates the reciprocal of the condition number (1-norm) of a real symmetric positive definite matrix \c A using the Cholesky factorization from POTRF.

    void POCON(const char UPLO, const OrdinalType n, const ScalarType* A, const OrdinalType lda, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
      { throw_error("POCON"); }

    //! Computes the solution to a real system of linear equations \c A*X=B, where \c A is a symmetric positive definite matrix and the \c nrhs solutions are returned in \c B.
    void POSV(const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
      { throw_error("POSV"); }

    //! Computes row and column scalings intended to equilibrate a symmetric positive definite matrix \c A and reduce its condition number (w.r.t. 2-norm).
    void POEQU(const OrdinalType n, const ScalarType* A, const OrdinalType lda, MagnitudeType* S, MagnitudeType* scond, MagnitudeType* amax, OrdinalType* info) const
      { throw_error("POEQU"); }

    //! Improves the computed solution to a system of linear equations when the coefficient matrix is symmetric positive definite, and provides error bounds and backward error estimates for the solution.
    void PORFS(const char UPLO, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
      { throw_error("PORFS"); }

    //! Uses the Cholesky factorization to compute the solution to a real system of linear equations \c A*X=B, where \c A is symmetric positive definite.  System can be equilibrated by POEQU and iteratively refined by PORFS, if requested.
    void POSVX(const char FACT, const char UPLO, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* AF, const OrdinalType ldaf, char EQUED, ScalarType* S, ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* rcond, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
      { throw_error("POSVX"); }
    //@}

    //! @name General Linear System Routines.
    //@{

    //! Solves an over/underdetermined real \c m by \c n linear system \c A using QR or LQ factorization of A.
    void GELS(const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("GELS"); }

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
    void GELSS(const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, MagnitudeType* S, const MagnitudeType rcond, OrdinalType* rank, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const
      { throw_error("GELSS"); }

    //! Legacy GELSS interface for real-valued ScalarType.
    void GELSS(const OrdinalType m, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* S, const ScalarType rcond, OrdinalType* rank, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("GELSS"); }

    //! Solves the linear equality-constrained least squares (LSE) problem where \c A is an \c m by \c n matrix,\c B is a \c p by \c n matrix \c C is a given \c m-vector, and D is a given \c p-vector.
    void GGLSE(const OrdinalType m, const OrdinalType n, const OrdinalType p, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* C, ScalarType* D, ScalarType* X, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("GGLSE"); }

    //! Computes a QR factorization of a general \c m by \c n matrix \c A.
    void GEQRF( const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("GEQRF"); }

    //! Computes an LU factorization of a general \c m by \c n matrix \c A using partial pivoting with row interchanges.
    void GETRF(const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const
      { throw_error("GETRF"); }

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B with a general \c n by \c n matrix \c A using the LU factorization computed by GETRF.
    void GETRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
      { throw_error("GETRS"); }

    //! Multiplies the \c m by \c n matrix \c A by the real scalar \c cto/cfrom.
    void LASCL(const char TYPE, const OrdinalType kl, const OrdinalType ku, const MagnitudeType cfrom, const MagnitudeType cto, const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* info) const
      { throw_error("LASCL"); }

    //! Computes a QR factorization with column pivoting of a matrix \c A: A*P = Q*R using Level 3 BLAS
    void
    GEQP3(const OrdinalType m,
          const OrdinalType n, ScalarType* A,
          const OrdinalType lda,
          OrdinalType *jpvt,
          ScalarType* TAU,
          ScalarType* WORK,
          const OrdinalType lwork,
          MagnitudeType* RWORK,
          OrdinalType* info ) const
      { throw_error("GEQP3"); }

    //! Apply a series of row interchanges to the matrix A.
    void
    LASWP (const OrdinalType N,
           ScalarType A[],
           const OrdinalType LDA,
           const OrdinalType K1,
           const OrdinalType K2,
           const OrdinalType IPIV[],
           const OrdinalType INCX) const
      { throw_error("LASWP"); }

    //! Computes an LU factorization of a general banded \c m by \c n matrix \c A using partial pivoting with row interchanges.
    void GBTRF(const OrdinalType m, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, OrdinalType* info) const
      { throw_error("GBTRF"); }

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B with a general banded \c n by \c n matrix \c A using the LU factorization computed by GBTRF.
    void GBTRS(const char TRANS, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
      { throw_error("GBTRS"); }

    //! Computes an LU factorization of a \c n by \c n tridiagonal matrix \c A using partial pivoting with row interchanges.
    void GTTRF(const OrdinalType n, ScalarType* dl, ScalarType* d, ScalarType* du, ScalarType* du2, OrdinalType* IPIV, OrdinalType* info) const
      { throw_error("GTTRF"); }

    //! Solves a system of linear equations \c A*X=B or \c A'*X=B or \c A^H*X=B with a tridiagonal matrix \c A using the LU factorization computed by GTTRF.
    void GTTRS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* dl, const ScalarType* d, const ScalarType* du, const ScalarType* du2, const OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
      { throw_error("GTTRS"); }

    //! Computes the inverse of a matrix \c A using the LU factorization computed by GETRF.
    void GETRI(const OrdinalType n, ScalarType* A, const OrdinalType lda, const OrdinalType* IPIV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("GETRI"); }

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
           OrdinalType* INFO) const
      { throw_error("LATRS"); }

    //! Estimates the reciprocal of the condition number of a general real matrix \c A, in either the 1-norm or the infinity-norm, using the LU factorization computed by GETRF.
    void GECON(const char NORM, const OrdinalType n, const ScalarType* A, const OrdinalType lda, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
      { throw_error("GECON"); }

    //! Estimates the reciprocal of the condition number of a general banded real matrix \c A, in either the 1-norm or the infinity-norm, using the LU factorization computed by GETRF.
    void GBCON(const char NORM, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, const ScalarType anorm, ScalarType* rcond, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
      { throw_error("GBCON"); }

    //! Returns the value of the one norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of an \c n by \c n band matrix \c A, with \c kl sub-diagonals and \c ku super-diagonals.
    typename ScalarTraits<ScalarType>::magnitudeType LANGB(const char NORM, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const ScalarType* A, const OrdinalType lda, MagnitudeType* WORK) const
      { throw_error("LANGB"); return MagnitudeType(0); }

    //! Computes the solution to a real system of linear equations \c A*X=B, where \c A is factored through GETRF and the \c nrhs solutions are computed through GETRS.
    void GESV(const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, OrdinalType* IPIV, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
      { throw_error("GESV"); }

    //! Computes row and column scalings intended to equilibrate an \c m by \c n matrix \c A and reduce its condition number.
    void GEEQU(const OrdinalType m, const OrdinalType n, const ScalarType* A, const OrdinalType lda, MagnitudeType* R, MagnitudeType* C, MagnitudeType* rowcond, MagnitudeType* colcond, MagnitudeType* amax, OrdinalType* info) const
      { throw_error("GEEQU"); }

    //! Improves the computed solution to a system of linear equations and provides error bounds and backward error estimates for the solution.  Use after GETRF/GETRS.
    void GERFS(const char TRANS, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, MagnitudeType* FERR, MagnitudeType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
      { throw_error("GERFS"); }

    //! Computes row and column scalings intended to equilibrate an \c m by \c n banded matrix \c A and reduce its condition number.
    void GBEQU(const OrdinalType m, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const ScalarType* A, const OrdinalType lda, MagnitudeType* R, MagnitudeType* C, MagnitudeType* rowcond, MagnitudeType* colcond, MagnitudeType* amax, OrdinalType* info) const
      { throw_error("GBEQU"); }

    //! Improves the computed solution to a banded system of linear equations and provides error bounds and backward error estimates for the solution.  Use after GBTRF/GBTRS.
    void GBRFS(const char TRANS, const OrdinalType n, const OrdinalType kl, const OrdinalType ku, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, const ScalarType* AF, const OrdinalType ldaf, const OrdinalType* IPIV, const ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
      { throw_error("GBRFS"); }

    //! Uses the LU factorization to compute the solution to a real system of linear equations \c A*X=B, returning error bounds on the solution and a condition estimate.
    void GESVX(const char FACT, const char TRANS, const OrdinalType n, const OrdinalType nrhs, ScalarType* A, const OrdinalType lda, ScalarType* AF, const OrdinalType ldaf, OrdinalType* IPIV, char EQUED, ScalarType* R, ScalarType* C, ScalarType* B, const OrdinalType ldb, ScalarType* X, const OrdinalType ldx, ScalarType* rcond, ScalarType* FERR, ScalarType* BERR, ScalarType* WORK, OrdinalType* IWORK, OrdinalType* info) const
      { throw_error("GESVX"); }

    /*! \brief Reduces a real symmetric matrix \c A to tridiagonal form by orthogonal similarity transformations.
        \note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    void SYTRD(const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* D, ScalarType* E, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("SYTRD"); }

    //! Reduces a real general matrix \c A to upper Hessenberg form by orthogonal similarity transformations.
    void GEHRD(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* A, const OrdinalType lda, ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("GEHRD"); }

    //! Solves a triangular linear system of the form \c A*X=B or \c A**T*X=B, where \c A is a triangular matrix.
    void TRTRS(const char UPLO, const char TRANS, const char DIAG, const OrdinalType n, const OrdinalType nrhs, const ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, OrdinalType* info) const
      { throw_error("TRTRS"); }

    //! Computes the inverse of an upper or lower triangular matrix \c A.
    void TRTRI(const char UPLO, const char DIAG, const OrdinalType n, const ScalarType* A, const OrdinalType lda, OrdinalType* info) const
      { throw_error("TRTRI"); }
    //@}

    //! @name Symmetric Eigenproblem Routines
    //@{
    /*! \brief Computes the eigenvalues and, optionally, eigenvectors of a symmetric \c n by \c n matrix \c A in packed storage.
        \note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    void SPEV(const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* AP, ScalarType* W, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, OrdinalType* info) const
      { throw_error("SPEV"); }

    /*! \brief Computes all the eigenvalues and, optionally, eigenvectors of a symmetric \c n by \c n matrix A.
        \note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    void SYEV(const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* W, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("SYEV"); }

    /*! \brief Computes all the eigenvalues and, optionally, eigenvectors of a symmetric \c n by \c n matrix pencil \c {A,B}, where \c A is symmetric and \c B is symmetric positive-definite.
        \note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    void SYGV(const OrdinalType itype, const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, ScalarType* W, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("SYGV"); }

    /*! \brief Computes all the eigenvalues and, optionally, eigenvectors of a Hermitian \c n by \c n matrix A.
        \note This method will call SYEV when ScalarType is \c float or \c double.
    */
    void HEEV(const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, MagnitudeType* W, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const
      { throw_error("HEEV"); }

    /*! \brief Computes all the eigenvalues and, optionally, eigenvectors of a generalized Hermitian-definite \c n by \c n matrix pencil \c {A,B}, where \c A is Hermitian and \c B is Hermitian positive-definite.
        \note This method will call SYGV when ScalarType is \c float or \c double.
    */
    void HEGV(const OrdinalType itype, const char JOBZ, const char UPLO, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, MagnitudeType* W, ScalarType* WORK, const OrdinalType lwork, MagnitudeType *RWORK, OrdinalType* info) const
      { throw_error("HEGV"); }

    //! Computes the eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal \c n by \c n matrix \c A using implicit QL/QR.  The eigenvectors can only be computed if \c A was reduced to tridiagonal form by SYTRD.
    void STEQR(const char COMPZ, const OrdinalType n, ScalarType* D, ScalarType* E, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, OrdinalType* info) const
      { throw_error("STEQR"); }
    //@}

    //! @name Non-Hermitian Eigenproblem Routines
    //@{
    //! Computes the eigenvalues of a real upper Hessenberg matrix \c H and, optionally, the matrices \c T and \c Z from the Schur decomposition, where T is an upper quasi-triangular matrix and Z contains the Schur vectors.
    void HSEQR(const char JOB, const char COMPZ, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* H, const OrdinalType ldh, ScalarType* WR, ScalarType* WI, ScalarType* Z, const OrdinalType ldz, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("HSEQR"); }

    /*! Computes for an \c n by \c n nonsymmetric matrix \c A, the eigenvalues, the Schur form \c T, and, optionally, the matrix of Schur vectors \c Z. When \c ScalarType is \c float or \c double, the real Schur form is computed.
       \note (This is the version used for \c float and \c double, where \c select requires two arguments to represent a complex eigenvalue.)
    */
    void GEES(const char JOBVS, const char SORT, OrdinalType (*ptr2func)(ScalarType*, ScalarType*), const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, ScalarType* WR, ScalarType* WI, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, OrdinalType* BWORK, OrdinalType* info) const
      { throw_error("GEES"); }

    /*! Computes for an \c n by \c n nonsymmetric matrix \c A, the eigenvalues, the Schur form \c T, and, optionally, the matrix of Schur vectors \c Z. When \c ScalarType is \c float or \c double, the real Schur form is computed.
       \note (This is the version used for \c std::complex<float> and \c std::complex<double>, where \c select requires one arguments to represent a complex eigenvalue.)
    */
    void GEES(const char JOBVS, const char SORT, OrdinalType (*ptr2func)(ScalarType*), const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, ScalarType* W, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* BWORK, OrdinalType* info) const
      { throw_error("GEES"); }

    /*! Computes for an \c n by \c n nonsymmetric matrix \c A, the eigenvalues, the Schur form \c T, and, optionally, the matrix of Schur vectors \c Z. When \c ScalarType is \c float or \c double, the real Schur form is computed.
       \note (This is the version used for any \c ScalarType, when the user doesn't want to enable the sorting functionality.)
    */
    void GEES(const char JOBVS, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType* sdim, MagnitudeType* WR, MagnitudeType* WI, ScalarType* VS, const OrdinalType ldvs, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* BWORK, OrdinalType* info) const
      { throw_error("GEES"); }

    /// \brief Computes for an \c n by \c n real nonsymmetric matrix \c A, the eigenvalues and, optionally, the left and/or right eigenvectors.
    ///
    /// Real and imaginary parts of the eigenvalues are returned in
    /// separate arrays, WR for real and WI for complex.  The RWORK
    /// array is only referenced if ScalarType is complex.
    void GEEV(const char JOBVL, const char JOBVR, const OrdinalType n, ScalarType* A, const OrdinalType lda, MagnitudeType* WR, MagnitudeType* WI, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const
      { throw_error("GEEV"); }

    /*! Computes for an \c n by \c n real nonsymmetric matrix \c A, the eigenvalues and, optionally, the left and/or right eigenvectors.
        Optionally, it can compute a balancing transformation to improve the conditioning of the eigenvalues and eigenvectors.
        \note (This is the function is only defined for \c ScalarType = \c float or \c double.)
    */
    void GEEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* WR, ScalarType* WI, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, OrdinalType* ilo, OrdinalType* ihi, MagnitudeType* SCALE, MagnitudeType* abnrm, MagnitudeType* RCONDE, MagnitudeType* RCONDV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* IWORK, OrdinalType* info) const
      { throw_error("GEEVX"); }

    /*! Computes for a pair of \c n by \c n nonsymmetric matrices (\c A,\c B) the generalized eigenvalues, and optionally, the left and/or right generalized eigenvectors.
        Optionally, it can compute a balancing transformation to improve the conditioning of the eigenvalues and eigenvectors.
        \note (This is the function is only defined for \c ScalarType = \c float or \c double.)
    */
    void GGEVX(const char BALANC, const char JOBVL, const char JOBVR, const char SENSE, const OrdinalType n, ScalarType* A, const OrdinalType lda, ScalarType* B, const OrdinalType ldb, MagnitudeType* ALPHAR, MagnitudeType* ALPHAI, ScalarType* BETA, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, OrdinalType* ilo, OrdinalType* ihi, MagnitudeType* lscale, MagnitudeType* rscale, MagnitudeType* abnrm, MagnitudeType* bbnrm, MagnitudeType* RCONDE, MagnitudeType* RCONDV, ScalarType* WORK, const OrdinalType lwork, OrdinalType* IWORK, OrdinalType* BWORK, OrdinalType* info) const
      { throw_error("GGEVX"); }

    /*! Computes for a pair of \c n by \c n nonsymmetric matrices (\c A,\c B) the generalized eigenvalues, and optionally, the left and/or right generalized eigenvectors.
       \note (This is the function is only defined for \c ScalarType = \c float or \c double.)
    */
    void GGEV(const char JOBVL, const char JOBVR, const OrdinalType n, ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb, MagnitudeType *ALPHAR, MagnitudeType *ALPHAI, ScalarType *BETA, ScalarType *VL, const OrdinalType ldvl, ScalarType *VR, const OrdinalType ldvr, ScalarType *WORK, const OrdinalType lwork, OrdinalType *info) const
      { throw_error("GGEV"); }


    /*! Reorders the real Schur factorization of a real matrix so that a selected cluster of eigenvalues appears in the leading diagonal blocks of the upper quasi-triangular matrix \c T, and the leading columns of \c Q form an orthonormal basis of the corresponding right invariant subspace.
       \note (This function is only defined for \c ScalarType = \c float or \c double.)
    */
  void TRSEN(const char JOB, const char COMPQ, const OrdinalType *SELECT, const OrdinalType n, ScalarType *T, const OrdinalType ldt, ScalarType *Q, const OrdinalType ldq, MagnitudeType *WR, MagnitudeType *WI, OrdinalType *M, ScalarType *S, MagnitudeType *SEP, ScalarType *WORK, const OrdinalType lwork, OrdinalType *IWORK, const OrdinalType liwork, OrdinalType *info ) const
    { throw_error("TRSEN"); }


    /*! Reorders the generalized real Schur decomposition of a real matrix pair (\c A, \c B), so that a selected cluster of eigenvalues appears in the leading diagonal blocks of the upper quasi-triangular matrix \c A and the upper triangular \c B.
       \note (This function is only defined for \c ScalarType = \c float or \c double.)
    */
  void TGSEN(const OrdinalType ijob, const OrdinalType wantq, const OrdinalType wantz, const OrdinalType *SELECT, const OrdinalType n, ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb, MagnitudeType *ALPHAR, MagnitudeType *ALPHAI, MagnitudeType *BETA, ScalarType *Q, const OrdinalType ldq, ScalarType *Z, const OrdinalType ldz, OrdinalType *M, MagnitudeType *PL, MagnitudeType *PR, MagnitudeType *DIF, ScalarType *WORK, const OrdinalType lwork, OrdinalType *IWORK, const OrdinalType liwork, OrdinalType *info ) const
    { throw_error("TGSEN"); }


    /*! Computes for a pair of \c n by \c n nonsymmetric matrices (\c A,\c B) the generalized eigenvalues, the generalized real Schur form (\c S,\c T), optionally, the left and/or right matrices of Schur vectors.
       \note (This is the function is only defined for \c ScalarType = \c float or \c double.)
    */
    void GGES(const char JOBVL, const char JOBVR, const char SORT, OrdinalType (*ptr2func)(ScalarType *, ScalarType *, ScalarType *), const OrdinalType n, ScalarType *A, const OrdinalType lda, ScalarType *B, const OrdinalType ldb, OrdinalType *sdim, MagnitudeType *ALPHAR, MagnitudeType *ALPHAI, MagnitudeType *BETA, ScalarType *VL, const OrdinalType ldvl, ScalarType *VR, const OrdinalType ldvr, ScalarType *WORK, const OrdinalType lwork, OrdinalType *BWORK, OrdinalType *info ) const
      { throw_error("GGES"); }

    //@}


    //! @name Singular Value Decompositon Routines
    //@{
    //! Computes the singular values (and optionally, vectors) of a real matrix \c A.
    void GESVD(const char JOBU, const char JOBVT, const OrdinalType m, const OrdinalType n, ScalarType* A, const OrdinalType lda, MagnitudeType* S, ScalarType* U, const OrdinalType ldu, ScalarType* V, const OrdinalType ldv, ScalarType* WORK, const OrdinalType lwork, MagnitudeType* RWORK, OrdinalType* info) const
      { throw_error("GESVD"); }
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
    void ORMQR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("ORMQR"); }

    /// \brief Apply Householder reflectors (complex case).
    ///
    /// Overwrite the general complex \c m by \c n matrix \c C with
    /// the product of \c Q and \c C, where Q is the product of \c k
    /// elementary (Householder) reflectors as returned by \c GEQRF.
    ///
    /// \note This method will call \c ORMQR when ScalarType is real.
    /// (Unitary real matrices are orthogonal.)
    void UNMQR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("UNMQR"); }

    /// \brief Compute explicit Q factor from QR factorization (GEQRF) (real case).
    ///
    /// Generate the \c m by \c n matrix Q with orthonormal columns
    /// corresponding to the first \c n columns of a product of \c k
    /// elementary reflectors of order \c m, as returned by \c GEQRF.
    ///
    /// \note This method is not defined when ScalarType is complex.
    /// Call \c UNGQR in that case.  ("OR" stands for "orthogonal";
    /// "UN" stands for "unitary.")
    void ORGQR(const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("ORGQR"); }

    /// \brief Compute explicit QR factor from QR factorization (GEQRF) (complex case).
    ///
    /// Generate the \c m by \c n matrix Q with orthonormal columns
    /// corresponding tothe first \c n columns of a product of \c k
    /// elementary reflectors of order \c m, as returned by \c GEQRF.
    ///
    /// \note This method will call \c ORGQR when ScalarType is real.
    /// (Unitary real matrices are orthogonal.)
    void UNGQR(const OrdinalType m, const OrdinalType n, const OrdinalType k, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("UNGQR"); }

    /*! \brief Generates a real orthogonal matrix \c Q which is the product of \c ihi-ilo elementary reflectors of order \c n, as returned by GEHRD.  On return \c Q is stored in \c A.
    \note This method is not defined when ScalarType is complex.
    */
    void ORGHR(const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("ORGHR"); }

    /*! \brief Overwrites the general real \c m by \c n matrix \c C with the product of \c C and \c Q, which is a product of \c ihi-ilo elementary reflectors, as returned by GEHRD.
    \note This method is not defined when ScalarType is complex.
    */
    void ORMHR(const char SIDE, const char TRANS, const OrdinalType m, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const ScalarType* A, const OrdinalType lda, const ScalarType* TAU, ScalarType* C, const OrdinalType ldc, ScalarType* WORK, const OrdinalType lwork, OrdinalType* info) const
      { throw_error("ORMHR"); }
    //@}

    //! @name Triangular Matrix Routines
    //@{

    /*! Computes some or all of the right and/or left eigenvectors of an upper triangular matrix \c T. If ScalarType is \c float or \c double, then the matrix is quasi-triangular and arugments \c RWORK is ignored.
    */
    void TREVC(const char SIDE, const char HOWMNY, OrdinalType* select, const OrdinalType n, const ScalarType* T, const OrdinalType ldt, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, ScalarType* WORK, OrdinalType* info) const
      { throw_error("TREVC"); }

    /*! Computes some or all of the right and/or left eigenvectors of an upper triangular matrix \c T. If ScalarType is \c float or \c double, then the matrix is quasi-triangular and arugments \c RWORK is ignored.
       \note (This is the version used for any \c ScalarType, when the user doesn't want to enable the selecting functionality, with HOWMNY='A'.)
    */
    void TREVC(const char SIDE, const OrdinalType n, const ScalarType* T, const OrdinalType ldt, ScalarType* VL, const OrdinalType ldvl, ScalarType* VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType* m, ScalarType* WORK, MagnitudeType* RWORK, OrdinalType* info) const
      { throw_error("TREVC"); }

    /*! Reorders the Schur factorization of a matrix \c T via unitary similarity transformations so that the diagonal element of \c T with row index \c ifst is moved to row \c ilst. If \c ScalarType is \c float or \c double, then \c T should be in real Schur form and the operation affects the diagonal block referenced by \c ifst.
      \note This method will ignore the WORK vector when ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    void TREXC(const char COMPQ, const OrdinalType n, ScalarType* T, const OrdinalType ldt, ScalarType* Q, const OrdinalType ldq, OrdinalType ifst, OrdinalType ilst, ScalarType* WORK, OrdinalType* info) const
      { throw_error("TREXC"); }

    /*! Computes some or all of the right and/or left eigenvectors of a pair of real matrices ( \c S, \c P ), where \c S is a quasi-triangular matrix and \c P is upper triangular.
       \note This method is only defined for \c ScalarType = \c float or \c double.
    */
    void TGEVC(const char SIDE, const char HOWMNY, const OrdinalType *SELECT, const OrdinalType n, ScalarType *S, const OrdinalType lds, ScalarType *P, const OrdinalType ldp, ScalarType *VL, const OrdinalType ldvl, ScalarType *VR, const OrdinalType ldvr, const OrdinalType mm, OrdinalType *M, ScalarType *WORK, OrdinalType *info) const
      { throw_error("TGEVC"); }


    //@}

    //! @name Rotation/Reflection generators
    //@{

    //! Gnerates a plane rotation that zeros out the second component of the input vector.
    void LARTG( const ScalarType f, const ScalarType g, MagnitudeType* c, ScalarType* s, ScalarType* r ) const
      { throw_error("LARTG"); }

    //! Generates an elementary reflector of order \c n that zeros out the last \c n-1 components of the input vector.
    void LARFG( const OrdinalType n, ScalarType* alpha, ScalarType* x, const OrdinalType incx, ScalarType* tau ) const
      { throw_error("LARFG"); }

    //@}

    //! @name Matrix Balancing Routines
    //@{

    //! Balances a general matrix A, through similarity transformations to make the rows and columns as close in norm as possible.
    void GEBAL(const char JOBZ, const OrdinalType n, ScalarType* A, const OrdinalType lda, OrdinalType ilo, OrdinalType ihi, MagnitudeType* scale, OrdinalType* info) const
      { throw_error("GEBAL"); }

    //! Forms the left or right eigenvectors of a general matrix that has been balanced by GEBAL by backward transformation of the computed eigenvectors \c V.
    void GEBAK(const char JOBZ, const char SIDE, const OrdinalType n, const OrdinalType ilo, const OrdinalType ihi, const MagnitudeType* scale , const OrdinalType m, ScalarType* V, const OrdinalType ldv, OrdinalType* info) const
      { throw_error("GEBAK"); }

    //@}

    //! @name Random number generators
    //@{
    //! Returns a random number from a uniform or normal distribution.
    ScalarType LARND( const OrdinalType idist, OrdinalType* seed ) const
      { throw_error("LARND"); return ScalarType(0); }

    //! Returns a vector of random numbers from a chosen distribution.
    void LARNV( const OrdinalType idist, OrdinalType* seed, const OrdinalType n, ScalarType* v ) const
      { throw_error("LARNV"); }
    //@}

    //! @name Machine Characteristics Routines.
    //@{
    /*! \brief Determines machine parameters for floating point characteristics.
        \note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    ScalarType LAMCH(const char CMACH) const
      { throw_error("LAMCH"); return ScalarType(0); }

    /*! \brief Chooses problem-dependent parameters for the local environment.
        \note This method should give parameters for good, but not optimal, performance on many currently
        available computers.
    */
    OrdinalType ILAENV( const OrdinalType ispec, const std::string& NAME, const std::string& OPTS, const OrdinalType N1 = -1, const OrdinalType N2 = -1, const OrdinalType N3 = -1, const OrdinalType N4 = -1 ) const
      { throw_error("ILAENV"); return OrdinalType(0); }
    //@}

    //! @name Miscellaneous Utilities.
    //@{
    /*! \brief Computes x^2 + y^2 safely, to avoid overflow.
        \note This method is not defined when the ScalarType is \c std::complex<float> or \c std::complex<double>.
    */
    ScalarType LAPY2(const ScalarType x, const ScalarType y) const
      { throw_error("LAPY2"); return ScalarType(0); }
    //@}

  private:

    void throw_error(const char *func) const {
      TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::logic_error,
        func << ":  Not implemented for Sacado::UQ::PCE scalar type!");
    }

  };

} // namespace Teuchos

#endif // _TEUCHOS_LAPACK__UQ_PCE_HPP_
