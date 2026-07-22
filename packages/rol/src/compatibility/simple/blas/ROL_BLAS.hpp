// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BLAS_H
#define ROL_BLAS_H

/** \class ROL::BLAS
  \brief Provides interface to BLAS
  */

/* A) Define PREFIX and ROL_fcd based on platform. */


// Following two macros were stolen from TEUCHOS_LAPACK.cpp
/* for INTEL_CXML, the second arg may need to be changed to 'one'.  If so
the appropriate declaration of one will need to be added back into
functions that include the macro:
*/
#ifdef CHAR_MACRO
#undef CHAR_MACRO
#endif
#if defined (INTEL_CXML)
#define CHAR_MACRO(char_var) &char_var, one
#else
#define CHAR_MACRO(char_var) &char_var
#endif

// Following macros were stolen from TEUCHOS_LAPACK_wrappers.hpp
#if defined(INTEL_CXML)
#  define PREFIX __stdcall
#  define ROL_fcd const char *, unsigned int
#elif defined(INTEL_MKL)
#  define PREFIX
#  define ROL_fcd const char *
#else
#  define PREFIX
#  define ROL_fcd const char *
#endif

// in progress - added for EQUED which is a modified char *, not const
#if defined(INTEL_CXML)
#  define PREFIX __stdcall
#  define ROL_nonconst_fcd char *, unsigned int // Need to evaluate unsigned int - CXML deprecated
#elif defined(INTEL_MKL)
#  define PREFIX
#  define ROL_nonconst_fcd char *
#else
#  define PREFIX
#  define ROL_nonconst_fcd char *
#endif

#define DROTG_F77   F77_BLAS_MANGLE(drotg,DROTG)
#define DROT_F77    F77_BLAS_MANGLE(drot,DROT)
#define DASUM_F77   F77_BLAS_MANGLE(dasum,DASUM)
#define DAXPY_F77   F77_BLAS_MANGLE(daxpy,DAXPY)
#define DCOPY_F77   F77_BLAS_MANGLE(dcopy,DCOPY)
#define DDOT_F77    F77_BLAS_MANGLE(ddot,DDOT)
#define DNRM2_F77   F77_BLAS_MANGLE(dnrm2,DNRM2)
#define DSCAL_F77   F77_BLAS_MANGLE(dscal,DSCAL)
#define IDAMAX_F77  F77_BLAS_MANGLE(idamax,IDAMAX)
#define DGEMV_F77   F77_BLAS_MANGLE(dgemv,DGEMV)
#define DGER_F77    F77_BLAS_MANGLE(dger,DGER)
#define DTRMV_F77   F77_BLAS_MANGLE(dtrmv,DTRMV)
#define DGEMM_F77   F77_BLAS_MANGLE(dgemm,DGEMM)
#define DSWAP_F77   F77_BLAS_MANGLE(dswap,DSWAP)
#define DSYMM_F77   F77_BLAS_MANGLE(dsymm,DSYMM)
#define DSYRK_F77   F77_BLAS_MANGLE(dsyrk,DSYRK)
#define DTRMM_F77   F77_BLAS_MANGLE(dtrmm,DTRMM)
#define DTRSM_F77   F77_BLAS_MANGLE(dtrsm,DTRSM)

namespace ROL { 
  extern "C" {
    void PREFIX DROTG_F77(double* da, double* db, double* c, double* s);
    void PREFIX DROT_F77(const int* n, double* dx, const int* incx, double* dy, const int* incy, double* c, double* s);
    double PREFIX DASUM_F77(const int* n, const double x[], const int* incx);
    void PREFIX DAXPY_F77(const int* n, const double* alpha, const double x[], const int* incx, double y[], const int* incy);
    void PREFIX DCOPY_F77(const int* n, const double *x, const int* incx, double *y, const int* incy);
    double PREFIX DDOT_F77(const int* n, const double x[], const int* incx, const double y[], const int* incy);
    double PREFIX DNRM2_F77(const int* n, const double x[], const int* incx);
    void PREFIX DSCAL_F77(const int* n, const double* alpha, double *x, const int* incx);
    void PREFIX DSWAP_F77(const int* const n, double* const x, const int* const incx,
                          double* const y, const int* const incy);
    int PREFIX IDAMAX_F77(const int* n, const double *x, const int* incx);
    void PREFIX DGEMV_F77(Teuchos_fcd, const int* m, const int* n, const double* alpha, const double A[], const int* lda,const double x[], const int* incx, const double* beta, double y[], const int* incy);
    void PREFIX DTRMV_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *n,const double *a, const int *lda, double *x, const int *incx);
    void PREFIX DGER_F77(const int *m, const int *n, const double *alpha, const double *x, const int *incx, const double *y,const int *incy, double *a, const int *lda);
    void PREFIX DGEMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *lda,
                    const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
    void PREFIX DSYMM_F77(Teuchos_fcd, Teuchos_fcd, const int *m, const int * n, const double *alpha, const double *a, const int *lda,
                    const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
    void PREFIX DSYRK_F77(Teuchos_fcd, Teuchos_fcd, const int *n, const int * k, const double *alpha, const double *a, const int *lda, const double *beta, double *c, const int *ldc);
    void PREFIX DTRMM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *m, const int *n, const double *alpha, const double *a, const int * lda, double *b, const int *ldb);
    void PREFIX DTRSM_F77(Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, Teuchos_fcd, const int *m, const int *n, const double *alpha, const double *a, const int *lda, double *b, const int *ldb);

  template<typename Index, typename Real>
  struct BLAS { 
    //! Computes a Givens plane rotation.
    void ROTG(ScalarType* da, ScalarType* db, rotg_c_type* c, ScalarType* s) const;

    //! Applies a Givens plane rotation.
    void ROT(const OrdinalType& n, ScalarType* dx, const OrdinalType& incx, ScalarType* dy, const OrdinalType& incy, MagnitudeType* c, ScalarType* s) const;

    //! Scale the vector \c x by the constant \c alpha.
    void SCAL(const OrdinalType& n, const ScalarType& alpha, ScalarType* x, const OrdinalType& incx) const;

    //! Copy the vector \c x to the vector \c y.
    void COPY(const OrdinalType& n, const ScalarType* x, const OrdinalType& incx, ScalarType* y, const OrdinalType& incy) const;

    //! Perform the operation: \c y \c <- \c y+alpha*x.
    template <typename alpha_type, typename x_type>
    void AXPY(const OrdinalType& n, const alpha_type alpha, const x_type* x, const OrdinalType& incx, ScalarType* y, const OrdinalType& incy) const;

    //! Sum the absolute values of the entries of \c x.
    typename ScalarTraits<ScalarType>::magnitudeType ASUM(const OrdinalType& n, const ScalarType* x, const OrdinalType& incx) const;

    //! Form the dot product of the vectors \c x and \c y.
    template <typename x_type, typename y_type>
    ScalarType DOT(const OrdinalType& n, const x_type* x, const OrdinalType& incx, const y_type* y, const OrdinalType& incy) const;

    //! Compute the 2-norm of the vector \c x.
    typename ScalarTraits<ScalarType>::magnitudeType NRM2(const OrdinalType& n, const ScalarType* x, const OrdinalType& incx) const;

    //! Return the index of the element of \c x with the maximum magnitude.
    OrdinalType IAMAX(const OrdinalType& n, const ScalarType* x, const OrdinalType& incx) const;
    //@}

    //! @name Level 2 BLAS Routines.
    //@{

    //! Performs the matrix-vector operation:  \c y \c <- \c alpha*A*x+beta*y or \c y \c <- \c alpha*A'*x+beta*y where \c A is a general \c m by \c n matrix.
    template <typename alpha_type, typename A_type, typename x_type, typename beta_type>
    void GEMV(ETransp trans, const OrdinalType& m, const OrdinalType& n, const alpha_type alpha, const A_type* A,
              const OrdinalType& lda, const x_type* x, const OrdinalType& incx, const beta_type beta, ScalarType* y, const OrdinalType& incy) const;

    //! Performs the matrix-vector operation:  \c x \c <- \c A*x or \c x \c <- \c A'*x where \c A is a unit/non-unit \c n by \c n upper/lower triangular matrix.
    template <typename A_type>
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const OrdinalType& n, const A_type* A,
              const OrdinalType& lda, ScalarType* x, const OrdinalType& incx) const;

    //! \brief Performs the rank 1 operation:  \c A \c <- \c alpha*x*y'+A.
    /// \note  For complex arithmetic, this routine performs [Z/C]GERU.
    template <typename alpha_type, typename x_type, typename y_type>
    void GER(const OrdinalType& m, const OrdinalType& n, const alpha_type alpha, const x_type* x, const OrdinalType& incx,
             const y_type* y, const OrdinalType& incy, ScalarType* A, const OrdinalType& lda) const;
    //@}

    //! @name Level 3 BLAS Routines.
    //@{

    /// \brief General matrix-matrix multiply.
    ///
    /// This computes C = alpha*op(A)*op(B) + beta*C.  op(X) here may
    /// be either X, the transpose of X, or the conjugate transpose of
    /// X.  op(A) has m rows and k columns, op(B) has k rows and n
    /// columns, and C has m rows and n columns.
    template <typename alpha_type, typename A_type, typename B_type, typename beta_type>
    void GEMM(ETransp transa, ETransp transb, const OrdinalType& m, const OrdinalType& n, const OrdinalType& k, const alpha_type alpha, const A_type* A, const OrdinalType& lda, const B_type* B, const OrdinalType& ldb, const beta_type beta, ScalarType* C, const OrdinalType& ldc) const;

    //! Swap the entries of x and y.
    void
    SWAP (const OrdinalType& n, ScalarType* const x, const OrdinalType& incx,
          ScalarType* const y, const OrdinalType& incy) const;

    //! Performs the matrix-matrix operation: \c C \c <- \c alpha*A*B+beta*C or \c C \c <- \c alpha*B*A+beta*C where \c A is an \c m by \c m or \c n by \c n symmetric matrix and \c B is a general matrix.
    template <typename alpha_type, typename A_type, typename B_type, typename beta_type>
    void SYMM(ESide side, EUplo uplo, const OrdinalType& m, const OrdinalType& n, const alpha_type alpha, const A_type* A, const OrdinalType& lda, const B_type* B, const OrdinalType& ldb, const beta_type beta, ScalarType* C, const OrdinalType& ldc) const;

    //! Performs the symmetric rank k operation: \c C <- \c alpha*A*A'+beta*C or \c C <- \c alpha*A'*A+beta*C, where \c alpha and \c beta are scalars, \c C is an \c n by \c n symmetric matrix and \c A is an \c n by \c k matrix in the first case or \c k by \c n matrix in the second case.
    template <typename alpha_type, typename A_type, typename beta_type>
    void SYRK(EUplo uplo, ETransp trans, const OrdinalType& n, const OrdinalType& k, const alpha_type alpha, const A_type* A, const OrdinalType& lda, const beta_type beta, ScalarType* C, const OrdinalType& ldc) const;

    //! Performs the matrix-matrix operation: \c B \c <- \c alpha*op(A)*B or \c B \c <- \c alpha*B*op(A) where \c op(A) is an unit/non-unit, upper/lower triangular matrix and \c B is a general matrix.
    template <typename alpha_type, typename A_type>
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType& m, const OrdinalType& n,
                const alpha_type alpha, const A_type* A, const OrdinalType& lda, ScalarType* B, const OrdinalType& ldb) const;

    //! Solves the matrix equations:  \c op(A)*X=alpha*B or \c X*op(A)=alpha*B where \c X and \c B are \c m by \c n matrices, \c A is a unit/non-unit, upper/lower triangular matrix and \c op(A) is \c A or \c A'.  The matrix \c X is overwritten on \c B.
    template <typename alpha_type, typename A_type>
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const OrdinalType& m, const OrdinalType& n,
                const alpha_type alpha, const A_type* A, const OrdinalType& lda, ScalarType* B, const OrdinalType& ldb) const;
    //@}
  };

  template<>
  struct BLAS<int, double> { 

    void ROTG(double* da, double* db, double* c, double* s) const
    { DROTG_F77(da, db, c, s); }
  
    void ROT(const int& n, double* dx, const int& incx, double* dy, const int& incy, double* c, double* s) const
    { DROT_F77(&n, dx, &incx, dy, &incy, c, s); }
  
    double ASUM(const int& n, const double* x, const int& incx) const
    { return DASUM_F77(&n, x, &incx); }
  
    void AXPY(const int& n, const double& alpha, const double* x, const int& incx, double* y, const int& incy) const
    { DAXPY_F77(&n, &alpha, x, &incx, y, &incy); }
  
    void COPY(const int& n, const double* x, const int& incx, double* y, const int& incy) const
    { DCOPY_F77(&n, x, &incx, y, &incy); }
  
    double DOT(const int& n, const double* x, const int& incx, const double* y, const int& incy) const
    {
      return DDOT_F77(&n, x, &incx, y, &incy);
    }
  
    int IAMAX(const int& n, const double* x, const int& incx) const
    { return IDAMAX_F77(&n, x, &incx); }
  
    double NRM2(const int& n, const double* x, const int& incx) const
    { return DNRM2_F77(&n, x, &incx); }
  
    void SCAL(const int& n, const double& alpha, double* x, const int& incx) const
    { DSCAL_F77(&n, &alpha, x, &incx); }
  
    void GEMV(ETransp trans, const int& m, const int& n, const double& alpha, const double* A, const int& lda, const double* x, const int& incx, const double& beta, double* y, const int& incy) const
    { DGEMV_F77(CHAR_MACRO(ETranspChar[trans]), &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy); }
  
    void GER(const int& m, const int& n, const double& alpha, const double* x, const int& incx, const double* y, const int& incy, double* A, const int& lda) const
    { DGER_F77(&m, &n, &alpha, x, &incx, y, &incy, A, &lda); }
  
    void TRMV(EUplo uplo, ETransp trans, EDiag diag, const int& n, const double* A, const int& lda, double* x, const int& incx) const
    { DTRMV_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), CHAR_MACRO(EDiagChar[diag]), &n, A, &lda, x, &incx); }
  
    void GEMM(ETransp transa, ETransp transb, const int& m, const int& n, const int& k, const double& alpha, const double* A, const int& lda, const double* B, const int& ldb, const double& beta, double* C, const int& ldc) const
    { DGEMM_F77(CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(ETranspChar[transb]), &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
    void SWAP(const int& n, double* const x, const int& incx, double* const y, const int& incy) const
    {
      DSWAP_F77 (&n, x, &incx, y, &incy);
    }
  
    void SYMM(ESide side, EUplo uplo, const int& m, const int& n, const double& alpha, const double* A, const int& lda, const double* B, const int& ldb, const double& beta, double* C, const int& ldc) const
    { DSYMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc); }
  
    void SYRK(EUplo uplo, ETransp trans, const int& n, const int& k, const double& alpha, const double* A, const int& lda, const double& beta, double* C, const int& ldc) const
    { DSYRK_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), &n, &k, &alpha, A, &lda, &beta, C, &ldc); }
  
    void HERK(EUplo uplo, ETransp trans, const int& n, const int& k, const double& alpha, const double* A, const int& lda, const double& beta, double* C, const int& ldc) const
    { DSYRK_F77(CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[trans]), &n, &k, &alpha, A, &lda, &beta, C, &ldc); }
  
    void TRMM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int& m, const int& n, const double& alpha, const double* A, const int& lda, double* B, const int& ldb) const
    { DTRMM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }
  
    void TRSM(ESide side, EUplo uplo, ETransp transa, EDiag diag, const int& m, const int& n, const double& alpha, const double* A, const int& lda, double* B, const int& ldb) const
    { DTRSM_F77(CHAR_MACRO(ESideChar[side]), CHAR_MACRO(EUploChar[uplo]), CHAR_MACRO(ETranspChar[transa]), CHAR_MACRO(EDiagChar[diag]), &m, &n, &alpha, A, &lda, B, &ldb); }
  };

}

#endif
