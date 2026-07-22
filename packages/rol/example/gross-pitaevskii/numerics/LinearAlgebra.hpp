// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include"Teuchos_LAPACK.hpp"
#include"ROL_StdVector.hpp"

#ifndef __LINEAR_ALGEBRA__
#define __LINEAR_ALGEBRA__

/** \brief Solve a tridiagonal system 
  @param[in]    lapack   pointer to LAPACK interface 
  @param[in]    a        subdiagonal band    ( length N-1 )
  @param[in]    b        diagonal band       ( length N   )
  @param[in]    c        superidiagonal band ( length N-1 )
  @param[in]    r        right-hand-side     ( length N   )
  @param[out]   x        solution vector     ( length N   )   */
template<class Real>
void trisolve(ROL::Ptr<Teuchos::LAPACK<int,Real> > lapack,
              const std::vector<Real>& a,
              const std::vector<Real>& b,
              const std::vector<Real>& c,
              const std::vector<Real>& r,
                    std::vector<Real>& x ) {

    const char TRANS = 'N'; 
    const int  N = b.size();
    const int  NRHS = 1;
    int info;
   
    std::vector<Real> dl(a);
    std::vector<Real> d(b);
    std::vector<Real> du(c);
    std::vector<Real> du2(N-2,0.0);
    std::vector<int>  ipiv(N);
   
    // Do matrix factorization, overwriting the LU bands 
    lapack->GTTRF(N,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&info);
   
    x = r;
    
    // Solve the system using the banded LU factors 
    lapack->GTTRS(TRANS,N,NRHS,&dl[0],&d[0],&du[0],&du2[0],&ipiv[0],&x[0],N,&info);
}


/** Solve a general system AX=B, with optional multiple (M) right-hand-sides using the LU factorization
  @param[in]    lapack   pointer to LAPACK interface 
  @param[in]    A        vector of column-stacked LHS matrix elements      ( length N*N ) 
  @param[in]    B        vector of column-stacked RHS matrix elements      ( length N*M ) 
  @param[in]    X        vector of column-stacked solution matrix elements ( length N*M ) 
 */
template<class Real>
void lusolve(ROL::Ptr<Teuchos::LAPACK<int,Real> > lapack,
             const std::vector<Real> &A,
             const std::vector<Real> &B,
                   std::vector<Real> &X) {

    const char TRANS = 'N';  
    const int  N2 = A.size();
    const int  N  = sqrt(N2);
    const int  M  = int(B.size()/N);
    const int  LDA = N;
    int info;
    std::vector<int> ipiv(N);
    std::vector<Real> PLU(A);
    X = B; 

    // Do matrix factorization
    lapack->GETRF(N,N,&PLU[0],LDA,&ipiv[0],&info);

    // Solve LU-factored system 
    lapack->GETRS(TRANS,N,M,&PLU[0],LDA,&ipiv[0],&X[0],LDA,&info);
}   
               


/** Solve a general system AX=B, with optional multiple (M) right-hand-sides using the Cholesky factorization
  @param[in]    lapack   pointer to LAPACK interface 
  @param[in]    A        vector of column-stacked LHS matrix elements      ( length N*N ) 
  @param[in]    B        vector of column-stacked RHS matrix elements      ( length N*M ) 
  @param[in]    X        vector of column-stacked solution matrix elements ( length N*M ) 
 */
template<class Real>
void cholsolve(ROL::Ptr<Teuchos::LAPACK<int,Real> > lapack,
               const std::vector<Real> &A,
               const std::vector<Real> &B,
                     std::vector<Real> &X) {

    const char UPLO = 'U';
    const int  N2   = A.size();
    const int  N    = sqrt(N2);
    const int  M    = int(B.size()/N);
    const int  LDA  = N;
    int info;
    std::vector<Real> C(A);
    X = B; 

    // Do Cholesky matrix factorization
    lapack->POTRF(UPLO,N,&C[0],LDA,&info);

    // Solve Cholesky-factored system 
    lapack->POTRS(UPLO,N,&C[0],LDA,&X[0],LDA,&info);
}   
 
#endif


