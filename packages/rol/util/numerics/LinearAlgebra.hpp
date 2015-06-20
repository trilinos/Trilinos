// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
void trisolve(Teuchos::LAPACK<int,Real> &lapack,
              const ROL::Vector<Real>& a,
              const ROL::Vector<Real>& b,
              const ROL::Vector<Real>& c,
              const ROL::Vector<Real>& r,
                    ROL::Vector<Real>& x ) {


    Teuchos::RCP<const std::vector<Real> > ap = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(a))).getVector();
    Teuchos::RCP<const std::vector<Real> > bp = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(b))).getVector();
    Teuchos::RCP<const std::vector<Real> > cp = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(c))).getVector();
    Teuchos::RCP<const std::vector<Real> > rp = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(r))).getVector();
    Teuchos::RCP<std::vector<Real> > xp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(x)).getVector()); 

    const char TRANS = 'N'; 
    const int  N = bp->size();
    const int  NRHS = 1;
    int info;
   
    Teuchos::RCP<std::vector<Real> > dlp   = Teuchos::rcp( new std::vector<Real>(*ap) );
    Teuchos::RCP<std::vector<Real> > dp    = Teuchos::rcp( new std::vector<Real>(*bp) );
    Teuchos::RCP<std::vector<Real> > dup   = Teuchos::rcp( new std::vector<Real>(*cp) );
    Teuchos::RCP<std::vector<Real> > du2p  = Teuchos::rcp( new std::vector<Real>(N-2,0.0) );
    Teuchos::RCP<std::vector<int>  > ipivp = Teuchos::rcp( new std::vector<int>(N,0.0) );
   
    // Do matrix factorization, overwriting the LU bands 
    lapack.GTTRF(N,&(*dlp)[0],&(*dp)[0],&(*dup)[0],&(*du2p)[0],&(*ipivp)[0],&info);
   
    std::copy(rp->begin(),rp->end(),xp->begin());
    
    // Solve the system using the banded LU factors 
    lapack.GTTRS(TRANS,N,NRHS,&(*dlp)[0],&(*dp)[0],&(*dup)[0],&(*du2p)[0],&(*ipivp)[0],&(*xp)[0],N,&info);
}


/** Solve a general system AX=B, with optional multiple (M) right-hand-sides using the LU factorization
  @param[in]    lapack   pointer to LAPACK interface 
  @param[in]    A        vector of column-stacked LHS matrix elements      ( length N*N ) 
  @param[in]    B        vector of column-stacked RHS matrix elements      ( length N*M ) 
  @param[in]    X        vector of column-stacked solution matrix elements ( length N*M ) 
 */
template<class Real>
void lusolve(Teuchos::LAPACK<int,Real> &lapack,
             const ROL::Vector<Real> &A,
             const ROL::Vector<Real> &B,
                   ROL::Vector<Real> &X) {

   Teuchos::RCP<const std::vector<Real> > Ap = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(A))).getVector();
    Teuchos::RCP<const std::vector<Real> > Bp = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(B))).getVector();
    Teuchos::RCP<std::vector<Real> > Xp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(X)).getVector()); 

    const char TRANS = 'N';  
    const int  N2 = Ap->size();
    const int  N  = sqrt(N2);
    const int  M  = int(Bp->size()/N);
    const int  LDA = N;
    int info;
    Teuchos::RCP<std::vector<int> >  ipivp = Teuchos::rcp(new std::vector<int> (N));
    Teuchos::RCP<std::vector<Real> > PLUp  = Teuchos::rcp(new std::vector<Real> (*Ap));

    std::copy(Bp->begin(),Bp->end(),Xp->begin());

    // Do matrix factorization
    lapack.GETRF(N,N,&(*PLUp)[0],LDA,&(*ipivp)[0],&info);

    // Solve LU-factored system 
    lapack.GETRS(TRANS,N,M,&(*PLUp)[0],LDA,&(*ipivp)[0],&(*Xp)[0],LDA,&info);
}   
               


/** Solve a general system AX=B, with optional multiple (M) right-hand-sides using the Cholesky factorization
  @param[in]    lapack   pointer to LAPACK interface 
  @param[in]    A        vector of column-stacked LHS matrix elements      ( length N*N ) 
  @param[in]    B        vector of column-stacked RHS matrix elements      ( length N*M ) 
  @param[in]    X        vector of column-stacked solution matrix elements ( length N*M ) 
 */
template<class Real>
void cholsolve(Teuchos::LAPACK<int,Real> &lapack,
               const ROL::Vector<Real> &A,
               const ROL::Vector<Real> &B,
                     ROL::Vector<Real> &X) {

    Teuchos::RCP<const std::vector<Real> > Ap = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(A))).getVector();
    Teuchos::RCP<const std::vector<Real> > Bp = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(B))).getVector();
    Teuchos::RCP<std::vector<Real> > Xp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(X)).getVector()); 

    const char UPLO = 'U';
    const int  N2   = Ap->size();
    const int  N    = sqrt(N2);
    const int  M    = int(Bp->size()/N);
    const int  LDA  = N;
    int info;

    Teuchos::RCP<std::vector<Real> > Cp = Teuchos::rcp(new std::vector<Real>(*Ap));

    std::copy(Bp->begin(),Bp->end(),Xp->begin());

    // Do Cholesky matrix factorization
    lapack.POTRF(UPLO,N,&(*Cp)[0],LDA,&info);

    // Solve Cholesky-factored system 
    lapack.POTRS(UPLO,N,&(*Cp)[0],LDA,&(*Xp)[0],LDA,&info);
}   
 
#endif


