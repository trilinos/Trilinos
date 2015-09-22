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

#include<cmath>
#include"ROL_StdVector.hpp"
#include"Teuchos_LAPACK.hpp"
#include"LinearAlgebra.hpp"
 

#ifndef __ORTHOGONAL_POLYNOMIALS__
#define __ORTHOGONAL_POLYNOMIALS__
           
/** \brief Generate the Jacobi polynomial recursion coeffcients \f$a_k,b_k\f$
    
    The Jacobi polynomials satisfy the recurrence relation
    \f[ P^{(\alpha,\beta)}_{k+1}(x) = (x-a_k)P_{k}^{(\alpha,\beta)}(x) - b_k P_{k-1}^{(\alpha,\beta)}(x) \f]
    and form an orthogonal basis on \f$[-1,1]\f$ with respect to the weight function 
    \f$w(x)=(1-x)^\alpha(1+x)^\beta\f$.

    @param[in]    alpha    is a parameter that defines the weight function
    @param[in]    beta     is a parameter that defines the weight function
    @param[out]   ap       is a vector of recursion coefficients
    @param[out]   bp       is a vector of recursion coefficients
  
    Adapted from the MATLAB code by Dirk Laurie and Walter Gautschi
    http://www.cs.purdue.edu/archives/2002/wxg/codes/r_jacobi.m  */
template<class Real>
void rec_jacobi( const double alpha, 
                 const double beta,
                 ROL::Vector<Real> &a,
                 ROL::Vector<Real> &b ) {

    Teuchos::RCP<std::vector<Real> > ap = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(a)).getVector()); 
    Teuchos::RCP<std::vector<Real> > bp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(b)).getVector()); 

    int N = ap->size();
    Real nu = (beta-alpha)/double(alpha+beta+2.0);
    Real mu  = pow(2.0,alpha+beta+1.0)*tgamma(alpha+1.0)*tgamma(beta+1.0)/tgamma(alpha+beta+2.0);
    Real nab;
    Real sqdif = pow(beta,2)-pow(alpha,2);
                              
    (*ap)[0] = nu;
    (*bp)[0] = mu;
    
    if(N>1){
     
        for(int n=1;n<N;++n) {
            nab = 2*n+alpha+beta;
            (*ap)[n] = sqdif/(nab*(nab+2));              
        }

        (*bp)[1] = 4.0*(alpha+1.0)*(beta+1.0)/(pow(alpha+beta+2.0,2)*(alpha+beta+3.0));

        if(N>2) {
            for(int n=2;n<N;++n) {
                nab = 2*n+alpha+beta;
                (*bp)[n] = 4.0*(n+alpha)*(n+beta)*n*(n+alpha+beta)/(nab*nab*(nab+1.0)*(nab-1.0));
            }
        }
    }    
}

/** \brief Construct the generalized Vandermonde matrix (in column stacked form)
           based upon the recurrence coefficients (a,b) on the grid x 

    @param[in]   a   vector of recursion coefficients
    @param[in]   b   vector of recursion coefficients
    @param[in]   x   vector of quadrature nodes
    @param[in]   V   column-stacked Vandermonde matrix
*/   

template<class Real>
void vandermonde( const ROL::Vector<Real> &a,
                  const ROL::Vector<Real> &b,
                  ROL::Vector<Real> &x,
                  ROL::Vector<Real> &V) {

    Teuchos::RCP<const std::vector<Real> > ap = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(a))).getVector();
    Teuchos::RCP<const std::vector<Real> > bp = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(b))).getVector();
    Teuchos::RCP<const std::vector<Real> > xp = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
    Teuchos::RCP<std::vector<Real> > Vp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(V)).getVector()); 

    const int n = ap->size();

    for(int i=0;i<n;++i) {
        // Zeroth polynomial is always 1
        (*Vp)[i]   = 1.0;
        // First polynomial is (x-a) 
        (*Vp)[i+n] = (*xp)[i] - (*ap)[0]; 
    }
    for(int j=1;j<n-1;++j) {
        for(int i=0;i<n;++i) { 
            (*Vp)[i+(j+1)*n] = ((*xp)[i] - (*ap)[j])*(*Vp)[i+j*n] - (*bp)[j]*(*Vp)[i+(j-1)*n];  
        }      
    }        
}


/** \brief Compute the Gauss quadrature nodes and weights for the polynomials generated by the 
           recurrence coefficients.
    
    @param[in]   lapack   pointer to the Teuchos::LAPACK interface
    @param[in]   a        vector of recursion coefficients
    @param[in]   b        vector of recursion coefficients
    @param[out]  x        vector of quadrature nodes
    @param[out]  w        vector of quadrature weights

    Adapted from the MATLAB code by Walter Gautschi
    http://www.cs.purdue.edu/archives/2002/wxg/codes/gauss.m  */
template<class Real>
void gauss( const Teuchos::LAPACK<int,Real> &lapack,
            const ROL::Vector<Real> &a,
            const ROL::Vector<Real> &b,
            ROL::Vector<Real> &x,
            ROL::Vector<Real> &w ) {
    int INFO;

    Teuchos::RCP<const std::vector<Real> > ap = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(a))).getVector();
    Teuchos::RCP<const std::vector<Real> > bp = 
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(b))).getVector();
    Teuchos::RCP<std::vector<Real> > xp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(x)).getVector()); 
    Teuchos::RCP<std::vector<Real> > wp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(w)).getVector()); 

    const int N = ap->size();  
    const int LDZ = N;
    const char COMPZ = 'I';

    Teuchos::RCP<std::vector<Real> > Dp = Teuchos::rcp(new std::vector<Real>(N,0.0));
    Teuchos::RCP<std::vector<Real> > Ep = Teuchos::rcp(new std::vector<Real>(N,0.0));
    Teuchos::RCP<std::vector<Real> > WORKp = Teuchos::rcp(new std::vector<Real>(4*N,0.0));

    // Column-stacked matrix of eigenvectors needed for weights
    Teuchos::RCP<std::vector<Real> > Zp = Teuchos::rcp(new std::vector<Real>(N*N,0));

    // D = a
    std::copy(ap->begin(),ap->end(),Dp->begin());
     
    for(int i=0;i<N-1;++i) {
        (*Ep)[i] = sqrt((*bp)[i+1]);  
    }

    // Eigenvalue Decomposition 
    lapack.STEQR(COMPZ,N,&(*Dp)[0],&(*Ep)[0],&(*Zp)[0],LDZ,&(*WORKp)[0],&INFO);

    for(int i=0;i<N;++i){
        (*xp)[i] = (*Dp)[i];
        (*wp)[i] = (*bp)[0]*pow((*Zp)[N*i],2);
    } 
}



/** \brief Modify the given recurrence coefficients so that the set of zeros of the maximal
           order polynomial include the two prescribed points
    
    @param[in]      lapack   pointer to the Teuchos::LAPACK interface
    @param[in]      xl1      location of one pre-assigned node
    @param[in]      xl2      location of another pre-assigned node
    @param[in/out]  ap       pointer to vector of recursion coefficients
    @param[in/out]  bp       pointer to vector of recursion coefficients

    Based on the section 7 of the paper 
    "Some modified matrix eigenvalue problems" 
    by Gene Golub, SIAM Review Vol 15, No. 2, April 1973, pp.318--334 */
template<class Real>
void rec_lobatto( Teuchos::LAPACK<int,Real> &lapack,
                  const double xl1, 
                  const double xl2,
                  ROL::Vector<Real> &a,
                  ROL::Vector<Real> &b ) {

    Teuchos::RCP<std::vector<Real> > ap = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(a)).getVector()); 
    Teuchos::RCP<std::vector<Real> > bp = 
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(b)).getVector()); 


    const int N = ap->size()-1;

    Teuchos::RCP<std::vector<Real> > amodp = Teuchos::rcp(new std::vector<Real> (N,0.0));
    Teuchos::RCP<std::vector<Real> > bmodp = Teuchos::rcp(new std::vector<Real> (N-1,0.0));
    Teuchos::RCP<std::vector<Real> > enp   = Teuchos::rcp(new std::vector<Real> (N,0.0));
    Teuchos::RCP<std::vector<Real> > gp    = Teuchos::rcp(new std::vector<Real> (N,0.0));

    // Nth canonical vector
    (*enp)[N-1] = 1.0;

    for(int i=0;i<N-1;++i) {
        (*bmodp)[i] = sqrt((*bp)[i+1]);
    }

    for(int i=0;i<N;++i) {
        (*amodp)[i] = (*ap)[i]-xl1;
    }
    
    ROL::StdVector<Real> amod(amodp);  
    ROL::StdVector<Real> bmod(bmodp);  
    ROL::StdVector<Real> en(enp);  
    ROL::StdVector<Real> g(gp);  

    trisolve(lapack,bmod,amod,bmod,en,g);         
    Real g1 = (*gp)[N-1];

    for(int i=0;i<N;++i) {
        (*amodp)[i] = (*ap)[i]-xl2;
    }
    

    trisolve(lapack,bmod,amod,bmod,en,g);         
    Real g2 = (*gp)[N-1];

    (*ap)[N] = (g1*xl2-g2*xl1)/(g1-g2);
    (*bp)[N] = (xl2-xl1)/(g1-g2);
    
}

 
#endif

             
