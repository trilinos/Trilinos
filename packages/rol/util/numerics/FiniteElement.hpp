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


#ifndef __FINITE_ELEMENT__
#define __FINITE_ELEMENT__

#include <string>
#include <iostream>

#include "Sacado.hpp"  // IMPORTANT: Sacado.hpp must precede ROL_StdVector.hpp
#include "ROL_StdVector.hpp"

#include "Lagrange.hpp"
#include "OrthogonalPolynomials.hpp"

#include "Teuchos_LAPACK.hpp"

template<class Real>
void printvec(std::vector<Real> v, const std::string &name ){
    int n = v.size();
    std::cout << "-------------------" << std::endl;
    std::cout << "Printing " << name << std::endl;
    for(int i=0;i<n;++i) {
        std::cout << v[i] << std::endl; 
    }
    std::cout << "-------------------" << std::endl;
}

using namespace ROL;

//! \brief Nodal basis functions on the reference element [-1,1]
template<class Real>
struct NodalBasis {

    // \param ni_ Number of interpolation points/basis functions
    const int ni_;

    // \param nq_ Number of quadrature points
    const int nq_;

    // \param xip_ pointer to vector of interpolation points
    Teuchos::RCP<std::vector<Real> > xip_;

    // \param xqp_ pointer to vector of quadrature points2
    Teuchos::RCP< std::vector<Real> > xqp_;

    // \param wqp_ pointer to vector of quadrature weights  
    Teuchos::RCP<std::vector<Real> > wqp_;

    // \param Lp_ pointer to vector containing Lagrange interpolant 
    Teuchos::RCP<std::vector<Real> > Lp_;
        
    // \param Dp_ pointer to vector containing derivative of Lagrange interpolant
    Teuchos::RCP<std::vector<Real> > Dp_; 
       
    NodalBasis(const int ni, const int nq);       

};


template<class Real>
NodalBasis<Real>::NodalBasis(const int ni, const int nq) : ni_(ni),nq_(nq), 
    xip_(Teuchos::rcp( new std::vector<Real>(ni,0) )),
    xqp_(Teuchos::rcp( new std::vector<Real>(nq,0) )),
    wqp_(Teuchos::rcp( new std::vector<Real>(nq,0) )),
    Lp_(Teuchos::rcp( new std::vector<Real>(ni*nq,0) )), 
    Dp_(Teuchos::rcp( new std::vector<Real>(ni*nq,0) )) {

    Teuchos::LAPACK<int,Real> lapack; 

    Teuchos::RCP<std::vector<Real> > wip = Teuchos::rcp( new std::vector<Real>(ni_,0) );
    StdVector<Real> xi(xip_);
    StdVector<Real> wi(wip);

    // ----[ Begin computing interpolation nodes ]---- //
    if(ni_<2) { // Require at least two points 
        TEUCHOS_TEST_FOR_EXCEPTION( (ni_<2), std::invalid_argument, ">>>  ERROR in FiniteElement Constructor: " 
                                                                    "Require number of interpolation points ni>=2");
    } 
    else if(ni==2) { // Linear Finite Elements
        (*xip_)[0] = -1.0;
        (*xip_)[1] = 1.0;
    }
    else {
  
        Teuchos::RCP<std::vector<Real> > aip = Teuchos::rcp( new std::vector<Real>(ni_,0) );
        Teuchos::RCP<std::vector<Real> > bip = Teuchos::rcp( new std::vector<Real>(ni_,0) );
 
        StdVector<Real> ai(aip);
        StdVector<Real> bi(bip);
     
        // Compute Legendre-Gauss recursion coefficients
        rec_jacobi(0,0,ai,bi);

        // Modify to Legendre-Gauss-Lobatto
        rec_lobatto(lapack,-1.0,1.0,ai,bi);

        // Compute interpolation nodes and weights (weights not used)
        gauss(lapack,ai,bi,xi,wi);
    }
    // ----[ End computing interpolation nodes ]---- //


    // ----[ Begin computing quadrature nodes and weights ]---- //
    Teuchos::RCP<std::vector<Real> > aqp = Teuchos::rcp( new std::vector<Real>(nq_,0) );
    Teuchos::RCP<std::vector<Real> > bqp = Teuchos::rcp( new std::vector<Real>(nq_,0) );

    StdVector<Real> aq(aqp);
    StdVector<Real> bq(bqp);

    StdVector<Real> xq(xqp_);
    StdVector<Real> wq(wqp_);

    // Compute Legendre-Gauss recursion coefficients
    rec_jacobi(0,0,aq,bq);

    gauss(lapack,aq,bq,xq,wq);
    // ----[ End computing quadrature nodes and weights ]---- //

    Lagrange<Real> lagrange(xi,xq); 
 
    // Pointer to canonical vector
    Teuchos::RCP<std::vector<Real> > ep = Teuchos::rcp(new std::vector<Real>(ni_,0) );

    // Allocate space to store an interpolant 
    Teuchos::RCP<std::vector<Real> > fp = Teuchos::rcp(new std::vector<Real>(nq_,0) );
    StdVector<Real> f(fp);

    // Evaluate Lagrange interpolants and their derivatives and store them 
    for(int i=0;i<ni_;++i) {
         
        lagrange.interpolant(i,f);
        std::copy(fp->begin(),fp->end(),Lp_->begin()+i*nq_);

        lagrange.derivative(i,f);
        std::copy(fp->begin(),fp->end(),Dp_->begin()+i*nq_);
    }

} // End FiniteElement Constructor



template<class Real,class Functor>
class VectorFunction {

    private:
        Real xl_;
        Real xr_;
        Real dx_;  
        bool deriv_; 
        Functor func_;
        Teuchos::RCP<NodalBasis<Real> > basisp_;
        int ni_;
        int nq_;  

    public:
        
        VectorFunction(Real xl, Real xr, bool deriv, Teuchos::RCP<NodalBasis<Real> > basisp);

        template<class ScalarT> 
        void evaluate(const Vector<ScalarT> &u, const Vector<ScalarT> &z, Vector<ScalarT> &f);

};


template<class Real,class Functor>
VectorFunction<Real,Functor>::VectorFunction(Real xl, Real xr, bool deriv, Teuchos::RCP<NodalBasis<Real> > basisp) : 
    xl_(xl), xr_(xr), dx_(xr-xl), deriv_(deriv), basisp_(basisp), ni_(basisp->ni_), nq_(basisp->nq_) {
}



template<class Real,class Functor>
template<class ScalarT>
void VectorFunction<Real,Functor>::evaluate(const Vector<ScalarT> &u, const Vector<ScalarT> &z, Vector<ScalarT> &f) {

    Teuchos::RCP<const std::vector<ScalarT> > up =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > zp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(z))).getVector();

    Teuchos::RCP<std::vector<ScalarT> > fp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (f)).getVector());

    // Need AD type of x for evaluating nonlinear function
    Teuchos::RCP<std::vector<ScalarT> > uqp   = Teuchos::rcp( new std::vector<ScalarT>(nq_,0) );
    Teuchos::RCP<std::vector<ScalarT> > u_xqp = Teuchos::rcp( new std::vector<ScalarT>(nq_,0) );
    Teuchos::RCP<std::vector<ScalarT> > zqp   = Teuchos::rcp( new std::vector<ScalarT>(nq_,0) );
    Teuchos::RCP<std::vector<ScalarT> > fqp   = Teuchos::rcp( new std::vector<ScalarT>(nq_,0) );

    // Interpolate z, u, u_x from onto the quadrature grid
    for(int i=0;i<ni_;++i) {
        for(int j=0;j<nq_;++j) {
          
            (*zqp)[j]   += (*basisp_->Lp_)[j+nq_*i]*(*zp)[i]; // y    
            (*u_xqp)[j] += (*basisp_->Dp_)[j+nq_*i]*(*up)[i]; // y'  
            (*uqp)[j]   += (*basisp_->Lp_)[j+nq_*i]*(*up)[i]; // u   
        }
    }
 
    // Evaluate a(u,u',z,x) on the quadrature grid
    for(int j=0;j<nq_;++j) {

        // Mapped x variable
        ScalarT xs = 0.5*(xl_*(1.0-(*basisp_->xqp_)[j])+xr_*(1.0+(*basisp_->xqp_)[j]));
        ScalarT u_xs = (*u_xqp)[j]/dx_;
        (*fqp)[j] = (func_)( (*uqp)[j], u_xs, (*zqp)[j], xs);      
    } 

    // Integrate against all first derivatives of interpolants
    if(deriv_) {
        for(int i=0;i<ni_;++i) {
            for(int j=0;j<nq_;++j) {
                (*fp)[i] += (*basisp_->wqp_)[j]*(*fqp)[j]*(*basisp_->Dp_)[j+nq_*i];   
            }
        } 
    }
    else { // Integrate against all interpolants 
        for(int i=0;i<ni_;++i) {
            for(int j=0;j<nq_;++j) {
                (*fp)[i] += (*basisp_->wqp_)[j]*(*fqp)[j]*(*basisp_->Lp_)[j+nq_*i]*dx_;   
            }
        } 
    }
}


#endif

