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




/** \brief Evaluate the following generic variable coefficient type term 
    \f[ f_j = \left(\varphi_j^{(i)}, a(y,y',u,x)\right),\quad i=0,1 \f]
    @param[in] y the simulation variable
    @param[in] u the optimization variable
    @param[in] deriv is a bool for whether to differentiate the test functions
    @param[in] xlr is a two element vector containing the locations of the left and right ends of an element
    @param[in] funcp the variable coefficient function which depends on \f$y,y',u,x\f$
    @param[in] basis is a structure containing terms needed for computing numerical inner products
    @param[out] f vector containing these inner products */
template<class Functor, class ScalarT, class Real>
void vectorFunction(const Vector<ScalarT> &y,
                    const Vector<ScalarT> &u,
                    const bool deriv, 
                    const Real xlr[2], 
                    Functor* funcp,
                    Teuchos::RCP<NodalBasis<Real> > basisp,  
                    Vector<ScalarT> &f){

    int ni = basisp->ni_;
    int nq = basisp->nq_;     

    Real dx = 0.5*(xlr[1]-xlr[0]); // relative length of element

    Teuchos::RCP<const std::vector<ScalarT> > yp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(y))).getVector();

    Teuchos::RCP<const std::vector<ScalarT> > up =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

    Teuchos::RCP<std::vector<ScalarT> > fp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (f)).getVector());

    // Need AD type of x for evaluating nonlinear function
    Teuchos::RCP<std::vector<ScalarT> > yqp   = Teuchos::rcp( new std::vector<ScalarT>(nq,0) );
    Teuchos::RCP<std::vector<ScalarT> > y_xqp = Teuchos::rcp( new std::vector<ScalarT>(nq,0) );
    Teuchos::RCP<std::vector<ScalarT> > uqp   = Teuchos::rcp( new std::vector<ScalarT>(nq,0) );
    Teuchos::RCP<std::vector<ScalarT> > fqp   = Teuchos::rcp( new std::vector<ScalarT>(nq,0) );

    // Interpolate y, u, y_x from onto the quadrature grid
    for(int i=0;i<ni;++i) {
        for(int j=0;j<nq;++j) {
          
            (*yqp)[j]   += (*basisp->Lp_)[j+nq*i]*(*yp)[i]; // y    
            (*y_xqp)[j] += (*basisp->Dp_)[j+nq*i]*(*yp)[i]; // y'  
            (*uqp)[j]   += (*basisp->Lp_)[j+nq*i]*(*up)[i]; // u   
        }
    }
 
    // Evaluate a(y,y',u,x) on the quadrature grid
    for(int j=0;j<nq;++j) {

        // Mapped x variable
        ScalarT xs = 0.5*(xlr[0]*(1-(*basisp->xqp_)[j])+xlr[1]*(1+(*basisp->xqp_)[j]));
        ScalarT y_xs = (*y_xqp)[j]/dx;
        (*fqp)[j] = (*funcp)( (*yqp)[j], y_xs, (*uqp)[j], xs);      
    } 

    // Integrate against all first derivatives of interpolants
    if(deriv) {
        for(int i=0;i<ni;++i) {
            for(int j=0;j<nq;++j) {
                (*fp)[i] += (*basisp->wqp_)[j]*(*fqp)[j]*(*basisp->Dp_)[j+nq*i];   
            }
        } 
    }
    else { // Integrate against all interpolants 
        for(int i=0;i<ni;++i) {
            for(int j=0;j<nq;++j) {
                (*fp)[i] += (*basisp->wqp_)[j]*(*fqp)[j]*(*basisp->Lp_)[j+nq*i]*dx;   
            }
        } 
    }

} // End vectorFunction




/** \brief Evaluate the action of the Jacobian block of a generic variable coefficient type term 
    \f[ f_j = \left(\varphi_j^{(i)}, a(y,y',u,x)\right),\quad i=0,1 \f]
    on a direction vector v.  
    @param[in] y the simulation variable
    @param[in] u the optimization variable
    @param[in] deriv is a bool for whether to differentiate the test functions
    @param[in] xlr is a two element vector containing the locations of the left and right ends of an element
    @param[in] funcp the variable coefficient function which depends on \f$y,y',u,x\f$
    @param[in] basisp is a pointer to a finite element object
    @param[in] v is a direction vector
    @param[in] blk determines whether to differentiate with respect to sim variable (0) or opt (1)
    @param[out] jvp is the jacobian block in the v direction */
template<class Functor,class ScalarT, class Real>
void applyJacobianBlock(const Vector<ScalarT> &y, 
                        const Vector<ScalarT> &u,
                        const bool deriv,
                        const Real xlr[2],
                        Functor *funcp,  
                        Teuchos::RCP<NodalBasis<Real> > basisp,
                        const Vector<ScalarT> &v,
                        const int blk,
                        Vector<ScalarT> &jv) {


   typedef Sacado::Fad::SFad<ScalarT,1> FadType;

   Teuchos::RCP<const std::vector<ScalarT> > yp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(y))).getVector();

   Teuchos::RCP<const std::vector<ScalarT> > up =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();

   Teuchos::RCP<const std::vector<ScalarT> > vp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();

    Teuchos::RCP<std::vector<ScalarT> > jvp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (jv)).getVector());

    int ni = yp->size();

    Teuchos::RCP<std::vector<FadType> > y_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    y_fad_rcp->reserve(ni);

    Teuchos::RCP<std::vector<FadType> > u_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    u_fad_rcp->reserve(ni);

    Teuchos::RCP<std::vector<FadType> > v_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    v_fad_rcp->reserve(ni);

    Teuchos::RCP<std::vector<FadType> > f_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    f_fad_rcp->reserve(ni);

    if(blk==0) { // Differentiate with respect to simulation variable
        for(int i=0; i<ni; ++i) {

             // Call constructors on each element 
             y_fad_rcp->push_back(FadType(1,0)); 
             u_fad_rcp->push_back(0);
             v_fad_rcp->push_back(0);
             f_fad_rcp->push_back(0);

             // Set values
             (*y_fad_rcp)[i].val() = (*yp)[i];
             (*u_fad_rcp)[i].val() = (*up)[i];
             (*v_fad_rcp)[i].val() = (*vp)[i];

             // Set directional derivative
             (*y_fad_rcp)[i].fastAccessDx(0) = (*v_fad_rcp)[i].val();
        }
    }
    else if(blk==1) { // Differentiate with respect to the optimization variable   
        for(int i=0; i<ni; ++i) {

             // Call constructors on each element 
             y_fad_rcp->push_back(0);
             u_fad_rcp->push_back(FadType(1,0)); 
             v_fad_rcp->push_back(0);
             f_fad_rcp->push_back(0);

             // Set values 
             (*y_fad_rcp)[i].val() = (*yp)[i]; 
             (*u_fad_rcp)[i].val() = (*up)[i];
             (*v_fad_rcp)[i].val() = (*vp)[i];
                 
             // Set direction derivative
             (*u_fad_rcp)[i].fastAccessDx(0) = (*v_fad_rcp)[i].val(); 
       }
    }
    else { // Undefined block index
        TEUCHOS_TEST_FOR_EXCEPTION( (blk>1), std::invalid_argument, 
            ">>>  ERROR in FiniteElement::applyJacobianBlock : " 
            "block index must be 0 or 1");
    }

    StdVector<FadType> y_fad(y_fad_rcp);
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> f_fad(f_fad_rcp);

    // Evaluate the vector
    vectorFunction(y_fad,u_fad,deriv,xlr,funcp,basisp,f_fad);

    // Compute its directional derivative
    for(int i=0; i<ni; ++i) {
        (*jvp)[i] = (*f_fad_rcp)[i].dx(0);
    }
} // End applyJacobianBlock



/** \brief Construct the Jacobian block of a generic variable coefficient type term 
    \f[ f_j = \left(\varphi_j^{(i)}, a(y,y',u,x)\right),\quad i=0,1 \f]
    @param[in] y the simulation variable
    @param[in] u the optimization variable
    @param[in] deriv is a bool for whether to differentiate the test functions
    @param[in] xlr is a two element vector containing the locations of the left and right ends of an element
    @param[in] funcp the variable coefficient function which depends on \f$y,y',u,x\f$
    @param[in] basisp is a pointer to a finite element object
    @param[in] blk determines whether to differentiate with respect to sim variable (0) or opt (1)
    @param[out] jac is the Jacobian block (column-stacked) */
template<class Functor, class ScalarT, class Real>
void buildJacobianBlock(const Vector<ScalarT> &y, 
                        const Vector<ScalarT> &u,
                        const bool deriv,
                        const Real xlr[2],
                        Functor *funcp,  
                        Teuchos::RCP<NodalBasis<Real> > basisp,
                        const int blk,
                        Vector<ScalarT> &jac){

   Teuchos::RCP<std::vector<ScalarT> > jacp = 
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (jac)).getVector());

   int ni  = basisp->ni_; 
   int ni2 = jacp->size();

   if(ni2!=ni*ni) { // Require at least two points 
        TEUCHOS_TEST_FOR_EXCEPTION( (ni2!=ni*ni), 
                                     std::invalid_argument, ">>>  ERROR in FiniteElement::buildloadVectorJacobian : " 
                                                            "loadVectorJacobian must have ni^2 elements");
   } 

   // Canonical vector
   Teuchos::RCP<std::vector<ScalarT> > e_fad_rcp = Teuchos::rcp(new std::vector<ScalarT>(ni,0)); 

   // Storage vector
   Teuchos::RCP<std::vector<ScalarT> > je_fad_rcp = Teuchos::rcp(new std::vector<ScalarT>(ni,0));

   StdVector<ScalarT> e_fad(e_fad_rcp);
   StdVector<ScalarT> je_fad(je_fad_rcp);

   for(int i=0;i<ni;++i) {
       if(i>0) {
           (*e_fad_rcp)[i-1] = 0.0;
       }

       (*e_fad_rcp)[i] = 1.0;

       applyJacobianBlock(y,u,deriv,xlr,funcp,basisp,e_fad,blk,je_fad);

       std::copy(je_fad_rcp->begin(),je_fad_rcp->end(),jacp->begin()+i*ni);
   }
} // End build JacobianBlock




/** \brief Apply a block of the constraint Hessian to vectors v1, v2 as in taking the second Gateaux derivative
    of the equality constraint. 
    \details \f[ \frac{\partial^2}{\partial s\partial t} c(\mathbf{x} + s \mathbf{v}_1 + t\mathbf{v}_2)\big|_{s,t=0} \f] 
    where \f$ \mathbf{x}=(\mathbf{y},\mathbf{u}) \f$ and we are either differentiating only in 
    optimization or simulation vector each time. That is to say that v1 and v2 are zero in the optimization
    block and nonzero in the simulation block or vice versa.  
    @param[in] y the simulation variable
    @param[in] u the optimization variable
    @param[in] deriv is a bool for whether to differentiate the test functions
    @param[in] xlr is a two element vector containing the locations of the left and right ends of an element
    @param[in] funcp the variable coefficient function which depends on \f$y,y',u,x\f$
    @param[in] basisp is a pointer to a finite element object
    @param[in] v is a direction vector
    @param[in] blk determines whether to differentiate with respect to sim variable (0) or opt (1)
*/

template<class Functor, class ScalarT, class Real>
void applyHessianBlock(const Vector<ScalarT> &y,
                       const Vector<ScalarT> &u,
                       const bool deriv,
                       const Real xlr[2],
                       Functor *funcp,
                       Teuchos::RCP<NodalBasis<Real> > basisp,
                       const Vector<ScalarT> &v1,
                       const Vector<ScalarT> &v2,
                       const int blk[2],
                       Vector<ScalarT> &hv) {

    typedef Sacado::Fad::SFad<Real,1> FadType;

    Teuchos::RCP<const std::vector<ScalarT> > yp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(y))).getVector();
    Teuchos::RCP<const std::vector<ScalarT> > up =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(u))).getVector();
    Teuchos::RCP<const std::vector<ScalarT> > v1p =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v1))).getVector();
    Teuchos::RCP<const std::vector<ScalarT> > v2p =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v2))).getVector();

    int ni = yp->size();

    Teuchos::RCP<std::vector<ScalarT> > hvp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (hv)).getVector());


    Teuchos::RCP<std::vector<FadType> > y_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    y_fad_rcp->reserve(ni);
    Teuchos::RCP<std::vector<FadType> > u_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    u_fad_rcp->reserve(ni);
    Teuchos::RCP<std::vector<FadType> > v1_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    v1_fad_rcp->reserve(ni);
    Teuchos::RCP<std::vector<FadType> > jv1_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    jv1_fad_rcp->reserve(ni);
    
    if(blk[1]==0) { // Second deriv wrt simulation variable
        for(int i=0; i<ni; ++i) {
            y_fad_rcp->push_back(FadType(1,(*yp)[i].val()));
            u_fad_rcp->push_back((*up)[i].val());
            v1_fad_rcp->push_back((*v1p)[i].val());
            jv1_fad_rcp->push_back(0);
        }
        for(int i=0; i<ni; ++i) {
            (*y_fad_rcp)[i].fastAccessDx(0) = (*v2p)[i].val();
        }
    }
    else if(blk[1]==1) { // Second deriv wrt optimization variable
        for(int i=0; i<ni; ++i) {
            y_fad_rcp->push_back((*yp)[i].val());
            u_fad_rcp->push_back(FadType(1,(*up)[i].val()));
            v1_fad_rcp->push_back((*v1p)[i].val());
            jv1_fad_rcp->push_back(0);
        }
        for(int i=0; i<ni; ++i) {
            (*u_fad_rcp)[i].fastAccessDx(0) = (*v2p)[i].val();
        }
    }
    else {
        TEUCHOS_TEST_FOR_EXCEPTION( (blk[1]>1), std::invalid_argument, 
            ">>>  ERROR in FiniteElement::applyHessianBlock : " 
            "block index must be 0 or 1");

    }

    StdVector<FadType> y_fad(y_fad_rcp);
    StdVector<FadType> u_fad(u_fad_rcp);
    StdVector<FadType> v1_fad(v1_fad_rcp);
    StdVector<FadType> jv1_fad(jv1_fad_rcp);

    applyJacobianBlock(y_fad,u_fad,deriv,xlr,funcp,basisp,v1_fad,blk[0],jv1_fad);
   
    for(int i=0; i<ni; ++i) {
        (*hvp)[i] = (*jv1_fad_rcp)[i].dx(0);
    }

} // End applyHessianBlock




#endif

