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

/** \file
    \brief An example combining ROL and Sacado to mimize the Zakharov 
           function. The gradient and the action of the Hessian on a given
           vector are computed by Sacado using automatic differentiation.    
                    
           This implementation is far from optimal as vectors of AD type
           are being created repeatedly. A more efficient implementation 
           would encapsulate the functions Zakharov, objgrad, and hessVec
           in an object so that the AD vectors can be instantiated once. 

    \author Created by G. von Winckel
**/

#include <iostream>
#include <iomanip>

#include "Sacado.hpp"

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_Algorithm.hpp"

#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


using namespace ROL;


/** \brief A Sacado-accessible version of the Zakharov function to differentiate
    \f[f(\mathbf{x}) = \mathbf{x}^\top\mathbf{x} 
                     + \frac{1}{4}(\mathbf{k}^\top \mathbf{x})^2 +
                       \frac{1}{16}(\mathbf{k}^\top \mathbf{x})^4 \f]
    Where \f$\mathbf{k}=(1,\cdots,n)\f$

    @param[in] x is the optimization vector 

    Returns the value of the objective function.  
*/ 
template <typename ScalarT>
ScalarT Zakharov(const Vector<ScalarT>& x) {
     Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    int n = xp->size();

    ScalarT xdotx = 0;
    ScalarT kdotx = 0;
    ScalarT J = 0;
   
    // Compute dot products 
    for(int i=0; i<n; ++i) {
        xdotx += pow((*xp)[i],2);       // (k,x)
        kdotx += double(i+1)*(*xp)[i];  // (x,x)
    }

    // Sum terms in objective function
    J = xdotx + pow(kdotx,2)/4.0 + pow(kdotx,4)/16.0;
    
    return J;
}



/** \brief A Sacado-accesible function which computes the gradient of the Zakharov function 
    using automatic differentiation.
    @param[in]  x is the optimization vector
    @param[out] J is the value of the objective function
    @param[out] g is the gradient of the objective function    
*/
template <typename ScalarT>
void objgrad(const Vector<ScalarT>& x, ScalarT& J, Vector<ScalarT>& g) {

    // Data type which supports automatic differentiation 
    typedef Sacado::Fad::DFad<ScalarT> FadType;

    // Get a pointer to the optimization vector
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    // Get a pointer to the gradient vector
    Teuchos::RCP<std::vector<ScalarT> > gp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (g)).getVector());
    
    int n = xp->size();
 
    // Create a vector of independent variables
    Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    x_fad_rcp->reserve(n);
   
    // Initialize constructor for each element
    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(n,i,(*xp)[i])); 
    }

    StdVector<FadType> x_fad(x_fad_rcp);

    // AD access to objective function
    FadType J_fad = Zakharov(x_fad);

    // Evaluate objective function
    J = J_fad.val();

    // Evaluate gradient
    for(int i=0; i<n; ++i) {
        (*gp)[i] = J_fad.dx(i);
    }
}

/** \brief A Sacado-accesible function which computes the action of the Hessian on a given direction vector
    @param[in]   x is the optimization vector
    @param[in]   v is a given direction vector 
    @param[out]  J is the value of the objective function
    @param[out]  g is the gradient of the objective function
    @param[out]  hv is the action of the Hessian the direction v 
*/
template <typename ScalarT>
void applyHessian(const Vector<ScalarT>& x, const Vector<ScalarT>& v,
             ScalarT& J, Vector<ScalarT>& g, Vector<ScalarT>& hv) {
  
    // Get a pointer to the optimization vector 
    Teuchos::RCP<const std::vector<ScalarT> > xp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(x))).getVector();

    // Get a pointer to the gradient vector
    Teuchos::RCP<std::vector<ScalarT> > gp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> > ((Teuchos::dyn_cast<StdVector<ScalarT> > (g)).getVector());

    // Get a pointer to the direction vector
    Teuchos::RCP<const std::vector<ScalarT> > vp =
        (Teuchos::dyn_cast<StdVector<ScalarT> >(const_cast<Vector<ScalarT> &>(v))).getVector();

    Teuchos::RCP<std::vector<ScalarT> > hvp =
        Teuchos::rcp_const_cast<std::vector<ScalarT> >((Teuchos::dyn_cast<StdVector<ScalarT> >(hv)).getVector());

    typedef Sacado::Fad::SFad<ScalarT,1> FadType;

    int n = xp->size();
   
    FadType J_fad; 

    // Create a vector of independent variables
    Teuchos::RCP<std::vector<FadType> > x_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    x_fad_rcp->reserve(n);
 
    // Allocate for gradient   
    Teuchos::RCP<std::vector<FadType> > g_fad_rcp = Teuchos::rcp( new std::vector<FadType> );
    g_fad_rcp->reserve(n);

    for(int i=0; i<n; ++i) {
        x_fad_rcp->push_back(FadType(1,(*xp)[i]));
    }

    // Set directional derivative    
    for(int i=0; i<n; ++i) {
        (*x_fad_rcp)[i].fastAccessDx(0) = (*vp)[i];
    }
    
    StdVector<FadType> x_fad(x_fad_rcp);
    StdVector<FadType> g_fad(g_fad_rcp);

    objgrad(x_fad,J_fad,g_fad);

    J = J_fad.val();

    for(int i=0; i<n; ++i) {
        (*gp)[i] = (*g_fad_rcp)[i].val();
        (*hvp)[i] = (*g_fad_rcp)[i].dx(0);            
    }
}




/** \brief Objective Class */
template<class Real>
class Objective_Sacado : public Objective<Real> {

    public:
    Objective_Sacado(void) {}
 
    /* Evaluate the objective function at x */
    Real value( const Vector<Real> &x, Real &tol ) {
      StdVector<Real> & ex =
        Teuchos::dyn_cast<StdVector<Real> >(const_cast <Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      Real J = Zakharov(x);

      return J;
    }

    /* Evaluate the gradient at x */
    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
        Teuchos::RCP<const std::vector<Real> > xp =
            (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
        Teuchos::RCP<std::vector<Real> > gp =
            Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());

        Real J;
        objgrad(x,J,g);

    }

    /* Compute the action of the Hessian evaluated at x on a vector v */
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
        Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

        int n = xp->size();

        Real J;

        Teuchos::RCP<std::vector<Real> > gp = Teuchos::rcp( new std::vector<Real>(n,0.0) );
        StdVector<Real> g(gp);

        applyHessian(x,v,J,g,hv);
    }
};

/** Objective Class Definition end **/

