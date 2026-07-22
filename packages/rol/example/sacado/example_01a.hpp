// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \brief An example equality constrained problem combining ROL and Sacado 
           This is the same problem as found in rol/examples/simple-eq-constr
           with the objective gradient, objective Hessian direction, constraint 
           Jacobian direction, constraint adjoint Jacobian direction, and
           constraint adjoint Hessian direction computed via automatic 
           differentiation with Sacado.  

    \author Created by G. von Winckel
**/



#include "ROL_StdVector.hpp"

using namespace ROL;

template<class Real>
class Zakharov {

    public:

        template<class ScalarT>
        ScalarT value(const Vector<ScalarT> &x, Real &tol );

};



/** \brief A Sacado-accessible version of the Zakharov function to differentiate
    \f[f(\mathbf{x}) = \mathbf{x}^\top\mathbf{x} 
                     + \frac{1}{4}(\mathbf{k}^\top \mathbf{x})^2 +
                       \frac{1}{16}(\mathbf{k}^\top \mathbf{x})^4 \f]
    Where \f$\mathbf{k}=(1,\cdots,n)\f$

    @param[in] x is the optimization vector 

    Returns the value of the objective function.  
*/ 
template<class Real>
template<class ScalarT>
ScalarT Zakharov<Real>::value(const Vector<ScalarT>& x, Real &tol) {

      

    ROL::Ptr<const std::vector<ScalarT> > xp = (dynamic_cast<const StdVector<ScalarT>&>(x)).getVector();

    int n = xp->size();

    ScalarT xdotx = 0;
    ScalarT kdotx = 0;
    ScalarT J = 0;
   
    // Compute dot products 
    for(int i=0; i<n; ++i) {
        xdotx += pow((*xp)[i],2);       // (k,x)
        kdotx += Real(i+1)*(*xp)[i];  // (x,x)
    }

    // Sum terms in objective function
    J = xdotx + pow(kdotx,2)/4.0 + pow(kdotx,4)/16.0;
    
    return J;
}


