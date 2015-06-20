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


