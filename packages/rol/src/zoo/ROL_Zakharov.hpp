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
    \brief  Contains definitions for the Zakharov function as evaluated using only the 
            ROL::Vector interface.
    \details This is a nice example not only because the gradient, Hessian, and inverse Hessian
             are easy to derive, but because they only require dot products, meaning this
             code can be used with any class that inherits from ROL::Vector.

    Objective function: 
    \f[f(\mathbf{x}) = \mathbf{x}^\top\mathbf{x} + \frac{1}{4}(\mathbf{k}^\top \mathbf{x})^2 +
                                                   \frac{1}{16}(\mathbf{k}^\top \mathbf{x})^4 \f]
    Where \f$\mathbf{k}=(1,\cdots,n)\f$
 
    Gradient:
    \f[
    g=\nabla f(\mathbf{x}) = 2\mathbf{x} + 
                           \frac{1}{4}\left(2(\mathbf{k}^\top\mathbf{x})+(\mathbf{k}^\top\mathbf{x})^3\right)\mathbf{k} 
    \f]

    Hessian: 
    \f[
    H=\nabla^2 f(\mathbf{x}) = 2 I + \frac{1}{4}[2+3(\mathbf{k}^\top\mathbf{x})^2]\mathbf{kk}^\top
    \f]
 
    The Hessian is a multiple of the identity plus a rank one symmetric 
    matrix, therefore the action of the inverse Hessian can be 
    performed using the Sherman-Morrison formula.

    \f[
    H^{-1}\mathbf{v} = \frac{1}{2}\mathbf{v}-\frac{(\mathbf{k}^\top\mathbf{v})}
                                             {\frac{16}{2+3(\mathbf{k}^\top\mathbf{x})^2}+2\mathbf{k^\top}\mathbf{k}}\mathbf{k}
    \f]

    \author Created by G. von Winckel
**/

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_ZAKHAROV_HPP
#define ROL_ZAKHAROV_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

namespace ROL {
namespace ZOO {

  /** \brief Zakharov function.
   */
  template<class Real>
  class Objective_Zakharov : public Objective<Real> {
  private:
      Teuchos::RCP<Vector<Real> > k_;  

  public:
    
    // Create using a ROL::Vector containing 1,2,3,...,n
    Objective_Zakharov(const Teuchos::RCP<Vector<Real> > k) : k_(k) {}

    Real value( const Vector<Real> &x, Real &tol ) {

        Real xdotx = x.dot(x); 
        Real kdotx = x.dot(*k_); 

        Real val = xdotx + pow(kdotx,2)/4.0 + pow(kdotx,4)/16.0;

        return val;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {

        Real kdotx = x.dot(*k_);
        Real coeff = 0.25*(2.0*kdotx+pow(kdotx,3.0));

        g.set(x);
        g.scale(2.0);
        g.axpy(coeff,*k_);
    }

#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

        Real kdotx = x.dot(*k_);
        Real kdotv = v.dot(*k_);
        Real coeff = 0.25*(2.0+3.0*pow(kdotx,2.0))*kdotv;

        hv.set(v);
        hv.scale(2.0);
        hv.axpy(coeff,*k_);
    }
#endif
    void invHessVec( Vector<Real> &ihv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    
        Real kdotv = v.dot(*k_);
        Real kdotx = x.dot(*k_);
        Real kdotk = (*k_).dot(*k_);
        Real coeff = -kdotv/(2.0*kdotk+16.0/(2.0+3.0*pow(kdotx,2.0)));
        
        ihv.set(v);
        ihv.scale(0.5);
        ihv.axpy(coeff,*k_); 
    }
};



template<class Real>
void getZakharov( Teuchos::RCP<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x ) {
    // Cast Initial Guess and Solution Vectors
    Teuchos::RCP<std::vector<Real> > x0p =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x0)).getVector());
    Teuchos::RCP<std::vector<Real> > xp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x)).getVector());
    int n = xp->size();
    // Resize Vectors
    n = 10;

    Teuchos::RCP<std::vector<Real> > k_rcp = Teuchos::rcp(new std::vector<Real>(n,0));
    for(int i=0;i<n;++i) {
        (*k_rcp)[i] = i+1.0;    
    }    
    Teuchos::RCP<Vector<Real> > k = Teuchos::rcp(new StdVector<Real>(k_rcp));
    x0p->resize(n);
    xp->resize(n);
    // Instantiate Objective Function
    obj = Teuchos::rcp( new Objective_Zakharov<Real>(k) );
    // Get Initial Guess
    (*x0p)[0] =  3.0;
    (*x0p)[1] =  3.0;
    // Get Solution
    (*xp)[0] = 0.0;
    (*xp)[1] = 0.0;
  }


}// End ZOO Namespace
}// End ROL Namespace

#endif
