// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include "ROL_TestProblem.hpp"
#include "ROL_StdVector.hpp"


namespace ROL {
namespace ZOO {

/** \brief Zakharov function.
 */
template<class Real>
class Objective_Zakharov : public Objective<Real> {
private:
    ROL::Ptr<Vector<Real> > k_;  

public:
  
  // Create using a ROL::Vector containing 1,2,3,...,n
  Objective_Zakharov(const ROL::Ptr<Vector<Real> > k) : k_(k) {}

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

  Real dirDeriv( const Vector<Real> &x, const Vector<Real> &d, Real &tol ) {

      Real kdotd = d.dot(*k_);
      Real kdotx = x.dot(*k_);
      Real xdotd = x.dot(d);
      
      Real coeff = 0.25*(2.0*kdotx+pow(kdotx,3.0));

      Real deriv = 2*xdotd + coeff*kdotd;

      return deriv;

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
class getZakharov : public TestProblem<Real> {
public:
  getZakharov(void) {}

  Ptr<Objective<Real>> getObjective(void) const {
    // Problem dimension
    int n = 10;
    // Instantiate Objective Function
    ROL::Ptr<std::vector<Real> > k_ptr = ROL::makePtr<std::vector<Real>>(n,0.0);
    for ( int i = 0; i < n; i++ ) {
      (*k_ptr)[i] = i+1.0;
    }
    ROL::Ptr<Vector<Real> > k = ROL::makePtr<StdVector<Real>>(k_ptr);
    return ROL::makePtr<Objective_Zakharov<Real>>(k);
  }

  Ptr<Vector<Real>> getInitialGuess(void) const {
    // Problem dimension
    int n = 10;
    // Get Initial Guess
    ROL::Ptr<std::vector<Real> > x0p = ROL::makePtr<std::vector<Real>>(n,3.0);
    return ROL::makePtr<StdVector<Real>>(x0p);
  }

  Ptr<Vector<Real>> getSolution(const int i = 0) const {
    // Problem dimension
    int n = 10;
    // Get Solution
    ROL::Ptr<std::vector<Real> > xp = ROL::makePtr<std::vector<Real>>(n,0.0);
    return ROL::makePtr<StdVector<Real>>(xp);
  }
};


}// End ZOO Namespace
}// End ROL Namespace

#endif
