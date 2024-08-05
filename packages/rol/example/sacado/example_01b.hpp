// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \brief Minimize Zakharov's function by combining ROL and Sacado AD.
           In this example the implementation of the ROL objective function
           uses three instances, Scalar, Fad and Fad Fad, of the
           FunctionZakharov class. 

    \author Created by Denis Ridzal.
**/



#include "ROL_StdVector.hpp"
#include "Sacado.hpp"

using namespace ROL;

template<class ScalarT>
class FunctionZakharov {

  

  public:

    ScalarT eval(const std::vector<ScalarT> &x);

};

/** \brief A Sacado-accessible version of the Zakharov function to differentiate
    \f[f(\mathbf{x}) = \mathbf{x}^\top\mathbf{x} 
                     + \frac{1}{4}(\mathbf{k}^\top \mathbf{x})^2 +
                       \frac{1}{16}(\mathbf{k}^\top \mathbf{x})^4 \f]
    Where \f$\mathbf{k}=(1,\cdots,n)\f$

    @param[in] x is the optimization vector 

    Returns the value of the objective function.  
*/ 
template<class ScalarT>
ScalarT FunctionZakharov<ScalarT>::eval(const std::vector<ScalarT> & x) {

    typedef typename std::vector<ScalarT>::size_type uint;

    ScalarT xdotx = 0;
    ScalarT kdotx = 0;
    ScalarT J = 0;
   
    // Compute dot products 
    for(uint i=0; i<x.size(); ++i) {
        xdotx += pow(x[i],2);       // (k,x)
        kdotx += ScalarT(i+1)*x[i];  // (x,x)
    }

    // Sum terms in objective function
    J = xdotx + pow(kdotx,2)/4.0 + pow(kdotx,4)/16.0;
    
    return J;
}


template<class Real>
class Zakharov_Sacado_Objective : public Objective<Real> {

  typedef std::vector<Real>  vector;
  typedef Vector<Real>       V;
  typedef StdVector<Real>    SV;

  typedef Sacado::Fad::DFad<Real>          GradType;
  typedef Sacado::Fad::SFad<Real,1>        DirDerivType;
  typedef Sacado::Fad::DFad<DirDerivType>  HessVecType;

  typedef typename vector::size_type uint;

  // In C++11, we could use template typedefs:
  // template <typename T>       using  GradTypeT = Sacado::Fad::DFad<T>;
  // typedef Sacado::Fad::SFad<Real,1>  DirDerivType;
  // typedef GradTypeT<Real>            GradType;
  // typedef GradTypeT<DirDerivType>    HessVecType;

  private:

    FunctionZakharov<Real>        zfunc_; 
    FunctionZakharov<GradType>    zfuncGrad_; 
    FunctionZakharov<HessVecType> zfuncHessVec_; 

    ROL::Ptr<const vector> getVector( const V& x ) {
      
      return dynamic_cast<const SV&>(x).getVector();
    }

    ROL::Ptr<vector> getVector( V& x ) {
      
      return dynamic_cast<SV&>(x).getVector();
    }

  public:

    Zakharov_Sacado_Objective() {}

    /* Evaluate the objective function at x */
    Real value( const Vector<Real> &x, Real &tol ) {
      
      ROL::Ptr<const vector> xp = getVector(x);      
      return zfunc_.eval(*xp);
    }

    /* Evaluate the gradient at x */
    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      
      ROL::Ptr<const vector> xp = getVector(x);
      ROL::Ptr<vector> gp = getVector(g);

      uint n = xp->size();

      std::vector<GradType> x_grad(n);

      for(uint i=0; i<n; ++i) {
        x_grad[i] = (*xp)[i]; // Set values x(i).
        x_grad[i].diff(i,n);  // Choose canonical directions.
      }

      GradType J_grad = zfuncGrad_.eval(x_grad);

      for(uint i=0; i<n; ++i) {
        (*gp)[i] = J_grad.dx(i);
      }

    }

    /* Compute the action of the Hessian evaluated at x on a vector v */
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {

      
      ROL::Ptr<vector> hvp = getVector(hv);
      ROL::Ptr<const vector> vp = getVector(v);
      ROL::Ptr<const vector> xp = getVector(x);

      uint n = xp->size();

      std::vector<HessVecType>  x_hessvec(n);

      for(uint i=0; i<n; ++i) {
        DirDerivType tmp(1,(*xp)[i]);  // Set values x(i).
        tmp.fastAccessDx(0)= (*vp)[i]; // Set direction values v(i).
        x_hessvec[i] = tmp;            // Use tmp to define hessvec-able x.
        x_hessvec[i].diff(i,n);        // Choose directions.
      }

      // Compute Hessian-vector product (and other currently irrelevant things).
      HessVecType J_hessvec = zfuncHessVec_.eval(x_hessvec);

      for(uint i=0; i<n; ++i) {
        (*hvp)[i]   =  (J_hessvec.dx(i)).fastAccessDx(0);
        // hessvec  =  get gradient      get dir deriv
      }
    }
};

