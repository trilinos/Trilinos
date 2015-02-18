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

    ScalarT xdotx = 0;
    ScalarT kdotx = 0;
    ScalarT J = 0;
   
    // Compute dot products 
    for(unsigned i=0; i<x.size(); ++i) {
        xdotx += pow(x[i],2);       // (k,x)
        kdotx += double(i+1)*x[i];  // (x,x)
    }

    // Sum terms in objective function
    J = xdotx + pow(kdotx,2)/4.0 + pow(kdotx,4)/16.0;
    
    return J;
}


template<class Real>
class Zakharov_Sacado_Objective : public Objective<Real> {

  typedef Sacado::Fad::DFad<Real>          GradType;
  typedef Sacado::Fad::SFad<Real,1>        DirDerivType;
  typedef Sacado::Fad::DFad<DirDerivType>  HessVecType;
  // In C++11, we could use template typedefs:
  // template <typename T>       using  GradTypeT = Sacado::Fad::DFad<T>;
  // typedef Sacado::Fad::SFad<Real,1>  DirDerivType;
  // typedef GradTypeT<Real>            GradType;
  // typedef GradTypeT<DirDerivType>    HessVecType;

  private:

    FunctionZakharov<Real>        zfunc_; 
    FunctionZakharov<GradType>    zfuncGrad_; 
    FunctionZakharov<HessVecType> zfuncHessVec_; 

  public:

    Zakharov_Sacado_Objective() {}

    /* Evaluate the objective function at x */
    Real value( const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<const StdVector<Real> >(x)).getVector();
      return zfunc_.eval(*xp);
    }

    /* Evaluate the gradient at x */
    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<std::vector<Real> > gp =
        (Teuchos::dyn_cast<StdVector<Real> >(g)).getVector();
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<const StdVector<Real> >(x)).getVector();

      unsigned n = xp->size();

      std::vector<GradType> x_grad(n);

      for(unsigned i=0; i<n; ++i) {
        x_grad[i] = (*xp)[i]; // Set values x(i).
        x_grad[i].diff(i,n);  // Choose canonical directions.
      }

      GradType J_grad = zfuncGrad_.eval(x_grad);

      for(unsigned i=0; i<n; ++i) {
        (*gp)[i] = J_grad.dx(i);
      }

    }

    /* Compute the action of the Hessian evaluated at x on a vector v */
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<std::vector<Real> > hvp =
        (Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<const StdVector<Real> >(v)).getVector();
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<const StdVector<Real> >(x)).getVector();

      unsigned n = xp->size();

      std::vector<HessVecType>  x_hessvec(n);

      for(unsigned i=0; i<n; ++i) {
        DirDerivType tmp(1,(*xp)[i]);  // Set values x(i).
        tmp.fastAccessDx(0)= (*vp)[i]; // Set direction values v(i).
        x_hessvec[i] = tmp;            // Use tmp to define hessvec-able x.
        x_hessvec[i].diff(i,n);        // Choose directions.
      }

      // Compute Hessian-vector product (and other currently irrelevant things).
      HessVecType J_hessvec = zfuncHessVec_.eval(x_hessvec);

      for(unsigned i=0; i<n; ++i) {
        (*hvp)[i]   =  (J_hessvec.dx(i)).fastAccessDx(0);
        // hessvec  =  get gradient      get dir deriv
      }
    }
};

