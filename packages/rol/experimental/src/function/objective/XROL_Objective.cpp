
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

#include "XROL.hpp"

namespace XROL {

/** \class XROL::Objective::Impl
    \brief Implementation of Objective's default finite difference approximations 
           and derivative/symmetry checks. This interface can only be seen by
           Objective 
*/

template<class XPrim, class XDual> 
class ObjectiveImpl {  
private:

  friend class Objective<XPrim,XDual>;

  using Obj     = Objective<XPrim,XDual>;
  using Real    = magnitude_t<XPrim>;
  using Scalar  = element_t<XPrim>;


  Real eps_;
  Real sqrt_eps_;
 
  ObjectiveImpl() : eps_(ROL::ROL_EPSILON<Real>()), 
                    sqrt_eps_(std::sqrt(ROL::ROL_EPSILON<Real>())) { }

  virtual ~ObjectiveImpl(){};

  auto dirDeriv( Obj& obj, const XPrim& x, const XPrim& d, Real& tol ) {
    auto ftol = sqrt_eps_;
    auto xd = clone(d);
    set(xd,x);
    axpy(xd,tol,d);
    return ( obj.value(*xd) - obj.value(x) )/tol;
  }

  void gradient( Obj& obj, XDual& g, const XPrim& x, Real& tol ) {
    zero(g);
    Real deriv = 0.0;
    Real h     = 0.0;   
    Real xi    = 0.0;
    for( index_t<XDual> i=0; i<dimension(g); ++i ) {
      auto exi = *basis(x,i);
      auto egi = *basis(g,i);
      xi     = std::abs(dot(x,exi));
      h      = ( xi < eps_ ? 1.0 : xi )*tol;
      deriv  = obj.dirDeriv(x,exi,h);
      axpy(g,deriv,egi);
    }
  }
 
  void hessVec( Obj& obj, XDual &hv, const XPrim &v, const XPrim &x, Real& tol ) {
    Real _zero(0), _one(1);
    if( norm(v) == _zero ) {
      zero(hv);
    }
    else {
      auto gtol = sqrt_eps_;
      auto h = std::max(_one,norm(x)/norm(v))*tol;
      auto g = *clone(hv);
      obj.gradient(g,x,gtol);

      auto xnew = *clone(x);
      set(xnew,x);
      axpy(xnew,h,v);
      obj.update(xnew);
 
      zero(hv);
      obj.gradient(hv,xnew,gtol);

      axpy(hv,-_one,g);
      scale(_one/h);
    }
  }

  auto checkGradient( Obj& obj, const XPrim &x, const XDual& g, const XPrim& d, 
                      std::ostream &os, Teuchos::ParameterList& parlist ) {

    using ROL::Finite_Difference_Arrays::shifts;
    using ROL::Finite_Difference_Arrays::weights;

    using Real = element_t<XPrim>;
    using uint = typename std::vector<Real>::size_type;

    auto objlist = parlist.sublist("Objective");
    auto gclist  = objlist.sublist("Gradient Check"); 
    
    bool printToStream = gclist.get("Print to Stream", true);

    uint numVals = 4;

    uint order    = static_cast<uint>(gclist.get("Finite Difference Order",1));
    uint numSteps = static_cast<uint>(gclist.get("Number of Steps",10));

    std::vector<Real> steps(numSteps);
    
    if( gclist.isParameter("Finite Difference Steps") ) { 
      steps = Teuchos::getArrayFromStringParameter<Real>(gclist,
                "Finite Difference Steps",numSteps).toVector();
    } else {
      double stepFactor = gclist.get("Step Reduction Factor",10);
      for( uint i=0; i<numSteps; ++i)
        steps[i] = std::pow(stepFactor,-i);
    }

    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real>> gCheck(numSteps, tmp);

    // Evaluate objective value at x. 
    obj.update(x);
    Real val = obj.value(x,sqrt_eps_);

    // Compute gradient at x.
    auto xnew = *clone(x);
    auto gtmp = *clone(g);

    //    obj.update(x); // Why do we need to update again ?
    obj.gradient(gtmp, x, sqrt_eps_);

    // Evaluate the dual
    dual(xnew,gtmp);

    Real dtg = dot(d,xnew);
    
    // Temporary vectors.

    for (uint i=0; i<numSteps; i++) {

      Real eta = steps[i];

      set(xnew,x);

      // Compute gradient, finite-difference gradient, and absolute error.
      gCheck[i][0] = eta;
      gCheck[i][1] = dtg;

      gCheck[i][2] = weights[order-1][0] * val;

      for(uint j=0; j<order; ++j) {
        // Evaluate at x <- x+eta*c_i*d.
        axpy(xnew,eta*shifts[order-1][j], d);

        // Only evaluate at shifts where the weight is nonzero  
        if( weights[order-1][j+1] != 0 ) {
          obj.update(xnew);
          gCheck[i][2] += weights[order-1][j+1] * obj.value(xnew,sqrt_eps_);
        }
      }

      gCheck[i][2] /= eta;

      gCheck[i][3] = std::abs(gCheck[i][2] - gCheck[i][1]);

      if (printToStream) {
        if (i==0) {
          os << std::right
             << std::setw(20) << "Step size"
             << std::setw(20) << "grad'*dir"
             << std::setw(20) << "FD approx"
             << std::setw(20) << "abs error"
             << "\n"
             << std::setw(20) << "---------"
             << std::setw(20) << "---------"
             << std::setw(20) << "---------"
             << std::setw(20) << "---------"
             << "\n";
        }
        os << std::scientific << std::setprecision(11) << std::right
           << std::setw(20) << gCheck[i][0]
           << std::setw(20) << gCheck[i][1]
           << std::setw(20) << gCheck[i][2]
           << std::setw(20) << gCheck[i][3]
           << "\n";
        }
      }
    
      return gCheck;   
  } // checkGradient


  auto checkHessVec( Obj& obj, const XPrim& x, const XDual& hv, const XPrim& v, 
                     std::ostream &os, Teuchos::ParameterList& parlist ) {

    using ROL::Finite_Difference_Arrays::shifts;
    using ROL::Finite_Difference_Arrays::weights;

    using Real = element_t<XPrim>;
    using uint = typename std::vector<Real>::size_type;

    auto objlist = parlist.sublist("Objective");
    auto hvlist  = objlist.sublist("Hessian Check"); 
    
    bool printToStream = hvlist.get("Print to Stream", true);

    int order    = hvlist.get("Finite Difference Order",1);
    uint numVals = 4;

    uint numSteps = static_cast<uint>(hvlist.get("Number of Steps",10));

    std::vector<Real> steps(numSteps);
    
    if( hvlist.isParameter("Finite Difference Steps") ) { 
      steps = Teuchos::getArrayFromStringParameter<Real>(hvlist,
                "Finite Difference Steps",numSteps).toVector();
    } else {
      double stepFactor = hvlist.get("Step Reduction Factor",10);
      for( uint i=0; i<numSteps; ++i)
        steps[i] = std::pow(stepFactor,-i);
    }

    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real>> hvCheck(numSteps, tmp);
  
    auto g = *clone(hv);
    obj.update(x);
    obj.gradient(g,x,sqrt_eps_);
    
    auto Hv = *clone(hv);
    obj.hessVec(Hv,v,x,sqrt_eps_);
    Real normHv = norm(Hv);

    auto gdif = *clone(hv);
    auto gnew = *clone(hv);
    auto xnew = *clone(x);

    for (uint i=0; i<numSteps; i++) {

      Real eta = steps[i]; 
  
      // Evaluate objective value at x+eta*d.
      set(xnew,x);
  
      set(gdif,g);
      scale(gdif,weights[order-1][0]);
  
      for(uint j=0; j<order; ++j) {
  
          // Evaluate at x <- x+eta*c_i*d.
          xnew->axpy(eta*shifts[order-1][j], v);
  
          // Only evaluate at shifts where the weight is nonzero  
          if( weights[order-1][j+1] != 0 ) {
              obj.update(xnew);
              obj.gradient(gnew, xnew, sqrt_eps_); 
              axpy(gdif,weights[order-1][j+1],gnew);
          }
         
      }
  
      scale(gdif,1.0/eta);    
  
      // Compute norms of hessvec, finite-difference hessvec, and error.
      hvCheck[i][0] = eta;
      hvCheck[i][1] = normHv;
      hvCheck[i][2] = gdif->norm();
      axpy(gdif,-1.0, *Hv);
      hvCheck[i][3] = gdif->norm();
  
      if (printToStream) {
        if (i==0) {
        os << std::right
           << std::setw(20) << "Step size"
           << std::setw(20) << "norm(Hess*vec)"
           << std::setw(20) << "norm(FD approx)"
           << std::setw(20) << "norm(abs error)"
           << "\n"
           << std::setw(20) << "---------"
           << std::setw(20) << "--------------"
           << std::setw(20) << "---------------"
           << std::setw(20) << "---------------"
           << "\n";
        }
        os << std::scientific << std::setprecision(11) << std::right
           << std::setw(20) << hvCheck[i][0]
           << std::setw(20) << hvCheck[i][1]
           << std::setw(20) << hvCheck[i][2]
           << std::setw(20) << hvCheck[i][3]
           << "\n";

      } // end if(printToStream)
  
    } // end for(uint i=0; i<numSteps; i++) 

    return hvCheck;

  } // checkHessVec


  auto checkHessSym( Obj& obj, const XPrim& x, const XDual& hv, const XPrim& v, const XPrim &w,
                     std::ostream &os, Teuchos::ParameterList& parlist ) {

    auto h     = *clone(hv);
    auto hdual = *clone(x);
    obj.hessVec(h,v,x,sqrt_eps_);

    dual(hdual,h);

    Real wHv = dot(v,hdual);

    obj.hessVec(h,w,x,sqrt_eps_);

    Real vHw = dot(w,hdual);


    std::vector<Real> hsymCheck(3, 0);

    hsymCheck[0] = wHv;
    hsymCheck[1] = vHw;
    hsymCheck[2] = std::abs(vHw-wHv);
    
    auto objlist = parlist.sublist("Objective");
    auto hvlist  = objlist.sublist("Hessian Check"); 
    
    bool printToStream = hvlist.get("Print to Stream", true);


    if (printToStream) {
      os << std::right
         << std::setw(20) << "<w, H(x)v>"
         << std::setw(20) << "<v, H(x)w>"
         << std::setw(20) << "abs error"
         << "\n";
      os << std::scientific << std::setprecision(11) << std::right
         << std::setw(20) << hsymCheck[0]
         << std::setw(20) << hsymCheck[1]
         << std::setw(20) << hsymCheck[2]
         << "\n";
    }
    return hsymCheck;

  } // checkHessSym

 

}; // class Objective::Impl

} // namespace XROL
