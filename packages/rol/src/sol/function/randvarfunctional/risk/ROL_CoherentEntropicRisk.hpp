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

#ifndef ROL_COHERENTEXPUTILITY_HPP
#define ROL_COHERENTEXPUTILITY_HPP

#include "ROL_RandVarFunctional.hpp"

/** @ingroup risk_group
    \class ROL::CoherentEntropicRisk
    \brief Provides the interface for the coherent entropic risk measure.

    The coherent entropic risk measure is
    \f[
       \mathcal{R}(X) = \inf_{\lambda > 0} \left\{
         \lambda \log\mathbb{E}\left[\exp\left(\frac{X}{\lambda}\right)\right]
         \right\}.
    \f]
    \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$\lambda\f$, then minimizes jointly for \f$(x_0,\lambda)\f$.
*/

namespace ROL {

template<class Real>
class CoherentEntropicRisk : public RandVarFunctional<Real> {
private:
  Real dval1_;
  Real dval2_;
  Real dval3_;

  using RandVarFunctional<Real>::val_;
  using RandVarFunctional<Real>::gv_;
  using RandVarFunctional<Real>::g_;
  using RandVarFunctional<Real>::hv_;
  using RandVarFunctional<Real>::dualVector_;

  using RandVarFunctional<Real>::weight_;

  using RandVarFunctional<Real>::computeValue;
  using RandVarFunctional<Real>::computeGradient;
  using RandVarFunctional<Real>::computeGradVec;
  using RandVarFunctional<Real>::computeHessVec;

public:
  CoherentEntropicRisk(void) : RandVarFunctional<Real>(),
    dval1_(0), dval2_(0), dval3_(0) {}

  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    Real zero(0);
    dval1_ = zero; dval2_ = zero; dval3_ = zero;
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    val_    += weight_ * std::exp(val/xstat[0]);
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real ev(0);
    sampler.sumAll(&val_,&ev,1);
    return xstat[0]*std::log(ev);
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real ev  = std::exp(val/xstat[0]);
    val_    += weight_ * ev;
    gv_     += weight_ * ev * val;
    computeGradient(*dualVector_,obj,x,tol);
    g_->axpy(weight_*ev,*dualVector_);
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    const Real one(1);
    // Perform sum over batches
    std::vector<Real> myval(2,0), val(2,0);
    myval[0] = val_;
    myval[1] = gv_;
    sampler.sumAll(&myval[0],&val[0],2);

    sampler.sumAll(*g_,g);
    g.scale(one/val[0]);
    gstat[0] = std::log(val[0]) - val[1]/(val[0]*xstat[0]);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real ev  = std::exp(val/xstat[0]);
    val_    += weight_ * ev;

    Real gv  = computeGradVec(*dualVector_,obj,v,x,tol);
    gv_     += weight_ * ev * gv;
    g_->axpy(weight_*ev,*dualVector_);
    hv_->axpy(weight_*ev*(gv-val*vstat[0]/xstat[0])/xstat[0],*dualVector_);

    dval1_  += weight_ * ev * val;
    dval2_  += weight_ * ev * val * val;
    dval3_  += weight_ * ev * val * gv;

    computeHessVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_*ev,*dualVector_);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    const Real one(1);
    std::vector<Real> myval(5,0), val(5,0);
    myval[0] = val_;
    myval[1] = gv_;
    myval[2] = dval1_;
    myval[3] = dval2_;
    myval[4] = dval3_;
    sampler.sumAll(&myval[0],&val[0],5);

    Real xs2 = xstat[0]*xstat[0];
    Real xs3 = xs2*xstat[0];
    Real v02 = val[0]*val[0];
    Real h11 = (val[3]*val[0] - val[2]*val[2])/(v02*xs3) * vstat[0];
    Real h12 = (val[1]*val[2] - val[4]*val[0])/(v02*xs2); 
    hvstat[0] = h11+h12;
    sampler.sumAll(*hv_,hv);
    hv.scale(one/val[0]);

    dualVector_->zero();
    sampler.sumAll(*g_,*dualVector_);
    hv.axpy((vstat[0]*val[2]/xs2-val[1]/xstat[0])/v02,*dualVector_);
  }
};

}

#endif
