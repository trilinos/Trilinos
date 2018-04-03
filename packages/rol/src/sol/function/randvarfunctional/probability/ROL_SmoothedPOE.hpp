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

#ifndef ROL_SMOOTHEDPOE_HPP
#define ROL_SMOOTHEDPOE_HPP

#include "ROL_RandVarFunctional.hpp"

/** @ingroup stochastic_group 
    \class ROL::SmoothedPOE
    \brief Provides the implementation of the smoothed probability of exceedance.

    Let \f$(\Omega,\mathcal{F},\mathbb{P})\f$ be a complete space.
    Here, \f$\Omega\f$ is the set of outcomes,
    \f$\mathcal{F}\subseteq 2^\Omega\f$ is a \f$\sigma\f$-algebra of events and
    \f$\mathbb{P}:\mathcal{F}\to[0,1]\f$ is a probability measure.  Moreover,
    let \f$\mathcal{X}\f$ be a class of random variables.

    ROL's SmoothedPOE class inherits from ROL::RandVarFunctional which is
    written in a way to exploit parallel sampling.
*/

namespace ROL {

template<class Real>
class SmoothedPOE : public RandVarFunctional<Real> {
private:
  Real threshold_;
  Real eps_;

  using RandVarFunctional<Real>::val_;
  using RandVarFunctional<Real>::gv_;
  using RandVarFunctional<Real>::g_;
  using RandVarFunctional<Real>::hv_;
  using RandVarFunctional<Real>::dualVector_;

  using RandVarFunctional<Real>::point_;
  using RandVarFunctional<Real>::weight_;

  using RandVarFunctional<Real>::computeValue;
  using RandVarFunctional<Real>::computeGradient;
  using RandVarFunctional<Real>::computeGradVec;
  using RandVarFunctional<Real>::computeHessVec;

  Real smoothHeaviside(const Real x, const int deriv = 0) const {
    const Real one(1), two(2);
    Real val(0);
    if (deriv == 0) {
      Real ex = std::exp(-two*x/eps_);
      val = one/(one+ex);
    }
    else if (deriv == 1) {
      Real ex = std::exp(-two*x/eps_);
      val = (two/eps_)*ex/std::pow(one+ex,2);
    }
    else if (deriv == 2) {
      Real ex = std::exp(two*x/eps_);
      val = std::pow(two/eps_,2)*ex*(one-ex)/std::pow(one+ex,3);
    }
    return val;
  }

public:
  SmoothedPOE(const Real threshold, const Real eps)
    : RandVarFunctional<Real>(),
      threshold_(threshold), eps_(eps) {}

  SmoothedPOE(Teuchos::ParameterList &parlist)
    : RandVarFunctional<Real>() {
    Teuchos::ParameterList &list = parlist.sublist("SOL").sublist("Probability").sublist("Smoothed POE");
    threshold_ = list.get<Real>("Threshold");
    eps_       = list.get<Real>("Smoothing Parameter");
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real sp  = smoothHeaviside(val-threshold_,0);
    if ( std::abs(sp) > ROL_EPSILON<Real>() ) {
      val_ += weight_*sp;
    }
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real spoe(0);
    sampler.sumAll(&val_,&spoe,1);
    return spoe;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real sp  = smoothHeaviside(val-threshold_,1);
    if ( std::abs(sp) > ROL_EPSILON<Real>() ) {
      computeGradient(*dualVector_,obj,x,tol);
      g_->axpy(weight_*sp,*dualVector_);
    }
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    sampler.sumAll(*g_,g);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real sp1 = smoothHeaviside(val-threshold_,1);
    Real sp2 = smoothHeaviside(val-threshold_,2);
    if ( std::abs(sp1) > ROL_EPSILON<Real>() ) {
      // Hessian only
      computeHessVec(*dualVector_,obj,v,x,tol);
      hv_->axpy(weight_*sp1,*dualVector_);
    }
    if ( std::abs(sp2) > ROL_EPSILON<Real>() ) {
      // Gradient only
      Real gv = computeGradVec(*dualVector_,obj,v,x,tol);
      hv_->axpy(weight_*sp2*gv,*dualVector_);
    }
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    sampler.sumAll(*hv_,hv);
  }
};

}

#endif
