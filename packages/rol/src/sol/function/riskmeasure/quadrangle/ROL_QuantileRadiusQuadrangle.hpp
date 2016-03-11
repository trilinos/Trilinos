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

#ifndef ROL_QUANTILERADIUSQUADRANGLE_HPP
#define ROL_QUANTILERADIUSQUADRANGLE_HPP

#include "ROL_RiskMeasure.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_RiskVector.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real>
class QuantileRadiusQuadrangle : public RiskMeasure<Real> {
private:
  Teuchos::RCP<PlusFunction<Real> > plusFunction_;

  Real prob_;
  Real coeff_;

  Teuchos::RCP<Vector<Real> > dualVector_;
  std::vector<Real> xvar_;
  std::vector<Real> vvar_;

  std::vector<Real> vec_;

  bool firstReset_;

public:

  QuantileRadiusQuadrangle( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>(), firstReset_(true) {
    Real zero(0), half(0.5), one(1);
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Quantile-Radius Quadrangle");
    // Grab probability and coefficient arrays
    prob_  = list.get("Confidence Level",half);
    coeff_ = list.get("Coefficient",one);
    // Check inputs
    TEUCHOS_TEST_FOR_EXCEPTION((prob_>one || prob_<zero), std::invalid_argument,
      ">>> ERROR (ROL::QuantileRadiusQuadrangle): Confidence level out of range!");
    TEUCHOS_TEST_FOR_EXCEPTION((coeff_<zero), std::invalid_argument,
      ">>> ERROR (ROL::QuantileRadiusQuadrangle): Coefficient is negative!");
    // Build (approximate) plus function
    plusFunction_ = Teuchos::rcp(new PlusFunction<Real>(list));
    // Initialize temporary storage
    xvar_.clear(); xvar_.resize(2,zero);
    vvar_.clear(); vvar_.resize(2,zero);
    vec_.clear();  vec_.resize(2,zero);
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    Real zero(0);
    RiskMeasure<Real>::reset(x0,x);
    for (int i = 0; i < 2; i++) {
      xvar_[i] = Teuchos::dyn_cast<const RiskVector<Real> >(
                 Teuchos::dyn_cast<const Vector<Real> >(x)).getStatistic(i);
      vec_[i]  = zero;
    }
    if ( firstReset_ ) {
      dualVector_ = (x0->dual()).clone();
      firstReset_ = false;
    }
    dualVector_->zero();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const RiskVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(v)).getVector());
    for (int i = 0; i < 2; i++) {
      vvar_[i] = Teuchos::dyn_cast<const RiskVector<Real> >(
                 Teuchos::dyn_cast<const Vector<Real> >(v)).getStatistic(i);
    }
  }

  void update(const Real val, const Real weight) {
    Real half(0.5), one(1);
    Real pf1 = plusFunction_->evaluate(val-xvar_[0],0);
    Real pf2 = plusFunction_->evaluate(-val-xvar_[1],0);
    RiskMeasure<Real>::val_ += weight*(val + half*coeff_/(one-prob_)*(pf1 + pf2));
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real half(0.5), one(1);
    Real pf1 = plusFunction_->evaluate(val-xvar_[0],1);
    Real pf2 = plusFunction_->evaluate(-val-xvar_[1],1);
    Real c   = half*weight*coeff_/(one-prob_);
    vec_[0] -= c*pf1;
    vec_[1] -= c*pf2;
    RiskMeasure<Real>::g_->axpy(weight + c * (pf1 - pf2),g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real half(0.5), one(1);
    Real pf11 = plusFunction_->evaluate(val-xvar_[0],1);
    Real pf12 = plusFunction_->evaluate(val-xvar_[0],2);
    Real pf21 = plusFunction_->evaluate(-val-xvar_[1],1);
    Real pf22 = plusFunction_->evaluate(-val-xvar_[1],2);
    Real c    = half*weight*coeff_/(one-prob_);
    vec_[0]  -= c*pf12*(gv-vvar_[0]);
    vec_[1]  -= c*pf22*(-gv-vvar_[1]);
    RiskMeasure<Real>::hv_->axpy(c*(pf12*(gv-vvar_[0]) + pf22*(-gv-vvar_[1])),g);
    RiskMeasure<Real>::hv_->axpy(weight + c * (pf11 - pf21),hv);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val  = RiskMeasure<Real>::val_, cvar(0), half(0.5);
    sampler.sumAll(&val,&cvar,1);
    cvar += half*coeff_*(xvar_[0] + xvar_[1]);
    return cvar;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    Real half(0.5);
    RiskVector<Real> &gs = Teuchos::dyn_cast<RiskVector<Real> >(Teuchos::dyn_cast<Vector<Real> >(g));
    std::vector<Real> var(2);
    sampler.sumAll(&vec_[0],&var[0],2);
    
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector_);
    var[0] += half*coeff_;
    var[1] += half*coeff_;
    gs.setStatistic(var);
    gs.setVector(*(Teuchos::rcp_dynamic_cast<Vector<Real> >(dualVector_))); 
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    RiskVector<Real> &hs = Teuchos::dyn_cast<RiskVector<Real> >(Teuchos::dyn_cast<Vector<Real> >(hv));
    std::vector<Real> var(2);
    sampler.sumAll(&vec_[0],&var[0],2);

    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector_);
    hs.setStatistic(var);
    hs.setVector(*(Teuchos::rcp_dynamic_cast<Vector<Real> >(dualVector_)));
  }
};

}

#endif
