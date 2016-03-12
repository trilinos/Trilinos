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

#ifndef ROL_HMCR_HPP
#define ROL_HMCR_HPP

#include "ROL_RiskMeasure.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_RiskVector.hpp"

namespace ROL {

template<class Real>
class HMCR : public RiskMeasure<Real> {
private:
  Teuchos::RCP<PlusFunction<Real> > plusFunction_;

  Real prob_;
  Real coeff_;

  unsigned order_;

  Teuchos::RCP<Vector<Real> > dualVector_;
  Teuchos::RCP<Vector<Real> > pHMdualVec0_;
  Teuchos::RCP<Vector<Real> > HMdualVec0_;
  Teuchos::RCP<Vector<Real> > pHMdualVec1_;
  Teuchos::RCP<Vector<Real> > HMdualVec1_;
  Teuchos::RCP<Vector<Real> > pHMdualVec2_;
  Teuchos::RCP<Vector<Real> > HMdualVec2_;
  Teuchos::RCP<Vector<Real> > pHMdualVec3_;
  Teuchos::RCP<Vector<Real> > HMdualVec3_;
  Real xvar_;
  Real vvar_;

  Real pnorm_;
  Real dpnorm_;
  Real dpnorm1_;
  Real pgv_;
  Real pgv1_;

  bool firstReset_;

public:

  HMCR( Real prob, Real coeff, unsigned order, Teuchos::RCP<PlusFunction<Real> > &pf )
    : RiskMeasure<Real>(), plusFunction_(pf), xvar_(0.0), vvar_(0.0),
      pnorm_(0.0), dpnorm_(0.0), dpnorm1_(0.0), pgv_(0.0), pgv1_(0.0),
      firstReset_(true) {
    Real zero(0), half(0.5), one(1);
    prob_  = ((prob  >= zero) ? ((prob  <= one) ? prob  : half) : half);
    coeff_ = ((coeff >= zero) ? ((coeff <= one) ? coeff : one) : one);
    order_ = ((order < 2) ? 2 : order);
  }

  HMCR( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>(), xvar_(0.0), vvar_(0.0),
      pnorm_(0.0), dpnorm_(0.0), dpnorm1_(0.0), pgv_(0.0), pgv1_(0.0),
      firstReset_(true) {
    Real zero(0), half(0.5), one(1);
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("HMCR");
    // Check HMCR inputs
    Real prob      = list.get("Confidence Level",half);
    prob_          = ((prob  >= zero) ? ((prob  <= one) ? prob  : half) : half);
    Real coeff     = list.get("Convex Combination Parameter",one);
    coeff_         = ((coeff >= zero) ? ((coeff <= one) ? coeff : one) : one);
    unsigned order = list.get("Order",2);
    order_         = ((order < 2) ? 2 : order);
    // Build (approximate) plus function
    plusFunction_  = Teuchos::rcp(new PlusFunction<Real>(list));
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    Real zero(0);
    RiskMeasure<Real>::reset(x0,x);
    xvar_ = Teuchos::dyn_cast<const RiskVector<Real> >(
              Teuchos::dyn_cast<const Vector<Real> >(x)).getStatistic();

    if ( firstReset_ ) {
      dualVector_  = (x0->dual()).clone();
      pHMdualVec0_ = (x0->dual()).clone();
      HMdualVec0_  = (x0->dual()).clone();
      pHMdualVec1_ = (x0->dual()).clone();
      HMdualVec1_  = (x0->dual()).clone();
      pHMdualVec2_ = (x0->dual()).clone();
      HMdualVec2_  = (x0->dual()).clone();
      pHMdualVec3_ = (x0->dual()).clone();
      HMdualVec3_  = (x0->dual()).clone();
      firstReset_  = false;
    }

    dualVector_->zero();
    pHMdualVec0_->zero(); HMdualVec0_->zero();

    pnorm_ = zero; dpnorm_ = zero;
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    Real zero(0);
    reset(x0,x);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(Teuchos::dyn_cast<const RiskVector<Real> >(
           Teuchos::dyn_cast<const Vector<Real> >(v)).getVector());
    vvar_ = Teuchos::dyn_cast<const RiskVector<Real> >(
              Teuchos::dyn_cast<const Vector<Real> >(v)).getStatistic();

    pHMdualVec1_->zero(); HMdualVec1_->zero();
    pHMdualVec2_->zero(); HMdualVec2_->zero();
    pHMdualVec3_->zero(); HMdualVec3_->zero();

    dpnorm1_ = zero; pgv_ = zero; pgv1_ = zero;
  }

  void update(const Real val, const Real weight) {
    // Expected value
    RiskMeasure<Real>::update(val,weight);
    // Higher moment
    Real pf = plusFunction_->evaluate(val-xvar_,0);
    pnorm_ += weight*std::pow(pf,(Real)order_);
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real one(1);
    // Expected value
    RiskMeasure<Real>::update(val,g,weight);
    // Higher moment
    Real pf0  = plusFunction_->evaluate(val-xvar_,0);
    Real pf1  = plusFunction_->evaluate(val-xvar_,1);

    Real rorder0 = (Real)order_;
    Real rorder1 = (Real)order_-one;

    Real pf0p0 = std::pow(pf0,rorder0);
    Real pf0p1 = std::pow(pf0,rorder1);

    pnorm_  += weight*pf0p0;
    dpnorm_ += weight*pf0p1*pf1;

    pHMdualVec0_->axpy(weight*pf0p1*pf1,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real one(1), two(2);
    // Expected value
    RiskMeasure<Real>::update(val,g,gv,hv,weight);
    // Higher moment
    Real pf0 = plusFunction_->evaluate(val-xvar_,0);
    Real pf1 = plusFunction_->evaluate(val-xvar_,1);
    Real pf2 = plusFunction_->evaluate(val-xvar_,2);

    Real rorder0 = (Real)order_;
    Real rorder1 = (Real)order_-one;
    Real rorder2 = (Real)order_-two;

    Real pf0p0 = std::pow(pf0,rorder0);
    Real pf0p1 = std::pow(pf0,rorder1);
    Real pf0p2 = std::pow(pf0,rorder2);

    Real coeff1 = pf0p1*pf1;
    Real coeff2 = rorder1*pf0p2*pf1*pf1 + pf0p1*pf2;

    pnorm_   += weight*pf0p0;
    dpnorm_  += weight*coeff1;
    dpnorm1_ += weight*coeff2;
    pgv_     += weight*coeff1*gv;
    pgv1_    += weight*coeff2*gv;

    pHMdualVec0_->axpy(weight*coeff1,g);
    pHMdualVec1_->axpy(weight*coeff2,g);
    pHMdualVec2_->axpy(weight*coeff1,hv);
    pHMdualVec3_->axpy(weight*coeff2*gv,g);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real one(1);
    std::vector<Real> val_in(2), val_out(2);
    val_in[0] = RiskMeasure<Real>::val_;
    val_in[1] = pnorm_;
    sampler.sumAll(&val_in[0],&val_out[0],2);
    return (one-coeff_)*val_out[0]
          + coeff_*(xvar_ + std::pow(val_out[1],one/(Real)order_)/(one-prob_));
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    Real zero(0), one(1);
    std::vector<Real> val_in(3), val_out(3);
    val_in[0] = RiskMeasure<Real>::val_;
    val_in[1] = pnorm_; val_in[2] = dpnorm_;

    sampler.sumAll(&val_in[0],&val_out[0],3);
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector_); 
    dualVector_->scale(one-coeff_);
    Real var = coeff_;

    if ( val_in[1] > zero ) {
      sampler.sumAll(*pHMdualVec0_,*HMdualVec0_);

      Real denom = (one-prob_)*std::pow(val_out[1],((Real)order_-one)/(Real)order_);

      var -= coeff_*((denom > zero) ? val_out[2]/denom : zero);

      dualVector_->axpy(coeff_/denom,*HMdualVec0_);
    }

    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setStatistic(var);
    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setVector(*dualVector_); 
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    Real zero(0), one(1);
    std::vector<Real> val_in(6), val_out(6);
    val_in[0] = RiskMeasure<Real>::val_; val_in[1] = pnorm_;
    val_in[2] = dpnorm_;                 val_in[3] = dpnorm1_;
    val_in[4] = pgv_;                    val_in[5] = pgv1_;

    sampler.sumAll(&val_in[0],&val_out[0],6);
    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector_);

    Real var = zero;
    dualVector_->scale(one-coeff_);

    if ( val_out[1] > zero ) {
      sampler.sumAll(*pHMdualVec0_,*HMdualVec0_); // E[pf^{p-1} pf' g]
      sampler.sumAll(*pHMdualVec1_,*HMdualVec1_); // E[{(p-1) pf^{p-2} pf' pf' + pf^{p-1} pf''} g]
      sampler.sumAll(*pHMdualVec2_,*HMdualVec2_); // E[pf^{p-1} pf' hv]
      sampler.sumAll(*pHMdualVec3_,*HMdualVec3_); // E[{(p-1) pf^{p-2} pf' pf' + pf^{p-1} pf''} g gv]

      Real rorder0 = (Real)order_;
      Real rorder1 = (Real)order_-one;
      Real rorder2 = (Real)(2*order_)-one;

      Real denom1 = (one-prob_)*std::pow(val_out[1],rorder1/rorder0);
      Real denom2 = (one-prob_)*std::pow(val_out[1],rorder2/rorder0);

      var = coeff_*((val_out[3]/denom1 - rorder1*val_out[2]*val_out[2]/denom2)*vvar_
                   -(val_out[5]/denom1 - rorder1*val_out[4]*val_out[2]/denom2));

      dualVector_->axpy(coeff_*(-vvar_/denom1),*HMdualVec1_);
      dualVector_->axpy(coeff_*(vvar_*rorder1*val_out[2]/denom2),*HMdualVec0_);
      dualVector_->axpy(coeff_/denom1,*HMdualVec3_);
      dualVector_->axpy(coeff_/denom1,*HMdualVec2_);
      dualVector_->axpy(coeff_*(-rorder1*val_out[4]/denom2),*HMdualVec0_);
    }

    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setStatistic(var);
    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setVector(*dualVector_); 
  }
};

}

#endif
