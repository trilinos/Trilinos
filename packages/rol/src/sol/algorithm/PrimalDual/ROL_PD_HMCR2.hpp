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

#ifndef ROL_PD_HMCR2_HPP
#define ROL_PD_HMCR2_HPP

#include "ROL_PD_RandVarFunctional.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
class PD_HMCR2 : public PD_RandVarFunctional<Real> {
private:
  Real beta_;
  Real mScalar1_, mScalar2_;

  Ptr<SampledScalar<Real>> values_;
  Ptr<SampledScalar<Real>> gradvecs_;
  Ptr<SampledVector<Real>> gradients_;
  Ptr<SampledVector<Real>> hessvecs_;

  using RandVarFunctional<Real>::val_;
  using RandVarFunctional<Real>::g_;
  using RandVarFunctional<Real>::gv_;
  using RandVarFunctional<Real>::hv_;
  using RandVarFunctional<Real>::dualVector_;

  using RandVarFunctional<Real>::point_;
  using RandVarFunctional<Real>::weight_;

  using RandVarFunctional<Real>::computeValue;
  using RandVarFunctional<Real>::computeGradient;
  using RandVarFunctional<Real>::computeGradVec;
  using RandVarFunctional<Real>::computeHessVec;

  using PD_RandVarFunctional<Real>::setValue;
  using PD_RandVarFunctional<Real>::setMultiplier;
  using PD_RandVarFunctional<Real>::getMultiplier;
  using PD_RandVarFunctional<Real>::getPenaltyParameter;

  void initializeStorage(void) {
    values_    = makePtr<SampledScalar<Real>>();
    gradvecs_  = makePtr<SampledScalar<Real>>();
    gradients_ = makePtr<SampledVector<Real>>();
    hessvecs_  = makePtr<SampledVector<Real>>();

    RandVarFunctional<Real>::setStorage(values_,gradients_);
    RandVarFunctional<Real>::setHessVecStorage(gradvecs_,hessvecs_);
  }

  void clear(void) {
    gradvecs_->update();
    hessvecs_->update();
  }

  void checkInputs(void) {
    Real zero(0), one(1);
    ROL_TEST_FOR_EXCEPTION((beta_ < zero) || (beta_ >= one), std::invalid_argument,
      ">>> ERROR (ROL::PD_HMCR2): Confidence parameter beta is out of range!");
    initializeStorage();
  }

public:
  PD_HMCR2(const Real beta)
    : PD_RandVarFunctional<Real>(), beta_(beta) {
    checkInputs();
  }

  void setStorage(const Ptr<SampledScalar<Real>> &value_storage,
                  const Ptr<SampledVector<Real>> &gradient_storage) {
    values_    = value_storage;
    gradients_ = gradient_storage;
    PD_RandVarFunctional<Real>::setStorage(values_,gradients_);
  }

  void setHessVecStorage(const Ptr<SampledScalar<Real>> &gradvec_storage,
                         const Ptr<SampledVector<Real>> &hessvec_storage) {
    gradvecs_ = gradvec_storage;
    hessvecs_ = hessvec_storage;
    PD_RandVarFunctional<Real>::setHessVecStorage(gradvecs_,hessvecs_);
  }

  Real computeDual(SampleGenerator<Real> &sampler) {
    const Real zero(0), two(2);
    Real val(0), lold(0), lnew(0), mdiff(0), gdiff(0), sum(0), gsum(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(val, sampler.getMyPoint(i));
      getMultiplier(lold, sampler.getMyPoint(i));
      lnew = std::max(zero, getPenaltyParameter()*val+lold);
      sum += sampler.getMyWeight(i) * std::pow(lnew,two);
    }
    sampler.sumAll(&sum,&gsum,1);
    gsum = std::sqrt(gsum);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(val, sampler.getMyPoint(i));
      getMultiplier(lold, sampler.getMyPoint(i));
      lnew = std::max(zero, getPenaltyParameter()*val+lold)/gsum;
      mdiff += sampler.getMyWeight(i) * std::pow(lnew-lold,2);
      setMultiplier(lnew, sampler.getMyPoint(i));
    }
    sampler.sumAll(&mdiff,&gdiff,1);
    gdiff = std::sqrt(gdiff);
    return gdiff;
  }
 
  void initialize(const Vector<Real> &x) {
    PD_RandVarFunctional<Real>::initialize(x);
    mScalar1_ = static_cast<Real>(0);
    mScalar2_ = static_cast<Real>(0);
    clear();
  }
  
  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    const Real zero(0), two(2);
    Real lam(0);
    getMultiplier(lam, point_);
    Real val = computeValue(obj, x, tol);
    Real arg = val - xstat[0];
    Real  pf = std::max(zero, arg + lam/getPenaltyParameter());
    val_    += weight_ * std::pow(pf,two);
    setValue(arg, point_);
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    const Real half(0.5), one(1);
    Real ev(0);
    sampler.sumAll(&val_, &ev, 1);
    Real norm = std::sqrt(ev);
    Real sig  = one/(one-beta_);
    Real val  = (norm <= sig/getPenaltyParameter()
                  ? half * getPenaltyParameter() * ev
                  : sig * (norm - sig*half/getPenaltyParameter()));
    return xstat[0] + val;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    const Real zero(0), two(2);
    Real lam(0);
    getMultiplier(lam, point_);
    Real val = computeValue(obj, x, tol);
    Real arg = val - xstat[0];
    Real  pf = std::max(zero, arg + lam/getPenaltyParameter());
    if ( pf > zero ) {
      val_ += weight_ * pf;
      gv_  += weight_ * std::pow(pf,two);
      computeGradient(*dualVector_, obj, x, tol);
      g_->axpy(weight_ * pf, *dualVector_);
    }
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    const Real one(1);
    std::vector<Real> mv = {val_, gv_};
    std::vector<Real> ev(2,0);
    sampler.sumAll(&mv[0], &ev[0], 2);
    Real norm = std::sqrt(ev[1]);
    Real sig  = one/(one-beta_);
    Real scal = (norm <= sig/getPenaltyParameter()
                  ? getPenaltyParameter()
                  : sig/norm);
    gstat[0] = one - scal * ev[0];
    sampler.sumAll(*g_, g);
    g.scale(scal);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    const Real zero(0), two(2);
    Real lam(0);
    getMultiplier(lam, point_);
    Real val = computeValue(obj, x, tol);
    Real arg = val - xstat[0];
    Real pf  = std::max(zero, arg + lam/getPenaltyParameter());
    if ( pf > zero ) {
      val_      += weight_ * std::pow(pf,two);
      mScalar1_ += weight_ * pf;

      Real gv = computeGradVec(*dualVector_, obj, v, x, tol);
      mScalar2_ += weight_ * pf * gv;
      gv_       += weight_ * (vstat[0] - gv);
      g_->axpy(weight_ * pf, *dualVector_);
      hv_->axpy(weight_ * (gv - vstat[0]), *dualVector_);
      computeHessVec(*dualVector_, obj, v, x, tol);
      hv_->axpy(weight_ * pf, *dualVector_);
    }
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    const Real one(1);
    std::vector<Real> mv = {val_, gv_, mScalar1_, mScalar2_};
    std::vector<Real> ev(4,0);
    sampler.sumAll(&mv[0],&ev[0],4);
    Real norm = std::sqrt(ev[0]);
    Real sig  = one/(one-beta_);
    Real scal = (norm <= sig/getPenaltyParameter()
                  ? getPenaltyParameter()
                  : sig/norm);
    hvstat[0] = scal * ev[1];
    sampler.sumAll(*hv_,hv);
    hv.scale(scal);
    if (norm > sig/getPenaltyParameter()) {
      Real norm3 = ev[0]*norm;
      hvstat[0] += sig/norm3 * (ev[3] - ev[2]*vstat[0]) * ev[2];
      dualVector_->zero();
      sampler.sumAll(*g_,*dualVector_);
      hv.axpy(sig/norm3 * (ev[2]*vstat[0] - ev[3]),*dualVector_);
    }
  }
};

}

#endif
