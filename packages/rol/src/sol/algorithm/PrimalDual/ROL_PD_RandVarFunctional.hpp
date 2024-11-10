// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PD_RANDVARFUNCTIONAL_HPP
#define ROL_PD_RANDVARFUNCTIONAL_HPP

#include "ROL_RandVarFunctional.hpp"

namespace ROL {

template<class Real>
class PD_RandVarFunctional : public RandVarFunctional<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  Real pen_;
  int update_;
  bool setData_;

  Ptr<ScalarController<Real>> values_;
  Ptr<ScalarController<Real>> multipliers_;
  Ptr<ScalarController<Real>> multipliers_new_;

protected:
  // Set value data at current parameter
  void setValue(const Real val, const std::vector<Real> &pt) {
    values_->set(val, pt);
  }

  // Get multiplier at current parameter
  void getMultiplier(Real &lam, const std::vector<Real> &pt) const {
    multipliers_->get(lam, pt);
  }

  void setMultiplier(Real &lam, const std::vector<Real> &pt) {
    multipliers_new_->set(lam, pt);
  }

  // Get penalty parameter
  Real getPenaltyParameter(void) const {
    return pen_;
  }

  // Smooth plus function approximation
  Real ppf(const Real x, const Real t, const Real r, const int deriv = 0) const {
    const Real zero(0), half(0.5), one(1), arg(r*x+t);
    Real val(0);
    if ( arg < zero ) {
      val = (deriv==0 ? -half*t*t/r : zero);
    }
    else if ( zero <= arg && arg <= one ) {
      val = (deriv==0 ? half*r*x*x+t*x
          : (deriv==1 ? arg : r));
    }
    else {
      val = (deriv==0 ? (arg-half*(t*t+one))/r
          : (deriv==1 ? one : zero));
    }
    return val;
  }

public:
  PD_RandVarFunctional(void)
    : RandVarFunctional<Real>(), pen_(1.0), update_(0), setData_(true) {
    values_          = makePtr<ScalarController<Real>>();
    multipliers_     = makePtr<ScalarController<Real>>();
    multipliers_new_ = makePtr<ScalarController<Real>>();
  }

  void setData(SampleGenerator<Real> &sampler, const Real pen, const Real lam = 0.0) {
    if (setData_) {
      pen_ = pen;
      for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
        multipliers_->set(lam, sampler.getMyPoint(i));
      }
      setData_ = false;
    }
  }

  virtual Real computeDual(SampleGenerator<Real> &sampler) {
    const Real zero(0), one(1);
    Real val(0), lold(0), lnew(0), mdiff(0), gdiff(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      values_->get(val, sampler.getMyPoint(i));
      multipliers_->get(lold, sampler.getMyPoint(i));
      if (update_ == 0) {
        //lnew = ppf(val, lold, pen_, 1);
        lnew = std::min(one, std::max(zero, pen_*val+lold));
      }
      else {
        lnew = (val < zero ? zero : one);
      }
      mdiff += sampler.getMyWeight(i) * std::pow(lnew-lold,2);
      multipliers_new_->set(lnew, sampler.getMyPoint(i));
    }
    sampler.sumAll(&mdiff,&gdiff,1);
    gdiff = std::sqrt(gdiff);
    return gdiff;
  }

  void updateDual(SampleGenerator<Real> &sampler) {
    Real lam(0);
    for (int i = sampler.start(); i < sampler.numMySamples(); ++i) {
      multipliers_new_->get(lam, sampler.getMyPoint(i));
      multipliers_->set(lam, sampler.getMyPoint(i));
    }
  }

  void updatePenalty(const Real pen) {
    pen_ = pen;
  }

  virtual void setStorage(const Ptr<ScalarController<Real>> &value_storage,
                          const Ptr<VectorController<Real>> &gradient_storage) {
    RandVarFunctional<Real>::setStorage(value_storage,gradient_storage);
  }

  virtual void setHessVecStorage(const Ptr<ScalarController<Real>> &gradvec_storage,
                                 const Ptr<VectorController<Real>> &hessvec_storage) {
    RandVarFunctional<Real>::setHessVecStorage(gradvec_storage,hessvec_storage);
  }
 
  virtual void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
  }
};

}

#endif
