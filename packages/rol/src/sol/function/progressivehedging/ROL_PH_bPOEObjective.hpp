// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PH_BPOEOBJECTIVE_H
#define PH_BPOEOBJECTIVE_H

#include "ROL_Objective.hpp"

/** @ingroup func_group
    \class ROL::PH_bPOEObjective
    \brief Provides the interface for the progressive hedging probability objective.

    ---
*/
namespace ROL {

template <class Real>
class PH_bPOEObjective : public Objective<Real> {
private:
  const Ptr<Objective<Real>> obj_;
  Real threshold_;
  Real order_;

  bool isValueComputed_;
  Real val_;

  bool isGradientInitialized_;
  bool isGradientComputed_;
  Ptr<Vector<Real>> g_;

  void getValue(const Vector<Real> &x, Real &tol) {
    if (!isValueComputed_) {
      val_ = obj_->value(x,tol);
      isValueComputed_ = true;
    }
  }

  void getGradient(const Vector<Real> &x, Real &tol) {
    if (!isGradientInitialized_) {
      g_ = x.dual().clone();
      isGradientInitialized_ = true;
    }
    if (!isGradientComputed_) {
      obj_->gradient(*g_,x,tol);
      isGradientComputed_ = true;
    }
  }

  Ptr<const Vector<Real>> getConstVector(const Vector<Real> &x) const {
    const RiskVector<Real> &xrv = dynamic_cast<const RiskVector<Real>&>(x);
    return xrv.getVector();
  }

  Ptr<Vector<Real>> getVector(Vector<Real> &x) const {
    RiskVector<Real> &xrv = dynamic_cast<RiskVector<Real>&>(x);
    return xrv.getVector();
  }

  Ptr<const std::vector<Real>> getConstStat(const Vector<Real> &x) const {
    const RiskVector<Real> &xrv = dynamic_cast<const RiskVector<Real>&>(x);
    Ptr<const std::vector<Real>> xstat = xrv.getStatistic();
    if (xstat == nullPtr) {
      xstat = makePtr<const std::vector<Real>>(0);
    }
    return xstat;
  }

  Ptr<std::vector<Real>> getStat(Vector<Real> &x) const {
    RiskVector<Real> &xrv = dynamic_cast<RiskVector<Real>&>(x);
    Ptr<std::vector<Real>> xstat = xrv.getStatistic();
    if (xstat == nullPtr) {
      xstat = makePtr<std::vector<Real>>(0);
    }
    return xstat;
  }

  // pth power of the positive part function
  Real pplus(const Real x, const int deriv = 0) const {
    const Real zero(0), one(1), two(2), three(3);
    Real val(0);
    if (x > zero) {
      if (deriv==0) {
        val = (order_==one ? x
               : std::pow(x,order_));
      }
      else if (deriv==1) {
        val = order_*(order_==one ? one
                      : (order_==two ? x
                         : std::pow(x,order_-one)));
      }
      else if (deriv==2) {
        val = order_*(order_-one)*(order_==one ? zero
                                   : (order_==two ? one
                                      : (order_==three ? x
                                         : std::pow(x,order_-two))));
      }
    }
    return val;
  }

  Real bPOEobjective(const Real t, const Real x, const int deriv = 0) const {
    const Real one(1);
    Real arg = t*(x-threshold_)+one;
    return pplus(arg,deriv);
  }

public:

  PH_bPOEObjective(const Ptr<Objective<Real>> &obj,
                   ParameterList              &parlist)
    : obj_(obj),
      isValueComputed_(false),
      isGradientInitialized_(false),
      isGradientComputed_(false) {
    ParameterList &list = parlist.sublist("SOL").sublist("Probability").sublist("bPOE");
    threshold_ = list.get<Real>("Threshold");
    order_     = list.get<Real>("Moment Order");
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    Ptr<const Vector<Real>> xvec = getConstVector(x);
    obj_->update(*xvec,flag,iter);
    isValueComputed_    = false;
    isGradientComputed_ = false;
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    Ptr<const Vector<Real>> xvec = getConstVector(x);
    Ptr<const std::vector<Real>> xstat = getConstStat(x);
    getValue(*xvec,tol);
    Real xt = (*xstat)[0];
    Real prob = bPOEobjective(xt,val_,0);
    return prob;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    Ptr<Vector<Real>>            gvec  = getVector(g);
    Ptr<std::vector<Real>>       gstat = getStat(g);
    Ptr<const Vector<Real>>      xvec  = getConstVector(x);
    Ptr<const std::vector<Real>> xstat = getConstStat(x);
    getValue(*xvec,tol);
    Real xt = (*xstat)[0], diff = val_-threshold_;
    Real prob = bPOEobjective(xt,val_,1);
    getGradient(*xvec,tol);
    gvec->set(*g_);
    gvec->scale(prob*xt);
    (*gstat)[0] = prob*diff;
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Ptr<Vector<Real>>            hvec  = getVector(hv);
    Ptr<std::vector<Real>>       hstat = getStat(hv);
    Ptr<const Vector<Real>>      vvec  = getConstVector(v);
    Ptr<const std::vector<Real>> vstat = getConstStat(v);
    Ptr<const Vector<Real>>      xvec  = getConstVector(x);
    Ptr<const std::vector<Real>> xstat = getConstStat(x);
    getValue(*xvec,tol);
    Real xt = (*xstat)[0], vt = (*vstat)[0], diff = val_-threshold_;
    Real prob1 = bPOEobjective(xt,val_,1);
    Real prob2 = bPOEobjective(xt,val_,2);
    getGradient(*xvec,tol);
    //Real gv   = vvec->dot(g_->dual());
    Real gv   = vvec->apply(*g_);
    obj_->hessVec(*hvec,*vvec,*xvec,tol);
    hvec->scale(prob1*xt);
    hvec->axpy(prob2*xt*(vt*diff+xt*gv)+vt*prob1,*g_);
    (*hstat)[0] = prob2*std::pow(diff,2)*vt+(prob2*diff*xt+prob1)*gv;
  }

  void setParameter(const std::vector<Real> &param) {
    obj_->setParameter(param);
    Objective<Real>::setParameter(param);
  }

};

}
#endif
