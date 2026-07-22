// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PH_REGRETOBJECTIVE_H
#define PH_REGRETOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_RegretMeasureFactory.hpp"

/** @ingroup func_group
    \class ROL::PH_RegretObjective
    \brief Provides the interface for the progressive hedging regret objective.

    ---
*/
namespace ROL {

template <class Real>
class PH_RegretObjective : public Objective<Real> {
private:
  const Ptr<Objective<Real>> obj_;
  Ptr<ExpectationQuad<Real>> quad_;

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

public:

  PH_RegretObjective(const Ptr<Objective<Real>> &obj,
                     ParameterList              &parlist)
    : obj_(obj),
      isValueComputed_(false),
      isGradientInitialized_(false),
      isGradientComputed_(false) {
    std::string regret = parlist.sublist("SOL").sublist("Regret Measure").get("Name","Mean Absolute Loss");
    ERegretMeasure ed = StringToERegretMeasure(regret);
    switch(ed) {
      case REGRETMEASURE_MEANABSOLUTELOSS:
             quad_ = makePtr<QuantileQuadrangle<Real>>(parlist);          break;
      case REGRETMEASURE_MOREAUYOSIDAMEANABSOLUTELOSS:
             quad_ = makePtr<MoreauYosidaCVaR<Real>>(parlist);            break;
      case REGRETMEASURE_GENMOREAUYOSIDAMEANABSOLUTELOSS:
             quad_ = makePtr<GenMoreauYosidaCVaR<Real>>(parlist);         break;
      case REGRETMEASURE_EXPONENTIAL:
             quad_ = makePtr<LogExponentialQuadrangle<Real>>(parlist);    break;
      case REGRETMEASURE_MEANL2:
             quad_ = makePtr<MeanVarianceQuadrangle<Real>>(parlist);      break;
      case REGRETMEASURE_TRUNCATEDMEAN:
             quad_ = makePtr<TruncatedMeanQuadrangle<Real>>(parlist);     break;
      case REGRETMEASURE_LOGQUANTILE:
             quad_ = makePtr<LogQuantileQuadrangle<Real>>(parlist);       break;
      case REGRETMEASURE_SMOOTHEDWORSTCASE:
             quad_ = makePtr<SmoothedWorstCaseQuadrangle<Real>>(parlist); break;
      default:
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               "Invalid regret measure type " << regret << "!");
    }
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
    isValueComputed_    = false;
    isGradientComputed_ = false;
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    getValue(x,tol);
    Real reg = quad_->regret(val_,0);
    return reg;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    getValue(x,tol);
    Real reg = quad_->regret(val_,1);
    getGradient(x,tol);
    g.set(*g_); g.scale(reg);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    getValue(x,tol);
    Real reg1 = quad_->regret(val_,1);
    Real reg2 = quad_->regret(val_,2);
    getGradient(x,tol);
    //Real gv   = v.dot(g_->dual());
    Real gv   = v.apply(*g_);
    obj_->hessVec(hv,v,x,tol);
    hv.scale(reg1); hv.axpy(reg2*gv,*g_);
  }

  void setParameter(const std::vector<Real> &param) {
    obj_->setParameter(param);
    Objective<Real>::setParameter(param);
  }

};

}
#endif
