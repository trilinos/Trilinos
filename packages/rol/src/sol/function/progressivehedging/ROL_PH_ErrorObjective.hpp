// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PH_ERROROBJECTIVE_H
#define PH_ERROROBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_ErrorMeasureFactory.hpp"

/** @ingroup func_group
    \class ROL::PH_ErrorObjective
    \brief Provides the interface for the progressive hedging error objective.

    ---
*/
namespace ROL {

template <class Real>
class PH_ErrorObjective : public Objective<Real> {
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

  PH_ErrorObjective(const Ptr<Objective<Real>> &obj,
                     ParameterList              &parlist)
    : obj_(obj),
      isValueComputed_(false),
      isGradientInitialized_(false),
      isGradientComputed_(false) {
    std::string risk = parlist.sublist("SOL").sublist("Error Measure").get("Name","Least Squares");
    EErrorMeasure ed = StringToEErrorMeasure(risk);
    switch(ed) {
      case ERRORMEASURE_MEANVARIANCEQUADRANGLE:
             quad_ = makePtr<MeanVarianceQuadrangle<Real>>(parlist);      break;
      case ERRORMEASURE_TRUNCATEDMEANQUADRANGLE:
             quad_ = makePtr<TruncatedMeanQuadrangle<Real>>(parlist);     break;
      case ERRORMEASURE_QUANTILEQUADRANGLE:
             quad_ = makePtr<QuantileQuadrangle<Real>>(parlist);          break;
      case ERRORMEASURE_MOREAUYOSIDACVAR:
             quad_ = makePtr<MoreauYosidaCVaR<Real>>(parlist);            break;
      case ERRORMEASURE_GENMOREAUYOSIDACVAR:
             quad_ = makePtr<GenMoreauYosidaCVaR<Real>>(parlist);         break;
      case ERRORMEASURE_LOGEXPONENTIALQUADRANGLE:
             quad_ = makePtr<LogExponentialQuadrangle<Real>>(parlist);    break;
      case ERRORMEASURE_LOGQUANTILEQUADRANGLE:
             quad_ = makePtr<LogQuantileQuadrangle<Real>>(parlist);       break;
      case ERRORMEASURE_SMOOTHEDWORSTCASEQUADRANGLE:
             quad_ = makePtr<SmoothedWorstCaseQuadrangle<Real>>(parlist); break;
      default:
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                               "Invalid error measure type " << risk << "!");
    }
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    obj_->update(x,flag,iter);
    isValueComputed_    = false;
    isGradientComputed_ = false;
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    getValue(x,tol);
    Real err = quad_->error(val_,0);
    return err;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    getValue(x,tol);
    Real err = quad_->error(val_,1);
    getGradient(x,tol);
    g.set(*g_); g.scale(err);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    getValue(x,tol);
    Real err1 = quad_->error(val_,1);
    Real err2 = quad_->error(val_,2);
    getGradient(x,tol);
    //Real gv   = v.dot(g_->dual());
    Real gv   = v.apply(*g_);
    obj_->hessVec(hv,v,x,tol);
    hv.scale(err1); hv.axpy(err2*gv,*g_);
  }

  void setParameter(const std::vector<Real> &param) {
    obj_->setParameter(param);
    Objective<Real>::setParameter(param);
  }

};

}
#endif
