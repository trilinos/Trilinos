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
//
#ifndef PH_RISKOBJECTIVE_H
#define PH_RISKOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_RiskMeasureFactory.hpp"

/** @ingroup func_group
    \class ROL::PH_RiskObjective
    \brief Provides the interface for the progressive hedging risk objective.

    ---
*/
namespace ROL {

template <class Real>
class PH_RiskObjective : public Objective<Real> {
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

public:

  PH_RiskObjective(const Ptr<Objective<Real>> &obj,
                        ParameterList              &parlist)
    : obj_(obj),
      isValueComputed_(false),
      isGradientInitialized_(false),
      isGradientComputed_(false) {
    std::string risk = parlist.sublist("SOL").sublist("Risk Measure").get("Name","CVaR");
    ERiskMeasure ed = StringToERiskMeasure(risk);
    switch(ed) {
      case RISKMEASURE_CVAR:
             quad_ = makePtr<QuantileQuadrangle<Real>>(parlist);          break;
      case RISKMEASURE_MOREAUYOSIDACVAR:
             quad_ = makePtr<MoreauYosidaCVaR<Real>>(parlist);            break;
      case RISKMEASURE_GENMOREAUYOSIDACVAR:
             quad_ = makePtr<GenMoreauYosidaCVaR<Real>>(parlist);         break;
      case RISKMEASURE_LOGEXPONENTIAL:
             quad_ = makePtr<LogExponentialQuadrangle<Real>>(parlist);    break;
      case RISKMEASURE_SAFETYMARGIN:
             quad_ = makePtr<MeanVarianceQuadrangle<Real>>(parlist);      break;
      case RISKMEASURE_TRUNCATEDMEAN:
             quad_ = makePtr<TruncatedMeanQuadrangle<Real>>(parlist);     break;
      case RISKMEASURE_LOGQUANTILE:
             quad_ = makePtr<LogQuantileQuadrangle<Real>>(parlist);       break;
      case RISKMEASURE_SMOOTHEDWORSTCASE:
             quad_ = makePtr<SmoothedWorstCaseQuadrangle<Real>>(parlist); break;
//      case RISKMEASURE_CHI2DIVERGENCE:
//             return makePtr<Chi2Divergence<Real>>(parlist);
//      case RISKMEASURE_KLDIVERGENCE:
//             return makePtr<KLDivergence<Real>>(parlist);
      default:
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
                                   "Invalid risk measure type " << risk << "!");
    }
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
    Real reg = quad_->regret(val_-(*xstat)[0],0);
    return (*xstat)[0] + reg;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    Ptr<Vector<Real>>            gvec  = getVector(g);
    Ptr<std::vector<Real>>       gstat = getStat(g);
    Ptr<const Vector<Real>>      xvec  = getConstVector(x);
    Ptr<const std::vector<Real>> xstat = getConstStat(x);
    getValue(*xvec,tol);
    Real reg = quad_->regret(val_-(*xstat)[0],1);
    getGradient(*xvec,tol);
    gvec->set(*g_); gvec->scale(reg);
    (*gstat)[0] = static_cast<Real>(1)-reg;
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Ptr<Vector<Real>>            hvec  = getVector(hv);
    Ptr<std::vector<Real>>       hstat = getStat(hv);
    Ptr<const Vector<Real>>      vvec  = getConstVector(v);
    Ptr<const std::vector<Real>> vstat = getConstStat(v);
    Ptr<const Vector<Real>>      xvec  = getConstVector(x);
    Ptr<const std::vector<Real>> xstat = getConstStat(x);
    getValue(*xvec,tol);
    Real reg1 = quad_->regret(val_-(*xstat)[0],1);
    Real reg2 = quad_->regret(val_-(*xstat)[0],2);
    getGradient(*xvec,tol);
    Real gv   = vvec->dot(g_->dual());
    obj_->hessVec(*hvec,*vvec,*xvec,tol);
    hvec->scale(reg1); hvec->axpy(reg2*(gv-(*vstat)[0]),*g_);
    (*hstat)[0] = reg2*((*vstat)[0]-gv);
  }

  void setParameter(const std::vector<Real> &param) {
    obj_->setParameter(param);
    Objective<Real>::setParameter(param);
  }

};

}
#endif
