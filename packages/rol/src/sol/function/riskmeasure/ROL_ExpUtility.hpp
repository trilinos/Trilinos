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

#ifndef ROL_EXPUTILITY_HPP
#define ROL_EXPUTILITY_HPP

#include "ROL_RiskMeasure.hpp"

namespace ROL {

/** @ingroup risk_group
    \class ROL::ExpUtility
    \brief Provides an interface for the entropic risk.

    The entropic risk measure (also called the exponential utility and the
    log-exponential risk measure) is
    \f[
       \mathcal{R}(X) = \lambda
       \log\mathbb{E}\left[\exp\left(\frac{X}{\lambda}\right)\right]
    \f]
    for \f$\lambda > 0\f$.  The entropic risk is convex, translation
    equivariant and monotonic.
*/

template<class Real>
class ExpUtility : public RiskMeasure<Real> {
private:
  Teuchos::RCP<Vector<Real> > scaledGradient_;
  Teuchos::RCP<Vector<Real> > dualVector1_;
  Teuchos::RCP<Vector<Real> > dualVector2_;
  bool firstReset_;

  Real coeff_;

  void checkInputs(void) const {
    Real zero(0);
    TEUCHOS_TEST_FOR_EXCEPTION((coeff_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::ExpUtility): Rate must be positive!");
  }

public:
  /** \brief Constructor.

      @param[in]     coeff    is the scale parameter \f$\lambda\f$
  */
  ExpUtility(const Real coeff = 1)
    : RiskMeasure<Real>(), firstReset_(true), coeff_(coeff) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measures"->"Exponential Utility"
      and withing the "Exponential Utility" sublist should have
      \li "Rate" (greater than 0).
  */
  ExpUtility(Teuchos::ParameterList &parlist)
    : RiskMeasure<Real>(), firstReset_(true) {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Exponential Utility");
    coeff_ = list.get<Real>("Rate");
    checkInputs();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    if ( firstReset_ ) {
      scaledGradient_ = (x0->dual()).clone();
      dualVector1_ = (x0->dual()).clone();
      dualVector2_ = (x0->dual()).clone();
      firstReset_ = false;
    }
    scaledGradient_->zero(); dualVector1_->zero(); dualVector2_->zero();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(
           Teuchos::dyn_cast<const RiskVector<Real> >(v).getVector());
  }

  void update(const Real val, const Real weight) {
    RiskMeasure<Real>::val_ += weight * std::exp(coeff_*val);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_, ev(0);
    sampler.sumAll(&val,&ev,1);
    return std::log(ev)/coeff_;
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real ev = std::exp(coeff_*val);
    RiskMeasure<Real>::val_ += weight * ev;
    RiskMeasure<Real>::g_->axpy(weight*ev,g);
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    Real val = RiskMeasure<Real>::val_, ev(0), one(1);
    sampler.sumAll(&val,&ev,1);

    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector1_);
    dualVector1_->scale(one/ev);

    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setVector(*dualVector1_);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
                      const Real weight) {
    Real ev = std::exp(coeff_*val);
    RiskMeasure<Real>::val_ += weight * ev;
    RiskMeasure<Real>::gv_  -= weight * ev * gv;
    RiskMeasure<Real>::g_->axpy(weight*ev,g);
    RiskMeasure<Real>::hv_->axpy(weight*ev,hv);
    scaledGradient_->axpy(weight*ev*gv,g);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    Real one(1);
    std::vector<Real> myval(2), val(2);
    myval[0] = RiskMeasure<Real>::val_;
    myval[1] = RiskMeasure<Real>::gv_;
    sampler.sumAll(&myval[0],&val[0],2);

    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector1_);

    sampler.sumAll(*scaledGradient_,*dualVector2_);
    dualVector1_->axpy(coeff_,*dualVector2_);
    dualVector1_->scale(one/val[0]);

    dualVector2_->zero();
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector2_);
    dualVector1_->axpy(coeff_*val[1]/(val[0]*val[0]),*dualVector2_);

    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setVector(*dualVector1_);
  }
};

}

#endif
