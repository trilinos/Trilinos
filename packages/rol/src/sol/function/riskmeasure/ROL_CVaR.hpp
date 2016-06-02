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

#ifndef ROL_CVAR_HPP
#define ROL_CVAR_HPP

#include "ROL_RiskMeasure.hpp"
#include "ROL_PlusFunction.hpp"
#include "ROL_RiskVector.hpp"

/** @ingroup risk_group
    \class ROL::CVaR
    \brief Provides an interface for a convex combination of the
           expected value and the conditional value-at-risk.

    The conditional value-at-risk (also called the average value-at-risk
    or the expected shortfall) with confidence level \f$0\le \beta < 1\f$
    is
    \f[
       \mathcal{R}(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \frac{1}{1-\beta} \mathbb{E}\left[(X-t)_+\right]
         \right\}
    \f]
    where \f$(x)_+ = \max\{0,x\}\f$.  If the distribution of \f$X\f$ is
    continuous, then \f$\mathcal{R}\f$ is the conditional expectation of
    \f$X\f$ exceeding the \f$\beta\f$-quantile of \f$X\f$ and the optimal
    \f$t\f$ is the \f$\beta\f$-quantile.
    Additionally, \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$t\f$, then minimizes jointly for \f$(x_0,t)\f$.

    When using derivative-based optimization, the user can provide a smooth
    approximation of \f$(\cdot)_+\f$ using the ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class CVaR : public RiskMeasure<Real> {
private:
  Teuchos::RCP<PlusFunction<Real> > plusFunction_;

  Real prob_;
  Real coeff_;

  Teuchos::RCP<Vector<Real> > dualVector_;
  Real xvar_;
  Real vvar_;

  bool firstReset_;

  void checkInputs(void) const {
    Real zero(0), one(1);
    TEUCHOS_TEST_FOR_EXCEPTION((prob_ <= zero) || (prob_ >= one), std::invalid_argument,
      ">>> ERROR (ROL::CVaR): Confidence level must be between 0 and 1!");
    TEUCHOS_TEST_FOR_EXCEPTION((coeff_ < zero) || (coeff_ > one), std::invalid_argument,
      ">>> ERROR (ROL::CVaR): Convex combination parameter must be positive!");
    TEUCHOS_TEST_FOR_EXCEPTION(plusFunction_ == Teuchos::null, std::invalid_argument,
      ">>> ERROR (ROL::CVaR): PlusFunction pointer is null!");
  }

public:

  /** \brief Constructor.

      @param[in]     prob    is the confidence level
      @param[in]     coeff   is the convex combination parameter (coeff=0
                             corresponds to the expected value whereas coeff=1
                             corresponds to the conditional value-at-risk)
      @param[in]     pf      is the plus function or an approximation
  */
  CVaR( const Real prob, const Real coeff,
        const Teuchos::RCP<PlusFunction<Real> > &pf )
    : RiskMeasure<Real>(), plusFunction_(pf), prob_(prob), coeff_(coeff),
      xvar_(0), vvar_(0), firstReset_(true) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"CVaR" and
      within the "CVaR" sublist should have the following parameters
      \li "Confidence Level" (between 0 and 1)
      \li "Convex Combination Parameter" (between 0 and 1)
      \li A sublist for plus function information.
  */
  CVaR( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>(), xvar_(0), vvar_(0), firstReset_(true) {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("CVaR");
    // Check CVaR inputs
    prob_  = list.get<Real>("Confidence Level");
    coeff_ = list.get<Real>("Convex Combination Parameter");
    // Build (approximate) plus function
    plusFunction_ = Teuchos::rcp(new PlusFunction<Real>(list));
    // Check Inputs
    checkInputs();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    xvar_ = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic();
    if ( firstReset_ ) {
      dualVector_ = (x0->dual()).clone();
      firstReset_ = false;
    }
    dualVector_->zero();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    const RiskVector<Real> &vr = Teuchos::dyn_cast<const RiskVector<Real> >(v);
    v0    = Teuchos::rcp_const_cast<Vector<Real> >(vr.getVector());
    vvar_ = vr.getStatistic();
  }

  void update(const Real val, const Real weight) {
    Real one(1);
    Real pf = plusFunction_->evaluate(val-xvar_,0);
    RiskMeasure<Real>::val_ += weight*((one-coeff_)*val + coeff_/(one-prob_)*pf);
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    Real one(1);
    Real pf = plusFunction_->evaluate(val-xvar_,1);
    RiskMeasure<Real>::val_ += weight*pf;
    Real c  = (one-coeff_) + coeff_/(one-prob_)*pf;
    RiskMeasure<Real>::g_->axpy(weight*c,g);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    Real one(1);
    Real pf1 = plusFunction_->evaluate(val-xvar_,1);
    Real pf2 = plusFunction_->evaluate(val-xvar_,2);
    RiskMeasure<Real>::val_ += weight*pf2*(vvar_-gv);
    Real c  = pf2*coeff_/(one-prob_)*(gv-vvar_);
    RiskMeasure<Real>::hv_->axpy(weight*c,g);
    c = (one-coeff_) + coeff_/(one-prob_)*pf1;
    RiskMeasure<Real>::hv_->axpy(weight*c,hv);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val  = RiskMeasure<Real>::val_, cvar(0);
    sampler.sumAll(&val,&cvar,1);
    cvar += coeff_*xvar_;
    return cvar;
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    RiskVector<Real> &gs = Teuchos::dyn_cast<RiskVector<Real> >(g);
    Real val = RiskMeasure<Real>::val_, var(0), one(1);
    sampler.sumAll(&val,&var,1);
    
    sampler.sumAll(*(RiskMeasure<Real>::g_),*dualVector_);
    var *= -coeff_/(one-prob_);
    var += coeff_;
    gs.setStatistic(var);
    gs.setVector(*dualVector_); 
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    RiskVector<Real> &hs = Teuchos::dyn_cast<RiskVector<Real> >(hv);
    Real val = RiskMeasure<Real>::val_, var(0), one(1);
    sampler.sumAll(&val,&var,1);

    sampler.sumAll(*(RiskMeasure<Real>::hv_),*dualVector_);
    var *= coeff_/(one-prob_);
    hs.setStatistic(var);
    hs.setVector(*dualVector_);
  }
};

}

#endif
