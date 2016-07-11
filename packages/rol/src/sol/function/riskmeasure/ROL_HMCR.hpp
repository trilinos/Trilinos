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

/** @ingroup risk_group
    \class ROL::HMCR
    \brief Provides an interface for a convex combination of the
           expected value and the higher moment coherent risk measure.

    The higher moment coherent risk measure of order \f$p\f$ with confidence
    level \f$0\le \beta < 1\f$ is
    \f[
       \mathcal{R}(X) = \inf_{t\in\mathbb{R}} \left\{
         t + \frac{1}{1-\beta} \mathbb{E}\left[(X-t)_+^p\right]^{1/p}
         \right\}
    \f]
    where \f$(x)_+ = \max\{0,x\}\f$.
    \f$\mathcal{R}\f$ is a law-invariant coherent risk measure.
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$t\f$, then minimizes jointly for \f$(x_0,t)\f$.

    The user can provide a smooth approximation of \f$(\cdot)_+\f$ using the
    ROL::PlusFunction class.
*/

namespace ROL {

template<class Real>
class HMCR : public RiskMeasure<Real> {
private:
  // Plus function (approximation)
  Teuchos::RCP<PlusFunction<Real> > plusFunction_;

  // User inputs
  Real prob_;
  Real lambda_;
  unsigned order_;

  // 1/(1-prob)
  Real coeff_;

  // Temporary vector storage
  Teuchos::RCP<Vector<Real> > mDualVector0_;
  Teuchos::RCP<Vector<Real> > gDualVector0_;
  Teuchos::RCP<Vector<Real> > mDualVector1_;
  Teuchos::RCP<Vector<Real> > gDualVector1_;

  // Statistic storage
  Real xvar_;
  Real vvar_;

  // Temporary scalar storage
  Real pnorm_;
  Real coeff0_;
  Real coeff1_;
  Real coeff2_;

  // Flag to initialized vector storage
  bool firstReset_;

  void checkInputs(void) const {
    const Real zero(0), one(1);
    TEUCHOS_TEST_FOR_EXCEPTION((prob_ <= zero) || (prob_ >= one), std::invalid_argument,
      ">>> ERROR (ROL::HMCR): Confidence level must be between 0 and 1!");
    TEUCHOS_TEST_FOR_EXCEPTION((lambda_ < zero) || (lambda_ > one), std::invalid_argument,
      ">>> ERROR (ROL::HMCR): Convex combination parameter must be positive!");
    TEUCHOS_TEST_FOR_EXCEPTION((order_ < 2), std::invalid_argument,
      ">>> ERROR (ROL::HMCR): Norm order is less than 2!");
    TEUCHOS_TEST_FOR_EXCEPTION(plusFunction_ == Teuchos::null, std::invalid_argument,
      ">>> ERROR (ROL::HMCR): PlusFunction pointer is null!");
  }

public:
  /** \brief Constructor.

      @param[in]     prob    is the confidence level
      @param[in]     lambda  is the convex combination parameter (lambda=0
                             corresponds to the expected value whereas lambda=1
                             corresponds to the higher moment coherent risk measure)
      @param[in]     order   is the order of higher moment coherent risk measure
      @param[in]     pf      is the plus function or an approximation
  */
  HMCR( const Real prob, const Real lambda, const unsigned order,
        const Teuchos::RCP<PlusFunction<Real> > &pf )
    : RiskMeasure<Real>(),
      plusFunction_(pf), prob_(prob), lambda_(lambda), order_(order),
      xvar_(0), vvar_(0), pnorm_(0), coeff0_(0), coeff1_(0), coeff2_(0),
      firstReset_(true) {
    checkInputs();
    const Real one(1);
    coeff_ = one/(one-prob_);
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"HMCR" and
      within the "HMCR" sublist should have the following parameters
      \li "Confidence Level" (between 0 and 1)
      \li "Convex Combination Parameter" (between 0 and 1)
      \li "Order" (unsigned integer)
      \li A sublist for plus function information.
  */
  HMCR( Teuchos::ParameterList &parlist )
    : RiskMeasure<Real>(),
      xvar_(0), vvar_(0), pnorm_(0), coeff0_(0), coeff1_(0), coeff2_(0),
      firstReset_(true) {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("HMCR");
    // Check HMCR inputs
    prob_  = list.get<Real>("Confidence Level");
    lambda_ = list.get<Real>("Convex Combination Parameter");
    order_ = (unsigned)list.get<int>("Order",2);
    // Build (approximate) plus function
    plusFunction_  = Teuchos::rcp(new PlusFunction<Real>(list));
    // Check inputs
    checkInputs();
    const Real one(1);
    coeff_ = one/(one-prob_);
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    RiskMeasure<Real>::reset(x0,x);
    xvar_ = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic();
    // Initialize additional vector storage
    if ( firstReset_ ) {
      mDualVector0_ = (x0->dual()).clone();
      gDualVector0_ = (x0->dual()).clone();
      mDualVector1_ = (x0->dual()).clone();
      gDualVector1_ = (x0->dual()).clone();
      firstReset_  = false;
    }
    // Zero temporary storage
    const Real zero(0);
    mDualVector0_->zero(); gDualVector0_->zero();
    pnorm_ = zero; coeff0_ = zero;
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    reset(x0,x);
    v0    = Teuchos::rcp_const_cast<Vector<Real> >(
            Teuchos::dyn_cast<const RiskVector<Real> >(v).getVector());
    vvar_ = Teuchos::dyn_cast<const RiskVector<Real> >(v).getStatistic();
    // Zero temporary storage
    const Real zero(0);
    mDualVector1_->zero(); gDualVector1_->zero();
    coeff1_ = zero; coeff2_ = zero;
  }

  void update(const Real val, const Real weight) {
    const Real rorder = static_cast<Real>(order_);
    // Expected value
    RiskMeasure<Real>::val_ += weight*val;
    // Higher moment
    Real pf = plusFunction_->evaluate(val-xvar_,0);
    pnorm_ += weight*std::pow(pf,rorder);
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    const Real one(1);
    const Real power = one/static_cast<Real>(order_);
    std::vector<Real> val_in(2), val_out(2);
    val_in[0] = RiskMeasure<Real>::val_;
    val_in[1] = pnorm_;
    sampler.sumAll(&val_in[0],&val_out[0],2);
    return (one-lambda_)*val_out[0]
          + lambda_*(xvar_ + coeff_*std::pow(val_out[1],power));
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    const Real one(1);
    const Real rorder0 = static_cast<Real>(order_);
    const Real rorder1 = rorder0 - one;
    // Expected value
    RiskMeasure<Real>::g_->axpy(weight,g);
    // Higher moment
    Real pf0 = plusFunction_->evaluate(val-xvar_,0);
    Real pf1 = plusFunction_->evaluate(val-xvar_,1);

    Real pf0p0 = std::pow(pf0,rorder0);
    Real pf0p1 = std::pow(pf0,rorder1);

    pnorm_  += weight*pf0p0;
    coeff0_ += weight*pf0p1*pf1;

    mDualVector0_->axpy(weight*pf0p1*pf1,g);
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    const Real zero(0), one(1);
    std::vector<Real> val_in(3), val_out(3);
    val_in[0] = RiskMeasure<Real>::val_;
    val_in[1] = pnorm_; val_in[2] = coeff0_;

    sampler.sumAll(&val_in[0],&val_out[0],3);
    sampler.sumAll(*(RiskMeasure<Real>::g_),*(RiskMeasure<Real>::dualVector_));
    RiskMeasure<Real>::dualVector_->scale(one-lambda_);
    Real var = lambda_;
    // If the higher moment term is positive then compute gradient
    if ( val_in[1] > zero ) {
      const Real rorder0 = static_cast<Real>(order_);
      const Real rorder1 = rorder0 - one;
      Real denom = std::pow(val_out[1],rorder1/rorder0);
      // Sum higher moment contribution
      sampler.sumAll(*mDualVector0_,*gDualVector0_);
      RiskMeasure<Real>::dualVector_->axpy(lambda_*coeff_/denom,*gDualVector0_);
      // Compute statistic gradient
      var -= lambda_*coeff_*((denom > zero) ? val_out[2]/denom : zero);
    }
    // Set gradients
    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setStatistic(var);
    (Teuchos::dyn_cast<RiskVector<Real> >(g)).setVector(*(RiskMeasure<Real>::dualVector_));
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    const Real one(1);
    const Real rorder0 = static_cast<Real>(order_);
    const Real rorder1 = rorder0-one;
    const Real rorder2 = rorder1-one;
    // Expected value
    RiskMeasure<Real>::hv_->axpy(weight,hv);
    // Higher moment
    Real pf0 = plusFunction_->evaluate(val-xvar_,0);
    Real pf1 = plusFunction_->evaluate(val-xvar_,1);
    Real pf2 = plusFunction_->evaluate(val-xvar_,2);

    Real pf0p0 = std::pow(pf0,rorder0);
    Real pf0p1 = std::pow(pf0,rorder1);
    Real pf0p2 = std::pow(pf0,rorder2);

    Real scale0 = (rorder1*pf0p2*pf1*pf1 + pf0p1*pf2)*(gv-vvar_);
    Real scale1 = pf0p1*pf1;

    pnorm_  += weight*pf0p0;
    coeff0_ += weight*scale0;
    coeff1_ += weight*scale1;
    coeff2_ += weight*rorder1*scale1*(vvar_-gv);

    mDualVector0_->axpy(weight*scale0,g);
    mDualVector0_->axpy(weight*scale1,hv);
    mDualVector1_->axpy(weight*scale1,g);
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    const Real zero(0), one(1);
    std::vector<Real> val_in(4), val_out(4);
    val_in[0] = pnorm_;  val_in[1] = coeff0_;
    val_in[2] = coeff1_; val_in[3] = coeff2_;

    sampler.sumAll(&val_in[0],&val_out[0],4);
    sampler.sumAll(*(RiskMeasure<Real>::hv_),*(RiskMeasure<Real>::dualVector_));

    Real var = zero;
    RiskMeasure<Real>::dualVector_->scale(one-lambda_);

    if ( val_out[0] > zero ) {
      const Real rorder0 = static_cast<Real>(order_);
      const Real rorder1 = rorder0-one;
      const Real rorder2 = rorder0 + rorder1;
      const Real coeff   = lambda_*coeff_;

      sampler.sumAll(*mDualVector0_,*gDualVector0_);
      sampler.sumAll(*mDualVector1_,*gDualVector1_);

      Real denom1 = std::pow(val_out[0],rorder1/rorder0);
      Real denom2 = std::pow(val_out[0],rorder2/rorder0);

      var = -coeff*(val_out[1]/denom1 + val_out[3]*val_out[2]/denom2);
      RiskMeasure<Real>::dualVector_->axpy(coeff/denom1,*gDualVector0_);
      RiskMeasure<Real>::dualVector_->axpy(coeff*val_out[3]/denom2,*gDualVector1_);
    }

    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setStatistic(var);
    (Teuchos::dyn_cast<RiskVector<Real> >(hv)).setVector(*(RiskMeasure<Real>::dualVector_));
  }
};

}

#endif
