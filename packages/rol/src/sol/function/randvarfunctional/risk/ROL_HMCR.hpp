// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_HMCR_HPP
#define ROL_HMCR_HPP

#include "ROL_RandVarFunctional.hpp"
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
class HMCR : public RandVarFunctional<Real> {
private:
  // User inputs
  ROL::Ptr<PlusFunction<Real> > plusFunction_;
  Real prob_;
  Real lambda_;
  unsigned order_;

  // 1/(1-prob)
  Real coeff_;

  // Temporary vector storage
  ROL::Ptr<Vector<Real> > mDualVector_;

  // Temporary scalar storage
  Real pnorm_;
  Real coeff0_;
  Real coeff1_;
  Real coeff2_;

  // Flag to initialized vector storage
  bool HMCR_firstReset_;

  using RandVarFunctional<Real>::val_;
  using RandVarFunctional<Real>::gv_;
  using RandVarFunctional<Real>::g_;
  using RandVarFunctional<Real>::hv_;
  using RandVarFunctional<Real>::dualVector_;

  using RandVarFunctional<Real>::point_;
  using RandVarFunctional<Real>::weight_;

  using RandVarFunctional<Real>::computeValue;
  using RandVarFunctional<Real>::computeGradient;
  using RandVarFunctional<Real>::computeGradVec;
  using RandVarFunctional<Real>::computeHessVec;

  void checkInputs(void) const {
    const Real zero(0), one(1);
    ROL_TEST_FOR_EXCEPTION((prob_ <= zero) || (prob_ >= one), std::invalid_argument,
      ">>> ERROR (ROL::HMCR): Confidence level must be between 0 and 1!");
    ROL_TEST_FOR_EXCEPTION((lambda_ < zero) || (lambda_ > one), std::invalid_argument,
      ">>> ERROR (ROL::HMCR): Convex combination parameter must be positive!");
    ROL_TEST_FOR_EXCEPTION((order_ < 2), std::invalid_argument,
      ">>> ERROR (ROL::HMCR): Norm order is less than 2!");
    ROL_TEST_FOR_EXCEPTION(plusFunction_ == ROL::nullPtr, std::invalid_argument,
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
        const ROL::Ptr<PlusFunction<Real> > &pf )
    : RandVarFunctional<Real>(),
      plusFunction_(pf), prob_(prob), lambda_(lambda), order_(order),
      pnorm_(0), coeff0_(0), coeff1_(0), coeff2_(0),
      HMCR_firstReset_(true) {
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
  HMCR( ROL::ParameterList &parlist )
    : RandVarFunctional<Real>(),
      pnorm_(0), coeff0_(0), coeff1_(0), coeff2_(0),
      HMCR_firstReset_(true) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("HMCR");
    // Check HMCR inputs
    prob_  = list.get<Real>("Confidence Level");
    lambda_ = list.get<Real>("Convex Combination Parameter");
    order_ = (unsigned)list.get<int>("Order",2);
    // Build (approximate) plus function
    plusFunction_ = ROL::makePtr<PlusFunction<Real> >(list);
    // Check inputs
    checkInputs();
    const Real one(1);
    coeff_ = one/(one-prob_);
  }

  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    // Initialize additional vector storage
    if ( HMCR_firstReset_ ) {
      mDualVector_ = x.dual().clone();
      HMCR_firstReset_  = false;
    }
    // Zero temporary storage
    const Real zero(0);
    mDualVector_->zero();
    pnorm_  = zero; coeff0_ = zero;
    coeff1_ = zero; coeff2_ = zero;
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    const Real rorder = static_cast<Real>(order_);
    // Expected value
    Real val = computeValue(obj,x,tol);
    val_    += weight_*val;
    // Higher moment
    Real pf  = plusFunction_->evaluate(val-xstat[0],0);
    pnorm_  += weight_*std::pow(pf,rorder);
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    const Real one(1);
    const Real power = one/static_cast<Real>(order_);
    std::vector<Real> val_in(2), val_out(2);
    val_in[0] = val_;
    val_in[1] = pnorm_;
    sampler.sumAll(&val_in[0],&val_out[0],2);
    return (one-lambda_)*val_out[0]
          + lambda_*(xstat[0] + coeff_*std::pow(val_out[1],power));
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    const Real one(1);
    const Real rorder0 = static_cast<Real>(order_);
    const Real rorder1 = rorder0 - one;
    // Expected value
    computeGradient(*dualVector_,obj,x,tol);
    g_->axpy(weight_,*dualVector_);
    // Higher moment
    Real val = computeValue(obj,x,tol);
    Real pf0 = plusFunction_->evaluate(val-xstat[0],0);
    Real pf1 = plusFunction_->evaluate(val-xstat[0],1);

    Real pf0p0 = std::pow(pf0,rorder0);
    Real pf0p1 = std::pow(pf0,rorder1);

    pnorm_  += weight_*pf0p0;
    coeff0_ += weight_*pf0p1*pf1;

    hv_->axpy(weight_*pf0p1*pf1,*dualVector_);
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    const Real zero(0), one(1);
    std::vector<Real> val_in(2), val_out(2);
    val_in[0] = pnorm_; val_in[1] = coeff0_;

    sampler.sumAll(&val_in[0],&val_out[0],2);
    sampler.sumAll(*g_,g);
    g.scale(one-lambda_);
    Real var = lambda_;
    // If the higher moment term is positive then compute gradient
    if ( val_out[0] > zero ) {
      const Real rorder0 = static_cast<Real>(order_);
      const Real rorder1 = rorder0 - one;
      Real denom = std::pow(val_out[0],rorder1/rorder0);
      // Sum higher moment contribution
      dualVector_->zero();
      sampler.sumAll(*hv_,*dualVector_);
      g.axpy(lambda_*coeff_/denom,*dualVector_);
      // Compute statistic gradient
      var -= lambda_*coeff_*((denom > zero) ? val_out[1]/denom : zero);
    }
    gstat[0] = var;
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    const Real one(1);
    const Real rorder0 = static_cast<Real>(order_);
    const Real rorder1 = rorder0-one;
    const Real rorder2 = rorder1-one;
    // Expected value
    computeHessVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_,*dualVector_);
    // Higher moment
    Real val = computeValue(obj,x,tol);
    Real pf0 = plusFunction_->evaluate(val-xstat[0],0);
    Real pf1 = plusFunction_->evaluate(val-xstat[0],1);
    Real pf2 = plusFunction_->evaluate(val-xstat[0],2);

    Real pf0p0 = std::pow(pf0,rorder0);
    Real pf0p1 = std::pow(pf0,rorder1);
    Real pf0p2 = std::pow(pf0,rorder2);

    Real scale1 = pf0p1*pf1;
    g_->axpy(weight_*scale1,*dualVector_);

    Real gv     = computeGradVec(*dualVector_,obj,v,x,tol);
    Real scale0 = (rorder1*pf0p2*pf1*pf1 + pf0p1*pf2)*(gv-vstat[0]);

    pnorm_  += weight_*pf0p0;
    coeff0_ += weight_*scale0;
    coeff1_ += weight_*scale1;
    coeff2_ += weight_*rorder1*scale1*(vstat[0]-gv);

    g_->axpy(weight_*scale0,*dualVector_);
    mDualVector_->axpy(weight_*scale1,*dualVector_);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    const Real zero(0), one(1);
    std::vector<Real> val_in(4), val_out(4);
    val_in[0] = pnorm_;  val_in[1] = coeff0_;
    val_in[2] = coeff1_; val_in[3] = coeff2_;

    sampler.sumAll(&val_in[0],&val_out[0],4);
    sampler.sumAll(*hv_,hv);

    Real var = zero;
    hv.scale(one-lambda_);

    if ( val_out[0] > zero ) {
      const Real rorder0 = static_cast<Real>(order_);
      const Real rorder1 = rorder0-one;
      const Real rorder2 = rorder0 + rorder1;
      const Real coeff   = lambda_*coeff_;

      Real denom1 = std::pow(val_out[0],rorder1/rorder0);
      Real denom2 = std::pow(val_out[0],rorder2/rorder0);

      dualVector_->zero();
      sampler.sumAll(*g_,*dualVector_);
      hv.axpy(coeff/denom1,*dualVector_);

      dualVector_->zero();
      sampler.sumAll(*mDualVector_,*dualVector_);
      hv.axpy(coeff*val_out[3]/denom2,*dualVector_);

      var = -coeff*(val_out[1]/denom1 + val_out[3]*val_out[2]/denom2);
    }
    hvstat[0] = var;
  }
};

}

#endif
