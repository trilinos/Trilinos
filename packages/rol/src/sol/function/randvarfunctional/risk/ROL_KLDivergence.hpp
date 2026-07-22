// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_KLDIVERGENCE_HPP
#define ROL_KLDIVERGENCE_HPP

#include "ROL_RandVarFunctional.hpp"

/** @ingroup risk_group
    \class ROL::KLDivergence
    \brief Provides an interface for the Kullback-Leibler distributionally robust
    expectation.

    This class defines a risk measure \f$\mathcal{R}\f$ which arises in distributionally
    robust stochastic programming.  \f$\mathcal{R}\f$ is given by
    \f[
       \mathcal{R}(X) = \sup_{\vartheta\in\mathfrak{A}}
           \mathbb{E}[\vartheta X]
    \f]
    where \f$\mathfrak{A}\f$ is called the ambiguity (or uncertainty) set and 
    is defined by a constraint on the Kullback-Leibler divergence, i.e.,
    \f[
       \mathfrak{A} = \{\vartheta\in\mathcal{X}^*\,:\,
         \mathbb{E}[\vartheta] = 1,\; \vartheta \ge 0,\;\text{and}\;
         \mathbb{E}[\vartheta\log(\vartheta)] \le \epsilon\}.
    \f]
    \f$\mathcal{R}\f$ is a law-invariant, coherent risk measure.  Moreover, by a
    duality argument, \f$\mathcal{R}\f$ can be reformulated as
    \f[
       \mathcal{R}(X) = \inf_{\lambda > 0}\left\{
             \lambda \epsilon + \lambda\mathbb{E}\left[\exp\left(
                \frac{X}{\lambda}\right)\right]\right\}.
    \f]
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$\lambda\f$, then minimizes jointly for \f$(x_0,\lambda)\f$.
*/

namespace ROL {

template<class Real>
class KLDivergence : public RandVarFunctional<Real> {
private:
  Real eps_;

  Real gval_;
  Real gvval_;
  Real hval_;
  ROL::Ptr<Vector<Real> > scaledGradient_;
  ROL::Ptr<Vector<Real> > scaledHessVec_;

  bool firstResetKLD_;

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
    Real zero(0);
    ROL_TEST_FOR_EXCEPTION((eps_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::KLDivergence): Threshold must be positive!");
  }

public:
  /** \brief Constructor.

      @param[in]     eps    is the tolerance for the KL divergence constraint
  */
  KLDivergence(const Real eps = 1.e-2)
    : RandVarFunctional<Real>(), eps_(eps), firstResetKLD_(true) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"KL Divergence" and
      within the "KL Divergence" sublist should have the following parameters
      \li "Threshold" (greater than 0)
  */
  KLDivergence(ROL::ParameterList &parlist)
    : RandVarFunctional<Real>(), firstResetKLD_(true) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("KL Divergence");
    eps_ = list.get<Real>("Threshold");
    checkInputs();
  }

  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    if ( firstResetKLD_ ) {
      scaledGradient_ = x.dual().clone();
      scaledHessVec_  = x.dual().clone();
      firstResetKLD_ = false;
    }
    const Real zero(0);
    gval_ = zero; gvval_ = zero; hval_ = zero;
    scaledGradient_->zero(); scaledHessVec_->zero();
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real ev  = exponential(val,xstat[0]*eps_);
    val_    += weight_ * ev;
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    if ( xstat[0] == static_cast<Real>(0) ) {
      return ROL_INF<Real>();
    }
    Real ev(0);
    sampler.sumAll(&val_,&ev,1);
    return (static_cast<Real>(1) + std::log(ev)/eps_)/xstat[0];
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real ev  = exponential(val,xstat[0]*eps_);
    val_    += weight_ * ev;
    gval_   += weight_ * ev * val;
    computeGradient(*dualVector_,obj,x,tol);
    g_->axpy(weight_*ev,*dualVector_);
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    std::vector<Real> local(2), global(2);
    local[0] = val_;
    local[1] = gval_;
    sampler.sumAll(&local[0],&global[0],2);
    Real ev = global[0], egval = global[1];

    sampler.sumAll(*g_,g);
    g.scale(static_cast<Real>(1)/ev);

    if ( xstat[0] == static_cast<Real>(0) ) {
      gstat[0] = ROL_INF<Real>();
    }
    else {
      gstat[0] = -((static_cast<Real>(1) + std::log(ev)/eps_)/xstat[0]
                 - egval/ev)/xstat[0];
    }
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    Real val = computeValue(obj,x,tol);
    Real ev  = exponential(val,xstat[0]*eps_);
    Real gv  = computeGradVec(*dualVector_,obj,v,x,tol);
    val_    += weight_ * ev;
    gv_     += weight_ * ev * gv;
    gval_   += weight_ * ev * val;
    gvval_  += weight_ * ev * val * gv;
    hval_   += weight_ * ev * val * val;
    g_->axpy(weight_*ev,*dualVector_);
    scaledGradient_->axpy(weight_*ev*gv,*dualVector_);
    scaledHessVec_->axpy(weight_*ev*val,*dualVector_);
    computeHessVec(*dualVector_,obj,v,x,tol);
    hv_->axpy(weight_*ev,*dualVector_);
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    std::vector<Real> local(5), global(5);
    local[0] = val_;
    local[1] = gv_;
    local[2] = gval_;
    local[3] = gvval_;
    local[4] = hval_;
    sampler.sumAll(&local[0],&global[0],5);
    Real ev     = global[0], egv   = global[1], egval = global[2];
    Real egvval = global[3], ehval = global[4];
    Real c0 = static_cast<Real>(1)/ev, c1 = c0*egval, c2 = c0*egv, c3 = eps_*c0;

    sampler.sumAll(*hv_,hv);
    dualVector_->zero();
    sampler.sumAll(*scaledGradient_,*dualVector_);
    hv.axpy(xstat[0]*eps_,*dualVector_);
    hv.scale(c0);

    dualVector_->zero();
    sampler.sumAll(*g_,*dualVector_);
    hv.axpy(-c3*(vstat[0]*c1 + xstat[0]*c2),*dualVector_);

    dualVector_->zero();
    sampler.sumAll(*scaledHessVec_,*dualVector_);
    hv.axpy(vstat[0]*c3,*dualVector_);

    if ( xstat[0] == static_cast<Real>(0) ) {
      hvstat[0] = ROL_INF<Real>();
    }
    else {
      Real xstat2 = static_cast<Real>(2)/(xstat[0]*xstat[0]);
      Real h11 = xstat2*((static_cast<Real>(1) + std::log(ev)/eps_)/xstat[0] - c1)
                 + (c3*ehval - eps_*c1*c1)/xstat[0];
      hvstat[0] = vstat[0] * h11 + (c3*egvval - eps_*c1*c2);
    }
  }

private:
  Real exponential(const Real arg1, const Real arg2) const {
    if ( arg1 < arg2 ) {
      return power(exponential(arg1),arg2);
    }
    else {
      return power(exponential(arg2),arg1);
    }
  }

  Real exponential(const Real arg) const {
    if ( arg >= std::log(ROL_INF<Real>()) ) {
      return ROL_INF<Real>();
    }
    else {
      return std::exp(arg);
    }
  }

  Real power(const Real arg, const Real pow) const {
    if ( arg >= std::pow(ROL_INF<Real>(),static_cast<Real>(1)/pow) ) {
      return ROL_INF<Real>();
    }
    else {
      return std::pow(arg,pow);
    }
  }
};

}

#endif
