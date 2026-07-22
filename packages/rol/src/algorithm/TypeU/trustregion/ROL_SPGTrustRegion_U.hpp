// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SPGTRUSTREGION_U_H
#define ROL_SPGTRUSTREGION_U_H

/** \class ROL::SPGTrustRegion_U
    \brief Provides interface for truncated CG trust-region subproblem solver.
*/

#include "ROL_TrustRegion_U.hpp"
#include "ROL_Types.hpp"

#include <deque>

namespace ROL { 

template<typename Real>
class SPGTrustRegion_U : public TrustRegion_U<Real> {
private:
  Ptr<Vector<Real>> dwa_, pwa_, pwa1_, gmod_, smin_;

  Real lambdaMin_;
  Real lambdaMax_;
  Real gamma_;
  int maxSize_;
  int maxit_;
  Real tol1_;
  Real tol2_;
  bool useMin_;
  bool useNMSP_;

public:

  // Constructor
  SPGTrustRegion_U( ParameterList &parlist ) {
    ParameterList &list = parlist.sublist("Step").sublist("Trust Region").sublist("SPG");
    // Spectral projected gradient parameters
    lambdaMin_ = list.sublist("Solver").get("Minimum Spectral Step Size",          1e-8);
    lambdaMax_ = list.sublist("Solver").get("Maximum Spectral Step Size",          1e8);
    gamma_     = list.sublist("Solver").get("Sufficient Decrease Tolerance",       1e-4);
    maxSize_   = list.sublist("Solver").get("Maximum Storage Size",                10);
    maxit_     = list.sublist("Solver").get("Iteration Limit",                     25);
    tol1_      = list.sublist("Solver").get("Absolute Tolerance",                  1e-4);
    tol2_      = list.sublist("Solver").get("Relative Tolerance",                  1e-2);
    useMin_    = list.sublist("Solver").get("Use Smallest Model Iterate",          true);
    useNMSP_   = list.sublist("Solver").get("Use Nonmonotone Search",              false);
  }

  void initialize(const Vector<Real> &x, const Vector<Real> &g) {
    pwa_  = x.clone();
    pwa1_ = x.clone();
    smin_ = x.clone();
    dwa_  = g.clone();
    gmod_ = g.clone();
  }

  void solve( Vector<Real>             &s,
              Real                     &snorm,
              Real                     &pRed,
              int                      &iflag,
              int                      &iter,
              const Real                del,
              TrustRegionModel_U<Real> &model ) {
    const Real zero(0), half(0.5), one(1), two(2), eps(std::sqrt(ROL_EPSILON<Real>()));
    Real tol(eps), alpha(1), sHs(0), alphaTmp(1), mmax(0), qmin(0), q(0);
    Real gnorm(0), ss(0), gs(0);
    std::deque<Real> mqueue; mqueue.push_back(0);
    gmod_->set(*model.getGradient());

    // Compute Cauchy point
    pwa1_->set(gmod_->dual());
    s.set(*pwa1_); s.scale(-one);
    model.hessVec(*dwa_,s,s,tol);
    gs = gmod_->apply(s);
    sHs = dwa_->apply(s);
    snorm = std::sqrt(std::abs(gs));
    alpha = -gs/sHs;
    if (alpha*snorm >= del || sHs <= zero) alpha = del/snorm;
    q = alpha*(gs+half*alpha*sHs);
    gmod_->axpy(alpha,*dwa_);
    s.scale(alpha);

    if (useNMSP_ && useMin_) { qmin = q; smin_->set(s);}

    // Compute initial projected gradient
    pwa1_->set(gmod_->dual());
    pwa_->set(s); pwa_->axpy(-one,*pwa1_);
    snorm = pwa_->norm();
    if (snorm > del) pwa_->scale(del/snorm);
    pwa_->axpy(-one,s);
    gnorm = pwa_->norm();
    if (gnorm == zero) {
      snorm = s.norm();
      pRed  = -q;
      return;
    }
    const Real gtol = std::min(tol1_,tol2_*gnorm);

    // Compute initial step
    Real lambda = std::max(lambdaMin_,std::min(one/gmod_->norm(),lambdaMax_));
    pwa_->set(s); pwa_->axpy(-lambda,*pwa1_);
    snorm = pwa_->norm();
    if (snorm > del) pwa_->scale(del/snorm);
    pwa_->axpy(-one,s);
    gs = gmod_->apply(*pwa_);
    ss = pwa_->dot(*pwa_);

    for (iter = 0; iter < maxit_; iter++) {
      // Evaluate model Hessian
      model.hessVec(*dwa_,*pwa_,s,tol);
      sHs = dwa_->apply(*pwa_);
      // Perform line search
      if (useNMSP_) { // Nonmonotone
        mmax     = *std::max_element(mqueue.begin(),mqueue.end());
        alphaTmp = (-(one-gamma_)*gs + std::sqrt(std::pow((one-gamma_)*gs,two)-two*sHs*(q-mmax)))/sHs;
      }
      else { // Exact
        alphaTmp = -gs/sHs;
      }
      alpha = (sHs > zero ? std::min(one,std::max(zero,alphaTmp)) : one);
      // Update model quantities
      q += alpha*(gs+half*alpha*sHs);
      gmod_->axpy(alpha,*dwa_);
      s.axpy(alpha,*pwa_);
      // Update nonmonotone line search information
      if (useNMSP_) {
        if (static_cast<int>(mqueue.size())==maxSize_) mqueue.pop_front();
        mqueue.push_back(q);
        if (useMin_ && q <= qmin) { qmin = q; smin_->set(s); }
      }
      // Compute Projected gradient norm
      pwa1_->set(gmod_->dual());
      pwa_->set(s); pwa_->axpy(-one,*pwa1_);
      snorm = pwa_->norm();
      if (snorm > del) pwa_->scale(del/snorm);
      pwa_->axpy(-one,s);
      gnorm = pwa_->norm();
      if (gnorm < gtol) break;
      // Compute new spectral step
      lambda = (sHs <= eps ? lambdaMax_ : std::max(lambdaMin_,std::min(ss/sHs,lambdaMax_)));
      pwa_->set(s); pwa_->axpy(-lambda,*pwa1_);
      snorm = pwa_->norm();
      if (snorm > del) pwa_->scale(del/snorm);
      pwa_->axpy(-one,s);
      gs = gmod_->apply(*pwa_);
      ss = pwa_->dot(*pwa_);
    }
    if (useNMSP_ && useMin_) { q = qmin; s.set(*smin_); }
    iflag = (iter==maxit_ ? 1 : 0);
    pRed = -q;
    snorm = s.norm();
  }
};

} // namespace ROL

#endif
