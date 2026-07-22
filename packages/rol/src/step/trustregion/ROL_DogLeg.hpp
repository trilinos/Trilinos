// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DOGLEG_H
#define ROL_DOGLEG_H

/** \class ROL::DogLeg
    \brief Provides interface for dog leg trust-region subproblem solver.
*/

#include "ROL_CauchyPoint.hpp"
#include "ROL_Types.hpp"

namespace ROL { 

template<class Real>
class DogLeg : public TrustRegion<Real> {
private:

  ROL::Ptr<CauchyPoint<Real> > cpt_;

  ROL::Ptr<Vector<Real> > s_;
  ROL::Ptr<Vector<Real> > Hp_;

  Real pRed_;

public:

  // Constructor
  DogLeg( ROL::ParameterList &parlist ) : TrustRegion<Real>(parlist), pRed_(0) {
    cpt_ = ROL::makePtr<CauchyPoint<Real>>(parlist);
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g) {
    TrustRegion<Real>::initialize(x,s,g);
    cpt_->initialize(x,s,g);
    s_  = s.clone();
    Hp_ = g.clone();
  }

  void run( Vector<Real>           &s,
            Real                   &snorm,
            int                    &iflag,
            int                    &iter,
            const Real              del,
            TrustRegionModel<Real> &model ) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    const Real zero(0), half(0.5), one(1), two(2);
    // Set s to be the (projected) gradient
    model.dualTransform(*Hp_,*model.getGradient());
    s.set(Hp_->dual());
    // Compute (quasi-)Newton step
    model.invHessVec(*s_,*Hp_,s,tol);
    Real sNnorm  = s_->norm();
    Real gsN     = -s_->dot(s);
    bool negCurv = (gsN > zero ? true : false);
    // Check if (quasi-)Newton step is feasible
    if ( negCurv ) {
      // Use Cauchy point
      cpt_->run(s,snorm,iflag,iter,del,model);
      pRed_ = cpt_->getPredictedReduction();
      iflag = 2;
    }
    else {
      // Approximately solve trust region subproblem using double dogleg curve
      if (sNnorm <= del) { // Use the (quasi-)Newton step
        s.set(*s_); 
        s.scale(-one);
        snorm = sNnorm;
        pRed_ = -half*gsN;
        iflag = 0;
      }
      else { // The (quasi-)Newton step is outside of trust region
        model.hessVec(*Hp_,s,s,tol);
        Real alpha  = zero;
        Real beta   = zero;
        Real gnorm  = s.norm();
        Real gnorm2 = gnorm*gnorm;
        Real gBg    = Hp_->dot(s.dual());
        Real gamma  = gnorm2/gBg;
        if ( gamma*gnorm >= del || gBg <= zero ) {
          // Use Cauchy point
          alpha = zero;
          beta  = del/gnorm;
          s.scale(-beta); 
          snorm = del;
          iflag = 2;
        }
        else {
          // Use a convex combination of Cauchy point and (quasi-)Newton step
          Real a = sNnorm*sNnorm + two*gamma*gsN + gamma*gamma*gnorm2;
          Real b = -gamma*gsN - gamma*gamma*gnorm2;
          Real c = gamma*gamma*gnorm2 - del*del;
          alpha  = (-b + sqrt(b*b - a*c))/a;
          beta   = gamma*(one-alpha);
          s.scale(-beta);
          s.axpy(-alpha,*s_);
          snorm = del;
          iflag = 1;
        }
        pRed_ = (alpha*(half*alpha-one)*gsN - half*beta*beta*gBg + beta*(one-alpha)*gnorm2);
      }
    }
    model.primalTransform(*s_,s);
    s.set(*s_);
    snorm = s.norm();
    // Update predicted reduction
    TrustRegion<Real>::setPredictedReduction(pRed_);
  }
};

}

#endif
