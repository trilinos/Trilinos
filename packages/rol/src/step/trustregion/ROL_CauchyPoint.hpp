// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TRUSTREGIONFACTORY_H
#include "ROL_TrustRegionFactory.hpp"
#else

#ifndef ROL_CAUCHYPOINT_H
#define ROL_CAUCHYPOINT_H

/** \class ROL::CauchyPoint
    \brief Provides interface for the Cauchy point trust-region subproblem solver.
*/

#include "ROL_TrustRegion.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Types.hpp"
#include "ROL_ParameterList.hpp"

namespace ROL {

template<class Real>
class CauchyPoint : public TrustRegion<Real> {
private:

  ROL::Ptr<Vector<Real> > g_;
  ROL::Ptr<Vector<Real> > p_;
  ROL::Ptr<Vector<Real> > Hp_;

  Real pRed_;
  Real eps_;
  Real alpha_;

  bool useCGTCP_;

public:

  // Constructor
  CauchyPoint( ROL::ParameterList &parlist )
    : TrustRegion<Real>(parlist), pRed_(0), alpha_(-1), useCGTCP_(false) {
    // Unravel Parameter List
    Real oe2(100);
    Real TRsafe = parlist.sublist("Step").sublist("Trust Region").get("Safeguard Size",oe2);
    eps_        = TRsafe*ROL_EPSILON<Real>();
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g) {
    TrustRegion<Real>::initialize(x,s,g);
    Hp_ = g.clone();
    p_  = s.clone();
//    if ( useCGTCP_ ) {
//      g_ = g.clone();
//      p_ = s.clone();
//    }
  }

  void run( Vector<Real> &s,
            Real &snorm,
            int &iflag,
            int &iter,
            const Real del,
            TrustRegionModel<Real> &model) {
    //if ( pObj.isConActivated() ) {
    //  if ( useCGTCP_ ) {
    //    cauchypoint_CGT( s, snorm, iflag, iter, del, model );
    //  }
    //  else {
    //    cauchypoint_M( s, snorm, iflag, iter, del, model );
    //  }
    //}
    //else {
    //  cauchypoint_unc( s, snorm, iflag, iter, del, model );
    //}
    cauchypoint_unc( s, snorm, iflag, iter, del, model );
    TrustRegion<Real>::setPredictedReduction(pRed_);
  }

private:
  void cauchypoint_unc( Vector<Real> &s,
                        Real &snorm,
                        int &iflag,
                        int &iter,
                        const Real del,
                        TrustRegionModel<Real> &model) {
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    // Set step to (projected) gradient
    model.dualTransform(*Hp_,*model.getGradient());
    s.set(Hp_->dual());
    // Apply (reduced) Hessian to (projected) gradient
    model.hessVec(*Hp_,s,s,tol);
    Real gBg   = Hp_->dot(s.dual());
    Real gnorm = s.dual().norm();
    Real gg    = gnorm*gnorm;
    Real alpha = del/gnorm;
    if ( gBg > ROL_EPSILON<Real>() ) {
      alpha = std::min(gg/gBg, del/gnorm);
    }

    s.scale(-alpha);
    model.primalTransform(*p_,s);
    s.set(*p_);
    snorm = s.norm(); //alpha*gnorm;
    iflag = 0;
    iter  = 0;
    pRed_ = alpha*(gg - static_cast<Real>(0.5)*alpha*gBg);
  }

//  void cauchypoint_M( Vector<Real>           &s,
//                      Real                   &snorm,
//                      int                    &iflag,
//                      int                    &iter,
//                      const Real              del,
//                      const Vector<Real>     &x,
//                      TrustRegionModel<Real> &model,
//                      BoundConstraint<Real>  &bnd) {
//    Real tol = std::sqrt(ROL_EPSILON<Real>()),
//    const Real zero(0), half(0.5), oe4(1.e4), two(2);
//    // Parameters
//    Real mu0(1.e-2), mu1(1), beta1(0), beta2(0);
//    bool decr = true;
//    bool stat = true;
//    // Initial step length
//    Real alpha    = (alpha_ > zero ? alpha_ : one);
//    Real alpha0   = alpha;
//    Real alphamax = oe4*alpha;
//    // Set step to (projected) gradient
//    s.zero();
//    model.gradient(*Hp_,s,tol);
//    s.set(Hp_->dual());
//    // Initial model value
//    s.scale(-alpha);
//    bnd.computeProjectedStep(s,x);
//    snorm = s.norm();
//    Real val  = model.value(s,tol);
//    Real val0 = val;
//
//    // Determine whether to increase or decrease alpha
//    if ( val > mu0 * gs || snorm > mu1 * del ) { 
//      beta1 = half; 
//      beta2 = half; 
//      decr  = true;
//    }
//    else {
//      beta1 = two;
//      beta2 = two;
//      decr  = false;
//    }
//
//    while ( stat ) {
//      // Update step length
//      alpha0 = alpha;
//      val0   = val;
//      alpha *= half*(beta1+beta2);
//  
//      // Update model value
//      s.set(grad.dual());
//      s.scale(-alpha);
//      pObj.computeProjectedStep(s,x);
//      snorm = s.norm();
//      pObj.hessVec(*Hp_,s,x,tol);
//      gs    = s.dot(grad.dual());
//      val   = gs + half*s.dot(Hp_->dual());
//
//      // Update termination criterion
//      if ( decr ) {
//        stat = ( val > mu0 * gs || snorm > mu1 * del );
//        if ( std::abs(val) < eps_ && std::abs(mu0 *gs) < eps_ ) {
//          stat = (snorm > mu1 * del);
//        }
//      }
//      else {
//        stat = !( val > mu0 * gs || snorm > mu1 * del );
//        if ( std::abs(val) < eps_ && std::abs(mu0 *gs) < eps_ ) {
//          stat = !(snorm > mu1 * del);
//        }
//        if ( alpha > alphamax ) {
//          stat = false;
//        }
//      } 
//    }
//    // Reset to last 'successful' step
//    val   = val0;
//    alpha = alpha0;
//    s.set(grad.dual());
//    s.scale(-alpha);
//    pObj.computeProjectedStep(s,x);
//    snorm = s.norm();
//    
//    alpha_ = alpha;
//    pRed_  = -val;
//  }
//
//  void cauchypoint_CGT( Vector<Real> &s, Real &snorm, Real &del, int &iflag, int &iter, const Vector<Real> &x,
//                        const Vector<Real> &grad, const Real &gnorm, ProjectedObjective<Real> &pObj ) {
//    Real tol = std::sqrt(ROL_EPSILON<Real>()), one(1), half(0.5), two(2);
//    bool tmax_flag = true;
//    int maxit      = 20;
//    Real t         = del/gnorm;
//    Real tmax(1.e10), tmin(0), gs(0), pgnorm(0);
//    Real c1(0.25), c2(0.75), c3(0.9), c4(0.25);
//    for ( int i = 0; i < maxit; i++ ) {
//      // Compute p = x + s = P(x - t*g)
//      p_->set(x);
//      p_->axpy(-t,grad.dual()); 
//      pObj.project(*p_);
//      // Compute s = p - x = P(x - t*g) - x
//      s.set(*p_);
//      s.axpy(-one,x);
//      snorm = s.norm();
//      // Evaluate Model
//      pObj.hessVec(*Hp_,s,x,tol);
//      gs = s.dot(grad.dual());
//      pRed_ = -gs - half*s.dot(Hp_->dual());
//
//      // Check Stopping Conditions
//      g_->set(grad);
//      pObj.pruneActive(*g_,grad,*p_); // Project gradient onto tangent cone at p
//      pgnorm = g_->norm();
//      if ( snorm > del || pRed_ < -c2*gs ) {
//        tmax = t;
//        tmax_flag = false;
//      }
//      else if ( snorm < c3*del && pRed_ > -c1*gs && pgnorm > c4*std::abs(gs)/del ) {
//        tmin = t;
//      } 
//      else {
//        break;
//      }
//   
//      // Update t
//      if ( tmax_flag ) {
//        t *= two;
//      }
//      else {
//        t = half*(tmax + tmin);
//      }
//    }
//  }
};

}

#endif
#endif
