// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_FLETCHEROBJECTIVEEDEF_H
#define ROL_FLETCHEROBJECTIVEEDEF_H

namespace ROL {

template<typename Real>
FletcherObjectiveE<Real>::FletcherObjectiveE(const ROL::Ptr<Objective<Real>> &obj,
                                             const ROL::Ptr<Constraint<Real>> &con,
                                             const Vector<Real> &xprim,
                                             const Vector<Real> &xdual,
                                             const Vector<Real> &cprim,
                                             const Vector<Real> &cdual,
                                             ROL::ParameterList &parlist)
  : FletcherObjectiveBase<Real>(obj, con, xprim, xdual, cprim, cdual, parlist) {
  Tv_    = xdual.clone();
  w_     = xdual.clone();
  wdual_ = xprim.clone();
  v_     = cdual.clone();
  wg_    = xdual.clone();
  vg_    = cdual.clone();

  xzeros_ = xdual.clone(); xzeros_->zero();
  czeros_ = cprim.clone(); czeros_->zero();
}

template<typename Real>
Real FletcherObjectiveE<Real>::value( const Vector<Real> &x, Real &tol ) {
  Real val(0);
  int key(0);
  bool isComputed = fPhi_->get(val,key);
  if( isComputed && multSolverError_*cnorm_ <= tol) {
    tol = multSolverError_*cnorm_;
  }
  else {
    // Reset tolerances
    Real origTol = tol;
    Real tol2 = origTol;
    // Compute penalty function value 
    Real fval = FletcherObjectiveBase<Real>::objValue(x, tol2); tol2 = origTol;
    multSolverError_ = origTol / (static_cast<Real>(2) * std::max(static_cast<Real>(1), cnorm_));
    FletcherObjectiveBase<Real>::computeMultipliers(*cdual_, *gLdual_, x, *xdual_, *cprim_, multSolverError_);
    tol =  multSolverError_*cnorm_;
    //val = fval - cprim_->dot(cdual_->dual());
    val = fval - cprim_->apply(*cdual_);
    if( quadPenaltyParameter_ > static_cast<Real>(0) )
      val += static_cast<Real>(0.5)*quadPenaltyParameter_*(cprim_->dot(*cprim_));
    // Store new penalty function value
    fPhi_->set(val,key);
  }
  return val;
}

template<typename Real>
void FletcherObjectiveE<Real>::gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
  int key(0);
  bool isComputed = gPhi_->get(g,key);
  if( isComputed && gradSolveError_ <= tol) {
    tol = gradSolveError_;
  }
  else {
    // Reset tolerances
    Real origTol = tol;
    Real tol2 = origTol;
    // Compute penalty function gradient 
    gradSolveError_ = origTol / static_cast<Real>(2);
    FletcherObjectiveBase<Real>::computeMultipliers(*cdual_, *gLdual_, x, *xdual_, *cprim_, gradSolveError_);
    gL_->set(gLdual_->dual());
    bool refine = isComputed;
    // gPhi = sum y_i H_i w + sigma w + sum v_i H_i gL - H w + gL
    solveAugmentedSystem( *wdual_, *vg_, *xzeros_, *cprim_, x, gradSolveError_, refine );
    gradSolveError_ += multSolverError_;
    tol = gradSolveError_;
    wg_->set(wdual_->dual());
    con_->applyAdjointHessian( g, *cdual_, *wdual_, x, tol2 ); tol2 = origTol;
    g.axpy( sigma_, *wg_ );
    obj_->hessVec( *Tv_, *wdual_, x, tol2 ); tol2 = origTol;
    g.axpy( static_cast<Real>(-1), *Tv_ );
    con_->applyAdjointHessian( *Tv_, *vg_, *gLdual_, x, tol2 ); tol2 = origTol;
    g.plus( *Tv_ );
    g.plus( *gL_ );
    if( quadPenaltyParameter_ > static_cast<Real>(0) ) {
      con_->applyAdjointJacobian( *Tv_, *cprim_, x, tol2 ); tol2 = origTol;
      g.axpy( quadPenaltyParameter_, *Tv_ );
    }
    gPhi_->set(g,key); 
  }
}

template<typename Real>
void FletcherObjectiveE<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
  // Reset tolerances
  Real origTol = tol;
  Real tol2 = origTol;
  int key(0);
  bool isComputed = y_->get(*cdual_,key);
  if( !isComputed || !useInexact_) {
    // hessVec tol is always set to ~1e-8. So if inexact linear system solves are used, then
    // the multipliers will always get re-evaluated to high precision, which defeats the purpose
    // of computing the objective/gradient approximately.
    // Thus if inexact linear system solves are used, we will not update the multiplier estimates
    // to high precision.
    FletcherObjectiveBase<Real>::computeMultipliers(*cdual_, *gLdual_, x, *xdual_, *cprim_, tol);
  }

  obj_->hessVec( hv, v, x, tol2 ); tol2 = origTol;
  con_->applyAdjointHessian( *Tv_, *cdual_, v, x, tol2 ); tol2 = origTol;
  hv.axpy(static_cast<Real>(-1), *Tv_ );

  tol2 = tol;
  solveAugmentedSystem( *w_, *v_, hv, *czeros_, x, tol2 ); tol2 = origTol;
  hv.scale( static_cast<Real>(-1) );
  hv.plus( *w_ );

  Tv_->set(v.dual());
  tol2 = tol;
  solveAugmentedSystem( *w_, *v_, *Tv_, *czeros_, x, tol2 ); tol2 = origTol;
  hv.axpy(static_cast<Real>(-2)*sigma_, *w_);

  wdual_->set(w_->dual());

  obj_->hessVec( *Tv_, *wdual_, x, tol2 ); tol2 = origTol;
  hv.plus( *Tv_ );
  con_->applyAdjointHessian( *Tv_, *cdual_, *wdual_, x, tol2 ); tol2 = origTol;
  hv.axpy( static_cast<Real>(-1), *Tv_ );

  hv.axpy( static_cast<Real>(2)*sigma_, v );

  if( quadPenaltyParameter_ > static_cast<Real>(0) ) {
    con_->applyJacobian( *b2_, v, x, tol2 ); tol2 = origTol;
    con_->applyAdjointJacobian( *Tv_, *b2_, x, tol2 ); tol2 = origTol;
    hv.axpy( quadPenaltyParameter_, *Tv_ );
    con_->applyAdjointHessian( *Tv_, *cprim_, v, x, tol2); tol2 = origTol;
    hv.axpy( -quadPenaltyParameter_, *Tv_ );
  }
}

template<typename Real>
void FletcherObjectiveE<Real>::solveAugmentedSystem(Vector<Real> &v1,
                                                    Vector<Real> &v2,
                                                    const Vector<Real> &b1,
                                                    const Vector<Real> &b2,
                                                    const Vector<Real> &x,
                                                    Real &tol,
                                                    bool refine) {
  // Ignore tol for now
  ROL::Ptr<LinearOperator<Real>>
    K = ROL::makePtr<AugSystem>(con_, makePtrFromRef(x), delta_);
  ROL::Ptr<LinearOperator<Real>>
    P = ROL::makePtr<AugSystemPrecond>(con_, makePtrFromRef(x), makePtrFromRef(b1));

  vv_->zero();
  bb_->zero();
  if( refine ) {
    // TODO: Make sure this tol is actually ok...
    Real origTol = tol;
    w1_->set(v1);
    w2_->set(v2);
    K->apply(*bb_, *ww_, tol); tol = origTol;
    bb_->scale(static_cast<Real>(-1));
  }
  b1_->plus(b1);
  b2_->plus(b2);

  // If inexact, change tolerance
  if( useInexact_ ) krylov_->resetAbsoluteTolerance(tol);

  //con_->solveAugmentedSystem(*v1_,*v2_,*b1_,*b2_,x,tol);
  flagKrylov_ = 0;
  tol = krylov_->run(*vv_,*K,*bb_,*P,iterKrylov_,flagKrylov_);

  if( refine ) {
    v1.plus(*v1_);
    v2.plus(*v2_);
  } else {
    v1.set(*v1_);
    v2.set(*v2_);
  }
}

} // namespace ROL

#endif
