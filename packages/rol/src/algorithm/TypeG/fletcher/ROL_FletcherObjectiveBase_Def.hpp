// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_FLETCHEROBJECTVEBASEDEF_H
#define ROL_FLETCHEROBJECTVEBASEDEF_H

namespace ROL {

template<typename Real>
FletcherObjectiveBase<Real>::FletcherObjectiveBase(const Ptr<Objective<Real>> &obj,
                                                   const Ptr<Constraint<Real>> &con,
                                                   const Vector<Real> &xprim,
                                                   const Vector<Real> &xdual,
                                                   const Vector<Real> &cprim,
                                                   const Vector<Real> &cdual,
                                                   ParameterList &parlist)
  : obj_(obj), con_(con), nfval_(0), ngval_(0), ncval_(0),
    fPhi_(makePtr<ScalarController<Real,int>>()),
    gPhi_(makePtr<VectorController<Real,int>>()),
    y_   (makePtr<VectorController<Real,int>>()),
    fval_(makePtr<ScalarController<Real,int>>()),
    g_   (makePtr<VectorController<Real,int>>()),
    c_   (makePtr<VectorController<Real,int>>()),
    multSolverError_(0), gradSolveError_(0),
    iterKrylov_(0), flagKrylov_(0) {
  gL_      = xdual.clone();
  gLdual_  = xprim.clone();
  scaledc_ = cprim.clone();
  xprim_   = xprim.clone();
  xdual_   = xdual.clone();
  cprim_   = cprim.clone();
  cdual_   = cdual.clone();

  v1_ = xprim.clone();
  v2_ = cdual.clone();
  vv_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>>>({v1_, v2_}));

  w1_ = xprim.clone();
  w2_ = cdual.clone();
  ww_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>>>({w1_, w2_}));

  b1_ = xdual.clone();
  b2_ = cprim.clone();
  bb_ = makePtr<PartitionedVector<Real>>(std::vector<Ptr<Vector<Real>>>({b1_, b2_}));

  ParameterList& sublist = parlist.sublist("Step").sublist("Fletcher");
  HessianApprox_         = sublist.get("Level of Hessian Approximation",           0);
  quadPenaltyParameter_  = sublist.get("Quadratic Penalty Parameter",              Real(0));
  useInexact_            = sublist.get("Inexact Solves",                           false);
  
  ROL::ParameterList krylovList;
  Real atol = static_cast<Real>(1e-12);
  Real rtol = static_cast<Real>(1e-2);
  krylovList.sublist("General").sublist("Krylov").set("Type", "GMRES");
  krylovList.sublist("General").sublist("Krylov").set("Absolute Tolerance", atol);
  krylovList.sublist("General").sublist("Krylov").set("Relative Tolerance", rtol);
  krylovList.sublist("General").sublist("Krylov").set("Iteration Limit", 200);
  krylov_ = KrylovFactory<Real>(krylovList);
}

template<typename Real>
void FletcherObjectiveBase<Real>::update( const Vector<Real> &x, UpdateType type, int iter ) {
  obj_->update(x,type,iter);
  con_->update(x,type,iter);
  fPhi_->objectiveUpdate(type);
  gPhi_->objectiveUpdate(type);
  y_->objectiveUpdate(type);
  fval_->objectiveUpdate(type);
  g_->objectiveUpdate(type);
  c_->objectiveUpdate(type);
}

// Accessors
template<typename Real>
Ptr<const Vector<Real>> FletcherObjectiveBase<Real>::getLagrangianGradient(const Vector<Real>& x) {
  // TODO: Figure out reasonable tolerance
  Real tol = static_cast<Real>(1e-12);
  computeMultipliers(*cdual_, *gLdual_, x, *xdual_, *cprim_, tol);
  gL_->set(gLdual_->dual());
  return gL_;
}

template<typename Real>
Ptr<const Vector<Real>> FletcherObjectiveBase<Real>::getConstraintVec(const Vector<Real>& x) {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  conValue(*cprim_, x, tol);  
  return cprim_;
}

template<typename Real>
Ptr<const Vector<Real>> FletcherObjectiveBase<Real>::getMultiplierVec(const Vector<Real>& x) {
  // TODO: Figure out reasonable tolerance
  Real tol = static_cast<Real>(1e-12);
  computeMultipliers(*cdual_, *gLdual_, x, *xdual_, *cprim_, tol);
  return cdual_;
}

template<typename Real>
Ptr<const Vector<Real>> FletcherObjectiveBase<Real>::getGradient(const Vector<Real>& x) {
  // TODO: Figure out reasonable tolerance
  Real tol = static_cast<Real>(1e-12);
  this->gradient(*xdual_, x, tol);
  return xdual_;
}

template<typename Real>
Real FletcherObjectiveBase<Real>::getObjectiveValue(const Vector<Real>& x) {
  Real tol = std::sqrt(ROL_EPSILON<Real>());
  return objValue(x, tol);
}

template<typename Real>
int FletcherObjectiveBase<Real>::getNumberFunctionEvaluations() const {
  return nfval_;
} 

template<typename Real>
int FletcherObjectiveBase<Real>::getNumberGradientEvaluations() const {
  return ngval_;
} 

template<typename Real>
int FletcherObjectiveBase<Real>::getNumberConstraintEvaluations() const {
  return ncval_;
}

template<typename Real>
void FletcherObjectiveBase<Real>::reset(Real sigma, Real delta) {
  sigma_ = sigma;
  delta_ = delta;
  fPhi_->reset(true);
  gPhi_->reset(true);
}

template<typename Real>
Real FletcherObjectiveBase<Real>::objValue(const Vector<Real>& x, Real &tol) {
  Real val(0);
  int key(0);
  bool isComputed = fval_->get(val,key);
  if( !isComputed ) {
    val = obj_->value(x,tol); nfval_++;
    fval_->set(val,key);
  }
  return val;
}

template<typename Real>
void FletcherObjectiveBase<Real>::objGrad(Vector<Real> &g, const Vector<Real>& x, Real &tol) {
  int key(0);
  bool isComputed = g_->get(g,key);
  if( !isComputed ) {
    obj_->gradient(g, x, tol); ngval_++;
    g_->set(g,key);
  }
}

template<typename Real>
void FletcherObjectiveBase<Real>::conValue(Vector<Real> &c, const Vector<Real>&x, Real &tol) {
  int key(0);
  bool isComputed = c_->get(c,key);
  if( !isComputed ) {
    con_->value(c, x, tol); ncval_++;
    c_->set(c,key);
  }
}

template<typename Real>
void FletcherObjectiveBase<Real>::computeMultipliers(Vector<Real> &y, Vector<Real> &gL, const Vector<Real> &x, Vector<Real> &g, Vector<Real> &c, Real tol) {
  int key(0);
  bool isComputed = y_->get(y,key);
  if (isComputed && multSolverError_ <= tol) return;
  if (!isComputed) {
    Real tol2 = tol;
    objGrad(g, x, tol2); tol2 = tol;
    conValue(c, x, tol2);
    scaledc_->set(c); scaledc_->scale(sigma_);
    cnorm_ = c.norm();
  }

  bool refine = isComputed;
  multSolverError_ = tol;
  solveAugmentedSystem(gL,y,g,*scaledc_,x,multSolverError_,refine);

  y_->set(y,key);
}

} // namespace ROL

#endif
