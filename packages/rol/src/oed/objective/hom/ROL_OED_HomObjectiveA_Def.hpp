// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEA_DEF_HPP
#define ROL_OED_HOMOBJECTIVEA_DEF_HPP

namespace ROL::OED::Hom {

template<typename Real>
void ObjectiveA<Real>::solveState(Vector<Real>& u, Vector<Real>& e, const Vector<Real>& z, int i) {
  std::vector<Real> param = {static_cast<Real>(i)};
  // Check if state has been computed.
  bool isComputed = storage_ ? stateStore_->get(u,i) : false;
  ts_->get(e,param);
  // Solve state equation if not done already.
  if (!isComputed || !storage_) {
    // Solve state equation.
    M_->applyInverse(u,e,z);
    // Store state.
    if (storage_) stateStore_->set(u,i);
  }
}

template<typename Real>
ObjectiveA<Real>::ObjectiveA( const Ptr<MomentOperator<Real>>& M,
                              const Ptr<const Vector<Real>>& theta,
                              const Ptr<TraceSampler<Real>>& ts,
                              const std::vector<Real>& weight,
                              bool storage)
  : M_(M), ts_(ts), weight_(weight),
    stateStore_(makePtr<VectorController<Real,int>>()),
    u_(theta->clone()), e_(theta->dual().clone()), s_(u_->clone()),
    storage_(storage) {}

template<typename Real>
void ObjectiveA<Real>::update(const Vector<Real>& z, UpdateType type, int iter) {
  M_->update(z,type,iter);
  stateStore_->objectiveUpdate(type);
}

template<typename Real>
Real ObjectiveA<Real>::value( const Vector<Real>& z, Real& tol ) {
  const int dim = weight_.size();
  Real val(0);
  for (int i = 0; i < dim; ++i) {
    // Solve state equation
    solveState(*u_,*e_,z,i);
    // Get objective function value
    val += weight_[i]*e_->apply(*u_);
  }
  return val;
}

template<typename Real>
void ObjectiveA<Real>::gradient( Vector<Real>& g, const Vector<Real>& z, Real& tol ) {
  if (g_==nullPtr) g_ = g.clone();
  g.zero();
  const int dim = weight_.size();
  for (int i = 0; i < dim; ++i) {
    std::vector<Real> param = {static_cast<Real>(i)};
    // Solve state equation
    solveState(*u_,*e_,z,i); // Sets both 
    // Build gradient
    M_->applySampleMatrices(*g_,*u_,*u_);
    g.axpy(-weight_[i],*g_);
  }
}

template<typename Real>
void ObjectiveA<Real>::hessVec( Vector<Real>& hv, const Vector<Real>& v, const Vector<Real>& z, Real& tol ) {
  if (g_==nullPtr) g_ = hv.clone();
  const int dim = weight_.size();
  hv.zero();
  for (int i = 0; i < dim; ++i) {
    std::vector<Real> param = {static_cast<Real>(i)};
    // Solve state equation
    solveState(*u_,*e_,z,i); // Sets both 
    // Solve state sensitivity equation
    M_->applyDeriv(*e_,*u_,v);
    M_->applyInverse(*s_,*e_,z);
    // Build hessVec
    M_->applySampleMatrices(*g_,*s_,*u_);
    hv.axpy(static_cast<Real>(2)*weight_[i],*g_);
  }
}

} // END ROL::OED::Hom Namespace

#endif
