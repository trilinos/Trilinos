// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HOMOBJECTIVEC_DEF_HPP
#define ROL_OED_HOMOBJECTIVEC_DEF_HPP

namespace ROL::OED::Hom {

template<typename Real>
void ObjectiveC<Real>::solveState(Vector<Real>& u, const Vector<Real>& z) {
  if (!isStateComputed_) {
    M_->applyInverse(u,*c_,z);
    isStateComputed_ = true;
  }
}

template<typename Real>
ObjectiveC<Real>::ObjectiveC(const Ptr<MomentOperator<Real>>& M,
                             const Ptr<const Vector<Real>>& c,
                             bool storage)
  : M_(M), c_(c), u_(c->dual().clone()), ucache_(u_->clone()),
    r_(c->clone()), s_(u_->clone()), storage_(storage),
    isStateComputed_(false), isStateCached_(false) {}

template<typename Real>
void ObjectiveC<Real>::update( const Vector<Real>& z, UpdateType type, int iter) {
  M_->update(z,type,iter);
  if (storage_) {
    switch(type) {
      case ROL::UpdateType::Initial: {
        isStateComputed_ = false;
        isStateCached_   = false;
	break;
      }
      case ROL::UpdateType::Trial: {
        isStateCached_   = isStateComputed_;
        isStateComputed_ = false;
	break;
      }
      case ROL::UpdateType::Accept: {
        if (isStateComputed_)
          ucache_->set(*u_);
	break;
      }
      case ROL::UpdateType::Revert: {
        if (isStateCached_) {
          u_->set(*ucache_);
          isStateComputed_ = true;
        }
        else
          isStateComputed_ = false;
	break;
      }
      case ROL::UpdateType::Temp: {
        isStateComputed_ = false;
	break;
      }
      default: {
        isStateComputed_ = false;
	break;
      }
    }
  }
  else {
    isStateComputed_ = false;
  }
}

template<typename Real>
Real ObjectiveC<Real>::value( const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  solveState(*u_,z);
  // Assemble objective value
  return c_->apply(*u_);
}

template<typename Real>
void ObjectiveC<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  solveState(*u_,z);
  // Assemble gradient
  M_->applySampleMatrices(g,*u_,*u_);
  g.scale(static_cast<Real>(-1));
}

template<typename Real>
void ObjectiveC<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  // Solve state equation
  solveState(*u_,z);
  // Solve state sensitivity equation
  M_->applyDeriv(*r_,*u_,v);
  M_->applyInverse(*s_,*r_,z);
  // Assemble Hessian application
  M_->applySampleMatrices(hv,*s_,*u_);
  hv.scale(static_cast<Real>(2));
}

} // END ROL::OED::Hom Namespace

#endif
