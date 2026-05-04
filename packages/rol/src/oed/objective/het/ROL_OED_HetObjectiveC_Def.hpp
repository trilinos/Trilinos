// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_HETOBJECTIVEC_DEF_HPP
#define ROL_OED_HETOBJECTIVEC_DEF_HPP

namespace ROL::OED::Het {

template<typename Real>
void ObjectiveC<Real>::solveState(Vector<Real>& u, Vector<Real>& Mu, const Vector<Real>& z) {
  if (!isStateComputed_) {
    M1_->applyInverse(u,*c_,z);
    M0_->apply(Mu,u,z);
    isStateComputed_ = true;
  }
}

template<typename Real>
void ObjectiveC<Real>::solveAdjoint(Vector<Real>& p, const Vector<Real>& Mu, const Vector<Real>& z) {
  if (!isAdjointComputed_) {
    M1_->applyInverse(p,Mu,z);
    p.scale(static_cast<Real>(-2));
    isAdjointComputed_ = true;
  }
}

template<typename Real>
ObjectiveC<Real>::ObjectiveC(const Ptr<MomentOperator<Real>>& M0,
                             const Ptr<MomentOperator<Real>>& M1,
                             const Ptr<const Vector<Real>>& c,
                             bool storage)
  : M0_(M0), M1_(M1), c_(c), u_(c->dual().clone()), ucache_(u_->clone()),
    Mu_(c->clone()), Mucache_(Mu_->clone()),
    p_(u_->clone()), pcache_(p_->clone()), r1_(c->clone()), r2_(c->clone()),
    s_(u_->clone()), q_(p_->clone()), storage_(storage),
    isStateComputed_(false), isStateCached_(false),
    isAdjointComputed_(false), isAdjointCached_(false) {}

template<typename Real>
void ObjectiveC<Real>::update( const Vector<Real>& z, UpdateType type, int iter) {
  M0_->update(z,type,iter);
  M1_->update(z,type,iter);
  if (storage_) {
    switch(type) {
      case ROL::UpdateType::Initial: {
        isStateComputed_   = false;
        isStateCached_     = false;
        isAdjointComputed_ = false;
        isAdjointCached_   = false;
	break;
      }
      case ROL::UpdateType::Trial: {
        isStateCached_     = isStateComputed_;
        isAdjointCached_   = isAdjointComputed_;
        isStateComputed_   = false;
        isAdjointComputed_ = false;
	break;
      }
      case ROL::UpdateType::Accept: {
        if (isStateComputed_) {
          ucache_->set(*u_);
          Mucache_->set(*Mu_);
        }
        if (isAdjointComputed_)
          pcache_->set(*p_);
	break;
      }
      case ROL::UpdateType::Revert: {
        if (isStateCached_) {
          u_->set(*ucache_);
	  Mu_->set(*Mucache_);
          isStateComputed_ = true;
        }
        else
          isStateComputed_ = false;
        if (isAdjointCached_) {
          p_->set(*pcache_);
          isAdjointComputed_ = true;
        }
        else
          isAdjointComputed_ = false;
	break;
      }
      case ROL::UpdateType::Temp: {
        isStateComputed_   = false;
        isAdjointComputed_ = false;
	break;
      }
      default: {
        isStateComputed_   = false;
        isAdjointComputed_ = false;
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
  solveState(*u_,*Mu_,z);
  // Assemble objective function value
  return Mu_->apply(*u_);
}

template<typename Real>
void ObjectiveC<Real>::gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = g.clone();  
  // Solve state equation
  solveState(*u_,*Mu_,z);
  // Solve adjoint equation
  solveAdjoint(*p_,*Mu_,z);
  // Adjoint gradient
  M0_->applySampleMatrices(*g_,*u_,*u_);
  M1_->applySampleMatrices(g,*u_,*p_);
  g.plus(*g_);
}

template<typename Real>
void ObjectiveC<Real>::hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
  if (g_ == nullPtr) g_ = hv.clone();  
  // Solve state equation
  solveState(*u_,*Mu_,z);
  // Solve adjoint equation
  solveAdjoint(*p_,*Mu_,z);
  // Solve state sensitivity equation
  M1_->applyDeriv(*r1_,*u_,v);
  M1_->applyInverse(*s_,*r1_,z);
  // Solve adjoint sensitivity equation
  const Real one(1), two(2);
  M0_->applyDeriv(*r1_,*u_,v);
  M1_->applyDeriv(*r2_,*p_,v);
  r2_->axpy(two,*r1_);
  M0_->apply(*r1_,*s_,z);
  r1_->scale(two);
  r1_->axpy(-one,*r2_);
  M1_->applyInverse(*q_,*r1_,z);
  // Assemble hessVec
  M1_->applySampleMatrices(hv,*u_,*q_);
  M0_->applySampleMatrices(*g_,*u_,*s_);
  hv.axpy(-two,*g_);
  M1_->applySampleMatrices(*g_,*p_,*s_);
  hv.axpy(-one,*g_);
}

} // End ROL::OED::Het Namespace

#endif
