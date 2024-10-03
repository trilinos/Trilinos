// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_BILINEARCONSTRAINT_DEF_HPP
#define ROL_OED_BILINEARCONSTRAINT_DEF_HPP

namespace ROL {
namespace OED {

template<typename Real>
void BilinearConstraint<Real>::computeG(Vector<Real> &g) {
  startTimer("computeG");
  if (type_ != "C") {
    std::vector<Real> param = Constraint_SimOpt<Real>::getParameter();
    if ((type_ == "I" && !useTrace_) || type_ == "R") {
      factors_->evaluate(g,param);
    }
    else if (type_ == "A" || type_ == "D" || (type_ == "I" && useTrace_)) {
      traceSampler_->get(g,param);
    }
    else {
      throw Exception::NotImplemented(">>> ROL::OED::BilinearConstraint::computeG : Optimality type not implemented!");
    }
  }
  stopTimer("computeG");
}

// C optimality
template<typename Real>
BilinearConstraint<Real>::BilinearConstraint(const Ptr<Factors<Real>>        &factors,
                   const Ptr<MomentOperator<Real>> &M,
                   const Ptr<Vector<Real>>    &c)
  : ProfiledClass<Real,std::string>("OED::BilinearConstraint"),
    factors_(factors), M_(M), type_("C"), g_(c->clone()),
    isPinit_(false) {
  g_->set(*c);
}

// A, D, I, and R optimality
template<typename Real>
BilinearConstraint<Real>::BilinearConstraint(const Ptr<Factors<Real>> &factors,
                   const Ptr<MomentOperator<Real>> &M,
                   const std::string               &type,
                   const Ptr<TraceSampler<Real>>   &traceSampler)
  : ProfiledClass<Real,std::string>("OED::BilinearConstraint"),
    factors_(factors), M_(M), traceSampler_(traceSampler),
    type_(type), g_(factors_->get(0)->clone()),
    isPinit_(false) {
  if (type_ != "A" && type_ != "D" && type_ != "I" && type_ != "R") {
    std::stringstream ss;
    ss << ">>> ROL::OED::BilinearConstraint : Wrong constructor for " << type_ << "-optimality!";
    throw Exception::NotImplemented(ss.str());
  }
  useTrace_ = (type_ == "I" && traceSampler_ == nullPtr ? false : true);
  if (type_ == "A" || type_ == "D") {
    if (traceSampler_ == nullPtr) {
      Ptr<Vector<Real>> u = g_->dual().clone();
      traceSampler_ = makePtr<TraceSampler<Real>>(u);
    }
  }
}

template<typename Real>
void BilinearConstraint<Real>::update_2(const Vector<Real> &z,
              UpdateType type,
              int iter) {
  M_->update(z,type,iter);
}

template<typename Real>
void BilinearConstraint<Real>::value(Vector<Real> &c,
           const Vector<Real> &u,
           const Vector<Real> &z,
           Real &tol) {
  startTimer("value");
  M_->apply(c,u,z);
  computeG(*g_);
  c.axpy(static_cast<Real>(-1),*g_);
  stopTimer("value");
}

template<typename Real>
void BilinearConstraint<Real>::solve(Vector<Real> &c,
           Vector<Real> &u, 
           const Vector<Real> &z,
           Real &tol) {
  startTimer("solve");
  computeG(*g_);
  M_->applyInverse(u,*g_,z);
  value(c,u,z,tol);
  stopTimer("solve");
}

template<typename Real>
void BilinearConstraint<Real>::applyJacobian_1(Vector<Real> &jv,
                     const Vector<Real> &v,
                     const Vector<Real> &u,
                     const Vector<Real> &z,
                     Real &tol) {
  startTimer("applyJacobian_1");
  M_->apply(jv,v,z);
  stopTimer("applyJacobian_1");
}

template<typename Real>
void BilinearConstraint<Real>::applyJacobian_2(Vector<Real> &jv,
                     const Vector<Real> &v,
                     const Vector<Real> &u,
                     const Vector<Real> &z,
                     Real &tol) {
  startTimer("applyJacobian_2");
  M_->applyDeriv(jv,u,v);
  stopTimer("applyJacobian_2");
}

template<typename Real>
void BilinearConstraint<Real>::applyInverseJacobian_1(Vector<Real> &ijv,
                            const Vector<Real> &v,
                            const Vector<Real> &u,
                            const Vector<Real> &z,
                            Real &tol) {
  startTimer("applyInverseJacobian_1");
  M_->applyInverse(ijv,v,z);
  stopTimer("applyInverseJacobian_1");
}

template<typename Real>
void BilinearConstraint<Real>::applyAdjointJacobian_1(Vector<Real> &ajv,
                            const Vector<Real> &v,
                            const Vector<Real> &u,
                            const Vector<Real> &z,
                            Real &tol) {
  startTimer("applyAdjointJacobian_1");
  applyJacobian_1(ajv,v,u,z,tol);
  stopTimer("applyAdjointJacobian_1");
}

template<typename Real>
void BilinearConstraint<Real>::applyAdjointJacobian_2(Vector<Real> &ajv,
                            const Vector<Real> &v,
                            const Vector<Real> &u,
                            const Vector<Real> &z,
                            Real &tol) {
  startTimer("applyAdjointJacobian_2");
  if (!isPinit_) {
    p_ = ajv.clone();
    isPinit_ = true;
  }
  M_->applySampleMatrices(ajv,u,v);
  stopTimer("applyAdjointJacobian_2");
}

template<typename Real>
void BilinearConstraint<Real>::applyInverseAdjointJacobian_1(Vector<Real> &iajv,
                                   const Vector<Real> &v,
                                   const Vector<Real> &u,
                                   const Vector<Real> &z,
                                   Real &tol) {
  startTimer("applyInverseAdjointJacobian_1");
  applyInverseJacobian_1(iajv,v,u,z,tol);
  stopTimer("applyInverseAdjointJacobian_1");
}

template<typename Real>
void BilinearConstraint<Real>::applyAdjointHessian_11(Vector<Real> &ahwv,
                            const Vector<Real> &w,
                            const Vector<Real> &v,
                            const Vector<Real> &u,
                            const Vector<Real> &z,
                            Real &tol) {
  startTimer("applyAdjointHessian_11");
  ahwv.zero();
  stopTimer("applyAdjointHessian_11");
}

template<typename Real>
void BilinearConstraint<Real>::applyAdjointHessian_12(Vector<Real> &ahwv,
                            const Vector<Real> &w,
                            const Vector<Real> &v,
                            const Vector<Real> &u,
                            const Vector<Real> &z,
                            Real &tol) {
  startTimer("applyAdjointHessian_12");
  applyAdjointJacobian_2(ahwv,v,w,z,tol);
  stopTimer("applyAdjointHessian_12");
}

template<typename Real>
void BilinearConstraint<Real>::applyAdjointHessian_21(Vector<Real> &ahwv,
                            const Vector<Real> &w,
                            const Vector<Real> &v,
                            const Vector<Real> &u,
                            const Vector<Real> &z,
                            Real &tol) {
  startTimer("applyAdjointHessian_21");
  applyJacobian_2(ahwv,v,w,z,tol);
  stopTimer("applyAdjointHessian_21");
}

template<typename Real>
void BilinearConstraint<Real>::applyAdjointHessian_22(Vector<Real> &ahwv,
                            const Vector<Real> &w,
                            const Vector<Real> &v,
                            const Vector<Real> &u,
                            const Vector<Real> &z,
                            Real &tol) {
  startTimer("applyAdjointHessian_22");
  ahwv.zero();
  stopTimer("applyAdjointHessian_22");
}

template<typename Real>
void BilinearConstraint<Real>::getFactor(Vector<Real> &F, int k) const {
  F.set(*factors_->get(k));
}

template<typename Real>
void BilinearConstraint<Real>::getFactor(Vector<Real> &F, const std::vector<Real> &param) const {
  factors_->evaluate(F,param);
}

template<typename Real>
Real BilinearConstraint<Real>::getNoise(int k) const {
  return M_->getNoise(k);
}

template<typename Real>
void BilinearConstraint<Real>::getTraceSample(Vector<Real> &g, const std::vector<Real> &param) const {
  traceSampler_->get(g,param);
}

// template<typename Real>
// void BilinearConstraint<Real>::sumAll(Real *in, Real *out, int size) const {
//   factors_->sumAll(in,out,size);
// }

template<typename Real>
Real BilinearConstraint<Real>::logDeterminant(const Vector<Real> &z) {
  return M_->logDeterminant(z);
}

} // End OED Namespace
} // End ROL Namespace

#endif
