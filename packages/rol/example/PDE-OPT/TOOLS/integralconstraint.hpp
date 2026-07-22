// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PDE_INTEGRALCONSTRAINT_HPP
#define PDE_INTEGRALCONSTRAINT_HPP

#include "ROL_Constraint_SimOpt.hpp"
#include "integralobjective.hpp"

template<class Real>
class IntegralConstraint : public ROL::Constraint_SimOpt<Real> {
private:
  const ROL::Ptr<QoI<Real> > qoi_;
  const ROL::Ptr<Assembler<Real> > assembler_;
  ROL::Ptr<IntegralObjective<Real> > obj_;
  ROL::Ptr<ROL::Vector<Real> > dualUVector_;
  ROL::Ptr<ROL::Vector<Real> > dualZVector_;
  bool isUvecInitialized_;
  bool isZvecInitialized_;

public:
  IntegralConstraint(const ROL::Ptr<QoI<Real> > &qoi,
                     const ROL::Ptr<Assembler<Real> > &assembler)
    : qoi_(qoi), assembler_(assembler),
      isUvecInitialized_(false), isZvecInitialized_(false) {
    obj_ = ROL::makePtr<IntegralObjective<Real>>(qoi,assembler);
  }

  void setParameter(const std::vector<Real> &param) {
    ROL::Constraint_SimOpt<Real>::setParameter(param);
    obj_->setParameter(param);
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > cp =
      (dynamic_cast<ROL::StdVector<Real>&>(c)).getVector();
    (*cp)[0] = obj_->value(u,z,tol);
  }

  void applyJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    if ( !isUvecInitialized_ ) {
      dualUVector_ = u.dual().clone();
      isUvecInitialized_ = true;
    }
    ROL::Ptr<std::vector<Real> > jvp =
      (dynamic_cast<ROL::StdVector<Real>&>(jv)).getVector();
    obj_->gradient_1(*dualUVector_,u,z,tol);
    //(*jvp)[0] = v.dot(dualUVector_->dual());
    (*jvp)[0] = v.apply(*dualUVector_);
  }

  void applyJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                 const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    if ( !isZvecInitialized_ ) {
      dualZVector_ = z.dual().clone();
      isZvecInitialized_ = true;
    }
    ROL::Ptr<std::vector<Real> > jvp =
      (dynamic_cast<ROL::StdVector<Real>&>(jv)).getVector();
    obj_->gradient_2(*dualZVector_,u,z,tol);
    //(*jvp)[0] = v.dot(dualZVector_->dual());
    (*jvp)[0] = v.apply(*dualZVector_);
  }

  void applyAdjointJacobian_1(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > vp =
      (dynamic_cast<const ROL::StdVector<Real>&>(v)).getVector();
    obj_->gradient_1(jv,u,z,tol);
    jv.scale((*vp)[0]);
  }

  void applyAdjointJacobian_2(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                        const ROL::Vector<Real> &u,  const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > vp =
      (dynamic_cast<const ROL::StdVector<Real>&>(v)).getVector();
    obj_->gradient_2(jv,u,z,tol);
    jv.scale((*vp)[0]);
  }

  void applyAdjointHessian_11( ROL::Vector<Real> &ahwv,
                         const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, 
                         const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > wp =
      (dynamic_cast<const ROL::StdVector<Real>&>(w)).getVector();
    obj_->hessVec_11(ahwv,v,u,z,tol);
    ahwv.scale((*wp)[0]);
  }

  void applyAdjointHessian_12( ROL::Vector<Real> &ahwv,
                         const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, 
                         const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > wp =
      (dynamic_cast<const ROL::StdVector<Real>&>(w)).getVector();
    obj_->hessVec_12(ahwv,v,u,z,tol);
    ahwv.scale((*wp)[0]);
  }

  void applyAdjointHessian_21( ROL::Vector<Real> &ahwv,
                         const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, 
                         const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > wp =
      (dynamic_cast<const ROL::StdVector<Real>&>(w)).getVector();
    obj_->hessVec_21(ahwv,v,u,z,tol);
    ahwv.scale((*wp)[0]);
  }

  void applyAdjointHessian_22( ROL::Vector<Real> &ahwv,
                         const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, 
                         const ROL::Vector<Real> &u, const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > wp =
      (dynamic_cast<const ROL::StdVector<Real>&>(w)).getVector();
    obj_->hessVec_22(ahwv,v,u,z,tol);
    ahwv.scale((*wp)[0]);
  }
}; // class IntegralConstraint


template<class Real>
class IntegralOptConstraint : public ROL::Constraint<Real> {
private:
  const ROL::Ptr<QoI<Real> > qoi_;
  const ROL::Ptr<Assembler<Real> > assembler_;
  ROL::Ptr<IntegralOptObjective<Real> > obj_;
  ROL::Ptr<ROL::Vector<Real> > dualZVector_;
  bool isZvecInitialized_;

public:
  IntegralOptConstraint(const ROL::Ptr<QoI<Real> > &qoi,
                        const ROL::Ptr<Assembler<Real> > &assembler)
    : qoi_(qoi), assembler_(assembler), isZvecInitialized_(false) {
    obj_ = ROL::makePtr<IntegralOptObjective<Real>>(qoi,assembler);
  }

  void setParameter(const std::vector<Real> &param) {
    ROL::Constraint<Real>::setParameter(param);
    obj_->setParameter(param);
  }

  void value(ROL::Vector<Real> &c, const ROL::Vector<Real> &z, Real &tol) {
    ROL::Ptr<std::vector<Real> > cp =
      (dynamic_cast<ROL::StdVector<Real>&>(c)).getVector();
    (*cp)[0] = obj_->value(z,tol);
  }

  void applyJacobian(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &z, Real &tol ) {
    if ( !isZvecInitialized_ ) {
      dualZVector_ = z.dual().clone();
      isZvecInitialized_ = true;
    }
    ROL::Ptr<std::vector<Real> > jvp =
      (dynamic_cast<ROL::StdVector<Real>&>(jv)).getVector();
    obj_->gradient(*dualZVector_,z,tol);
    //(*jvp)[0] = v.dot(dualZVector_->dual());
    (*jvp)[0] = v.apply(*dualZVector_);
  }

  void applyAdjointJacobian(ROL::Vector<Real> &jv, const ROL::Vector<Real> &v,
                            const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > vp =
      (dynamic_cast<const ROL::StdVector<Real>&>(v)).getVector();
    obj_->gradient(jv,z,tol);
    jv.scale((*vp)[0]);
  }

  void applyAdjointHessian(ROL::Vector<Real> &ahwv,
                           const ROL::Vector<Real> &w, const ROL::Vector<Real> &v, 
                           const ROL::Vector<Real> &z, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > wp =
      (dynamic_cast<const ROL::StdVector<Real>&>(w)).getVector();
    obj_->hessVec(ahwv,v,z,tol);
    ahwv.scale((*wp)[0]);
  }

}; // class IntegralOptConstraint

#endif
