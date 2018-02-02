// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
    (*jvp)[0] = v.dot(dualUVector_->dual());
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
    (*jvp)[0] = v.dot(dualZVector_->dual());
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
    (*jvp)[0] = v.dot(dualZVector_->dual());
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
