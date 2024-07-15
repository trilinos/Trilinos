// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

// ROL_Types contains predefined constants and objects
#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Solver.hpp"
#include "ROL_PEBBL_IntegerConstraint.hpp"

template<class Real>
class Constraint_SimpleBinary : public ROL::Constraint<Real> {
private:
  const Real budget_;

  ROL::Ptr<std::vector<Real>> getVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<const std::vector<Real>> getConstVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

public:
  Constraint_SimpleBinary(int budget = 4) : budget_(budget) {}

  void value(ROL::Vector<Real> &c,
       const ROL::Vector<Real> &x, 
             Real &tol) {
    ROL::Ptr<std::vector<Real>>       cp = getVector(c);
    ROL::Ptr<const std::vector<Real>> xp = getConstVector(x);
    (*cp)[0] = -static_cast<Real>(budget_);
    for (const auto & y : *xp) (*cp)[0] += y;
  }

  void applyJacobian(ROL::Vector<Real> &jv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &x, 
                     Real &tol) {
    ROL::Ptr<std::vector<Real>>      jvp = getVector(jv);
    ROL::Ptr<const std::vector<Real>> vp = getConstVector(v);
    (*jvp)[0] = static_cast<Real>(0);
    for (const auto & y : *vp) (*jvp)[0] += y;
  }

  void applyAdjointJacobian(ROL::Vector<Real> &ajv,
                      const ROL::Vector<Real> &v,
                      const ROL::Vector<Real> &x, 
                            Real &tol) {
    ROL::Ptr<std::vector<Real>>      jvp = getVector(ajv);
    ROL::Ptr<const std::vector<Real>> vp = getConstVector(v);
    jvp->assign(jvp->size(),(*vp)[0]);
  }

  void applyAdjointHessian(ROL::Vector<Real> &ahwv,
                     const ROL::Vector<Real> &w,
                     const ROL::Vector<Real> &v,
                     const ROL::Vector<Real> &x,
                           Real &tol) {
    ahwv.zero();
  }
};

template<class Real>
class Objective_SimpleBinary : public ROL::Objective<Real> {
private:
  const std::vector<Real> alpha_, beta_;

  ROL::Ptr<std::vector<Real>> getVector(ROL::Vector<Real> &x) const {
    return dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
  }

  ROL::Ptr<const std::vector<Real>> getConstVector(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
  }

public:
  Objective_SimpleBinary(const std::vector<Real> &alpha)
    : alpha_(alpha), beta_(std::vector<Real>(alpha.size(),0)) {}
  Objective_SimpleBinary(const std::vector<Real> &alpha, const std::vector<Real> &beta)
    : alpha_(alpha), beta_(beta) {}

  Real value(const ROL::Vector<Real> &x, Real &tol) {
    ROL::Ptr<const std::vector<Real>> xp = getConstVector(x);
    const Real half(0.5);
    Real val(0);
    const int dim(xp->size());
    for (int i = 0; i < dim; ++i)
      val += half * alpha_[i] * (*xp)[i] * (*xp)[i] - beta_[i] * (*xp)[i];
    return val;
  }

  void gradient(ROL::Vector<Real> &g,
          const ROL::Vector<Real> &x,
                Real &tol) {
    ROL::Ptr<std::vector<Real>>       gp = getVector(g);
    ROL::Ptr<const std::vector<Real>> xp = getConstVector(x);
    const int dim(xp->size());
    for (int i = 0; i < dim; ++i)
      (*gp)[i] = alpha_[i] * (*xp)[i] - beta_[i];
  }

  void hessVec(ROL::Vector<Real> &hv,
         const ROL::Vector<Real> &v,
         const ROL::Vector<Real> &x,
               Real &tol ) {
    ROL::Ptr<std::vector<Real>>      hvp = getVector(hv);
    ROL::Ptr<const std::vector<Real>> vp = getConstVector(v);
    const int dim(vp->size());
    for (int i = 0; i < dim; ++i)
      (*hvp)[i] = alpha_[i] * (*vp)[i];
  }
};
