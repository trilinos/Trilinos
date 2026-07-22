// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEBINDICATOROBJECTIVE_H
#define ROL_TYPEBINDICATOROBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_PolyhedralProjectionFactory.hpp"

/** @ingroup func_group
    \class ROL::TypeBIndicatorObjective
    \brief Provides the interface to evaluate the indicator function of linear constraints.

        ---
*/


namespace ROL {

template<typename Real>
class TypeBIndicatorObjective : public Objective<Real> {
private:
  const Ptr<PolyhedralProjection<Real>> proj_;
  const Ptr<Vector<Real>> res_;
  bool isInit_;
  Real tol_;
 
public:

  TypeBIndicatorObjective(const Ptr<BoundConstraint<Real>> &bnd)
    : proj_(makePtr<PolyhedralProjection<Real>>(bnd)),
      isInit_(true), tol_(0) {}

  TypeBIndicatorObjective(const Vector<Real>               &xprim,
                          const Vector<Real>               &xdual,
                          const Ptr<BoundConstraint<Real>> &bnd,
                          const Ptr<Constraint<Real>>      &con,
                          const Vector<Real>               &mul,
                          const Vector<Real>               &res,
                          ParameterList                    &list)
    : proj_(PolyhedralProjectionFactory<Real>(xprim,xdual,bnd,con,mul,res,list)),
      res_(res.clone()), isInit_(false) {}

  TypeBIndicatorObjective(const Ptr<PolyhedralProjection<Real>> &proj)
    : proj_(proj), res_(proj->getResidual()->clone()), isInit_(false) {}

  void initialize(const Vector<Real> &x) {
    if (!isInit_) {
      auto xz = x.clone(); xz->zero();
      Real tol(std::sqrt(ROL_EPSILON<Real>()));
      tol_ = static_cast<Real>(1e-2)*tol;
      proj_->getLinearConstraint()->value(*res_,*xz,tol);
      Real rnorm = res_->norm();
      if (rnorm > ROL_EPSILON<Real>()) tol_ *= rnorm;
      isInit_ = true;
    }
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    initialize(x);
    const Real zero(0);
    bool isBndFeasible = proj_->getBoundConstraint()->isFeasible(x); 
    bool isConFeasible = true;
    if (res_ != nullPtr) {
      proj_->getLinearConstraint()->value(*res_,x,tol);
      if (res_->norm() > tol_) isConFeasible = false;
    }
    return (isBndFeasible && isConFeasible) ? zero : ROL_INF<Real>();
  }

  void prox( Vector<Real> &Pv, const Vector<Real> &v, Real t, Real &tol){
    Pv.set(v); proj_->project(Pv);
  }
}; // class TypeBIndicatorObjective

} // namespace ROL

#endif
