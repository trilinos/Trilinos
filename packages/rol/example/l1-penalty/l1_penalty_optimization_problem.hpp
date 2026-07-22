// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#ifndef L1_OPTIMIZATION_PROBLEM_HPP
#define L1_OPTIMIZATION_PROBLEM_HPP

#include "ROL_Bounds.hpp"
#include "ROL_BoundConstraint_Partitioned.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "l1_penalty_objective.hpp"
#include "l1_penalty_constraint.hpp"
#include "identity_constraint.hpp"

namespace ROL {

template<typename Real>
Ptr<OptimizationProblem<Real>>
make_l1_penalty_problem( const Ptr<Objective<Real>>&  obj,
                         const Ptr<Vector<Real>>&     x,
                         const Ptr<Constraint<Real>>& con,
                         const Ptr<Vector<Real>>&     emul,
                               Real                   gamma ) {
  auto r  = emul->clone();
  auto s  = emul->clone();
  auto lb = emul->clone();
  auto ub = emul->clone();
 
  r->setScalar(1.0);
  s->setScalar(1.0);
  lb->zero();
  ub->setScalar(ROL_INF<Real>());
  
  std::vector<Ptr<Vector<Real>>> xrs = {x,r,s};
  auto xp = makePtr<PartitionedVector<Real>>( xrs );

  auto x_bnd  = makePtr<BoundConstraint<Real>>(*x);
  auto rs_bnd = makePtr<Bounds<Real>>(lb,ub);

  x_bnd->deactivate();
  rs_bnd->activateLower();
  rs_bnd->deactivateUpper();

  auto bnds = std::vector<Ptr<BoundConstraint<Real>>>{{x_bnd,rs_bnd,rs_bnd}};
  auto bnd = makePtr<BoundConstraint_Partitioned<Real>>(bnds,xrs);

  auto lp_obj = makePtr<L1PenaltyObjective<Real>>(obj,*emul,gamma);
  auto lp_con = makePtr<L1PenaltyConstraint<Real>>(con);

  return makePtr<OptimizationProblem<Real>>(lp_obj,xp,bnd,lp_con,emul);
}

template<typename Real>
Ptr<OptimizationProblem<Real>>
make_l1_penalty_problem( const Ptr<Objective<Real>>& obj,
                         const Ptr<Vector<Real>>&    x,
                               Real                  gamma ) {
  Ptr<Constraint<Real>> con = makePtr<IdentityConstraint<Real>>();
  return make_l1_penalty_problem(obj,x,con,x,gamma);
}

} // namespace ROL


#endif // L1_OPTIMIZATION_PROBLEM_HPP

