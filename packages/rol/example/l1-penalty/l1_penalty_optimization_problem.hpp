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

