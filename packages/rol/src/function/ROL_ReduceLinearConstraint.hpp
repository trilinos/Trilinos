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

#ifndef ROL_REDUCE_LINEAR_CONSTRAINT_H
#define ROL_REDUCE_LINEAR_CONSTRAINT_H

#include "ROL_AffineTransformObjective.hpp"
#include "ROL_AffineTransformConstraint.hpp"
#include "ROL_NullSpaceOperator.hpp"
#include "ROL_RangeSpaceOperator.hpp"

/** @ingroup func_group
    \class ROL::ReduceLinearConstraint
    \brief Performs null-space transformation for reducible linear equality
           constraints.

    ---
*/

namespace ROL {

template<typename Real>
class ReduceLinearConstraint {
private:
  const Ptr<Constraint<Real>> lcon_;
  const Ptr<Vector<Real>>     x_;
  const Ptr<VectorController<Real>>  storage_;
  const Ptr<NullSpaceOperator<Real>> nsop_;

public:
  virtual ~ReduceLinearConstraint(void) {}

  ReduceLinearConstraint(const Ptr<Constraint<Real>>   &lcon,
                         const Ptr<Vector<Real>>       &x,
                         const Ptr<const Vector<Real>> &c);

  Ptr<Objective<Real>> transform(const Ptr<Objective<Real>> &obj) const;
  Ptr<Constraint<Real>> transform(const Ptr<Constraint<Real>> &con) const;
  Ptr<Constraint<Real>> getLinearConstraint(void) const;
  Ptr<const Vector<Real>> getFeasibleVector(void) const;
  void project(Vector<Real> &x, const Vector<Real> &y) const;
  void project(const Ptr<Vector<Real>> &x, const Ptr<const Vector<Real>> &y) const;

private:
  void feasible(const Ptr<const Vector<Real>> &c);

}; // class ReduceLinearConstraint

} // namespace ROL

#include "ROL_ReduceLinearConstraint_Def.hpp"

#endif
