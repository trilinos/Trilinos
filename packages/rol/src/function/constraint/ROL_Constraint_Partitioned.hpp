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

#ifndef ROL_CONSTRAINT_PARTITIONED_H
#define ROL_CONSTRAINT_PARTITIONED_H

#include "ROL_PartitionedVector.hpp"
#include "ROL_Constraint.hpp"

namespace ROL {

/** @ingroup func_group
 *  \class ROL::Constraint_Partitioned
 *  \brief Has both inequality and equality constraints.
 *        Treat inequality constraint as equality with slack variable
 */

template<typename Real>
class Constraint_Partitioned : public Constraint<Real> {
private:
  std::vector<Ptr<Constraint<Real>>> cvec_;
  std::vector<bool> isInequality_; // Label whether cvec_[i] is inequality
  const int offset_;               // Offset for slack variables
  Ptr<Vector<Real>> scratch_;      // Scratch vector for intermediate computation
  int  ncval_;                     // Number of constraint evaluations
  bool initialized_;               // Is scratch vector initialized?

public:
  Constraint_Partitioned(const std::vector<Ptr<Constraint<Real>>> &cvec,
                         bool isInequality = false,
                         int offset = 0);
  Constraint_Partitioned(const std::vector<Ptr<Constraint<Real>>> &cvec,
                         std::vector<bool> isInequality,
                         int offset = 0);

  int getNumberConstraintEvaluations(void) const;
  Ptr<Constraint<Real>> get(int ind = 0) const;

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  void value( Vector<Real> &c, const Vector<Real> &x, Real &tol ) override;
  void applyJacobian( Vector<Real> &jv,
                      const Vector<Real> &v,
                      const Vector<Real> &x,
                      Real &tol ) override;
  using Constraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian( Vector<Real> &ajv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol ) override;
  void applyAdjointHessian( Vector<Real> &ahuv,
                            const Vector<Real> &u,
                            const Vector<Real> &v,
                            const Vector<Real> &x,
                            Real &tol ) override;
  virtual void applyPreconditioner(Vector<Real> &pv,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   const Vector<Real> &g,
                                   Real &tol) override;

// Definitions for parametrized (stochastic) equality constraints
public:
  void setParameter(const std::vector<Real> &param) override;

private:
  Vector<Real>& getOpt( Vector<Real> &xs ) const;
  const Vector<Real>& getOpt( const Vector<Real> &xs ) const;
  Vector<Real>& getSlack( Vector<Real> &xs, int ind ) const;
  const Vector<Real>& getSlack( const Vector<Real> &xs, int ind ) const;

}; // class Constraint_Partitioned

} // namespace ROL

#include "ROL_Constraint_Partitioned_Def.hpp"

#endif
