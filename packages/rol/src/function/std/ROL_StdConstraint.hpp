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

#ifndef ROL_STDEQUALITY_CONSTRAINT_H
#define ROL_STDEQUALITY_CONSTRAINT_H

#include "ROL_Constraint.hpp"
#include "ROL_StdVector.hpp"

/** @ingroup func_group
    \class ROL::StdConstraint
    \brief Defines the equality constraint operator interface for StdVectors

*/

namespace ROL {

template<typename Real>
class StdConstraint : public virtual Constraint<Real> {
public:
  virtual ~StdConstraint() {}

  using Constraint<Real>::update;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  virtual void update( const std::vector<Real> &x, bool flag = true, int iter = -1 ) {}  
  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  virtual void update( const std::vector<Real> &x, UpdateType type, int iter = -1 ) {}

  using Constraint<Real>::value;
  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) override;
  virtual void value( std::vector<Real> &c, const std::vector<Real> &x, Real &tol ) = 0;

  using Constraint<Real>::applyJacobian;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, 
                             const Vector<Real> &x, Real &tol) override;
  virtual void applyJacobian( std::vector<Real> &jv, const std::vector<Real> &v, 
                              const std::vector<Real> &x, Real &tol );

  using Constraint<Real>::applyAdjointJacobian;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v,
                                    const Vector<Real> &x, Real &tol) override;
   virtual void applyAdjointJacobian( std::vector<Real> &ajv, const std::vector<Real> &v, 
                                      const std::vector<Real> &x, Real &tol );

  using Constraint<Real>::applyAdjointHessian;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                           const Vector<Real> &x, Real &tol) override;
  virtual void applyAdjointHessian( std::vector<Real> &ahuv, const std::vector<Real> &u,
                                    const std::vector<Real> &v, const std::vector<Real> &x,
                                    Real &tol );

  using Constraint<Real>::solveAugmentedSystem;
  std::vector<Real> solveAugmentedSystem(Vector<Real> &v1, Vector<Real> &v2,
                                         const Vector<Real> &b1, const Vector<Real> &b2,
                                         const Vector<Real> &x, Real &tol) override;
  virtual std::vector<Real> solveAugmentedSystem( std::vector<Real> &v1, std::vector<Real> &v2,
                                                  const std::vector<Real> &b1, const std::vector<Real> &b2,
                                                  const std::vector<Real> &x, Real tol );

  using Constraint<Real>::applyPreconditioner;
  void applyPreconditioner(Vector<Real> &pv, const Vector<Real> &v, const Vector<Real> &x,
                           const Vector<Real> &g, Real &tol) override;
  virtual void applyPreconditioner( std::vector<Real> &pv, const std::vector<Real> &v,
                                    const std::vector<Real> &x, const std::vector<Real> &g, Real &tol );

}; // class StdConstraint

} // namespace ROL

#include "ROL_StdConstraint_Def.hpp"

#endif
