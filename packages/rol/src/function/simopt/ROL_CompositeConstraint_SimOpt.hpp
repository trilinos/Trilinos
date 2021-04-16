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

#ifndef ROL_COMPOSITECONSTRAINT_SIMOPT_H
#define ROL_COMPOSITECONSTRAINT_SIMOPT_H

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_VectorController.hpp"

/** @ingroup func_group
    \class ROL::CompositeConstraint_SimOpt
    \brief Defines a composite equality constraint operator interface for
           simulation-based optimization.

    This equality constraint interface inherits from ROL_Constraint_SimOpt, for the
    use case when \f$\mathcal{X}=\mathcal{U}\times\mathcal{Z}\f$ where \f$\mathcal{U}\f$ and 
    \f$\mathcal{Z}\f$ are Banach spaces.  \f$\mathcal{U}\f$ denotes the "simulation space"
    and \f$\mathcal{Z}\f$ denotes the "optimization space" (of designs, controls, parameters).
    The simulation-based constraints are of the form
    \f[
      c(u,S(z)) = 0
    \f]
    where \f$S(z)\f$ solves the reducible constraint
    \f[
       c_0(S(z),z) = 0.
    \f]

    ---
*/

namespace ROL {

template<typename Real>
class CompositeConstraint_SimOpt : public Constraint_SimOpt<Real> {
private:
  // Constraints
  const ROL::Ptr<Constraint_SimOpt<Real>> conVal_, conRed_;
  // Additional vector storage for solve
  ROL::Ptr<Vector<Real>> Sz_, primRed_, dualRed_, primZ_, dualZ_, dualZ1_, primU_;
  // State storage through VectorController interface
  ROL::Ptr<VectorController<Real>> stateStore_;
  // Update information
  bool updateFlag_, newUpdate_;
  int updateIter_;
  UpdateType updateType_;
  // Boolean variables
  const bool storage_, isConRedParametrized_;

public:
  CompositeConstraint_SimOpt(const ROL::Ptr<Constraint_SimOpt<Real>> &conVal,
                             const ROL::Ptr<Constraint_SimOpt<Real>> &conRed,
                             const Vector<Real> &cVal,
                             const Vector<Real> &cRed,
                             const Vector<Real> &u,
                             const Vector<Real> &Sz,
                             const Vector<Real> &z,
                             bool storage = true,
                             bool isConRedParametrized = false);

  void update(const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1) override;
  void update_1(const Vector<Real> &u, bool flag = true, int iter = -1) override;
  void update_2(const Vector<Real> &z, bool flag = true, int iter = -1) override;
  void update(const Vector<Real> &u, const Vector<Real> &z, UpdateType type, int iter = -1) override;
  void update_1(const Vector<Real> &u, UpdateType type, int iter = -1) override;
  void update_2(const Vector<Real> &z, UpdateType type, int iter = -1) override;
  void solve_update(const Vector<Real> &u, const Vector<Real> &z, UpdateType type, int iter = -1) override;
  void value(Vector<Real> &c, const Vector<Real> &u, const Vector<Real> &z, Real &tol) override;
  void solve(Vector<Real> &c, Vector<Real> &u, const Vector<Real> &z, Real &tol) override;
  void applyJacobian_1(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u,
                       const Vector<Real> &z, Real &tol) override;
  void applyJacobian_2(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u,
                       const Vector<Real> &z, Real &tol) override; 
  void applyInverseJacobian_1(Vector<Real> &ijv, const Vector<Real> &v, const Vector<Real> &u,
                              const Vector<Real> &z, Real &tol) override;
  void applyAdjointJacobian_1(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &u,
                              const Vector<Real> &z, Real &tol) override;
  void applyAdjointJacobian_2(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &u,
                              const Vector<Real> &z, Real &tol) override;
  void applyInverseAdjointJacobian_1(Vector<Real> &ijv, const Vector<Real> &v, const Vector<Real> &u,
                                     const Vector<Real> &z, Real &tol) override;
  void applyAdjointHessian_11(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                              const Vector<Real> &u, const Vector<Real> &z, Real &tol) override;
  void applyAdjointHessian_12(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                              const Vector<Real> &u, const Vector<Real> &z, Real &tol) override;
  void applyAdjointHessian_21(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                              const Vector<Real> &u, const Vector<Real> &z, Real &tol) override;
  void applyAdjointHessian_22(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                              const Vector<Real> &u, const Vector<Real> &z, Real &tol) override;

// Definitions for parametrized (stochastic) equality constraints
public:
  void setParameter(const std::vector<Real> &param) override;

private:
  void solveConRed(Vector<Real> &Sz, const Vector<Real> &z, Real &tol);
  void applySens(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &Sz, const Vector<Real> &z, Real &tol);
  void applyAdjointSens(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &Sz, const Vector<Real> &z, Real &tol);
}; // class CompositeConstraint_SimOpt

} // namespace ROL

#include "ROL_CompositeConstraint_SimOpt_Def.hpp"

#endif
