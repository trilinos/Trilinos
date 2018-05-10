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

#ifndef ROL_COMPOSITE_EQUALITY_CONSTRAINT_SIMOPT_H
#define ROL_COMPOSITE_EQUALITY_CONSTRAINT_SIMOPT_H

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_SimController.hpp"

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

template <class Real>
class CompositeConstraint_SimOpt : public Constraint_SimOpt<Real> {
private:
  // Constraints
  const ROL::Ptr<Constraint_SimOpt<Real> > conVal_;
  const ROL::Ptr<Constraint_SimOpt<Real> > conRed_;
  // Additional vector storage for solve
  ROL::Ptr<Vector<Real> > Sz_;
  ROL::Ptr<Vector<Real> > primRed_;
  ROL::Ptr<Vector<Real> > dualRed_;
  ROL::Ptr<Vector<Real> > primZ_;
  ROL::Ptr<Vector<Real> > dualZ_;
  ROL::Ptr<Vector<Real> > dualZ1_;
  // State storage through SimController interface
  ROL::Ptr<SimController<Real> > stateStore_;
  // Update information
  bool updateFlag_;
  int updateIter_;
  // Boolean variables
  const bool storage_, isConRedParametrized_;

  void solveConRed(const Vector<Real> &z, Real &tol) {
    std::vector<Real> param = Constraint_SimOpt<Real>::getParameter();
    // Check if state has been computed.
    bool isComputed = false;
    if (storage_) {
      isComputed = stateStore_->get(*Sz_,param);
    }
    // Solve state equation if not done already.
    if (!isComputed || !storage_) {
      // Update equality constraint with new Opt variable.
      conRed_->update_2(z,updateFlag_,updateIter_);
      // Solve state equation.
      conRed_->solve(*primRed_,*Sz_,z,tol);
      // Update equality constraint with new Sim variable.
      conRed_->update_1(*Sz_,updateFlag_,updateIter_);
      // Update equality constraint.
      conRed_->update(*Sz_, z, updateFlag_, updateIter_);
      // Store state.
      if (storage_) {
        stateStore_->set(*Sz_,param);
      }
    }
  }

  void applySens(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &z, Real &tol) { 
    // Solve reducible constraint
    solveConRed(z, tol);
    // Solve linearization of reducible constraint in direction v
    conRed_->applyJacobian_2(*primRed_, v, *Sz_, z, tol);
    conRed_->applyInverseJacobian_1(jv, *primRed_, *Sz_, z, tol);
    jv.scale(static_cast<Real>(-1));
  }

  void applyAdjointSens(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
    // Solve reducible constraint
    solveConRed(z, tol);
    // Solve adjoint of linearized reducible constraint
    conRed_->applyInverseAdjointJacobian_1(*dualRed_, v, *Sz_, z, tol);
    conRed_->applyAdjointJacobian_2(ajv, *dualRed_, *Sz_, z, tol);
    ajv.scale(static_cast<Real>(-1));
  }

public:
  CompositeConstraint_SimOpt(const ROL::Ptr<Constraint_SimOpt<Real> > &conVal,
                             const ROL::Ptr<Constraint_SimOpt<Real> > &conRed,
                             const Vector<Real> &cVal, const Vector<Real> &cRed,
                             const Vector<Real> &u, const Vector<Real> &Sz, const Vector<Real> &z,
                             const bool storage = true, const bool isConRedParametrized = false)
    : Constraint_SimOpt<Real>(), conVal_(conVal), conRed_(conRed),
      updateFlag_(true), updateIter_(0), storage_(storage),
      isConRedParametrized_(isConRedParametrized) {
    Sz_      = Sz.clone();
    primRed_ = cRed.clone();
    dualRed_ = cRed.dual().clone();
    primZ_   = z.clone();
    dualZ_   = z.dual().clone();
    dualZ1_  = z.dual().clone();
    stateStore_ = ROL::makePtr<SimController<Real>>();
  }

  void update(const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1 ) {
    // Update this
    update_2(z, flag, iter);
    update_1(u, flag, iter);
  }

  void update_1( const Vector<Real> &u, bool flag = true, int iter = -1 ) {
    conVal_->update_1(u, flag, iter);
    // Update constraints with solution to reducible constraint
    conVal_->update(u, *Sz_, flag, iter);
  }

  void update_2( const Vector<Real> &z, bool flag = true, int iter = -1 ) {
    //conRed_->update_2(z, flag, iter);
    // Solve reducible constraint
    updateFlag_ = flag;
    updateIter_ = iter;
    Real ctol = std::sqrt(ROL_EPSILON<Real>());
    stateStore_->equalityConstraintUpdate(true);
    solveConRed(z, ctol);
  }

  void value(Vector<Real> &c, const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    solveConRed(z, tol);
    conVal_->value(c, u, *Sz_, tol);
  }

  void solve(Vector<Real> &c, Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    solveConRed(z, tol);
    conVal_->solve(c, u, *Sz_, tol);
  }

  void applyJacobian_1(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u,
                       const Vector<Real> &z, Real &tol) {
    solveConRed(z, tol);
    conVal_->applyJacobian_1(jv, v, u, *Sz_, tol);
  }

  void applyJacobian_2(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &u,
                       const Vector<Real> &z, Real &tol) { 
    applySens(*primZ_, v, z, tol);
    conVal_->applyJacobian_2(jv, *primZ_, u, *Sz_, tol);
  }

  void applyInverseJacobian_1(Vector<Real> &ijv, const Vector<Real> &v, const Vector<Real> &u,
                              const Vector<Real> &z, Real &tol) {
    solveConRed(z, tol);
    conVal_->applyInverseJacobian_1(ijv, v, u, *Sz_, tol);
  }

  void applyAdjointJacobian_1(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &u,
                              const Vector<Real> &z, Real &tol) {
    solveConRed(z, tol);
    conVal_->applyAdjointJacobian_1(ajv, v, u, *Sz_, tol);
  }

  void applyAdjointJacobian_2(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &u,
                              const Vector<Real> &z, Real &tol) {
    solveConRed(z, tol);
    conVal_->applyAdjointJacobian_2(*dualZ_, v, u, *Sz_, tol);
    applyAdjointSens(ajv, *dualZ_, z, tol);
  }

  void applyInverseAdjointJacobian_1(Vector<Real> &ijv, const Vector<Real> &v, const Vector<Real> &u,
                                     const Vector<Real> &z, Real &tol) {
    solveConRed(z, tol);
    conVal_->applyInverseAdjointJacobian_1(ijv, v, u, *Sz_, tol);
  }

  void applyAdjointHessian_11(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                              const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    solveConRed(z, tol);
    conVal_->applyAdjointHessian_11(ahwv, w, v, u, z, tol);
  }

  void applyAdjointHessian_12(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                              const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    solveConRed(z, tol);
    conVal_->applyAdjointHessian_12(*dualZ_, w, v, u, *Sz_, tol);
    applyAdjointSens(ahwv, *dualZ_, z, tol);
  }

  void applyAdjointHessian_21(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                              const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    applySens(*primZ_, v, z, tol);
    conVal_->applyAdjointHessian_21(ahwv, w, *primZ_, u, *Sz_, tol);
  }

  void applyAdjointHessian_22(Vector<Real> &ahwv, const Vector<Real> &w, const Vector<Real> &v,
                              const Vector<Real> &u, const Vector<Real> &z, Real &tol) {
    ahwv.zero();
    applySens(*primZ_, v, z, tol);

    conVal_->applyAdjointJacobian_2(*dualZ_, w, u, *Sz_, tol);
    conRed_->applyInverseAdjointJacobian_1(*dualRed_, *dualZ_, *Sz_, z, tol);
    conRed_->applyAdjointHessian_22(*dualZ_, *dualRed_, v, *Sz_, z, tol);
    ahwv.axpy(static_cast<Real>(-1), *dualZ_);
    conRed_->applyAdjointHessian_12(*dualZ_, *dualRed_, *primZ_, *Sz_, z, tol);
    ahwv.axpy(static_cast<Real>(-1), *dualZ_);

    conRed_->applyAdjointHessian_11(*dualZ1_, *dualRed_, *primZ_, *Sz_, z, tol);
    conRed_->applyAdjointHessian_21(*dualZ_, *dualRed_, v, *Sz_, z, tol);
    dualZ1_->plus(*dualZ_); 
    dualZ1_->scale(static_cast<Real>(-1));
    
    conVal_->applyAdjointHessian_22(*dualZ_, w, *primZ_, u, *Sz_, tol);
    dualZ1_->plus(*dualZ_); 

    applyAdjointSens(*dualZ_, *dualZ1_, z, tol);
    ahwv.plus(*dualZ_);
  }

// Definitions for parametrized (stochastic) equality constraints
public:
  void setParameter(const std::vector<Real> &param) {
    conVal_->setParameter(param);
    if (isConRedParametrized_) {
      conRed_->setParameter(param);
      Constraint_SimOpt<Real>::setParameter(param);
    }
  }
}; // class CompositeConstraint_SimOpt

} // namespace ROL

#endif
