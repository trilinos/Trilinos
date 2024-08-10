// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_REDUCED_CONSTRAINT_SIMOPT_H
#define ROL_REDUCED_CONSTRAINT_SIMOPT_H

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_VectorController.hpp"

namespace ROL {

template <class Real>
class Reduced_Constraint_SimOpt : public Constraint<Real> {
private:
  const ROL::Ptr<Constraint_SimOpt<Real>> conVal_;
  const ROL::Ptr<Constraint_SimOpt<Real>> conRed_;
  const ROL::Ptr<VectorController<Real>>  stateStore_;
  ROL::Ptr<VectorController<Real>>        adjointStore_;

  // Primal vectors
  ROL::Ptr<Vector<Real>> state_;
  ROL::Ptr<Vector<Real>> adjoint_;
  ROL::Ptr<Vector<Real>> residual_;
  ROL::Ptr<Vector<Real>> state_sens_;
  ROL::Ptr<Vector<Real>> adjoint_sens_;

  // Dual vectors
  ROL::Ptr<Vector<Real>> dualstate_;
  ROL::Ptr<Vector<Real>> dualstate1_;
  ROL::Ptr<Vector<Real>> dualadjoint_;
  ROL::Ptr<Vector<Real>> dualcontrol_;
  ROL::Ptr<Vector<Real>> dualresidual_;

  const bool storage_;
  const bool useFDhessVec_;

  bool updateFlag_;
  int  updateIter_;

  void solve_state_equation(const Vector<Real> &z, Real &tol) {
    // Check if state has been computed.
    bool isComputed = false;
    if (storage_) {
      isComputed = stateStore_->get(*state_,Constraint<Real>::getParameter());
    }
    // Solve state equation if not done already.
    if (!isComputed || !storage_) {
      // Update equality constraint with new Opt variable.
      conRed_->update_2(z,updateFlag_,updateIter_);
      // Solve state equation.
      conRed_->solve(*dualadjoint_,*state_,z,tol);
      // Update equality constraint with new Sim variable.
      conRed_->update_1(*state_,updateFlag_,updateIter_);
      // Update full objective function.
      conVal_->update(*state_,z,updateFlag_,updateIter_);
      // Store state.
      if (storage_) {
        stateStore_->set(*state_,Constraint<Real>::getParameter());
      }
    }
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation, solve 
             the adjoint equation 
             \f$c_u(u,z)^*\lambda + c_u(u,z)^*w = 0\f$ for \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_equation(const Vector<Real> &w, const Vector<Real> &z, Real &tol) { 
    // Check if adjoint has been computed.
    bool isComputed = false;
    if (storage_) {
      isComputed = adjointStore_->get(*adjoint_,Constraint<Real>::getParameter());
    }
    // Solve adjoint equation if not done already.
    if (!isComputed || !storage_) {
      // Evaluate the full gradient wrt u
      conVal_->applyAdjointJacobian_1(*dualstate_,w,*state_,z,tol);
      // Solve adjoint equation
      conRed_->applyInverseAdjointJacobian_1(*adjoint_,*dualstate_,*state_,z,tol);
      adjoint_->scale(static_cast<Real>(-1));
      // Store adjoint
      if (storage_) {
        adjointStore_->set(*adjoint_,Constraint<Real>::getParameter());
      }
    }
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation and 
             a direction \f$v\in\mathcal{Z}\f$, solve the state senstivity equation 
             \f$c_u(u,z)s + c_z(u,z)v = 0\f$ for \f$s=u_z(z)v\in\mathcal{U}\f$.
  */
  void solve_state_sensitivity(const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
    // Solve state sensitivity equation
    conRed_->applyJacobian_2(*dualadjoint_,v,*state_,z,tol);
    dualadjoint_->scale(static_cast<Real>(-1));
    conRed_->applyInverseJacobian_1(*state_sens_,*dualadjoint_,*state_,z,tol);
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$, the adjoint variable 
             \f$\lambda\in\mathcal{C}^*\f$, and a direction \f$v\in\mathcal{Z}\f$, solve the 
             adjoint sensitvity equation 
             \f$c_u(u,z)^*p + J_{uu}(u,z)s + J_{uz}(u,z)v + c_{uu}(u,z)(\cdot,s)^*\lambda 
                            + c_{zu}(u,z)(\cdot,v)^*\lambda = 0\f$
             for \f$p = \lambda_z(u(z),z)v\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_sensitivity(const Vector<Real> &w, const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
    // Evaluate full hessVec in the direction (s,v)
    conVal_->applyAdjointHessian_11(*dualstate_,w,*state_sens_,*state_,z,tol);
    conVal_->applyAdjointHessian_12(*dualstate1_,w,v,*state_,z,tol);
    dualstate_->plus(*dualstate1_);
    // Apply adjoint Hessian of constraint
    conRed_->applyAdjointHessian_11(*dualstate1_,*adjoint_,*state_sens_,*state_,z,tol);
    dualstate_->plus(*dualstate1_);
    conRed_->applyAdjointHessian_21(*dualstate1_,*adjoint_,v,*state_,z,tol);
    dualstate_->plus(*dualstate1_);
    // Solve adjoint sensitivity equation
    dualstate_->scale(static_cast<Real>(-1));
    conRed_->applyInverseAdjointJacobian_1(*adjoint_sens_,*dualstate_,*state_,z,tol);
  }

public:
  /** \brief Constructor.

      @param[in] obj          is a pointer to a SimOpt objective function.
      @param[in] con          is a pointer to a SimOpt equality constraint.
      @param[in] stateStore   is a pointer to a VectorController object.
      @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
      @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
      @param[in] storage      is a flag whether or not to store computed states and adjoints.
      @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
  */
  Reduced_Constraint_SimOpt(
      const ROL::Ptr<Constraint_SimOpt<Real> > &conVal, 
      const ROL::Ptr<Constraint_SimOpt<Real> > &conRed, 
      const ROL::Ptr<VectorController<Real> > &stateStore,
      const ROL::Ptr<Vector<Real> > &state, 
      const ROL::Ptr<Vector<Real> > &control, 
      const ROL::Ptr<Vector<Real> > &adjoint,
      const ROL::Ptr<Vector<Real> > &residual,
      const bool storage = true,
      const bool useFDhessVec = false) 
    : conVal_(conVal), conRed_(conRed), stateStore_(stateStore),
      storage_(storage), useFDhessVec_(useFDhessVec),
      updateFlag_(true), updateIter_(0) {
    adjointStore_ = ROL::makePtr<VectorController<Real>>();
    state_        = state->clone();
    adjoint_      = adjoint->clone();
    residual_     = residual->clone();
    state_sens_   = state->clone();
    adjoint_sens_ = adjoint->clone();
    dualstate_    = state->dual().clone();
    dualstate1_   = state->dual().clone();
    dualadjoint_  = adjoint->dual().clone();
    dualcontrol_  = control->dual().clone();
    dualresidual_ = residual->dual().clone();
  }

  /** \brief Secondary, general constructor for use with dual optimization vector spaces
             where the user does not define the dual() method.

      @param[in] obj          is a pointer to a SimOpt objective function.
      @param[in] con          is a pointer to a SimOpt equality constraint.
      @param[in] stateStore   is a pointer to a VectorController object.
      @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
      @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
      @param[in] dualstate    is a pointer to a dual state space vector, \f$\mathcal{U}^*\f$.
      @param[in] dualadjoint  is a pointer to a constraint space vector, \f$\mathcal{C}\f$.
      @param[in] storage      is a flag whether or not to store computed states and adjoints.
      @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
  */
  Reduced_Constraint_SimOpt(
      const ROL::Ptr<Constraint_SimOpt<Real> > &conVal, 
      const ROL::Ptr<Constraint_SimOpt<Real> > &conRed,
      const ROL::Ptr<VectorController<Real> > &stateStore, 
      const ROL::Ptr<Vector<Real> > &state,
      const ROL::Ptr<Vector<Real> > &control, 
      const ROL::Ptr<Vector<Real> > &adjoint,
      const ROL::Ptr<Vector<Real> > &residual,
      const ROL::Ptr<Vector<Real> > &dualstate,
      const ROL::Ptr<Vector<Real> > &dualcontrol, 
      const ROL::Ptr<Vector<Real> > &dualadjoint,
      const ROL::Ptr<Vector<Real> > &dualresidual,
      const bool storage = true,
      const bool useFDhessVec = false)
    : conVal_(conVal), conRed_(conRed), stateStore_(stateStore),
      storage_(storage), useFDhessVec_(useFDhessVec),
      updateFlag_(true), updateIter_(0) {
    adjointStore_ = ROL::makePtr<VectorController<Real>>();
    state_        = state->clone();
    adjoint_      = adjoint->clone();
    residual_     = residual->clone();
    state_sens_   = state->clone();
    adjoint_sens_ = adjoint->clone();
    dualstate_    = dualstate->clone();
    dualstate1_   = dualstate->clone();
    dualadjoint_  = dualadjoint->clone();
    dualcontrol_  = dualcontrol->clone();
    dualresidual_ = dualresidual->clone();
  }

  /** \brief Update the SimOpt objective function and equality constraint.
  */
  void update( const Vector<Real> &z, bool flag = true, int iter = -1 ) {
    updateFlag_ = flag;
    updateIter_ = iter;
    stateStore_->constraintUpdate(true);
    adjointStore_->constraintUpdate(flag);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the equality constraint 
             \f$\widehat{c}(z) = c(u(z),z)\f$ where 
             \f$u=u(z)\in\mathcal{U}\f$ solves \f$e(u,z) = 0\f$.
  */
  void value( Vector<Real> &c, const Vector<Real> &z, Real &tol ) {
    // Solve state equation
    solve_state_equation(z,tol);
    // Get constraint value
    conVal_->value(c,*state_,z,tol);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, apply the Jacobian to a vector
             \f$\widehat{c}'(z)v = c_u(u,z)s + c_z(u,z)v\f$ where 
             \f$s=s(u,z,v)\in\mathcal{U}^*\f$ solves 
             \f$e_u(u,z)s+e_z(u,z)v = 0\f$.
  */
  void applyJacobian( Vector<Real> &jv, const Vector<Real> &v,
                      const Vector<Real> &z, Real &tol ) {
    // Solve state equation.
    solve_state_equation(z,tol);
    // Solve state sensitivity equation.
    solve_state_sensitivity(v,z,tol);
    // Apply Sim Jacobian to state sensitivity.
    conVal_->applyJacobian_1(*residual_,*state_sens_,*state_,z,tol);
    // Apply Opt Jacobian to vector.
    conVal_->applyJacobian_2(jv,v,*state_,z,tol);
    jv.plus(*residual_);
  }

  void applyAdjointJacobian( Vector<Real> &ajw, const Vector<Real> &w,
                             const Vector<Real> &z, Real &tol ) {
    // Solve state equation
    solve_state_equation(z,tol);
    // Solve adjoint equation
    solve_adjoint_equation(w,z,tol);
    // Evaluate the full gradient wrt z
    conVal_->applyAdjointJacobian_2(*dualcontrol_,w,*state_,z,tol);
    // Build gradient
    conRed_->applyAdjointJacobian_2(ajw,*adjoint_,*state_,z,tol);
    ajw.plus(*dualcontrol_);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the Hessian of the objective function 
             \f$\nabla^2\widehat{J}(z)\f$ in the direction \f$v\in\mathcal{Z}\f$.
  */
  void applyAdjointHessian( Vector<Real> &ahwv, const Vector<Real> &w,
                            const Vector<Real> &v, const Vector<Real> &z,
                            Real &tol ) {
    if ( useFDhessVec_ ) {
      Constraint<Real>::applyAdjointHessian(ahwv,w,v,z,tol);
    }
    else {
      // Solve state equation
      solve_state_equation(z,tol);
      // Solve adjoint equation
      solve_adjoint_equation(w,z,tol);
      // Solve state sensitivity equation
      solve_state_sensitivity(v,z,tol);
      // Solve adjoint sensitivity equation
      solve_adjoint_sensitivity(w,v,z,tol);
      // Build hessVec
      conRed_->applyAdjointJacobian_2(ahwv,*adjoint_sens_,*state_,z,tol);
      conVal_->applyAdjointHessian_21(*dualcontrol_,w,*state_sens_,*state_,z,tol);
      ahwv.plus(*dualcontrol_);
      conVal_->applyAdjointHessian_22(*dualcontrol_,w,v,*state_,z,tol);
      ahwv.plus(*dualcontrol_);
      conRed_->applyAdjointHessian_12(*dualcontrol_,*adjoint_,*state_sens_,*state_,z,tol);
      ahwv.plus(*dualcontrol_);
      conRed_->applyAdjointHessian_22(*dualcontrol_,*adjoint_,v,*state_,z,tol);
      ahwv.plus(*dualcontrol_);
    }
  }

// For parametrized (stochastic) objective functions and constraints
public:
  void setParameter(const std::vector<Real> &param) {
    Constraint<Real>::setParameter(param);
    conVal_->setParameter(param);
    conRed_->setParameter(param);
  }
}; // class Reduced_Constraint_SimOpt

} // namespace ROL

#endif
