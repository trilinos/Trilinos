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


#ifndef ROL_REDUCED_OBJECTIVE_SIMOPT_H
#define ROL_REDUCED_OBJECTIVE_SIMOPT_H

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_SimController.hpp"

namespace ROL {

template <class Real>
class Reduced_Objective_SimOpt : public Objective<Real> {
private:
  const ROL::Ptr<Objective_SimOpt<Real> > obj_;          
  const ROL::Ptr<Constraint_SimOpt<Real> > con_; 
  ROL::Ptr<SimController<Real> > stateStore_;
  ROL::Ptr<SimController<Real> > adjointStore_;

  // Primal vectors
  ROL::Ptr<Vector<Real> > state_;                              
  ROL::Ptr<Vector<Real> > adjoint_;                            
  ROL::Ptr<Vector<Real> > state_sens_;                              
  ROL::Ptr<Vector<Real> > adjoint_sens_;                            

  // Dual vectors
  ROL::Ptr<Vector<Real> > dualstate_;
  ROL::Ptr<Vector<Real> > dualstate1_;
  ROL::Ptr<Vector<Real> > dualadjoint_;
  ROL::Ptr<Vector<Real> > dualcontrol_;

  const bool storage_;             
  const bool useFDhessVec_;
  
  bool updateFlag_;
  int  updateIter_;

  void solve_state_equation(const Vector<Real> &z, Real &tol) { 
    // Check if state has been computed.
    bool isComputed = false;
    if (storage_) {
      isComputed = stateStore_->get(*state_,Objective<Real>::getParameter());
    }
    // Solve state equation if not done already.
    if (!isComputed || !storage_) {
      // Update equality constraint with new Opt variable.
      con_->update_2(z,updateFlag_,updateIter_);
      // Solve state equation.
      con_->solve(*dualadjoint_,*state_,z,tol);
      // Update equality constraint with new Sim variable.
      con_->update_1(*state_,updateFlag_,updateIter_);
      // Update full objective function.
      obj_->update(*state_,z,updateFlag_,updateIter_);
      // Store state.
      if (storage_) {
        stateStore_->set(*state_,Objective<Real>::getParameter());
      }
    }
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation, solve 
             the adjoint equation 
             \f$c_u(u,z)^*\lambda + J_u(u,z) = 0\f$ for \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_equation(const Vector<Real> &z, Real &tol) { 
    // Check if adjoint has been computed.
    bool isComputed = false;
    if (storage_) {
      isComputed = adjointStore_->get(*adjoint_,Objective<Real>::getParameter());
    }
    // Solve adjoint equation if not done already.
    if (!isComputed || !storage_) {
      // Evaluate the full gradient wrt u
      obj_->gradient_1(*dualstate_,*state_,z,tol);
      // Solve adjoint equation
      con_->applyInverseAdjointJacobian_1(*adjoint_,*dualstate_,*state_,z,tol);
      adjoint_->scale(static_cast<Real>(-1));
      // Store adjoint
      if (storage_) {
        adjointStore_->set(*adjoint_,Objective<Real>::getParameter());
      }
    }
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation and 
             a direction \f$v\in\mathcal{Z}\f$, solve the state senstivity equation 
             \f$c_u(u,z)s + c_z(u,z)v = 0\f$ for \f$s=u_z(z)v\in\mathcal{U}\f$.
  */
  void solve_state_sensitivity(const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
    // Solve state sensitivity equation
    con_->applyJacobian_2(*dualadjoint_,v,*state_,z,tol);
    dualadjoint_->scale(static_cast<Real>(-1));
    con_->applyInverseJacobian_1(*state_sens_,*dualadjoint_,*state_,z,tol);
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$, the adjoint variable 
             \f$\lambda\in\mathcal{C}^*\f$, and a direction \f$v\in\mathcal{Z}\f$, solve the 
             adjoint sensitvity equation 
             \f$c_u(u,z)^*p + J_{uu}(u,z)s + J_{uz}(u,z)v + c_{uu}(u,z)(\cdot,s)^*\lambda 
                            + c_{zu}(u,z)(\cdot,v)^*\lambda = 0\f$
             for \f$p = \lambda_z(u(z),z)v\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_sensitivity(const Vector<Real> &v, const Vector<Real> &z, Real &tol) {
    // Evaluate full hessVec in the direction (s,v)
    obj_->hessVec_11(*dualstate_,*state_sens_,*state_,z,tol);
    obj_->hessVec_12(*dualstate1_,v,*state_,z,tol);
    dualstate_->plus(*dualstate1_);
    // Apply adjoint Hessian of constraint
    con_->applyAdjointHessian_11(*dualstate1_,*adjoint_,*state_sens_,*state_,z,tol);
    dualstate_->plus(*dualstate1_);
    con_->applyAdjointHessian_21(*dualstate1_,*adjoint_,v,*state_,z,tol);
    dualstate_->plus(*dualstate1_);
    // Solve adjoint sensitivity equation
    dualstate_->scale(static_cast<Real>(-1));
    con_->applyInverseAdjointJacobian_1(*adjoint_sens_,*dualstate_,*state_,z,tol);
  }

public:
  /** \brief Constructor.

      @param[in] obj          is a pointer to a SimOpt objective function.
      @param[in] con          is a pointer to a SimOpt equality constraint.
      @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
      @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
      @param[in] storage      is a flag whether or not to store computed states and adjoints.
      @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
  */
  Reduced_Objective_SimOpt(
      const ROL::Ptr<Objective_SimOpt<Real> > &obj, 
      const ROL::Ptr<Constraint_SimOpt<Real> > &con, 
      const ROL::Ptr<Vector<Real> > &state, 
      const ROL::Ptr<Vector<Real> > &control, 
      const ROL::Ptr<Vector<Real> > &adjoint,
      const bool storage = true,
      const bool useFDhessVec = false) 
    : obj_(obj), con_(con),
      storage_(storage), useFDhessVec_(useFDhessVec),
      updateFlag_(true), updateIter_(0) {
    stateStore_   = ROL::makePtr<SimController<Real>>();
    adjointStore_ = ROL::makePtr<SimController<Real>>();
    state_        = state->clone();
    adjoint_      = adjoint->clone();
    state_sens_   = state->clone();
    adjoint_sens_ = adjoint->clone();
    dualstate_    = state->dual().clone();
    dualstate1_   = state->dual().clone();
    dualadjoint_  = adjoint->dual().clone();
    dualcontrol_  = control->dual().clone();
  }

  /** \brief Secondary, general constructor for use with dual optimization vector spaces
             where the user does not define the dual() method.

      @param[in] obj          is a pointer to a SimOpt objective function.
      @param[in] con          is a pointer to a SimOpt equality constraint.
      @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
      @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
      @param[in] dualstate    is a pointer to a dual state space vector, \f$\mathcal{U}^*\f$.
      @param[in] dualadjoint  is a pointer to a constraint space vector, \f$\mathcal{C}\f$.
      @param[in] storage      is a flag whether or not to store computed states and adjoints.
      @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
  */
  Reduced_Objective_SimOpt(
      const ROL::Ptr<Objective_SimOpt<Real> > &obj,
      const ROL::Ptr<Constraint_SimOpt<Real> > &con,
      const ROL::Ptr<Vector<Real> > &state,
      const ROL::Ptr<Vector<Real> > &control, 
      const ROL::Ptr<Vector<Real> > &adjoint,
      const ROL::Ptr<Vector<Real> > &dualstate,
      const ROL::Ptr<Vector<Real> > &dualcontrol, 
      const ROL::Ptr<Vector<Real> > &dualadjoint,
      const bool storage = true,
      const bool useFDhessVec = false)
    : obj_(obj), con_(con),
      storage_(storage), useFDhessVec_(useFDhessVec),
      updateFlag_(true), updateIter_(0) {
    stateStore_   = ROL::makePtr<SimController<Real>>();
    adjointStore_ = ROL::makePtr<SimController<Real>>();
    state_        = state->clone();
    adjoint_      = adjoint->clone();
    state_sens_   = state->clone();
    adjoint_sens_ = adjoint->clone();
    dualstate_    = dualstate->clone();
    dualstate1_   = dualstate->clone();
    dualadjoint_  = dualadjoint->clone();
    dualcontrol_  = dualcontrol->clone();
  }

  /** \brief Constructor.

      @param[in] obj          is a pointer to a SimOpt objective function.
      @param[in] con          is a pointer to a SimOpt equality constraint.
      @param[in] stateStore   is a pointer to a SimController object.
      @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
      @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
      @param[in] storage      is a flag whether or not to store computed states and adjoints.
      @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
  */
  Reduced_Objective_SimOpt(
      const ROL::Ptr<Objective_SimOpt<Real> > &obj, 
      const ROL::Ptr<Constraint_SimOpt<Real> > &con, 
      const ROL::Ptr<SimController<Real> > &stateStore, 
      const ROL::Ptr<Vector<Real> > &state, 
      const ROL::Ptr<Vector<Real> > &control, 
      const ROL::Ptr<Vector<Real> > &adjoint,
      const bool storage = true,
      const bool useFDhessVec = false) 
    : obj_(obj), con_(con), stateStore_(stateStore),
      storage_(storage), useFDhessVec_(useFDhessVec),
      updateFlag_(true), updateIter_(0) {
    adjointStore_ = ROL::makePtr<SimController<Real>>();
    state_        = state->clone();
    adjoint_      = adjoint->clone();
    state_sens_   = state->clone();
    adjoint_sens_ = adjoint->clone();
    dualstate_    = state->dual().clone();
    dualstate1_   = state->dual().clone();
    dualadjoint_  = adjoint->dual().clone();
    dualcontrol_  = control->dual().clone();
  }

  /** \brief Secondary, general constructor for use with dual optimization vector spaces
             where the user does not define the dual() method.

      @param[in] obj          is a pointer to a SimOpt objective function.
      @param[in] con          is a pointer to a SimOpt equality constraint.
      @param[in] stateStore   is a pointer to a SimController object.
      @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
      @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
      @param[in] dualstate    is a pointer to a dual state space vector, \f$\mathcal{U}^*\f$.
      @param[in] dualadjoint  is a pointer to a constraint space vector, \f$\mathcal{C}\f$.
      @param[in] storage      is a flag whether or not to store computed states and adjoints.
      @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
  */
  Reduced_Objective_SimOpt(
      const ROL::Ptr<Objective_SimOpt<Real> > &obj,
      const ROL::Ptr<Constraint_SimOpt<Real> > &con,
      const ROL::Ptr<SimController<Real> > &stateStore, 
      const ROL::Ptr<Vector<Real> > &state,
      const ROL::Ptr<Vector<Real> > &control, 
      const ROL::Ptr<Vector<Real> > &adjoint,
      const ROL::Ptr<Vector<Real> > &dualstate,
      const ROL::Ptr<Vector<Real> > &dualcontrol, 
      const ROL::Ptr<Vector<Real> > &dualadjoint,
      const bool storage = true,
      const bool useFDhessVec = false)
    : obj_(obj), con_(con), stateStore_(stateStore),
      storage_(storage), useFDhessVec_(useFDhessVec),
      updateFlag_(true), updateIter_(0) {
    adjointStore_ = ROL::makePtr<SimController<Real>>();
    state_        = state->clone();
    adjoint_      = adjoint->clone();
    state_sens_   = state->clone();
    adjoint_sens_ = adjoint->clone();
    dualstate_    = dualstate->clone();
    dualstate1_   = dualstate->clone();
    dualadjoint_  = dualadjoint->clone();
    dualcontrol_  = dualcontrol->clone();
  }

  /** \brief Update the SimOpt objective function and equality constraint.
  */
  void update( const Vector<Real> &z, bool flag = true, int iter = -1 ) {
    updateFlag_ = flag;
    updateIter_ = iter;
    stateStore_->objectiveUpdate(true);
    adjointStore_->objectiveUpdate(flag);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the objective function 
             \f$\widehat{J}(z) = J(u(z),z)\f$ where 
             \f$u=u(z)\in\mathcal{U}\f$ solves \f$e(u,z) = 0\f$.
  */
  Real value( const Vector<Real> &z, Real &tol ) {
    // Solve state equation
    solve_state_equation(z,tol);
    // Get objective function value
    return obj_->value(*state_,z,tol);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the gradient of the objective function 
             \f$\nabla\widehat{J}(z) = J_z(z) + c_z(u,z)^*\lambda\f$ where 
             \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$ solves 
             \f$e_u(u,z)^*\lambda+J_u(u,z) = 0\f$.
  */
  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
    // Solve state equation
    solve_state_equation(z,tol);
    // Solve adjoint equation
    solve_adjoint_equation(z,tol);
    // Evaluate the full gradient wrt z
    obj_->gradient_2(*dualcontrol_,*state_,z,tol);
    // Build gradient
    con_->applyAdjointJacobian_2(g,*adjoint_,*state_,z,tol);
    g.plus(*dualcontrol_);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the Hessian of the objective function 
             \f$\nabla^2\widehat{J}(z)\f$ in the direction \f$v\in\mathcal{Z}\f$.
  */
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
    if ( useFDhessVec_ ) {
      Objective<Real>::hessVec(hv,v,z,tol);
    }
    else {
      // Solve state equation
      solve_state_equation(z,tol);
      // Solve adjoint equation
      solve_adjoint_equation(z,tol);
      // Solve state sensitivity equation
      solve_state_sensitivity(v,z,tol);
      // Solve adjoint sensitivity equation
      solve_adjoint_sensitivity(v,z,tol);
      // Build hessVec
      con_->applyAdjointJacobian_2(hv,*adjoint_sens_,*state_,z,tol);
      obj_->hessVec_21(*dualcontrol_,*state_sens_,*state_,z,tol);
      hv.plus(*dualcontrol_);
      obj_->hessVec_22(*dualcontrol_,v,*state_,z,tol);
      hv.plus(*dualcontrol_);
      con_->applyAdjointHessian_12(*dualcontrol_,*adjoint_,*state_sens_,*state_,z,tol);
      hv.plus(*dualcontrol_);
      con_->applyAdjointHessian_22(*dualcontrol_,*adjoint_,v,*state_,z,tol);
      hv.plus(*dualcontrol_);
    }
  }

  /** \brief Apply a reduced Hessian preconditioner.
  */
  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
    Pv.set(v.dual());
  }

// For parametrized (stochastic) objective functions and constraints
public:
  void setParameter(const std::vector<Real> &param) {
    Objective<Real>::setParameter(param);
    con_->setParameter(param);
    obj_->setParameter(param);
  }
}; // class Reduced_Objective_SimOpt

} // namespace ROL

#endif
