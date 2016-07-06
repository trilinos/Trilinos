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


#ifndef ROL_REDUCED_PARAMETRIZEDOBJECTIVE_SIMOPT_H
#define ROL_REDUCED_PARAMETRIZEDOBJECTIVE_SIMOPT_H

#include "ROL_ParametrizedObjective_SimOpt.hpp"
#include "ROL_ParametrizedEqualityConstraint_SimOpt.hpp"
#include "ROL_Vector_SimOpt.hpp"

namespace ROL {

template <class Real>
class Reduced_ParametrizedObjective_SimOpt : public ParametrizedObjective<Real> {
private:
  Teuchos::RCP<ParametrizedObjective_SimOpt<Real> > obj_;          
  Teuchos::RCP<ParametrizedEqualityConstraint_SimOpt<Real> > con_; 

  // Primal vectors
  Teuchos::RCP<Vector<Real> > state_;                              
  Teuchos::RCP<Vector<Real> > state_sens_;                              
  Teuchos::RCP<Vector<Real> > adjoint_;                            
  Teuchos::RCP<Vector<Real> > adjoint_sens_;                            

  // Dual vectors
  Teuchos::RCP<Vector<Real> > dualstate_;
  Teuchos::RCP<Vector<Real> > dualstate1_;
  Teuchos::RCP<Vector<Real> > dualadjoint_;
  Teuchos::RCP<Vector<Real> > dualcontrol_;

  // Vector storage
  std::map<std::vector<Real>,Teuchos::RCP<Vector<Real> > > state_storage_; 
  std::map<std::vector<Real>,Teuchos::RCP<Vector<Real> > > adjoint_storage_;

  bool storage_;             
  bool useFDhessVec_;
  bool is_initialized_;

  void solve_state_equation(const Vector<Real> &x, Real &tol, bool flag = true, int iter = -1) { 
    // Solve state equation if not done already
    if ( state_storage_.count(this->getParameter()) && storage_ ) {
      state_->set(*state_storage_[this->getParameter()]);
    }
    else {
      // Update equality constraint
      con_->update_2(x,flag,iter);
      // Solve state
      con_->solve(*dualadjoint_,*state_,x,tol);
      // Update equality constraint
      con_->update_1(*state_,flag,iter);
      // Update full objective function
      obj_->update(*state_,x,flag,iter);
      // Store state
      if ( storage_ ) {
        state_storage_.insert(
          std::pair<std::vector<Real>,Teuchos::RCP<Vector<Real> > >(
            this->getParameter(),state_->clone()));
        state_storage_[this->getParameter()]->set(*state_);
      }
    }
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation, solve 
             the adjoint equation 
             \f$c_u(u,z)^*\lambda + J_u(u,z) = 0\f$ for \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_equation(const Vector<Real> &x, Real &tol) { 
    // Solve adjoint equation if not done already
    if ( adjoint_storage_.count(this->getParameter()) && storage_ ) {
      adjoint_->set(*adjoint_storage_[this->getParameter()]);
    }
    else {
      // Evaluate the full gradient wrt u
      obj_->gradient_1(*dualstate_,*state_,x,tol);
      // Solve adjoint equation
      con_->applyInverseAdjointJacobian_1(*adjoint_,*dualstate_,*state_,x,tol);
      adjoint_->scale(-1.0);
      // Store adjoint
      if ( storage_ ) {
        adjoint_storage_.insert(
          std::pair<std::vector<Real>,Teuchos::RCP<Vector<Real> > >(
            this->getParameter(),adjoint_->clone()));
        adjoint_storage_[this->getParameter()]->set(*adjoint_);
      }
    }
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation and 
             a direction \f$v\in\mathcal{Z}\f$, solve the state senstivity equation 
             \f$c_u(u,z)s + c_z(u,z)v = 0\f$ for \f$s=u_z(z)v\in\mathcal{U}\f$.
  */
  void solve_state_sensitivity(const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    // Solve state sensitivity equation
    con_->applyJacobian_2(*dualadjoint_,v,*state_,x,tol);
    dualadjoint_->scale(-1.0);
    con_->applyInverseJacobian_1(*state_sens_,*dualadjoint_,*state_,x,tol);
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$, the adjoint variable 
             \f$\lambda\in\mathcal{C}^*\f$, and a direction \f$v\in\mathcal{Z}\f$, solve the 
             adjoint sensitvity equation 
             \f$c_u(u,z)^*p + J_{uu}(u,z)s + J_{uz}(u,z)v + c_{uu}(u,z)(\cdot,s)^*\lambda 
                            + c_{zu}(u,z)(\cdot,v)^*\lambda = 0\f$
             for \f$p = \lambda_z(u(z),z)v\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_sensitivity(const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    // Evaluate full hessVec in the direction (s,v)
    obj_->hessVec_11(*dualstate_,*state_sens_,*state_,x,tol);
    obj_->hessVec_12(*dualstate1_,v,*state_,x,tol);
    dualstate_->plus(*dualstate1_);
    // Apply adjoint Hessian of constraint
    con_->applyAdjointHessian_11(*dualstate1_,*adjoint_,*state_sens_,*state_,x,tol);
    dualstate_->plus(*dualstate1_);
    con_->applyAdjointHessian_21(*dualstate1_,*adjoint_,v,*state_,x,tol);
    dualstate_->plus(*dualstate1_);
    // Solve adjoint sensitivity equation
    dualstate_->scale(-1.0);
    con_->applyInverseAdjointJacobian_1(*adjoint_sens_,*dualstate_,*state_,x,tol);
  }

public:
  /** \brief Constructor.

      @param[in] obj     is a pointer to a SimOpt objective function.
      @param[in] con     is a pointer to a SimOpt equality constraint.
      @param[in] state   is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] adjoint is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
      @param[in] storage is a flag whether or not to store computed states and adjoints.
  */
  Reduced_ParametrizedObjective_SimOpt(Teuchos::RCP<ParametrizedObjective_SimOpt<Real> > &obj, 
                                       Teuchos::RCP<ParametrizedEqualityConstraint_SimOpt<Real> > &con, 
                                       Teuchos::RCP<Vector<Real> > &state, 
                                       Teuchos::RCP<Vector<Real> > &adjoint,
                                       bool storage = true, bool useFDhessVec = false) 
    : obj_(obj), con_(con),
      state_(state), state_sens_(state->clone()),
      adjoint_(adjoint), adjoint_sens_(adjoint->clone()),
      dualstate_(state->dual().clone()), dualstate1_(state->dual().clone()),
      dualadjoint_(adjoint->dual().clone()), dualcontrol_(Teuchos::null),
      storage_(storage), useFDhessVec_(useFDhessVec), is_initialized_(false) {
    state_storage_.clear();
    adjoint_storage_.clear();
  }

  /** \brief Secondary, general constructor for use with dual optimization vector spaces
             where the user does not define the dual() method.

      @param[in] obj          is a pointer to a SimOpt objective function.
      @param[in] con          is a pointer to a SimOpt equality constraint.
      @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
      @param[in] dualstate    is a pointer to a dual state space vector, \f$\mathcal{U}^*\f$.
      @param[in] dualadjoint  is a pointer to a constraint space vector, \f$\mathcal{C}\f$.
      @param[in] storage      is a flag whether or not to store computed states and adjoints.
      @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
  */
  Reduced_ParametrizedObjective_SimOpt(Teuchos::RCP<ParametrizedObjective_SimOpt<Real> > &obj,
                                       Teuchos::RCP<ParametrizedEqualityConstraint_SimOpt<Real> > &con,
                                       Teuchos::RCP<Vector<Real> > &state,
                                       Teuchos::RCP<Vector<Real> > &adjoint,
                                       Teuchos::RCP<Vector<Real> > &dualstate,
                                       Teuchos::RCP<Vector<Real> > &dualadjoint,
                                       bool storage = true, bool useFDhessVec = false)
    : obj_(obj), con_(con),
      state_(state), state_sens_(state->clone()),
      adjoint_(adjoint), adjoint_sens_(adjoint->clone()),
      dualstate_(dualstate), dualstate1_(dualstate->clone()),
      dualadjoint_(dualadjoint->clone()), dualcontrol_(Teuchos::null),
      storage_(storage), useFDhessVec_(useFDhessVec), is_initialized_(false) {
    state_storage_.clear();
    adjoint_storage_.clear();
  }


  void setParameter(const std::vector<Real> &param) {
    ParametrizedObjective<Real>::setParameter(param);
    con_->setParameter(param);
    obj_->setParameter(param);
  }

  /** \brief Update the SimOpt objective function and equality constraint.
  */
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    // Reset storage flags
    state_storage_.clear();
    if ( flag ) {
      adjoint_storage_.clear();
    }
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the objective function 
             \f$\widehat{J}(z) = J(u(z),z)\f$ where 
             \f$u=u(z)\in\mathcal{U}\f$ solves \f$e(u,z) = 0\f$.
  */
  Real value( const Vector<Real> &x, Real &tol ) {
    // Solve state equation
    solve_state_equation(x,tol);
    // Get objective function value
    return obj_->value(*state_,x,tol);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the gradient of the objective function 
             \f$\nabla\widehat{J}(z) = J_z(z) + c_z(u,z)^*\lambda\f$ where 
             \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$ solves 
             \f$e_u(u,z)^*\lambda+J_u(u,z) = 0\f$.
  */
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    if (!is_initialized_) {
      dualcontrol_ = g.clone();
      is_initialized_ = true;
    }
    // Solve state equation
    solve_state_equation(x,tol);
    // Solve adjoint equation
    solve_adjoint_equation(x,tol);
    // Evaluate the full gradient wrt z
    obj_->gradient_2(*dualcontrol_,*state_,x,tol);
    // Build gradient
    con_->applyAdjointJacobian_2(g,*adjoint_,*state_,x,tol);
    g.plus(*dualcontrol_);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the Hessian of the objective function 
             \f$\nabla^2\widehat{J}(z)\f$ in the direction \f$v\in\mathcal{Z}\f$.
  */
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    if ( useFDhessVec_ ) {
      ParametrizedObjective<Real>::hessVec(hv,v,x,tol);
    }
    else {
      if (!is_initialized_) {
        dualcontrol_ = hv.clone();
        is_initialized_ = true;
      }
      // Solve state equation
      solve_state_equation(x,tol);
      // Solve adjoint equation
      solve_adjoint_equation(x,tol);
      // Solve state sensitivity equation
      solve_state_sensitivity(v,x,tol);
      // Solve adjoint sensitivity equation
      solve_adjoint_sensitivity(v,x,tol);
      // Build hessVec
      con_->applyAdjointJacobian_2(hv,*adjoint_sens_,*state_,x,tol);
      obj_->hessVec_21(*dualcontrol_,*state_sens_,*state_,x,tol);
      hv.plus(*dualcontrol_);
      obj_->hessVec_22(*dualcontrol_,v,*state_,x,tol);
      hv.plus(*dualcontrol_);
      con_->applyAdjointHessian_12(*dualcontrol_,*adjoint_,*state_sens_,*state_,x,tol);
      hv.plus(*dualcontrol_);
      con_->applyAdjointHessian_22(*dualcontrol_,*adjoint_,v,*state_,x,tol);
      hv.plus(*dualcontrol_);
    }
  }

  /** \brief Apply a reduced Hessian preconditioner.
  */
  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Pv.set(v.dual());
  }

}; // class Reduced_Objective_SimOpt

} // namespace ROL

#endif
