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
  Teuchos::RCP<Vector<Real> > state_;                              
  Teuchos::RCP<Vector<Real> > adjoint_;                            
  std::map<std::vector<Real>,Teuchos::RCP<Vector<Real> > > state_storage_; 
  std::map<std::vector<Real>,Teuchos::RCP<Vector<Real> > > adjoint_storage_;
  Teuchos::RCP<const Vector<Real> > dualstate_;         ///< Dual state vector, used for cloning only.
  Teuchos::RCP<const Vector<Real> > dualadjoint_;       ///< Dual adjoint vector, used for cloning only.

  bool storage_;             

  bool useFDhessVec_;

  void solve_state_equation(const Vector<Real> &x, Real &tol, bool flag = true, int iter = -1) { 
    // Solve state equation if not done already
    if ( state_storage_.count(this->getParameter()) ) {
      state_->set(*state_storage_[this->getParameter()]);
    }
    else {
      con_->solve(*state_,x,tol);
      // Update full objective function
      obj_->update(*state_,x,flag,iter);
      // Update equality constraint
      con_->update(*state_,x,flag,iter);
      // Store state
      Teuchos::RCP<Vector<Real> > tmp = state_->clone();
      state_storage_.insert(std::pair<std::vector<Real>,Teuchos::RCP<Vector<Real> > >(this->getParameter(),tmp));
      state_storage_[this->getParameter()]->set(*state_);
    }
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation, solve 
             the adjoint equation 
             \f$c_u(u,z)^*\lambda + J_u(u,z) = 0\f$ for \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_equation(const Vector<Real> &x, Real &tol) { 
    // Solve state equation if not done already
    solve_state_equation(x,tol);
    // Solve adjoint equation if not done already
    if ( adjoint_storage_.count(this->getParameter()) ) {
      adjoint_->set(*adjoint_storage_[this->getParameter()]);
    }
    else {
      // Evaluate the full gradient wrt u
      Teuchos::RCP<Vector<Real> > gu = dualstate_->clone();
      obj_->gradient_1(*gu,*state_,x,tol);
      // Solve adjoint equation
      con_->applyInverseAdjointJacobian_1(*adjoint_,*gu,*state_,x,tol);
      adjoint_->scale(-1.0);
      // Store adjoint
      Teuchos::RCP<Vector<Real> > tmp = adjoint_->clone();
      adjoint_storage_.insert(std::pair<std::vector<Real>,Teuchos::RCP<Vector<Real> > >(this->getParameter(),tmp));
      adjoint_storage_[this->getParameter()]->set(*adjoint_);
    }
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation and 
             a direction \f$v\in\mathcal{Z}\f$, solve the state senstivity equation 
             \f$c_u(u,z)s + c_z(u,z)v = 0\f$ for \f$s=u_z(z)v\in\mathcal{U}\f$.
  */
  void solve_state_sensitivity(Vector<Real> &s, const Vector<Real> &v, 
                               const Vector<Real> &x, Real &tol) {
    // Solve state equation if not done already
    solve_state_equation(x,tol);
    // Solve state sensitivity equation
    Teuchos::RCP<Vector<Real> > Bv = dualadjoint_->clone();
    con_->applyJacobian_2(*Bv,v,*state_,x,tol);
    Bv->scale(-1.0);
    con_->applyInverseJacobian_1(s,*Bv,*state_,x,tol);
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$, the adjoint variable 
             \f$\lambda\in\mathcal{C}^*\f$, and a direction \f$v\in\mathcal{Z}\f$, solve the 
             adjoint sensitvity equation 
             \f$c_u(u,z)^*p + J_{uu}(u,z)s + J_{uz}(u,z)v + c_{uu}(u,z)(\cdot,s)^*\lambda 
                            + c_{zu}(u,z)(\cdot,v)^*\lambda = 0\f$
             for \f$p = \lambda_z(u(z),z)v\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_sensitivity(Vector<Real> &p, const Vector<Real> &s,
                                 const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
    // Solve state equation if not done already
    solve_state_equation(x,tol);
    // Solve adjoint equation if not done already
    solve_adjoint_equation(x,tol);
    // Evaluate full hessVec in the direction (s,v)
    Teuchos::RCP<Vector<Real> > hv11 = dualstate_->clone();
    obj_->hessVec_11(*hv11,s,*state_,x,tol);
    Teuchos::RCP<Vector<Real> > hv12 = dualstate_->clone();
    obj_->hessVec_12(*hv12,v,*state_,x,tol);
    // Apply adjoint Hessian of constraint
    Teuchos::RCP<Vector<Real> > hc11 = dualstate_->clone();
    con_->applyAdjointHessian_11(*hc11,*adjoint_,s,*state_,x,tol);
    Teuchos::RCP<Vector<Real> > hc21 = dualstate_->clone();
    con_->applyAdjointHessian_21(*hc21,*adjoint_,v,*state_,x,tol);
    // Solve adjoint sensitivity equation
    Teuchos::RCP<Vector<Real> > r = dualstate_->clone();
    r->set(*hv11);
    r->plus(*hv12);
    r->plus(*hc11);
    r->plus(*hc21);
    r->scale(-1.0);
    con_->applyInverseAdjointJacobian_1(p,*r,*state_,x,tol);
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
    : obj_(obj), con_(con), state_(state), adjoint_(adjoint), storage_(storage), useFDhessVec_(useFDhessVec) {
    state_storage_.clear();
    adjoint_storage_.clear();
    dualstate_   = Teuchos::rcpFromRef(state_->dual());
    dualadjoint_ = Teuchos::rcpFromRef(adjoint_->dual());
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
  Reduced_ParametrizedObjective_SimOpt(Teuchos::RCP<Objective_SimOpt<Real> > &obj,
                                       Teuchos::RCP<EqualityConstraint_SimOpt<Real> > &con,
                                       Teuchos::RCP<Vector<Real> > &state,
                                       Teuchos::RCP<Vector<Real> > &adjoint,
                                       Teuchos::RCP<Vector<Real> > &dualstate,
                                       Teuchos::RCP<Vector<Real> > &dualadjoint,
                                       bool storage = true, bool useFDhessVec = false)
    : obj_(obj), con_(con),
      state_(state), adjoint_(adjoint), dualstate_(dualstate), dualadjoint_(dualadjoint),
      storage_(storage), useFDhessVec_(useFDhessVec) {
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
    // Solve state equation
    solve_state_equation(x,tol);
    // Solve adjoint equation
    solve_adjoint_equation(x,tol);
    // Evaluate the full gradient wrt z
    Teuchos::RCP<Vector<Real> > gz = x.clone();
    obj_->gradient_2(*gz,*state_,x,tol);
    // Build gradient
    con_->applyAdjointJacobian_2(g,*adjoint_,*state_,x,tol);
    g.plus(*gz);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the Hessian of the objective function 
             \f$\nabla^2\widehat{J}(z)\f$ in the direction \f$v\in\mathcal{Z}\f$.
  */
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    if ( useFDhessVec_ ) {
      ParametrizedObjective<Real>::hessVec(hv,v,x,tol);
    }
    else {
      // Solve state equation
      solve_state_equation(x,tol);
      // Solve adjoint equation
      solve_adjoint_equation(x,tol);
      // Solve state sensitivity equation
      Teuchos::RCP<Vector<Real> > s = state_->clone();
      solve_state_sensitivity(*s,v,x,tol);
      // Solve adjoint sensitivity equation
      Teuchos::RCP<Vector<Real> > p = adjoint_->clone();
      solve_adjoint_sensitivity(*p,*s,v,x,tol);
      // Build hessVec
      con_->applyAdjointJacobian_2(hv,*p,*state_,x,tol);
      Teuchos::RCP<Vector<Real> > tmp = x.clone();
      obj_->hessVec_21(*tmp,*s,*state_,x,tol);
      hv.plus(*tmp);
      obj_->hessVec_22(*tmp,v,*state_,x,tol);
      hv.plus(*tmp);
      con_->applyAdjointHessian_12(*tmp,*adjoint_,*s,*state_,x,tol);
      hv.plus(*tmp);
      con_->applyAdjointHessian_22(*tmp,*adjoint_,v,*state_,x,tol);
      hv.plus(*tmp);
    }
  }

  /** \brief Apply a reduced Hessian preconditioner.
  */
  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Pv.set(v);
  }

}; // class Reduced_Objective_SimOpt

} // namespace ROL

#endif
