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
#include "ROL_EqualityConstraint_SimOpt.hpp"
#include "ROL_Vector_SimOpt.hpp"

/** @ingroup func_group
    \class ROL::Reduced_Objective_SimOpt
    \brief Provides the interface to evaluate simulation-based reduced objective functions.

           The reduced simulation-based objective function is 
           \f$\widehat{J}(z) = J(u(z),z)\f$ where \f$u(z)=u\f$ solves \f$c(u,z) = 0\f$.
*/


namespace ROL {

template <class Real>
class Reduced_Objective_SimOpt : public Objective<Real> {
private:
  Teuchos::RCP<Objective_SimOpt<Real> > obj_;           ///< SimOpt objective function. 
  Teuchos::RCP<EqualityConstraint_SimOpt<Real> > con_;  ///< SimOpt equality constraint.
  Teuchos::RCP<Vector<Real> > state_;                   ///< Storage for the state variable.
  Teuchos::RCP<Vector<Real> > adjoint_;                 ///< Storage for the adjoint variable.

  bool storage_;                                        ///< Flag whether or not to store computed quantities.
  bool is_state_computed_;                              ///< Flag whether or not to store the state variable.
  bool is_adjoint_computed_;                            ///< Flag whether or not to store the adjoint variable.

  /** \brief Given \f$z\in\mathcal{Z}\f$, solve the state equation \f$c(u,z) = 0\f$ for 
      \f$u=u(z)\in\mathcal{U}\f$.
  */
  void solve_state_equation(const ROL::Vector<Real> &x, Real &tol, bool flag = true, int iter = -1) { 
    // Solve state equation if not done already
    if (!this->is_state_computed_ || !this->storage_) {
      (this->con_)->solve(*(this->state_),x,tol);
      // Update full objective function
      (this->obj_)->update(*(this->state_),x,flag,iter);
      // Update equality constraint
      (this->con_)->update(*(this->state_),x,flag,iter);
      // Reset storage flags
      this->is_state_computed_ = true;
    }
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation, solve 
             the adjoint equation 
             \f$c_u(u,z)^*\lambda + J_u(u,z) = 0\f$ for \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_equation(const ROL::Vector<Real> &x, Real &tol) { 
    // Solve state equation if not done already
    this->solve_state_equation(x,tol);
    // Solve adjoint equation if not done already
    if(!this->is_adjoint_computed_ || !this->storage_) {
      // Evaluate the full gradient wrt u
      Teuchos::RCP<Vector<Real> > gu = (this->state_)->clone();
      (this->obj_)->gradient_1(*gu,*(this->state_),x,tol);
      // Solve adjoint equation
      (this->con_)->applyInverseAdjointJacobian_1(*(this->adjoint_),*gu,*(this->state_),x,tol);
      (this->adjoint_)->scale(-1.0);
      // Reset storage flags
      this->is_adjoint_computed_ = true;
    }
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation and 
             a direction \f$v\in\mathcal{Z}\f$, solve the state senstivity equation 
             \f$c_u(u,z)s + c_z(u,z)v = 0\f$ for \f$s=u_z(z)v\in\mathcal{U}\f$.
  */
  void solve_state_sensitivity(ROL::Vector<Real> &s, const ROL::Vector<Real> &v, 
                               const ROL::Vector<Real> &x, Real &tol) {
    // Solve state equation if not done already
    this->solve_state_equation(x,tol);
    // Solve state sensitivity equation
    Teuchos::RCP<Vector<Real> > Bv = (this->state_)->clone();
    (this->con_)->applyJacobian_2(*Bv,v,*(this->state_),x,tol);
    Bv->scale(-1.0);
    (this->con_)->applyInverseJacobian_1(s,*Bv,*(this->state_),x,tol);
  }

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$, the adjoint variable 
             \f$\lambda\in\mathcal{C}^*\f$, and a direction \f$v\in\mathcal{Z}\f$, solve the 
             adjoint sensitvity equation 
             \f$c_u(u,z)^*p + J_{uu}(u,z)s + J_{uz}(u,z)v + c_{uu}(u,z)(\cdot,s)^*\lambda 
                            + c_{zu}(u,z)(\cdot,v)^*\lambda = 0\f$
             for \f$p = \lambda_z(u(z),z)v\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_sensitivity(ROL::Vector<Real> &p, const ROL::Vector<Real> &s,
                                 const ROL::Vector<Real> &v, const ROL::Vector<Real> &x, Real &tol) {
    // Solve state equation if not done already
    this->solve_state_equation(x,tol);
    // Solve adjoint equation if not done already
    this->solve_adjoint_equation(x,tol);
    // Evaluate full hessVec in the direction (s,v)
    Teuchos::RCP<Vector<Real> > hv11 = (this->state_)->clone();
    (this->obj_)->hessVec_11(*hv11,s,*(this->state_),x,tol);
    Teuchos::RCP<Vector<Real> > hv12 = (this->state_)->clone();
    (this->obj_)->hessVec_12(*hv12,v,*(this->state_),x,tol);
    // Apply adjoint Hessian of constraint
    Teuchos::RCP<Vector<Real> > hc11 = (this->state_)->clone();
    (this->con_)->applyAdjointHessian_11(*hc11,*(this->adjoint_),s,*(this->state_),x,tol);
    Teuchos::RCP<Vector<Real> > hc21 = (this->state_)->clone();
    (this->con_)->applyAdjointHessian_21(*hc21,*(this->adjoint_),v,*(this->state_),x,tol);
    // Solve adjoint sensitivity equation
    Teuchos::RCP<Vector<Real> > r = (this->state_)->clone();
    r->set(*hv11);
    r->plus(*hv12);
    r->plus(*hc11);
    r->plus(*hc21);
    r->scale(-1.0);
    (this->con_)->applyInverseAdjointJacobian_1(p,*r,*(this->state_),x,tol);
  }

public:
  /** \brief Constructor.

      @param[in] obj     is a pointer to a SimOpt objective function.
      @param[in] con     is a pointer to a SimOpt equality constraint.
      @param[in] state   is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] adjoint is a pointer to a dual constraint space vector, \f$\mathcal{C}^*\f$.
      @param[in] storage is a flag whether or not to store computed states and adjoints.
  */
  Reduced_Objective_SimOpt(Teuchos::RCP<Objective_SimOpt<Real> > &obj, 
                           Teuchos::RCP<EqualityConstraint_SimOpt<Real> > &con, 
                           Teuchos::RCP<Vector<Real> > &state, 
                           Teuchos::RCP<Vector<Real> > &adjoint,
                           bool storage = true) 
    : obj_(obj), con_(con), state_(state), adjoint_(adjoint), storage_(storage) {
    is_state_computed_   = false;
    is_adjoint_computed_ = false;
  }

  /** \brief Update the SimOpt objective function and equality constraint.
  */
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    // Reset storage flags
    this->is_state_computed_ = false;
    this->is_adjoint_computed_ = false;
    // Solve state equation
    Real tol = std::sqrt(ROL_EPSILON);
    this->solve_state_equation(x,tol,flag,iter);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the objective function 
             \f$\widehat{J}(z) = J(u(z),z)\f$ where 
             \f$u=u(z)\in\mathcal{U}\f$ solves \f$e(u,z) = 0\f$.
  */
  Real value( const Vector<Real> &x, Real &tol ) {
    // Solve state equation
    this->solve_state_equation(x,tol);
    // Get objective function value
    return (this->obj_)->value(*(this->state_),x,tol);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the gradient of the objective function 
             \f$\nabla\widehat{J}(z) = J_z(z) + c_z(u,z)^*\lambda\f$ where 
             \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$ solves 
             \f$e_u(u,z)^*\lambda+J_u(u,z) = 0\f$.
  */
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    // Solve state equation
    this->solve_state_equation(x,tol);
    // Solve adjoint equation
    this->solve_adjoint_equation(x,tol);
    // Evaluate the full gradient wrt z
    Teuchos::RCP<Vector<Real> > gz = x.clone();
    (this->obj_)->gradient_2(*gz,*(this->state_),x,tol);
    // Build gradient
    (this->con_)->applyAdjointJacobian_2(g,*(this->adjoint_),*(this->state_),x,tol);
    g.plus(*gz);
  }

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the Hessian of the objective function 
             \f$\nabla^2\widehat{J}(z)\f$ in the direction \f$v\in\mathcal{Z}\f$.
  */
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    // Solve state equation
    this->solve_state_equation(x,tol);
    // Solve adjoint equation
    this->solve_adjoint_equation(x,tol);
    // Solve state sensitivity equation
    Teuchos::RCP<Vector<Real> > s = (this->state_)->clone();
    this->solve_state_sensitivity(*s,v,x,tol);
    // Solve adjoint sensitivity equation
    Teuchos::RCP<Vector<Real> > p = (this->adjoint_)->clone();
    this->solve_adjoint_sensitivity(*p,*s,v,x,tol);
    // Build hessVec
    (this->con_)->applyAdjointJacobian_2(hv,*p,*(this->state_),x,tol);
    Teuchos::RCP<Vector<Real> > tmp = x.clone();
    (this->obj_)->hessVec_21(*tmp,*s,*(this->state_),x,tol);
    hv.plus(*tmp);
    (this->obj_)->hessVec_22(*tmp,v,*(this->state_),x,tol);
    hv.plus(*tmp);
    (this->con_)->applyAdjointHessian_12(*tmp,*(this->adjoint_),*s,*(this->state_),x,tol);
    hv.plus(*tmp);
    (this->con_)->applyAdjointHessian_22(*tmp,*(this->adjoint_),v,*(this->state_),x,tol);
    hv.plus(*tmp);
  }

  /** \brief Apply a reduced Hessian preconditioner.
  */
  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Pv.set(v);
  }

}; // class ROL::Reduced_Objective_SimOpt

} // namespace ROL

#endif
