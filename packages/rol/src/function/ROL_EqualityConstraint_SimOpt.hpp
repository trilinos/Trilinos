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

#ifndef ROL_EQUALITY_CONSTRAINT_SIMOPT_H
#define ROL_EQUALITY_CONSTRAINT_SIMOPT_H

#include "ROL_EqualityConstraint.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_Types.hpp"
#include <iostream>

/** @ingroup func_group
    \class ROL::EqualityConstraint_SimOpt
    \brief Defines the equality constraint operator interface for simulation-based optimization.

    This equality constraint interface inherits from ROL_EqualityConstraint, for the
    use case when \f$\mathcal{X}=\mathcal{U}\times\mathcal{Z}\f$ where \f$\mathcal{U}\f$ and 
    \f$\mathcal{Z}\f$ are Banach spaces.  \f$\mathcal{U}\f$ denotes the "simulation space"
    and \f$\mathcal{Z}\f$ denotes the "optimization space" (of designs, controls, parameters).
    The simulation-based constraints are of the form
    \f[
      c(u,z) = 0 \,.
    \f]
    The basic operator interface, to be implemented by the user, requires:
    \li #value -- constraint evaluation.
    \li #applyJacobian_1        -- action of the partial constraint Jacobian --derivatives are 
                                   with respect to the first component \f$\mathcal{U}\f$;
    \li #applyJacobian_2        -- action of the partial constraint Jacobian --derivatives are 
                                   with respect to the second component \f$\mathcal{Z}\f$;
    \li #applyAdjointJacobian_1 -- action of the adjoint of the partial constraint Jacobian --derivatives are 
                                   with respect to the first component \f$\mathcal{U}\f$;
    \li #applyAdjointJacobian_2 -- action of the adjoint of the partial constraint Jacobian --derivatives are 
                                   with respect to the second component \f$\mathcal{Z}\f$;

    The user may also overload:
    \li #applyAdjointHessian_11  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the first component only;
    \li #applyAdjointHessian_12  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the first and second components;
    \li #applyAdjointHessian_21  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the second and first components;
    \li #applyAdjointHessian_22  -- action of the adjoint of the partial constraint Hessian --derivatives 
                                    are with respect to the second component only;
    \li #solveAugmentedSystem -- solution of the augmented system --the default is an iterative
                                 scheme based on the action of the Jacobian and its adjoint.
    \li #applyPreconditioner  -- action of a constraint preconditioner --the default is null-op.

    ---
*/


namespace ROL {

template <class Real>
class EqualityConstraint_SimOpt : public EqualityConstraint<Real> {
public:
  /** \brief Update constraint functions.  
                x is the optimization variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &u, const Vector<Real> &z, bool flag = true, int iter = -1 ) {}


  /** \brief Evaluate the constraint operator \f$c:\mathcal{U}\times\mathcal{Z} \rightarrow \mathcal{C}\f$
             at \f$(u,z)\f$.

             @param[out]      c   is the result of evaluating the constraint operator at @b \f$(u,z)\f$; a constraint-space vector
             @param[in]       u   is the constraint argument; a simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{c} = c(u,z)\f$,
             where \f$\mathsf{c} \in \mathcal{C}\f$, \f$\mathsf{u} \in \mathcal{U}\f$, and $\f$\mathsf{z} \in\mathcal{Z}\f$.

             ---
  */
  virtual void value(Vector<Real> &c,
                     const Vector<Real> &u,
                     const Vector<Real> &z,
                     Real &tol) = 0;

  /** \brief Given \f$z\f$, solve \f$c(u,z)=0\f$ for \f$u\f$.

             @param[out]      u   is the solution vector; a simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             ---
  */
  virtual void solve(Vector<Real> &u, 
                     const Vector<Real> &z,
                     Real &tol) {
    u.zero();
  }

  /** \brief Apply the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_u(u,z) \in L(\mathcal{U}, \mathcal{C})\f$,
             to the vector \f$v\f$.

             @param[out]      jv  is the result of applying the constraint Jacobian to @b v at @b \f$(u,z)\f$; a constraint-space vector
             @param[in]       v   is a simulation-space vector
             @param[in]       u   is the constraint argument; an simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{jv} = c_u(u,z)v\f$, where
             \f$v \in \mathcal{U}\f$, \f$\mathsf{jv} \in \mathcal{C}\f$.

             ---
  */
  virtual void applyJacobian_1(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) {
    Real ctol = std::sqrt(ROL_EPSILON);
    // Compute step length
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Update state vector to u + hv
    Teuchos::RCP<Vector<Real> > unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    // Compute new constraint value
    update(*unew,z);
    value(jv,*unew,z,ctol);
    // Compute current constraint value
    Teuchos::RCP<Vector<Real> > cold = jv.clone();
    update(u,z);
    value(*cold,u,z,ctol);
    // Compute Newton quotient
    jv.axpy(-1.0,*cold);
    jv.scale(1.0/h);
  }


  /** \brief Apply the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_z(u,z) \in L(\mathcal{Z}, \mathcal{C})\f$,
             to the vector \f$v\f$.

             @param[out]      jv  is the result of applying the constraint Jacobian to @b v at @b \f$(u,z)\f$; a constraint-space vector
             @param[in]       v   is an optimization-space vector
             @param[in]       u   is the constraint argument; a simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{jv} = c_z(u,z)v\f$, where
             \f$v \in \mathcal{Z}\f$, \f$\mathsf{jv} \in \mathcal{C}\f$. 

             ---
  */
  virtual void applyJacobian_2(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) { 
    Real ctol = std::sqrt(ROL_EPSILON);
    // Compute step length
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Update state vector to u + hv
    Teuchos::RCP<Vector<Real> > znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    // Compute new constraint value
    update(u,*znew);
    value(jv,u,*znew,ctol);
    // Compute current constraint value
    Teuchos::RCP<Vector<Real> > cold = jv.clone();
    update(u,z);
    value(*cold,u,z,ctol);
    // Compute Newton quotient
    jv.axpy(-1.0,*cold);
    jv.scale(1.0/h);
  }

  /** \brief Apply the inverse partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_u(u,z)^{-1} \in L(\mathcal{C}, \mathcal{U})\f$,
             to the vector \f$v\f$.

             @param[out]      ijv is the result of applying the inverse constraint Jacobian to @b v at @b \f$(u,z)\f$; a simulation-space vector
             @param[in]       v   is a constraint-space vector
             @param[in]       u   is the constraint argument; a simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ijv} = c_u(u,z)^{-1}v\f$, where
             \f$v \in \mathcal{C}\f$, \f$\mathsf{ijv} \in \mathcal{U}\f$.

             ---
  */
  virtual void applyInverseJacobian_1(Vector<Real> &ijv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    ijv.zero();
  }

  /** \brief Apply the adjoint of the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_u(u,z)^* \in L(\mathcal{C}^*, \mathcal{U}^*)\f$,
             to the vector \f$v\f$.  This is the primary interface.

             @param[out]      ajv    is the result of applying the adjoint of the constraint Jacobian to @b v at @b (u,z); a dual simulation-space vector
             @param[in]       v      is a dual constraint-space vector
             @param[in]       u      is the constraint argument; a simulation-space vector
             @param[in]       z      is the constraint argument; an optimization-space vector
             @param[in,out]   tol    is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ajv} = c_u(u,z)^*v\f$, where
             \f$v \in \mathcal{C}^*\f$, \f$\mathsf{ajv} \in \mathcal{U}^*\f$.

             ---
  */
  virtual void applyAdjointJacobian_1(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    applyAdjointJacobian_1(ajv, v, u, z, v.dual(), tol);
  }


  /** \brief Apply the adjoint of the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_u(u,z)^* \in L(\mathcal{C}^*, \mathcal{U}^*)\f$,
             to the vector \f$v\f$.  This is the secondary interface, for use
             with dual spaces where the user does not define the dual() operation.

             @param[out]      ajv    is the result of applying the adjoint of the constraint Jacobian to @b v at @b (u,z); a dual simulation-space vector
             @param[in]       v      is a dual constraint-space vector
             @param[in]       u      is the constraint argument; a simulation-space vector
             @param[in]       z      is the constraint argument; an optimization-space vector
             @param[in]       dualv  is a vector used for temporary variables; a constraint-space vector
             @param[in,out]   tol    is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ajv} = c_u(u,z)^*v\f$, where
             \f$v \in \mathcal{C}^*\f$, \f$\mathsf{ajv} \in \mathcal{U}^*\f$.

             ---
  */
  virtual void applyAdjointJacobian_1(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      const Vector<Real> &dualv,
                                      Real &tol) {
    Real ctol = std::sqrt(ROL_EPSILON);
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    Teuchos::RCP<Vector<Real> > cold = dualv.clone();
    Teuchos::RCP<Vector<Real> > cnew = dualv.clone();
    update(u,z);
    value(*cold,u,z,ctol);
    Teuchos::RCP<Vector<Real> > unew = u.clone();
    ajv.zero();
    for (int i = 0; i < u.dimension(); i++) {
      unew->set(u);
      unew->axpy(h,*(u.basis(i)));
      update(*unew,z);
      value(*cnew,*unew,z,ctol);
      cnew->axpy(-1.0,*cold);
      cnew->scale(1.0/h);
      ajv.axpy(cnew->dot(v),*((u.dual()).basis(i)));
    }
    update(u,z);
  }


  /** \brief Apply the adjoint of the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_z(u,z)^* \in L(\mathcal{C}^*, \mathcal{Z}^*)\f$,
             to vector \f$v\f$.  This is the primary interface.

             @param[out]      ajv    is the result of applying the adjoint of the constraint Jacobian to @b v at @b \f$(u,z)\f$; a dual optimization-space vector
             @param[in]       v      is a dual constraint-space vector
             @param[in]       u      is the constraint argument; a simulation-space vector
             @param[in]       z      is the constraint argument; an optimization-space vector
             @param[in,out]   tol    is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ajv} = c_z(u,z)^*v\f$, where
             \f$v \in \mathcal{C}^*\f$, \f$\mathsf{ajv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointJacobian_2(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    applyAdjointJacobian_2(ajv, v, u, z, v.dual(), tol);
  }


  /** \brief Apply the adjoint of the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_z(u,z)^* \in L(\mathcal{C}^*, \mathcal{Z}^*)\f$,
             to vector \f$v\f$.  This is the secondary interface, for use
             with dual spaces where the user does not define the dual() operation.

             @param[out]      ajv    is the result of applying the adjoint of the constraint Jacobian to @b v at @b \f$(u,z)\f$; a dual optimization-space vector
             @param[in]       v      is a dual constraint-space vector
             @param[in]       u      is the constraint argument; a simulation-space vector
             @param[in]       z      is the constraint argument; an optimization-space vector
             @param[in]       dualv  is a vector used for temporary variables; a constraint-space vector
             @param[in,out]   tol    is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ajv} = c_z(u,z)^*v\f$, where
             \f$v \in \mathcal{C}^*\f$, \f$\mathsf{ajv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointJacobian_2(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      const Vector<Real> &dualv,
                                      Real &tol) {
    Real ctol = std::sqrt(ROL_EPSILON);
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    Teuchos::RCP<Vector<Real> > cold = dualv.clone();
    Teuchos::RCP<Vector<Real> > cnew = dualv.clone();
    update(u,z);
    value(*cold,u,z,ctol);
    Teuchos::RCP<Vector<Real> > znew = z.clone();
    ajv.zero();
    for (int i = 0; i < z.dimension(); i++) {
      znew->set(z);
      znew->axpy(h,*(z.basis(i)));
      update(u,*znew);
      value(*cnew,u,*znew,ctol);
      cnew->axpy(-1.0,*cold);
      cnew->scale(1.0/h);
      ajv.axpy(cnew->dot(v),*((z.dual()).basis(i)));
    }
    update(u,z);
  }

  /** \brief Apply the inverse of the adjoint of the partial constraint Jacobian at \f$(u,z)\f$, 
             \f$c_u(u,z)^{-*} \in L(\mathcal{U}^*, \mathcal{C}^*)\f$,
             to the vector \f$v\f$.

             @param[out]      iajv is the result of applying the inverse adjoint of the constraint Jacobian to @b v at @b (u,z); a dual constraint-space vector
             @param[in]       v   is a dual simulation-space vector
             @param[in]       u   is the constraint argument; a simulation-space vector
             @param[in]       z   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{iajv} = c_u(u,z)^{-*}v\f$, where
             \f$v \in \mathcal{U}^*\f$, \f$\mathsf{iajv} \in \mathcal{C}^*\f$.

             ---
  */
  virtual void applyInverseAdjointJacobian_1(Vector<Real> &iajv,
                                             const Vector<Real> &v,
                                             const Vector<Real> &u,
                                             const Vector<Real> &z,
                                             Real &tol) {
    iajv.zero();
  };

  /** \brief Apply the adjoint of the partial constraint Hessian at \f$(u,z)\f$,
             \f$c_{uu}(u,z)^* \in L(L(\mathcal{C}^*, \mathcal{U}^*), \mathcal{U}^*)\f$,
             to vector \f$v\f$ in direction \f$w\f$.

             @param[out]      ahwv is the result of applying the adjoint of the constraint Hessian to @b v at @b \f$(u,z)\f$ in direction @b w; a dual simulation-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a dual simulation-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$ \mathsf{ahwv} = c_{uu}(u,z)^*(w,v) \f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{U}^*\f$, and \f$\mathsf{ahwv} \in \mathcal{U}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_11(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON);
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new state
    Teuchos::RCP<Vector<Real> > unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    update(*unew,z);
    applyAdjointJacobian_1(ahwv,w,*unew,z,jtol);
    // Evaluate Jacobian at old state
    Teuchos::RCP<Vector<Real> > jv = ahwv.clone();
    update(u,z);
    applyAdjointJacobian_1(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
  }


  /** \brief Apply the adjoint of the partial constraint Hessian at \f$(u,z)\f$,
             \f$c_{uz}(u,z)^* \in L(L(\mathcal{C}^*, \mathcal{U}^*), \mathcal{Z}^*)\f$,
             to vector \f$v\f$ in direction \f$w\f$.

             @param[out]      ahwv is the result of applying the adjoint of the constraint Hessian to @b v at @b \f$(u,z)\f$ in direction @b w; a dual optimization-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a dual simulation-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$ \mathsf{ahwv} = c_{uz}(u,z)^*(w,v) \f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{U}^*\f$, and \f$\mathsf{ahwv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_12(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON);
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new state
    Teuchos::RCP<Vector<Real> > unew = u.clone();
    unew->set(u);
    unew->axpy(h,v);
    update(*unew,z);
    applyAdjointJacobian_2(ahwv,w,*unew,z,jtol);
    // Evaluate Jacobian at old state
    Teuchos::RCP<Vector<Real> > jv = ahwv.clone();
    update(u,z);
    applyAdjointJacobian_2(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
  }


  /** \brief Apply the adjoint of the partial constraint Hessian at \f$(u,z)\f$,
             \f$c_{zu}(u,z)^* \in L(L(\mathcal{C}^*, \mathcal{Z}^*), \mathcal{U}^*)\f$,
             to vector \f$v\f$ in direction \f$w\f$.

             @param[out]      ahwv is the result of applying the adjoint of the constraint Hessian to @b v at @b \f$(u,z)\f$ in direction @b w; a dual simulation-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a dual optimization-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$ \mathsf{ahwv} = c_{zu}(u,z)^*(w,v) \f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{Z}^*\f$, and \f$\mathsf{ahwv} \in \mathcal{U}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_21(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON);
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new control
    Teuchos::RCP<Vector<Real> > znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    update(u,*znew);
    applyAdjointJacobian_1(ahwv,w,u,*znew,jtol);
    // Evaluate Jacobian at old control
    Teuchos::RCP<Vector<Real> > jv = ahwv.clone();
    update(u,z);
    applyAdjointJacobian_1(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
  }

  /** \brief Apply the adjoint of the partial constraint Hessian at \f$(u,z)\f$,
             \f$c_{zz}(u,z)^* \in L(L(\mathcal{C}^*, \mathcal{Z}^*), \mathcal{Z}^*)\f$,
             to vector \f$v\f$ in direction \f$w\f$.

             @param[out]      ahwv is the result of applying the adjoint of the constraint Hessian to @b v at @b \f$(u,z)\f$ in direction @b w; a dual optimization-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a dual optimization-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$ \mathsf{ahwv} = c_{zz}(u,z)^*(w,v) \f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{Z}^*\f$, and \f$\mathsf{ahwv} \in \mathcal{Z}^*\f$. 

             ---
  */
  virtual void applyAdjointHessian_22(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) {
    Real jtol = std::sqrt(ROL_EPSILON);
    // Compute step size
    Real h = tol;
    if (v.norm() > std::sqrt(ROL_EPSILON)) {
      h = std::max(1.0,u.norm()/v.norm())*tol;
    }
    // Evaluate Jacobian at new control
    Teuchos::RCP<Vector<Real> > znew = z.clone();
    znew->set(z);
    znew->axpy(h,v);
    update(u,*znew);
    applyAdjointJacobian_2(ahwv,w,u,*znew,jtol);
    // Evaluate Jacobian at old control
    Teuchos::RCP<Vector<Real> > jv = ahwv.clone();
    update(u,z);
    applyAdjointJacobian_2(*jv,w,u,z,jtol);
    // Compute Newton quotient
    ahwv.axpy(-1.0,*jv);
    ahwv.scale(1.0/h);
}

  /** \brief Approximately solves the <em> augmented system </em>
             \f[
                 \begin{pmatrix}
                   I     & c'(x)^* \\
                   c'(x) & 0
                 \end{pmatrix}
                 \begin{pmatrix}
                   v_{1} \\
                   v_{2}
                 \end{pmatrix}
                 =
                 \begin{pmatrix}
                   b_{1} \\
                   b_{2}
                 \end{pmatrix}
             \f]
             where \f$v_{1} \in \mathcal{X}\f$, \f$v_{2} \in \mathcal{C}^*\f$,
             \f$b_{1} \in \mathcal{X}^*\f$, \f$b_{2} \in \mathcal{C}\f$,
             \f$I : \mathcal{X} \rightarrow \mathcal{X}^*\f$ is an identity
             operator, and \f$0 : \mathcal{C}^* \rightarrow \mathcal{C}\f$
             is a zero operator.

             @param[out]      v1  is the optimization-space component of the result
             @param[out]      v2  is the dual constraint-space component of the result
             @param[in]       b1  is the dual optimization-space component of the right-hand side
             @param[in]       b2  is the constraint-space component of the right-hand side
             @param[in]       x   is the constraint argument; an optimization-space vector
             @param[in,out]   tol is the nominal relative residual tolerance

             On return, \f$ [\mathsf{v1} \,\, \mathsf{v2}] \f$ approximately
             solves the augmented system, where the size of the residual is
             governed by special stopping conditions. \n\n
             The default implementation is the preconditioned generalized
             minimal residual (GMRES) method, which enables the use of
             nonsymmetric preconditioners.

             ---
  */
  virtual std::vector<Real> solveAugmentedSystem(Vector<Real> &v1,
                                                 Vector<Real> &v2,
                                                 const Vector<Real> &b1,
                                                 const Vector<Real> &b2,
                                                 const Vector<Real> &x,
                                                 Real &tol) {
    return EqualityConstraint<Real>::solveAugmentedSystem(v1,v2,b1,b2,x,tol);
  }


  /** \brief Apply a constraint preconditioner at \f$x\f$, \f$P(x) \in L(\mathcal{C}, \mathcal{C})\f$,
             to vector \f$v\f$.  In general, this preconditioner satisfies the following relationship:
             \f[
               c'(x) c'(x)^* P(x) v \approx v \,.
             \f]
             It is used by the #solveAugmentedSystem method.

             @param[out]      pv  is the result of applying the constraint preconditioner to @b v at @b x; a constraint-space vector
             @param[in]       v   is a constraint-space vector
             @param[in]       x   is the preconditioner argument; an optimization-space vector
             @param[in,out]   tol is a tolerance for inexact evaluations

             On return, \f$\mathsf{pv} = P(x)v\f$, where
             \f$v \in \mathcal{C}\f$, \f$\mathsf{pv} \in \mathcal{C}\f$. \n\n
             The default implementation is a null-op.

             ---
  */
  virtual void applyPreconditioner(Vector<Real> &pv,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   Real &tol) {
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    Teuchos::RCP<ROL::Vector<Real> > ijv = (xs.get_1())->clone();
    applyInverseJacobian_1(*ijv, v, *(xs.get_1()), *(xs.get_2()), tol);
    applyInverseAdjointJacobian_1(pv, ijv->dual(), *(xs.get_1()), *(xs.get_2()), tol);
  }


  EqualityConstraint_SimOpt(void) : EqualityConstraint<Real>() {}

  /** \brief Update constraint functions.  
                x is the optimization variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    update(*(xs.get_1()),*(xs.get_2()),flag,iter);
  }

  /** \brief Check if the vector, v, is feasible
  */
  virtual bool isFeasible( const Vector<Real> &v ) { return true; }

  virtual void value(Vector<Real> &c,
                     const Vector<Real> &x,
                     Real &tol) {
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    value(c,*(xs.get_1()),*(xs.get_2()),tol);
  }


  virtual void applyJacobian(Vector<Real> &jv,
                             const Vector<Real> &v,
                             const Vector<Real> &x,
                             Real &tol) { 
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    const Vector_SimOpt<Real> &vs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(v));
    applyJacobian_1(jv,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    Teuchos::RCP<Vector<Real> > jv2 = jv.clone();
    applyJacobian_2(*jv2,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    jv.plus(*jv2);
  }


  virtual void applyAdjointJacobian(Vector<Real> &ajv,
                                    const Vector<Real> &v,
                                    const Vector<Real> &x,
                                    Real &tol) { 
    Vector_SimOpt<Real> &ajvs = Teuchos::dyn_cast<Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<Vector<Real> >(ajv));
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    Teuchos::RCP<Vector<Real> > ajv1 = (ajvs.get_1())->clone();
    applyAdjointJacobian_1(*ajv1,v,*(xs.get_1()),*(xs.get_2()),tol);
    ajvs.set_1(*ajv1);
    Teuchos::RCP<Vector<Real> > ajv2 = (ajvs.get_2())->clone();
    applyAdjointJacobian_2(*ajv2,v,*(xs.get_1()),*(xs.get_2()),tol);
    ajvs.set_2(*ajv2);
  }


  virtual void applyAdjointHessian(Vector<Real> &ahwv,
                                   const Vector<Real> &w,
                                   const Vector<Real> &v,
                                   const Vector<Real> &x,
                                   Real &tol) {
    Vector_SimOpt<Real> &ahwvs = Teuchos::dyn_cast<Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<Vector<Real> >(ahwv));
    const Vector_SimOpt<Real> &xs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(x));
    const Vector_SimOpt<Real> &vs = Teuchos::dyn_cast<const Vector_SimOpt<Real> >(
      Teuchos::dyn_cast<const Vector<Real> >(v));
    // Block-row 1
    Teuchos::RCP<Vector<Real> > C11 = (ahwvs.get_1())->clone();
    Teuchos::RCP<Vector<Real> > C21 = (ahwvs.get_1())->clone();
    applyAdjointHessian_11(*C11,w,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    applyAdjointHessian_21(*C21,w,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    C11->plus(*C21);
    ahwvs.set_1(*C11); 
    // Block-row 2
    Teuchos::RCP<Vector<Real> > C12 = (ahwvs.get_2())->clone();
    Teuchos::RCP<Vector<Real> > C22 = (ahwvs.get_2())->clone();
    applyAdjointHessian_12(*C12,w,*(vs.get_1()),*(xs.get_1()),*(xs.get_2()),tol);
    applyAdjointHessian_22(*C22,w,*(vs.get_2()),*(xs.get_1()),*(xs.get_2()),tol);
    C22->plus(*C12);
    ahwvs.set_2(*C22); 
  }



  virtual Real checkSolve(const ROL::Vector<Real> &u, 
                          const ROL::Vector<Real> &z, 
                          const ROL::Vector<Real> &c,
                          const bool printToStream = true,
                          std::ostream & outStream = std::cout) {
    // Solve equality constraint for u. 
    Real tol = ROL_EPSILON;
    Teuchos::RCP<ROL::Vector<Real> > s = u.clone();
    solve(*s,z,tol);
    // Evaluate equality constraint residual at (u,z).
    Teuchos::RCP<ROL::Vector<Real> > cs = c.clone();
    value(*cs,*s,z,tol);
    // Output norm of residual.
    Real cnorm = cs->norm();
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt solve at feasible (u,z): \n  ||c(u,z)|| = " << cnorm << "\n";
      outStream << hist.str();
    }
    return cnorm;
  }


  /** \brief Check the consistency of the Jacobian and its adjoint.
             This is the primary interface.

             @param[out]      w              is a dual constraint-space vector
             @param[in]       v              is a simulation-space vector
             @param[in]       u              is the constraint argument; a simulation-space vector
             @param[in]       z              is the constraint argument; an optimization-space vector
             @param[in]       printToStream  is is a flag that turns on/off output
             @param[in]       outStream      is the output stream

             ---
  */
  virtual Real checkAdjointConsistencyJacobian_1(const Vector<Real> &w, 
                                                 const Vector<Real> &v, 
                                                 const Vector<Real> &u,
                                                 const Vector<Real> &z,
                                                 const bool printToStream = true,
                                                 std::ostream & outStream = std::cout) {
    return checkAdjointConsistencyJacobian_1(w, v, u, z, w.dual(), v.dual(), printToStream, outStream);
  }


  /** \brief Check the consistency of the Jacobian and its adjoint.
             This is the secondary interface, for use with dual spaces where
             the user does not define the dual() operation.

             @param[out]      w              is a dual constraint-space vector
             @param[in]       v              is a simulation-space vector
             @param[in]       u              is the constraint argument; a simulation-space vector
             @param[in]       z              is the constraint argument; an optimization-space vector
             @param[in]       dualw          is a constraint-space vector 
             @param[in]       dualv          is a dual simulation-space vector
             @param[in]       printToStream  is is a flag that turns on/off output
             @param[in]       outStream      is the output stream

             ---
  */
  virtual Real checkAdjointConsistencyJacobian_1(const Vector<Real> &w, 
                                                 const Vector<Real> &v, 
                                                 const Vector<Real> &u,
                                                 const Vector<Real> &z,
                                                 const Vector<Real> &dualw,
                                                 const Vector<Real> &dualv,
                                                 const bool printToStream = true,
                                                 std::ostream & outStream = std::cout) {
    Real tol = ROL_EPSILON;
    Teuchos::RCP<Vector<Real> > Jv = dualw.clone();
    applyJacobian_1(*Jv,v,u,z,tol);
    Real wJv = w.dot(Jv->dual());
    Teuchos::RCP<Vector<Real> > Jw = dualv.clone();
    applyAdjointJacobian_1(*Jw,w,u,z,tol);
    Real vJw = v.dot(Jw->dual());
    Real diff = std::abs(wJv-vJw);
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of Jacobian_1 and its adjoint: \n  |<w,Jv> - <adj(J)w,v>| = " 
           << diff << "\n";
      hist << "  |<w,Jv>|               = " << std::abs(wJv) << "\n";
      hist << "  Relative Error         = " << diff / (std::abs(wJv)+ROL_UNDERFLOW) << "\n";
      outStream << hist.str();
    }
    return diff;
  }


  /** \brief Check the consistency of the Jacobian and its adjoint.
             This is the primary interface.

             @param[out]      w              is a dual constraint-space vector
             @param[in]       v              is an optimization-space vector
             @param[in]       u              is the constraint argument; a simulation-space vector
             @param[in]       z              is the constraint argument; an optimization-space vector
             @param[in]       printToStream  is is a flag that turns on/off output
             @param[in]       outStream      is the output stream

             ---
  */
  virtual Real checkAdjointConsistencyJacobian_2(const Vector<Real> &w, 
                                                 const Vector<Real> &v, 
                                                 const Vector<Real> &u,
                                                 const Vector<Real> &z,
                                                 const bool printToStream = true,
                                                 std::ostream & outStream = std::cout) {
    return checkAdjointConsistencyJacobian_2(w, v, u, z, w.dual(), v.dual(), printToStream, outStream);
  }

  /** \brief Check the consistency of the Jacobian and its adjoint.
             This is the secondary interface, for use with dual spaces where
             the user does not define the dual() operation.

             @param[out]      w              is a dual constraint-space vector
             @param[in]       v              is an optimization-space vector
             @param[in]       u              is the constraint argument; a simulation-space vector
             @param[in]       z              is the constraint argument; an optimization-space vector
             @param[in]       dualw          is a constraint-space vector 
             @param[in]       dualv          is a dual optimization-space vector
             @param[in]       printToStream  is is a flag that turns on/off output
             @param[in]       outStream      is the output stream

             ---
  */
  virtual Real checkAdjointConsistencyJacobian_2(const Vector<Real> &w, 
                                                 const Vector<Real> &v, 
                                                 const Vector<Real> &u,
                                                 const Vector<Real> &z,
                                                 const Vector<Real> &dualw,
                                                 const Vector<Real> &dualv,
                                                 const bool printToStream = true,
                                                 std::ostream & outStream = std::cout) {
    Real tol = ROL_EPSILON;
    Teuchos::RCP<Vector<Real> > Jv = dualw.clone();
    applyJacobian_2(*Jv,v,u,z,tol);
    Real wJv = w.dot(Jv->dual());
    Teuchos::RCP<Vector<Real> > Jw = dualv.clone();
    applyAdjointJacobian_2(*Jw,w,u,z,tol);
    Real vJw = v.dot(Jw->dual());
    Real diff = std::abs(wJv-vJw);
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of Jacobian_2 and its adjoint: \n  |<w,Jv> - <adj(J)w,v>| = "
           << diff << "\n";
      hist << "  |<w,Jv>|               = " << std::abs(wJv) << "\n";
      hist << "  Relative Error         = " << diff / (std::abs(wJv)+ROL_UNDERFLOW) << "\n";
      outStream << hist.str();
    }
    return diff;
  }

  virtual Real checkInverseJacobian_1(const Vector<Real> &jv, 
                                      const Vector<Real> &v, 
                                      const Vector<Real> &u, 
                                      const Vector<Real> &z, 
                                      const bool printToStream = true,
                                      std::ostream & outStream = std::cout) {
    Real tol = ROL_EPSILON;
    Teuchos::RCP<Vector<Real> > Jv = jv.clone();
    applyJacobian_1(*Jv,v,u,z,tol);
    Teuchos::RCP<Vector<Real> > iJJv = u.clone();
    applyInverseJacobian_1(*iJJv,*Jv,u,z,tol);
    Teuchos::RCP<Vector<Real> > diff = v.clone();
    diff->set(v);
    diff->axpy(-1.0,*iJJv);
    Real dnorm = diff->norm();
    Real vnorm = v.norm();
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of inverse Jacobian_1: \n  ||v-inv(J)Jv|| = " 
           << dnorm << "\n";
      hist << "  ||v||          = " << vnorm << "\n";
      hist << "  Relative Error = " << dnorm / (vnorm+ROL_UNDERFLOW) << "\n";
      outStream << hist.str();
    }
    return dnorm;
  }

  virtual Real checkInverseAdjointJacobian_1(const Vector<Real> &jv, 
                                             const Vector<Real> &v, 
                                             const Vector<Real> &u, 
                                             const Vector<Real> &z, 
                                             const bool printToStream = true,
                                             std::ostream & outStream = std::cout) {
    Real tol = ROL_EPSILON;
    Teuchos::RCP<Vector<Real> > Jv = jv.clone();
    applyAdjointJacobian_1(*Jv,v,u,z,tol);
    Teuchos::RCP<Vector<Real> > iJJv = v.clone();
    applyInverseAdjointJacobian_1(*iJJv,*Jv,u,z,tol);
    Teuchos::RCP<Vector<Real> > diff = v.clone();
    diff->set(v);
    diff->axpy(-1.0,*iJJv);
    Real dnorm = diff->norm();
    Real vnorm = v.norm();
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt consistency of inverse adjoint Jacobian_1: \n  ||v-inv(adj(J))adj(J)v|| = "
           << dnorm << "\n";
      hist << "  ||v||                   = " << vnorm << "\n";
      hist << "  Relative Error          = " << dnorm / (vnorm+ROL_UNDERFLOW) << "\n";
      outStream << hist.str();
    }
    return dnorm;
  }

}; // class EqualityConstraint_SimOpt

} // namespace ROL

#endif
