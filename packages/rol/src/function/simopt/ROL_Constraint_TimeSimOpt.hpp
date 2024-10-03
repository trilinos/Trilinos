// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINT_TIMESIMOPT_H
#define ROL_CONSTRAINT_TIMESIMOPT_H

#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_VectorWorkspace.hpp"

/** @ingroup func_group
    \class ROL::Constraint_TimeSimOpt
    \brief Defines the time dependent constraint operator interface for simulation-based optimization.

    This constraint interface inherits from ROL_Constraint_SimOpt. Though the interface
    takes two simulation space vectors from spaces
    \f$\mathcal{U_o}\times\mathcal{U_n}\f$. The space \f$\mathcal{U_o}\f$ is ``old'' information
    that accounts for the initial condition on the time interval. The space \f$\mathcal{U_n}\f$ is the
    ``new'' variables that can be determined by satisfying constraints in the form
    \f[
      c(u_o,u_n,z) = 0 \,.
    \f]
    where \f$u_0 \in \mathcal{U_o},\; u_n\in\mathcal{U_n},\f$ and \f$z\in\mathcal{Z}\f$. In this way
    this constraint defines a sequence of state variables.
    The basic operator interface, to be implemented by the user, requires:
    \li #value -- constraint evaluation.
    \li #applyJacobian_1_old         -- action of the partial constraint Jacobian --derivatives are 
                                        with respect to the first component \f$\mathcal{U_o}\f$;
    \li #applyJacobian_1_new         -- action of the partial constraint Jacobian --derivatives are 
                                        with respect to the first component \f$\mathcal{U_n}\f$;
    \li #applyJacobian_2             -- action of the partial constraint Jacobian --derivatives are 
                                        with respect to the second component \f$\mathcal{Z}\f$;
    \li #applyAdjointJacobian_1_old  -- action of the adjoint of the partial constraint Jacobian --derivatives are 
                                        with respect to the first component \f$\mathcal{U_o}\f$;
    \li #applyAdjointJacobian_1_new  -- action of the adjoint of the partial constraint Jacobian --derivatives are 
                                        with respect to the first component \f$\mathcal{U_n}\f$;
    \li #applyAdjointJacobian_2_time -- action of the adjoint of the partial constraint Jacobian --derivatives are 
                                        with respect to the second component \f$\mathcal{Z}\f$; (note the time
                                        suffix here is to prevent collisions with the parent class)

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
class Constraint_TimeSimOpt : public Constraint_SimOpt<Real> {
private:

  // Get the end point of the time intervals vector
  Vector<Real> & getNewVector(Vector<Real> & x) const
  { 
    PartitionedVector<Real> & xpv = dynamic_cast<PartitionedVector<Real>&>(x);
    return *xpv.get(1);
  }

  const Vector<Real> & getNewVector(const Vector<Real> & x) const
  { 
    const PartitionedVector<Real> & xpv = dynamic_cast<const PartitionedVector<Real>&>(x);
    return *xpv.get(1);
  }
 
  // Get the start point of the time intervals vector
  Vector<Real> & getOldVector(Vector<Real> & x) const
  { 
    PartitionedVector<Real> & xpv = dynamic_cast<PartitionedVector<Real>&>(x);
    return *xpv.get(0);
  }

  const Vector<Real> & getOldVector(const Vector<Real> & x) const
  { 
    const PartitionedVector<Real> & xpv = dynamic_cast<const PartitionedVector<Real>&>(x);
    return *xpv.get(0);
  }

  mutable VectorWorkspace<Real> workspace_;

protected:

  VectorWorkspace<Real>& getVectorWorkspace() const { return workspace_; }

public:
  Constraint_TimeSimOpt()
    : Constraint_SimOpt<Real>()
  { }

  // Interface functions (to be overloaded)

  /** \brief Update constraint functions.  
                u_old Is the state from the end of the previous time step.
                u_new Is the state from the end of this time step.
                z Is the control variable
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update( const Vector<Real> & u_old,
                       const Vector<Real> & u_new,
                       const Vector<Real> &z, 
                       bool flag = true, int iter = -1 ) {
    update_1_old(u_old,flag,iter);
    update_1_new(u_new,flag,iter);
    update_2(z,flag,iter);  
  }

  /** \brief Update constraint functions with respect to Sim variable.  
                u_old is the state variable
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_1_old( const Vector<Real> &u_old, bool flag = true, int iter = -1 ) {}

  /** \brief Update constraint functions with respect to Sim variable.  
                u_new is the state variable
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_1_new( const Vector<Real> &u_new, bool flag = true, int iter = -1 ) {}

  /** \brief Update constraint functions with respect to Opt variable.
                z is the control variable, 
                flag = true if optimization variable is changed,
                iter is the outer algorithm iterations count.
  */
  virtual void update_2( const Vector<Real> &z, bool flag = true, int iter = -1 ) override {}

  /** \brief Evaluate the constraint operator \f$c:\mathcal{U_o}\times\mathcal{U_n}\times\mathcal{Z} \rightarrow \mathcal{C}\f$
             at \f$(u,z)\f$.

             @param[out]      c      is the result of evaluating the constraint operator at @b \f$(u,z)\f$; a constraint-space vector
             @param[in]       u_old  is the constraint argument; a simulation-space vector from the previous interval
             @param[in]       u_new  is the constraint argument; a simulation-space vector from the current interval
             @param[in]       z      is the constraint argument; an optimization-space vector
             @param[in,out]   tol    is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{c} = c(u,z)\f$,
             where \f$\mathsf{c} \in \mathcal{C}\f$, \f$\mathsf{u_o} \in \mathcal{U_o}\f$, 
                                                     \f$\mathsf{u_n} \in \mathcal{U_n}\f$, 
                                                     and $\f$\mathsf{z} \in\mathcal{Z}\f$.

             ---
  */
  virtual void value(Vector<Real> &c,
                     const Vector<Real> &u_old,
                     const Vector<Real> &u_new,
                     const Vector<Real> &z,
                     Real &tol) = 0;

  virtual void solve(Vector<Real> &c,
                     const Vector<Real> &u_old,
                     Vector<Real> &u_new,
                     const Vector<Real> &z,
                     Real &tol) = 0;

  virtual void applyJacobian_1_old(Vector<Real> &jv,
                                   const Vector<Real> &v_old,
                                   const Vector<Real> &u_old, const Vector<Real> &u_new,
                                   const Vector<Real> &z,
                                   Real &tol) = 0;

  virtual void applyJacobian_1_new(Vector<Real> &jv,
                                   const Vector<Real> &v_new,
                                   const Vector<Real> &u_old, const Vector<Real> &u_new,
                                   const Vector<Real> &z,
                                   Real &tol) = 0;

  virtual void applyInverseJacobian_1_new(Vector<Real> &ijv,
                                          const Vector<Real> &v_new,
                                          const Vector<Real> &u_old, const Vector<Real> &u_new,
                                          const Vector<Real> &z,
                                          Real &tol) = 0;


  virtual void applyJacobian_2(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u_old, const Vector<Real> &u_new,
                               const Vector<Real> &z,
                               Real &tol) = 0;

  virtual void applyAdjointJacobian_1_old(Vector<Real> &ajv_old,
                                      const Vector<Real> &dualv,
                                      const Vector<Real> &u_old, const Vector<Real> &u_new,
                                      const Vector<Real> &z,
                                      Real &tol) = 0;

  virtual void applyAdjointJacobian_1_new(Vector<Real> &ajv_new,
                                      const Vector<Real> &dualv,
                                      const Vector<Real> &u_old, const Vector<Real> &u_new,
                                      const Vector<Real> &z,
                                      Real &tol) = 0;

  virtual void applyInverseAdjointJacobian_1_new(Vector<Real> &iajv,
                                                 const Vector<Real> &v_new,
                                                 const Vector<Real> &u_old, const Vector<Real> &u_new,
                                                 const Vector<Real> &z,
                                                 Real &tol) = 0;

  virtual void applyAdjointJacobian_2_time(Vector<Real> &ajv,
                                      const Vector<Real> &dualv,
                                      const Vector<Real> &u_old, const Vector<Real> &u_new,
                                      const Vector<Real> &z,
                                      Real &tol) = 0;

  virtual void applyAdjointHessian_11_old(Vector<Real> &ahwv_old,
                                          const Vector<Real> &w,
                                          const Vector<Real> &v_new,
                                          const Vector<Real> &u_old, const Vector<Real> &u_new,
                                          const Vector<Real> &z,
                                          Real &tol) = 0;

  virtual void applyAdjointHessian_11_new(Vector<Real> &ahwv_new,
                                          const Vector<Real> &w,
                                          const Vector<Real> &v_new,
                                          const Vector<Real> &u_old, const Vector<Real> &u_new,
                                          const Vector<Real> &z,
                                          Real &tol) = 0;

  // Functions from SimOpt that are overriden
  ///////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void update( const Vector<Real> & u,
                       const Vector<Real> & z, 
                       bool flag = true, int iter = -1 ) override {
    update(getOldVector(u),
           getNewVector(u),
           z,
           flag,iter);  
  }

  virtual void value(Vector<Real> &c,
                     const Vector<Real> &u,
                     const Vector<Real> &z,
                     Real &tol) override {

    value(c,
          getOldVector(u),
          getNewVector(u),
          z,
          tol); 
  } 

  virtual void solve(Vector<Real> &c,
                     Vector<Real> &u, 
                     const Vector<Real> &z,
                     Real &tol) override {
    solve(c,
          getOldVector(u),
          getNewVector(u),
          z,
          tol);  
  }

  virtual void applyJacobian_1(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) override {
    const Vector<Real> & v_old = getOldVector(v);
    const Vector<Real> & v_new = getNewVector(v);
    const Vector<Real> & u_old = getOldVector(u);
    const Vector<Real> & u_new = getNewVector(u);

    // evaluate derivative against "old" time variable
    applyJacobian_1_old(jv,v_old,
                           u_old,u_new,
                           z,
                           tol);

    ROL::Ptr<Vector<Real> > jv_new = workspace_.clone(jv);

    // evaluate derivative against "new" time variable
    applyJacobian_1_new(*jv_new,v_new,
                           u_old,u_new,
                           z,
                           tol);

    jv.axpy(1.0,*jv_new);
  }

  virtual void applyJacobian_2(Vector<Real> &jv,
                               const Vector<Real> &v,
                               const Vector<Real> &u,
                               const Vector<Real> &z,
                               Real &tol) override { 
    const Vector<Real> & u_old = getOldVector(u);
    const Vector<Real> & u_new = getNewVector(u);

    // evaluate derivative against "old" time variable
    applyJacobian_2(jv,v,
                    u_old,u_new,
                    z,
                    tol);
  }

  virtual void applyInverseJacobian_1(Vector<Real> &ijv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override final {
    ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
      "The method applyInverseJacobian_1 is used but not implemented!\n");
  }

  virtual void applyAdjointJacobian_1(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override {
    Vector<Real> & ajv_old = getOldVector(ajv);
    Vector<Real> & ajv_new = getNewVector(ajv);
    const Vector<Real> & u_old = getOldVector(u);
    const Vector<Real> & u_new = getNewVector(u);
    
    applyAdjointJacobian_1_old(ajv_old,v,u_old,u_new,z,tol);
    applyAdjointJacobian_1_new(ajv_new,v,u_old,u_new,z,tol);
  }

  virtual void applyAdjointJacobian_2(Vector<Real> &ajv,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override {
    const Vector<Real> & u_old = getOldVector(u);
    const Vector<Real> & u_new = getNewVector(u);
    
    applyAdjointJacobian_2_time(ajv,v,u_old,u_new,z,tol);
  }

  virtual void applyInverseAdjointJacobian_1(Vector<Real> &iajv,
                                             const Vector<Real> &v,
                                             const Vector<Real> &u,
                                             const Vector<Real> &z,
                                             Real &tol) override final {
    ROL_TEST_FOR_EXCEPTION(true, std::logic_error,
      "The method applyInverseAdjointJacobian_1 is used but not implemented!\n");
  };

  /** \brief Apply the simulation-space derivative of the adjoint of the constraint
             simulation-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the
             direction \f$v\f$, according to \f$v\mapsto c_{uu}(u,z)(v,\cdot)^*w\f$.

             @param[out]      ahwv is the result of applying the simulation-space derivative of the adjoint of the constraint simulation-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual simulation-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a simulation-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ahwv} = c_{uu}(u,z)(v,\cdot)^*w\f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{U}\f$, and
             \f$\mathsf{ahwv} \in \mathcal{U}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_11(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override
  {
    Vector<Real> & ahwv_old = getOldVector(ahwv);
    Vector<Real> & ahwv_new = getNewVector(ahwv);
    const Vector<Real> & v_old = getOldVector(v);
    const Vector<Real> & v_new = getNewVector(v);
    const Vector<Real> & u_old = getOldVector(u);
    const Vector<Real> & u_new = getNewVector(u);

    // this implicitly assumes that there is no cross coupling
    // between the old state and the new state. Is that true? For 
    // simple (Euler, Theta method) integrators yes.
    applyAdjointHessian_11_old(ahwv_old,w,v_old,u_old,u_new,z,tol);
    applyAdjointHessian_11_new(ahwv_new,w,v_new,u_old,u_new,z,tol);
  }

  /** \brief Apply the optimization-space derivative of the adjoint of the constraint
             simulation-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the
             direction \f$v\f$, according to \f$v\mapsto c_{uz}(u,z)(v,\cdot)^*w\f$.

             @param[out]      ahwv is the result of applying the optimization-space derivative of the adjoint of the constraint simulation-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual optimization-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a simulation-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ahwv} = c_{uz}(u,z)(v,\cdot)^*w\f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{U}\f$, and
             \f$\mathsf{ahwv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_12(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override
  {
    ahwv.zero();
  }

  /** \brief Apply the simulation-space derivative of the adjoint of the constraint 
             optimization-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the 
             direction \f$v\f$, according to \f$v\mapsto c_{zu}(u,z)(v,\cdot)^*w\f$. 
 
             @param[out]      ahwv is the result of applying the simulation-space derivative of the adjoint of the constraint optimization-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual simulation-space vector 
             @param[in]       w    is the direction vector; a dual constraint-space vector 
             @param[in]       v    is a optimization-space vector 
             @param[in]       u    is the constraint argument; a simulation-space vector 
             @param[in]       z    is the constraint argument; an optimization-space vector 
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused 
 
             On return, \f$\mathsf{ahwv} = c_{zu}(u,z)(v,\cdot)^*w\f$, where 
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{Z}\f$, and 
             \f$\mathsf{ahwv} \in \mathcal{U}^*\f$. 
 
             --- 
  */ 
  virtual void applyAdjointHessian_21(Vector<Real> &ahwv, 
                                      const Vector<Real> &w, 
                                      const Vector<Real> &v, 
                                      const Vector<Real> &u, 
                                      const Vector<Real> &z, 
                                      Real &tol) override { 
    ahwv.zero();
  }

  /** \brief Apply the optimization-space derivative of the adjoint of the constraint
             optimization-space Jacobian at \f$(u,z)\f$ to the vector \f$w\f$ in the
             direction \f$v\f$, according to \f$v\mapsto c_{zz}(u,z)(v,\cdot)^*w\f$.

             @param[out]      ahwv is the result of applying the optimization-space derivative of the adjoint of the constraint optimization-space Jacobian at @b \f$(u,z)\f$ to the vector @b \f$w\f$ in direction @b \f$w\f$; a dual optimization-space vector
             @param[in]       w    is the direction vector; a dual constraint-space vector
             @param[in]       v    is a optimization-space vector
             @param[in]       u    is the constraint argument; a simulation-space vector
             @param[in]       z    is the constraint argument; an optimization-space vector
             @param[in,out]   tol  is a tolerance for inexact evaluations; currently unused

             On return, \f$\mathsf{ahwv} = c_{zz}(u,z)(v,\cdot)^*w\f$, where
             \f$w \in \mathcal{C}^*\f$, \f$v \in \mathcal{Z}\f$, and
             \f$\mathsf{ahwv} \in \mathcal{Z}^*\f$.

             ---
  */
  virtual void applyAdjointHessian_22(Vector<Real> &ahwv,
                                      const Vector<Real> &w,
                                      const Vector<Real> &v,
                                      const Vector<Real> &u,
                                      const Vector<Real> &z,
                                      Real &tol) override {
    ahwv.zero();
  }

  // We override the check solve routine because we are abusing SimOpt
  virtual Real checkSolve(const ROL::Vector<Real> &u,
                          const ROL::Vector<Real> &z,
                          const ROL::Vector<Real> &c,
                          const bool printToStream = true,
                          std::ostream & outStream = std::cout) override {
    // Solve constraint for u. 
    Real tol = ROL_EPSILON<Real>();
    ROL::Ptr<ROL::Vector<Real> > r = workspace_.clone(c);
    ROL::Ptr<ROL::Vector<Real> > s = workspace_.clone(u);
    s->set(u);
    solve(*r,*s,z,tol);
    // Evaluate constraint residual at (u,z).
    ROL::Ptr<ROL::Vector<Real> > cs = workspace_.clone(c);
    update(*s,z);
    value(*cs,*s,z,tol);
    // Output norm of residual.
    Real rnorm = r->norm();
    Real cnorm = cs->norm();
    if ( printToStream ) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\nTest SimOpt solve at feasible (u,z):\n";
      hist << "  Solver Residual = " << rnorm << "\n";
      hist << "       ||c(u,z)|| = " << cnorm << "\n";
      outStream << hist.str();
    }
    return cnorm;
  }
  
  // Verify that ||v-Jinv*J*v|| < tol
  virtual Real checkInverseJacobian_1_new( const ROL::Vector<Real> &c,
                                           const ROL::Vector<Real> &u_new,
                                           const ROL::Vector<Real> &u_old,
                                           const ROL::Vector<Real> &z,
                                           const ROL::Vector<Real> &v_new,
                                           const bool printToStream = true,
                                           std::ostream & outStream = std::cout) {
     Real tol = ROL_EPSILON<Real>();
     auto Jv   = workspace_.clone(c);
     update( u_new, u_old, z );
     applyJacobian_1_new( *Jv, v_new, u_old, u_new, z, tol );
     auto iJJv = workspace_.clone(u_new);
     update( u_new, u_old, z );
     applyInverseJacobian_1_new( *iJJv, *Jv, u_old, u_new, z, tol );
     auto diff = workspace_.clone(v_new);
     diff->set(v_new);
     diff->axpy(-1.0,*iJJv);
     Real dnorm = diff->norm();
     Real vnorm = v_new.norm();
     if ( printToStream ) {
       std::stringstream hist;
       hist << std::scientific << std::setprecision(8);
       hist << "\nTest TimeSimOpt consistency of inverse Jacobian_1_new: \n  ||v-inv(J)Jv|| = " 
            << dnorm << "\n";
       hist << "  ||v||          = " << vnorm << "\n";
       hist << "  Relative Error = " << dnorm / (vnorm+ROL_UNDERFLOW<Real>()) << "\n";
       outStream << hist.str();
     }
     return dnorm;
   }

  virtual Real checkInverseAdjointJacobian_1_new( const ROL::Vector<Real> &c,
                                                  const ROL::Vector<Real> &u_new,
                                                  const ROL::Vector<Real> &u_old,
                                                  const ROL::Vector<Real> &z,
                                                  const ROL::Vector<Real> &v_new,
                                                  const bool printToStream = true,
                                                  std::ostream & outStream = std::cout) {
     Real tol = ROL_EPSILON<Real>();
     auto Jv   = workspace_.clone(c);
     update( u_new, u_old, z );
     applyAdjointJacobian_1_new( *Jv, v_new, u_old, u_new, z, tol );
     auto iJJv = workspace_.clone(u_new);
     update( u_new, u_old, z );
     applyInverseAdjointJacobian_1_new( *iJJv, *Jv, u_old, u_new, z, tol );
     auto diff = workspace_.clone(v_new);
     diff->set(v_new);
     diff->axpy(-1.0,*iJJv);
     Real dnorm = diff->norm();
     Real vnorm = v_new.norm();
     if ( printToStream ) {
       std::stringstream hist;
       hist << std::scientific << std::setprecision(8);
       hist << "\nTest TimeSimOpt consistency of inverse adjoint Jacobian_1_new: \n  ||v-inv(adj(J))adj(J)v|| = " 
            << dnorm << "\n";
       hist << "  ||v||          = " << vnorm << "\n";
       hist << "  Relative Error = " << dnorm / (vnorm+ROL_UNDERFLOW<Real>()) << "\n";
       outStream << hist.str();
     }
     return dnorm;
   }

  std::vector<std::vector<Real> > checkApplyJacobian_1_new(const Vector<Real> &u_new,
                                                       const Vector<Real> &u_old,
                                                       const Vector<Real> &z,
                                                       const Vector<Real> &v,
                                                       const Vector<Real> &jv,
                                                       const bool printToStream = true,
                                                       std::ostream & outStream = std::cout,
                                                       const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                       const int order = 1) {
    std::vector<Real> steps(numSteps);
    for(int i=0;i<numSteps;++i) {
      steps[i] = pow(10,-i);
    }
   
    return checkApplyJacobian_1_new(u_new,u_old,z,v,jv,steps,printToStream,outStream,order);
  }
  
  std::vector<std::vector<Real> > checkApplyJacobian_1_new(const Vector<Real> &u_new, 
                                                         const Vector<Real> &u_old, 
                                                         const Vector<Real> &z,
                                                         const Vector<Real> &v,
                                                         const Vector<Real> &jv,
                                                         const std::vector<Real> &steps, 
                                                         const bool printToStream = true,
                                                         std::ostream & outStream = std::cout,
                                                         const int order = 1) {
 
    ROL_TEST_FOR_EXCEPTION( order<1 || order>4, std::invalid_argument, 
                                "Error: finite difference order must be 1,2,3, or 4" );
 
    Real one(1.0);
 
    using Finite_Difference_Arrays::shifts;
    using Finite_Difference_Arrays::weights;
 
    Real tol = std::sqrt(ROL_EPSILON<Real>());
 
    int numSteps = steps.size();
    int numVals = 4;
    std::vector<Real> tmp(numVals);
    std::vector<std::vector<Real> > jvCheck(numSteps, tmp);
 
    // Save the format state of the original outStream.
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);
 
    // Compute constraint value at x.
    ROL::Ptr<Vector<Real> > c = workspace_.clone(jv);
    this->update(u_new, u_old, z);
    this->value(*c, u_new, u_old, z, tol);
 
    // Compute (Jacobian at x) times (vector v).
    ROL::Ptr<Vector<Real> > Jv = workspace_.clone(jv);
    this->applyJacobian_1_new(*Jv, v, u_new, u_old, z, tol);
    Real normJv = Jv->norm();
 
    // Temporary vectors.
    ROL::Ptr<Vector<Real> > cdif = workspace_.clone(jv);
    ROL::Ptr<Vector<Real> > cnew = workspace_.clone(jv);
    ROL::Ptr<Vector<Real> > u_2  = workspace_.clone(u_new);
 
    for (int i=0; i<numSteps; i++) {
 
      Real eta = steps[i];
 
      u_2->set(u_new);
 
      cdif->set(*c);
      cdif->scale(weights[order-1][0]);
 
      for(int j=0; j<order; ++j) {
 
         u_2->axpy(eta*shifts[order-1][j], v);

         if( weights[order-1][j+1] != 0 ) {
             this->update(*u_2,u_old,z);
             this->value(*cnew,*u_2,u_old,z,tol);
             cdif->axpy(weights[order-1][j+1],*cnew);
         }

      }
 
      cdif->scale(one/eta);
 
      // Compute norms of Jacobian-vector products, finite-difference approximations, and error.
      jvCheck[i][0] = eta;
      jvCheck[i][1] = normJv;
      jvCheck[i][2] = cdif->norm();
      cdif->axpy(-one, *Jv);
      jvCheck[i][3] = cdif->norm();
 
      if (printToStream) {
        std::stringstream hist;
        if (i==0) {
        hist << std::right
             << std::setw(20) << "Step size"
             << std::setw(20) << "norm(Jac*vec)"
             << std::setw(20) << "norm(FD approx)"
             << std::setw(20) << "norm(abs error)"
             << "\n"
             << std::setw(20) << "---------"
             << std::setw(20) << "-------------"
             << std::setw(20) << "---------------"
             << std::setw(20) << "---------------"
             << "\n";
        }
        hist << std::scientific << std::setprecision(11) << std::right
             << std::setw(20) << jvCheck[i][0]
             << std::setw(20) << jvCheck[i][1]
             << std::setw(20) << jvCheck[i][2]
             << std::setw(20) << jvCheck[i][3]
             << "\n";
        outStream << hist.str();
      }
 
    }
 
    // Reset format state of outStream.
    outStream.copyfmt(oldFormatState);
 
    return jvCheck;
  } // checkApplyJacobian_1_new


}; // class Constraint_SimOpt

} // namespace ROL

#endif
