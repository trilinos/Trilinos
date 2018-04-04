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

#ifndef ROL_CONSTRAINT_TIMESIMOPT_H
#define ROL_CONSTRAINT_TIMESIMOPT_H

#include "ROL_Constraint_SimOpt.hpp"

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

    ROL::Ptr<Vector<Real> > jv_new = jv.clone();

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
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
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
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
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
    ROL::Ptr<ROL::Vector<Real> > r = c.clone();
    ROL::Ptr<ROL::Vector<Real> > s = u.clone();
    s->set(u);
    solve(*r,*s,z,tol);
    // Evaluate constraint residual at (u,z).
    ROL::Ptr<ROL::Vector<Real> > cs = c.clone();
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
  

}; // class Constraint_SimOpt

} // namespace ROL

#endif
