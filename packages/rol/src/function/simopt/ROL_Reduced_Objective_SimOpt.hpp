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
#include "ROL_VectorController.hpp"

namespace ROL {

template<typename Real>
class Reduced_Objective_SimOpt : public Objective<Real> {
private:
  const Ptr<Objective_SimOpt<Real>> obj_;
  const Ptr<Constraint_SimOpt<Real>> con_;
  Ptr<VectorController<Real>> stateStore_;
  Ptr<VectorController<Real>> adjointStore_;

  // Primal vectors
  Ptr<Vector<Real>> state_;
  Ptr<Vector<Real>> adjoint_;
  Ptr<Vector<Real>> state_sens_;
  Ptr<Vector<Real>> adjoint_sens_;

  // Dual vectors
  Ptr<Vector<Real>> dualstate_;
  Ptr<Vector<Real>> dualstate1_;
  Ptr<Vector<Real>> dualadjoint_;
  Ptr<Vector<Real>> dualcontrol_;

  const bool storage_;
  const bool useFDhessVec_;

  unsigned nupda_, nvalu_, ngrad_, nhess_, nprec_;
  unsigned nstat_, nadjo_, nssen_, nasen_;
  
  bool updateFlag_;
  int  updateIter_;
  UpdateType updateType_;
  bool newUpdate_;

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
      const Ptr<Objective_SimOpt<Real>> &obj, 
      const Ptr<Constraint_SimOpt<Real>> &con, 
      const Ptr<Vector<Real>> &state, 
      const Ptr<Vector<Real>> &control, 
      const Ptr<Vector<Real>> &adjoint,
      const bool storage = true,
      const bool useFDhessVec = false);

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
      const Ptr<Objective_SimOpt<Real>> &obj,
      const Ptr<Constraint_SimOpt<Real>> &con,
      const Ptr<Vector<Real>> &state,
      const Ptr<Vector<Real>> &control, 
      const Ptr<Vector<Real>> &adjoint,
      const Ptr<Vector<Real>> &dualstate,
      const Ptr<Vector<Real>> &dualcontrol, 
      const Ptr<Vector<Real>> &dualadjoint,
      const bool storage = true,
      const bool useFDhessVec = false);

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
  Reduced_Objective_SimOpt(
      const Ptr<Objective_SimOpt<Real>> &obj, 
      const Ptr<Constraint_SimOpt<Real>> &con, 
      const Ptr<VectorController<Real>> &stateStore, 
      const Ptr<Vector<Real>> &state, 
      const Ptr<Vector<Real>> &control, 
      const Ptr<Vector<Real>> &adjoint,
      const bool storage = true,
      const bool useFDhessVec = false);

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
  Reduced_Objective_SimOpt(
      const Ptr<Objective_SimOpt<Real>> &obj,
      const Ptr<Constraint_SimOpt<Real>> &con,
      const Ptr<VectorController<Real>> &stateStore, 
      const Ptr<Vector<Real>> &state,
      const Ptr<Vector<Real>> &control, 
      const Ptr<Vector<Real>> &adjoint,
      const Ptr<Vector<Real>> &dualstate,
      const Ptr<Vector<Real>> &dualcontrol, 
      const Ptr<Vector<Real>> &dualadjoint,
      const bool storage = true,
      const bool useFDhessVec = false);

  /** \brief Update the SimOpt objective function and equality constraint.
  */
  void update( const Vector<Real> &z, bool flag = true, int iter = -1 ) override;
  void update( const Vector<Real> &z, UpdateType type, int iter = -1 ) override;

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the objective function 
             \f$\widehat{J}(z) = J(u(z),z)\f$ where 
             \f$u=u(z)\in\mathcal{U}\f$ solves \f$e(u,z) = 0\f$.
  */
  Real value( const Vector<Real> &z, Real &tol ) override;

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the gradient of the objective function 
             \f$\nabla\widehat{J}(z) = J_z(z) + c_z(u,z)^*\lambda\f$ where 
             \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$ solves 
             \f$e_u(u,z)^*\lambda+J_u(u,z) = 0\f$.
  */
  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) override;

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the Hessian of the objective function 
             \f$\nabla^2\widehat{J}(z)\f$ in the direction \f$v\in\mathcal{Z}\f$.
  */
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) override;

  /** \brief Apply a reduced Hessian preconditioner.
  */
  virtual void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) override;

  /** Write summary to stream.
  */
  void summarize(std::ostream &stream) const;

  /** Reset summary data.
  */
  void reset();

// For parametrized (stochastic) objective functions and constraints
public:
  void setParameter(const std::vector<Real> &param) override;

private:
  void solve_state_equation(const Vector<Real> &z, Real &tol); 

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation, solve 
             the adjoint equation 
             \f$c_u(u,z)^*\lambda + J_u(u,z) = 0\f$ for \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_equation(const Vector<Real> &z, Real &tol);

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation and 
             a direction \f$v\in\mathcal{Z}\f$, solve the state senstivity equation 
             \f$c_u(u,z)s + c_z(u,z)v = 0\f$ for \f$s=u_z(z)v\in\mathcal{U}\f$.
  */
  void solve_state_sensitivity(const Vector<Real> &v, const Vector<Real> &z, Real &tol);

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$, the adjoint variable 
             \f$\lambda\in\mathcal{C}^*\f$, and a direction \f$v\in\mathcal{Z}\f$, solve the 
             adjoint sensitvity equation 
             \f$c_u(u,z)^*p + J_{uu}(u,z)s + J_{uz}(u,z)v + c_{uu}(u,z)(\cdot,s)^*\lambda 
                            + c_{zu}(u,z)(\cdot,v)^*\lambda = 0\f$
             for \f$p = \lambda_z(u(z),z)v\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_sensitivity(const Vector<Real> &v, const Vector<Real> &z, Real &tol);

}; // class Reduced_Objective_SimOpt

} // namespace ROL

#include "ROL_Reduced_Objective_SimOpt_Def.hpp"

#endif
