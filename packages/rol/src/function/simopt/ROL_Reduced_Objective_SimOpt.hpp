// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_REDUCED_OBJECTIVE_SIMOPT_H
#define ROL_REDUCED_OBJECTIVE_SIMOPT_H

#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"
#include "ROL_VectorController.hpp"
#include "ROL_BatchManager.hpp"

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
  bool isUpdated_;

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
  void summarize(std::ostream &stream, const Ptr<BatchManager<Real>> &bman = nullPtr) const;

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
