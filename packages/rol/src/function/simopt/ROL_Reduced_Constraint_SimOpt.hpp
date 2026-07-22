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
#include "ROL_BatchManager.hpp"

namespace ROL {

template <class Real>
class Reduced_Constraint_SimOpt : public Constraint<Real> {
private:
  const ROL::Ptr<Constraint_SimOpt<Real>> conVal_, conRed_;
  const ROL::Ptr<VectorController<Real>>  stateStore_, adjointStore_;

  // Primal vectors
  const ROL::Ptr<Vector<Real>> state_, adjoint_, residual_;
  const ROL::Ptr<Vector<Real>> state_sens_, adjoint_sens_;

  // Dual vectors
  const ROL::Ptr<Vector<Real>> dualstate_, dualstate1_, dualadjoint_;
  const ROL::Ptr<Vector<Real>> dualcontrol_, dualresidual_;

  const bool storage_;
  const bool useFDhessVec_;

  unsigned nupda_, nvalu_, njaco_, najac_, nhess_;
  unsigned nstat_, nadjo_, nssen_, nasen_;

  bool updateFlag_;
  int  updateIter_;
  UpdateType updateType_;
  bool newUpdate_;
  bool isUpdated_;

  void solve_state_equation(const Vector<Real> &z, Real &tol);

  /** \brief Given \f$(u,z)\in\mathcal{U}\times\mathcal{Z}\f$ which solves the state equation, solve 
             the adjoint equation 
             \f$c_u(u,z)^*\lambda + c_u(u,z)^*w = 0\f$ for \f$\lambda=\lambda(u,z)\in\mathcal{C}^*\f$.
  */
  void solve_adjoint_equation(const Vector<Real> &w, const Vector<Real> &z, Real &tol);

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
  void solve_adjoint_sensitivity(const Vector<Real> &w, const Vector<Real> &v, const Vector<Real> &z, Real &tol);

public:
  /** \brief Constructor.

      @param[in] conVal       is a pointer to a SimOpt constraint, to be evaluated.
      @param[in] conRed       is a pointer to a SimOpt constraint, to be reduced.
      @param[in] stateStore   is a pointer to a VectorController object.
      @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
      @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}_{\text{red}}^*\f$.
      @param[in] residual     is a pointer to a primal constraint space vector, \f$\mathcal{C}_{\text{val}}\f$.
      @param[in] storage      is a flag whether or not to store computed states and adjoints.
      @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
  */
  Reduced_Constraint_SimOpt(
      const ROL::Ptr<Constraint_SimOpt<Real>> &conVal,
      const ROL::Ptr<Constraint_SimOpt<Real>> &conRed,
      const ROL::Ptr<VectorController<Real>> &stateStore,
      const ROL::Ptr<Vector<Real>> &state,
      const ROL::Ptr<Vector<Real>> &control,
      const ROL::Ptr<Vector<Real>> &adjoint,
      const ROL::Ptr<Vector<Real>> &residual,
      bool storage = true,
      bool useFDhessVec = false);

  /** \brief Secondary, general constructor for use with dual optimization vector spaces
             where the user does not define the dual() method.

      @param[in] conVal       is a pointer to a SimOpt constraint, to be evaluated.
      @param[in] conRed       is a pointer to a SimOpt constraint, to be reduced.
      @param[in] stateStore   is a pointer to a VectorController object.
      @param[in] state        is a pointer to a state space vector, \f$\mathcal{U}\f$.
      @param[in] control      is a pointer to a optimization space vector, \f$\mathcal{Z}\f$.
      @param[in] adjoint      is a pointer to a dual constraint space vector, \f$\mathcal{C}_{\text{red}}^*\f$.
      @param[in] residual     is a pointer to a primal constraint space vector, \f$\mathcal{C}_{\text{val}}\f$.
      @param[in] dualstate    is a pointer to a dual state space vector, \f$\mathcal{U}^*\f$.
      @param[in] dualadjoint  is a pointer to a constraint space vector, \f$\mathcal{C}_{\text{red}}\f$.
      @param[in] dualresidual is a pointer to a dual constraint space vector, \f$\mathcal{C}_{\text{val}}^*\f$.
      @param[in] storage      is a flag whether or not to store computed states and adjoints.
      @param[in] useFDhessVec is a flag whether or not to use a finite-difference Hessian approximation.
  */
  Reduced_Constraint_SimOpt(
      const ROL::Ptr<Constraint_SimOpt<Real>> &conVal,
      const ROL::Ptr<Constraint_SimOpt<Real>> &conRed,
      const ROL::Ptr<VectorController<Real>> &stateStore,
      const ROL::Ptr<Vector<Real>> &state,
      const ROL::Ptr<Vector<Real>> &control,
      const ROL::Ptr<Vector<Real>> &adjoint,
      const ROL::Ptr<Vector<Real>> &residual,
      const ROL::Ptr<Vector<Real>> &dualstate,
      const ROL::Ptr<Vector<Real>> &dualcontrol,
      const ROL::Ptr<Vector<Real>> &dualadjoint,
      const ROL::Ptr<Vector<Real>> &dualresidual,
      bool storage = true,
      bool useFDhessVec = false);

  void summarize(std::ostream &stream, const Ptr<BatchManager<Real>> &bman = nullPtr) const;

  void reset();

  /** \brief Update the SimOpt objective function and equality constraint.
  */
  void update( const Vector<Real> &z, bool flag = true, int iter = -1 );
  void update( const Vector<Real> &z, UpdateType type, int iter = -1 );

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the equality constraint 
             \f$\widehat{c}(z) = c(u(z),z)\f$ where 
             \f$u=u(z)\in\mathcal{U}\f$ solves \f$e(u,z) = 0\f$.
  */
  void value( Vector<Real> &c, const Vector<Real> &z, Real &tol );

  /** \brief Given \f$z\in\mathcal{Z}\f$, apply the Jacobian to a vector
             \f$\widehat{c}'(z)v = c_u(u,z)s + c_z(u,z)v\f$ where 
             \f$s=s(u,z,v)\in\mathcal{U}^*\f$ solves 
             \f$e_u(u,z)s+e_z(u,z)v = 0\f$.
  */
  void applyJacobian( Vector<Real> &jv, const Vector<Real> &v,
                      const Vector<Real> &z, Real &tol );

  void applyAdjointJacobian( Vector<Real> &ajw, const Vector<Real> &w,
                             const Vector<Real> &z, Real &tol );

  /** \brief Given \f$z\in\mathcal{Z}\f$, evaluate the Hessian of the objective function 
             \f$\nabla^2\widehat{J}(z)\f$ in the direction \f$v\in\mathcal{Z}\f$.
  */
  void applyAdjointHessian( Vector<Real> &ahwv, const Vector<Real> &w,
                            const Vector<Real> &v, const Vector<Real> &z,
                            Real &tol );

// For parametrized (stochastic) objective functions and constraints
public:
  void setParameter(const std::vector<Real> &param) {
    Constraint<Real>::setParameter(param);
    conVal_->setParameter(param);
    conRed_->setParameter(param);
  }
}; // class Reduced_Constraint_SimOpt

} // namespace ROL

#include "ROL_Reduced_Constraint_SimOpt_Def.hpp"

#endif
