// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_D_HET_OBJECTIVE_HPP
#define ROL_OED_D_HET_OBJECTIVE_HPP

#include "ROL_OED_HetObjectiveBase.hpp"

/*! \file  OED_D_HetObjective.hpp
    \brief Implements the D-optimal criterion.

    Given a covariance operator \f$C(p)\f$ that depends
    on the current design \f$p\f$, this objective function
    implements the logarithm of the determinant of \f$C(p)\f$,
    i.e.,
    \f[
       J(p) = \log(\operatorname{det}(C(p))).
    \f]
*/

namespace ROL {
namespace OED {
namespace Het {

template<typename Real>
class D_Objective : public ObjectiveBase<Real,std::vector<Real>> {
private:
  const Ptr<BilinearConstraint<Real>> con0_;
  const Ptr<VectorController<Real,int>> ustore_;
  const Ptr<Vector<Real>> u_, g_;
  const int dim_;
  Ptr<Vector<Real>> p_;
  Real det0_, det1_;
  bool isDetComputed_;

  using ObjectiveBase<Real,std::vector<Real>>::setConstraint;
  using ObjectiveBase<Real,std::vector<Real>>::setStorage;
  using ObjectiveBase<Real,std::vector<Real>>::setUpdate;
  using ObjectiveBase<Real,std::vector<Real>>::initialize;
  using ObjectiveBase<Real,std::vector<Real>>::getConstraint;
  using ObjectiveBase<Real,std::vector<Real>>::getState;
  using ObjectiveBase<Real,std::vector<Real>>::solve_state_equation;
  using ObjectiveBase<Real,std::vector<Real>>::getStateSens;
  using ObjectiveBase<Real,std::vector<Real>>::solve_state_sensitivity;

  void solve(int key, const Vector<Real> &z, Real &tol);
  void solveSens(const Vector<Real> &v, const Vector<Real> &z, Real &tol);

public:
  D_Objective( const Ptr<BilinearConstraint<Real>> &con0, // Objective Moment Operator
               const Ptr<BilinearConstraint<Real>> &con1, // Constraint Moment Operator
               const Ptr<Vector<Real>>           &state,
               bool storage = true);

  void update(const Vector<Real> &z, UpdateType type, int iter = -1) override;
  Real value( const Vector<Real> &z, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) override;
};

} // END Het Namespace
} // END OED Namespace
} // END ROL Namespace

#include "ROL_OED_D_HetObjective_Def.hpp"

#endif
