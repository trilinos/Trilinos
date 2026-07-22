// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_I_HET_OBJECTIVE_HPP
#define ROL_OED_I_HET_OBJECTIVE_HPP

#include "ROL_OED_HetObjectiveBase.hpp"

namespace ROL {
namespace OED {
namespace Het {

template<typename Real>
class I_Objective : public ObjectiveBase<Real,std::vector<Real>> {
private:
  Ptr<Vector<Real>> dp_;

  using ObjectiveBase<Real,std::vector<Real>>::setConstraint;
  using ObjectiveBase<Real,std::vector<Real>>::setObjective;
  using ObjectiveBase<Real,std::vector<Real>>::setStorage;
  using ObjectiveBase<Real,std::vector<Real>>::initialize;
  using ObjectiveBase<Real,std::vector<Real>>::getConstraint;
  using ObjectiveBase<Real,std::vector<Real>>::getObjective;
  using ObjectiveBase<Real,std::vector<Real>>::getState;
  using ObjectiveBase<Real,std::vector<Real>>::getAdjoint;
  using ObjectiveBase<Real,std::vector<Real>>::getStateSens;
  using ObjectiveBase<Real,std::vector<Real>>::getAdjointSens;
  using ObjectiveBase<Real,std::vector<Real>>::solve_state_equation;
  using ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_equation;
  using ObjectiveBase<Real,std::vector<Real>>::solve_state_sensitivity;
  using ObjectiveBase<Real,std::vector<Real>>::solve_adjoint_sensitivity;

public:
  I_Objective( const Ptr<BilinearConstraint<Real>>   &con,
               const Ptr<QuadraticObjective<Real>> &obj,
               const Ptr<Vector<Real>>        &state,
               bool storage = true);

  Real value( const Vector<Real> &z, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) override;
};

} // END Het Namespace
} // END OED Namespace
} // END ROL Namespace

#include "ROL_OED_I_HetObjective_Def.hpp"

#endif
