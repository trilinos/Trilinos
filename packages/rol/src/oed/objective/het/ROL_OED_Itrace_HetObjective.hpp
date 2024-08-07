// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_ITRACE_HET_OBJECTIVE_HPP
#define ROL_OED_ITRACE_HET_OBJECTIVE_HPP

#include "ROL_OED_HetObjectiveBase.hpp"
#include "ROL_SampleGenerator.hpp"

namespace ROL {
namespace OED {
namespace Het {

template<typename Real>
class Itrace_Objective : public ObjectiveBase<Real,std::vector<Real>> {
private:
  const Ptr<BilinearConstraint<Real>> M_;
  const Ptr<SampleGenerator<Real>> sampler_;
  const std::vector<Real> weight_;
  const Ptr<Vector<Real>> state_, adjoint1_, adjoint2_,
                                    sens1_, sens2_, sadj1_, sadj2_, rhs_, rhs1_;
  const Ptr<VectorController<Real,int>> stateStore_, adjoint1Store_, adjoint2Store_;
  std::vector<Ptr<Vector<Real>>> b_;
  Ptr<Vector<Real>> Xa_;
  const bool storage_;

  using ObjectiveBase<Real,std::vector<Real>>::setConstraint;
  using ObjectiveBase<Real,std::vector<Real>>::setObjective;
  using ObjectiveBase<Real,std::vector<Real>>::setStorage;
  using ObjectiveBase<Real,std::vector<Real>>::initialize;
  using ObjectiveBase<Real,std::vector<Real>>::getConstraint;
  using ObjectiveBase<Real,std::vector<Real>>::getObjective;
  using ObjectiveBase<Real,std::vector<Real>>::getState;
  using ObjectiveBase<Real,std::vector<Real>>::getStateSens;
  using ObjectiveBase<Real,std::vector<Real>>::solve_state_equation;
  using ObjectiveBase<Real,std::vector<Real>>::solve_state_sensitivity;

  void solveStateEquation(Vector<Real> &state, const Vector<Real> &u,
                          const Vector<Real> &z, int i, Real &tol);
  void solveAdjointEquation(Vector<Real> &adjoint, const Vector<Real> &state,
                            const Vector<Real> &u, const Vector<Real> &z,
                            VectorController<Real,int> &store,
                            int i, Real &tol);
  void solveStateSensitivityEquation(Vector<Real> &sens,
                               const Vector<Real> &v, const Vector<Real> &u,
                               const Vector<Real> &z, Real &tol);
  void solveAdjointSensitivityEquation(Vector<Real> &sens,
                                 const Vector<Real> &v,
                                 const Vector<Real> &p, const Vector<Real> &s,
                                 const Vector<Real> &u, const Vector<Real> &z, Real &tol);

public:
  Itrace_Objective( const Ptr<BilinearConstraint<Real>>     &con,
                    const Ptr<QuadraticObjective<Real>>   &obj,
                    const Ptr<Vector<Real>>               &state,
                    const Ptr<SampleGenerator<Real>>      &sampler,
                    const std::vector<Real>               &weight,
                    bool storage = true);

  void update( const Vector<Real> &z,
               UpdateType type,
               int iter = -1 ) override;
  Real value( const Vector<Real> &z, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) override;
};

} // END Het Namespace
} // END OED Namespace
} // END ROL Namespace

#include "ROL_OED_Itrace_HetObjective_Def.hpp"

#endif
