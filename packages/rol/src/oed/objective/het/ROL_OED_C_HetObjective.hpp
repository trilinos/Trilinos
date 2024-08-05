// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_C_HET_OBJECTIVE_HPP
#define ROL_OED_C_HET_OBJECTIVE_HPP

#include "ROL_OED_HetObjectiveBase.hpp"

namespace ROL {
namespace OED {
namespace Het {

template<typename Real>
class C_Objective : public ObjectiveBase<Real,std::vector<Real>> {
private:
  std::vector<Real> param_;
  Ptr<Vector<Real>> dp_;

public:
  C_Objective( const Ptr<BilinearConstraint<Real>>   &con,
               const Ptr<QuadraticObjective<Real>> &obj,
               const Ptr<Vector<Real>>             &state,
               bool storage = true);

  Real value( const Vector<Real> &z, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) override;
};

} // END Het Namespace
} // END OED Namespace
} // END ROL Namespace

#include "ROL_OED_C_HetObjective_Def.hpp"

#endif
