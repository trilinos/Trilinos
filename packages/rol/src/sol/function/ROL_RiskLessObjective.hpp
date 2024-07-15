// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RISKLESSOBJECTIVE_HPP
#define ROL_RISKLESSOBJECTIVE_HPP

#include "ROL_RiskVector.hpp"
#include "ROL_Objective.hpp"

namespace ROL {

template<typename Real>
class RiskLessObjective : public Objective<Real> {
private:
  const Ptr<Objective<Real>> obj_;

public:
  RiskLessObjective(const Ptr<Objective<Real>> &obj);

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v,
                const Vector<Real> &x, Real &tol ) override;
  void precond( Vector<Real> &Pv, const Vector<Real> &v,
                const Vector<Real> &x, Real &tol ) override;
  void setParameter(const std::vector<Real> &param) override;
};

}

#include "ROL_RiskLessObjective_Def.hpp"

#endif
