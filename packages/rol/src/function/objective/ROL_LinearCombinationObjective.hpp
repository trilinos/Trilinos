// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINEARCOMBINATIONOBJECTIVE_H
#define ROL_LINEARCOMBINATIONOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Ptr.hpp"

namespace ROL {

template<typename Real>
class LinearCombinationObjective : public Objective<Real> {
private:
  const std::vector<Ptr<Objective<Real>>> obj_;
  std::vector<Real> weights_;
  size_t size_;

  Ptr<Vector<Real>> xdual_;
  bool initialized_;

public:
  LinearCombinationObjective(const std::vector<Ptr<Objective<Real>>> &obj);
  LinearCombinationObjective(const std::vector<Real> &weights,
                             const std::vector<Ptr<Objective<Real>>> &obj);

  void update(const Vector<Real> &x, UpdateType type, int iter = -1) override;
  void update(const Vector<Real> &x, bool flag = true, int iter = -1) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void setParameter(const std::vector<Real> &param) override;

}; // class LinearCombinationObjective

} // namespace ROL

#include "ROL_LinearCombinationObjective_Def.hpp"

#endif
