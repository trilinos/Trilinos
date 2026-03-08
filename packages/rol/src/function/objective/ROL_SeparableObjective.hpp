// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SEPARABLEOBJECTIVE_H
#define ROL_SEPARABLEOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Ptr.hpp"
#include "ROL_PartitionedVector.hpp"
#include <vector>

namespace ROL {

template<typename Real>
class SeparableObjective : public Objective<Real> {
private:
  const std::vector<Ptr<Objective<Real>>> obj_;
  size_t size_;

public:
  SeparableObjective(const std::vector<Ptr<Objective<Real>>> &obj);
  const Ptr<Objective<Real>> get(unsigned i) const;
  void update(const Vector<Real> &x, UpdateType type, int iter = -1) override;
  void update(const Vector<Real> &x, bool flag = true, int iter = -1) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void precond( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void setParameter(const std::vector<Real> &param) override;

}; // class SeparableObjective

} // namespace ROL

#include "ROL_SeparableObjective_Def.hpp"

#endif
