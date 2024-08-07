// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANVALUECONSTRAINT_HPP
#define ROL_MEANVALUECONSTRAINT_HPP

#include "ROL_Constraint.hpp"
#include "ROL_SampleGenerator.hpp"

namespace ROL {

template<class Real>
class MeanValueConstraint : public Constraint<Real> {
private:
  const Ptr<Constraint<Real>> con_;

public:
  MeanValueConstraint( const Ptr<Constraint<Real>>      &con,
                       const Ptr<SampleGenerator<Real>> &sampler );

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void value(Vector<Real> &c, const Vector<Real> &x, Real &tol ) override;
  void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v, const Vector<Real> &x, Real &tol) override;

private:
  std::vector<Real> computeSampleMean(const Ptr<SampleGenerator<Real>> &sampler) const;
};

}

#include "ROL_MeanValueConstraint_Def.hpp"

#endif
