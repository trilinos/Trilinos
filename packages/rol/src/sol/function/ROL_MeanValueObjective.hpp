// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MEANVALUEOBJECTIVE_HPP
#define ROL_MEANVALUEOBJECTIVE_HPP

#include "ROL_Objective.hpp"
#include "ROL_SampleGenerator.hpp"

namespace ROL {

template<typename Real>
class MeanValueObjective : public Objective<Real> {
private:
  const Ptr<Objective<Real>> obj_;

public:
  MeanValueObjective( const Ptr<Objective<Real>> &obj,
                      const Ptr<SampleGenerator<Real>> &sampler);

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

private:
  std::vector<Real> computeSampleMean(const Ptr<SampleGenerator<Real>> &sampler) const;
};

}

#include "ROL_MeanValueObjective_Def.hpp"

#endif
