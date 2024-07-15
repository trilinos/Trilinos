// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for std::vector bound constraints.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_STDBOUNDCONSTRAINT_HPP
#define ROL_STDBOUNDCONSTRAINT_HPP

#include "ROL_StdVector.hpp"
#include "ROL_BoundConstraint.hpp"

namespace ROL {

template<class Real>
class StdBoundConstraint : public BoundConstraint<Real> {
private:
  int dim_;

  std::vector<Real> x_lo_;
  std::vector<Real> x_up_;

  const Real scale_;
  const Real feasTol_;
  
  using BoundConstraint<Real>::lower_;
  using BoundConstraint<Real>::upper_;

  Real min_diff_;

  inline Real buildC(int i) const {
    const Real zeta(0.5), kappa(1);
    return std::min(zeta*(x_up_[i] - x_lo_[i]), kappa);
  }

  inline Real sgn(Real x) const {
    const Real zero(0), one(1);
    return x > zero ? one : (x < zero ? -one : zero);
  }

  void buildScalingFunction(Vector<Real> &d, const Vector<Real> &x, const Vector<Real> &g) const;

public:
  StdBoundConstraint(std::vector<Real> &x,
                     bool isLower = false,
                     Real scale = Real(1),
                     const Real feasTol = std::sqrt(ROL_EPSILON<Real>()));

  StdBoundConstraint(std::vector<Real> &l,
                     std::vector<Real> &u,
                     Real scale = Real(1),
                     const Real feasTol = std::sqrt(ROL_EPSILON<Real>()));

  void project( Vector<Real> &x ) override;

  void projectInterior( Vector<Real> &x ) override;

  void pruneUpperActive(Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0)) override;

  void pruneUpperActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) override;

  void pruneLowerActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0)) override;

  void pruneLowerActive(Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0)) override;

  bool isFeasible( const Vector<Real> &v ) override;

  void applyInverseScalingFunction( Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const override;

  void applyScalingFunctionJacobian(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const override;
};

}// End ROL Namespace

#include "ROL_StdBoundConstraint_Def.hpp"

#endif
