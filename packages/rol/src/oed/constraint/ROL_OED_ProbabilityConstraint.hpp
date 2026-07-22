// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_PROBABILITY_CONSTRAINT_HPP
#define ROL_OED_PROBABILITY_CONSTRAINT_HPP

#include "ROL_Constraint.hpp"
#include "ROL_ProbabilityVector.hpp"

namespace ROL {
namespace OED {

template<typename Real>
class ProbabilityConstraint : public Constraint<Real> {
private:
  const bool useScale_;
  Real scale_;

  /***************************************************************************/
  /* Begin Accessor Functions                                                */
  /***************************************************************************/
  std::vector<Real>& getData(Vector<Real> &x) const;
  const std::vector<Real>& getConstData(const Vector<Real> &x) const;
  void sumAll(Real *input, Real *output, int size, const Vector<Real> &x) const;
  /***************************************************************************/
  /* End Accessor Functions                                                  */
  /***************************************************************************/

public:
  ProbabilityConstraint(const Vector<Real> &p,
                        bool useScale = true,
                        Real scale = Real(-1));

  void value(Vector<Real> &c,const Vector<Real> &x,Real &tol) override;
  void applyJacobian(Vector<Real> &jv,const Vector<Real> &v,
                     const Vector<Real> &x,Real &tol) override;
  void applyAdjointJacobian(Vector<Real> &ajv,const Vector<Real> &v,
                            const Vector<Real> &x,Real &tol) override;
  void applyAdjointHessian(Vector<Real> &ahuv,const Vector<Real> &u,
                           const Vector<Real> &v,const Vector<Real> &x,Real &tol) override;
};

} // End OED Namespace
} // End ROL Namespace

#include "ROL_OED_ProbabilityConstraint_Def.hpp"

#endif
