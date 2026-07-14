// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_BINARYCONSTRAINT_HPP
#define ROL_OED_BINARYCONSTRAINT_HPP

#include "ROL_StdConstraint.hpp"
#include "ROL_OED_ProfiledClass.hpp"

namespace ROL::OED {

template<typename Real>
class BinaryConstraint : public StdConstraint<Real>, public ProfiledClass<Real,std::string> {
private:
  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

public:
  BinaryConstraint() : ProfiledClass<Real,std::string>("OED::BinaryConstraint") {}

  void value(std::vector<Real> &c,const std::vector<Real> &x,Real &tol) override;
  void applyJacobian(std::vector<Real> &jv,const std::vector<Real> &v,
                     const std::vector<Real> &x,Real &tol) override;
  void applyAdjointJacobian(std::vector<Real> &ajv,const std::vector<Real> &v,
                            const std::vector<Real> &x,Real &tol) override;
  void applyAdjointHessian(std::vector<Real> &ahwv,const std::vector<Real> &w,
                           const std::vector<Real> &v,
                           const std::vector<Real> &x,Real &tol) override;
}; // class BinaryConstraint

} // End ROL::OED Namespace

#include "ROL_OED_BinaryConstraint_Def.hpp"

#endif
