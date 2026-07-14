// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_BINARYPENALTYTRANSFORMATION_HPP
#define ROL_OED_BINARYPENALTYTRANSFORMATION_HPP

#include "ROL_StdConstraint.hpp"
#include "ROL_OED_ProfiledClass.hpp"

namespace ROL::OED {

template<typename Real>
class BinaryPenaltyTransformation : public StdConstraint<Real>, public ProfiledClass<Real,std::string> {
private:
  Real pmin_, ppow_;
  const int type_;

  using ProfiledClass<Real,std::string>::startTimer;
  using ProfiledClass<Real,std::string>::stopTimer;

public:
  BinaryPenaltyTransformation(ParameterList& parlist);

  void value(std::vector<Real> &c,const std::vector<Real> &x,Real &tol) override;
  void applyJacobian(std::vector<Real> &jv,const std::vector<Real> &v,
                     const std::vector<Real> &x,Real &tol) override;
  void applyAdjointJacobian(std::vector<Real> &ajv,const std::vector<Real> &v,
                            const std::vector<Real> &x,Real &tol) override;
  void applyAdjointHessian(std::vector<Real> &ahwv,const std::vector<Real> &w,
                           const std::vector<Real> &v,
                           const std::vector<Real> &x,Real &tol) override;
}; // class BinaryPenaltyTransformation

} // End ROL::OED Namespace

#include "ROL_OED_BinaryPenaltyTransformation_Def.hpp"

#endif
