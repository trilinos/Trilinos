// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ADVDIFF_INTEGERTRANSFORMATION_H
#define ROL_ADVDIFF_INTEGERTRANSFORMATION_H

#include "ROL_PEBBL_StdIntegerTransformation.hpp"
#include "../../TOOLS/pdevector.hpp"

template <class Real>
class AdvDiffIntegerTransformation : public ROL::PEBBL::StdIntegerTransformation<Real> {
private:
  ROL::Ptr<ROL::StdVector<Real>> getParameter(ROL::Vector<Real> &x) const {
    return dynamic_cast<PDE_OptVector<Real>&>(x).getParameter();
  }

public:
  AdvDiffIntegerTransformation(void)
    : ROL::PEBBL::StdIntegerTransformation<Real>() {}

  AdvDiffIntegerTransformation(const AdvDiffIntegerTransformation &T)
    : ROL::PEBBL::StdIntegerTransformation<Real>(T) {}

  void fixValues(ROL::Vector<Real> &c, bool zero = false) const {
    ROL::PEBBL::StdIntegerTransformation<Real>::fixValues(*getParameter(c),zero);
  }

}; // class AdvDiffIntegerTransformation

#endif
