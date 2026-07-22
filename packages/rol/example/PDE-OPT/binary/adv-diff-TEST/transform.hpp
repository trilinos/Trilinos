// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ADVDIFFTEST_INTEGERTRANSFORMATION_H
#define ROL_ADVDIFFTEST_INTEGERTRANSFORMATION_H

#include "ROL_PEBBL_StdIntegerTransformation.hpp"
#include "ROL_PEBBL_TpetraIntegerTransformation.hpp"
#include "../../TOOLS/pdevector.hpp"

template <class Real>
class StdAdvDiffIntegerTransformation : public ROL::PEBBL::StdIntegerTransformation<Real> {
private:
  ROL::Ptr<ROL::StdVector<Real>> getParameter(ROL::Vector<Real> &x) const {
    try {
      return ROL::makePtrFromRef(dynamic_cast<ROL::StdVector<Real>&>(x));
    }
    catch (std::exception &e) {
      return dynamic_cast<PDE_OptVector<Real>&>(x).getParameter();
    }
  }

public:
  StdAdvDiffIntegerTransformation(void)
    : ROL::PEBBL::StdIntegerTransformation<Real>() {}

  StdAdvDiffIntegerTransformation(const StdAdvDiffIntegerTransformation &T)
    : ROL::PEBBL::StdIntegerTransformation<Real>(T) {}

  void fixValues(ROL::Vector<Real> &c, bool zero = false) const {
    ROL::PEBBL::StdIntegerTransformation<Real>::fixValues(*getParameter(c),zero);
  }

}; // class StdAdvDiffIntegerTransformation

template <class Real>
class TpetraAdvDiffIntegerTransformation : public ROL::PEBBL::TpetraIntegerTransformation<Real> {
private:
  ROL::Ptr<ROL::TpetraMultiVector<Real>> getData(ROL::Vector<Real> &x) const {
    try {
      return ROL::makePtrFromRef(dynamic_cast<ROL::TpetraMultiVector<Real>&>(x));
    }
    catch (std::exception &e) {
      return dynamic_cast<PDE_OptVector<Real>&>(x).getField();
    }
  }

public:
  TpetraAdvDiffIntegerTransformation(void)
    : ROL::PEBBL::TpetraIntegerTransformation<Real>() {}

  TpetraAdvDiffIntegerTransformation(const TpetraAdvDiffIntegerTransformation &T)
    : ROL::PEBBL::TpetraIntegerTransformation<Real>(T) {}

  void fixValues(ROL::Vector<Real> &c, bool zero = false) const {
    ROL::PEBBL::TpetraIntegerTransformation<Real>::fixValues(*getData(c),zero);
  }

}; // class TpetraAdvDiffIntegerTransformation

#endif
