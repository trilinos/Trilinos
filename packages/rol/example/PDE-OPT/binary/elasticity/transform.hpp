// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_MULTIMAT_TRANSFORM_H
#define ROL_PDEOPT_MULTIMAT_TRANSFORM_H

#include "ROL_PEBBL_TpetraIntegerTransformation.hpp"
#include "../../TOOLS/pdevector.hpp"

template <class Real>
class TpetraMultiMatIntegerTransformation : public ROL::PEBBL::TpetraIntegerTransformation<Real> {
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
  TpetraMultiMatIntegerTransformation(void)
    : ROL::PEBBL::TpetraIntegerTransformation<Real>() {}

  TpetraMultiMatIntegerTransformation(const TpetraMultiMatIntegerTransformation &T)
    : ROL::PEBBL::TpetraIntegerTransformation<Real>(T) {}

  void fixValues(ROL::Vector<Real> &c, bool zero = false) const {
    ROL::PEBBL::TpetraIntegerTransformation<Real>::fixValues(*getData(c),zero);
  }

}; // class TpetraMultiMatIntegerTransformation

#endif
