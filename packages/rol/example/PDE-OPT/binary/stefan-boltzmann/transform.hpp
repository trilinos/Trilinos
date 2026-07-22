// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STEFANBOLTZMANN_INTEGERTRANSFORMATION_H
#define ROL_STEFANBOLTZMANN_INTEGERTRANSFORMATION_H

#include "ROL_PEBBL_TpetraIntegerTransformation.hpp"
#include "ROL_PEBBL_StdIntegerTransformation.hpp"
#include "../../TOOLS/pdevector.hpp"

template <class Real>
class TpetraStefanBoltzmannIntegerTransformation : public ROL::PEBBL::TpetraIntegerTransformation<Real> {
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
  TpetraStefanBoltzmannIntegerTransformation(void)
    : ROL::PEBBL::TpetraIntegerTransformation<Real>() {}

  TpetraStefanBoltzmannIntegerTransformation(const TpetraStefanBoltzmannIntegerTransformation &T)
    : ROL::PEBBL::TpetraIntegerTransformation<Real>(T) {}

  void fixValues(ROL::Vector<Real> &c, bool zero = false) const {
    ROL::PEBBL::TpetraIntegerTransformation<Real>::fixValues(*getData(c),zero);
  }

}; // class TpetraStefanBoltzmannIntegerTransformation

#endif
