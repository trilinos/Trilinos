// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STEFANBOLTZMANN_BRANCHHELPER_H
#define ROL_STEFANBOLTZMANN_BRANCHHELPER_H

#include "ROL_PEBBL_TpetraBranchHelper.hpp"
#include "ROL_PEBBL_StdBranchHelper.hpp"
#include "transform.hpp"

template <class Real>
class TpetraStefanBoltzmannBranchHelper : public ROL::PEBBL::TpetraBranchHelper<Real> {
private:
  ROL::Ptr<const ROL::TpetraMultiVector<Real>> getData(const ROL::Vector<Real> &x) const {
    try {
      return ROL::makePtrFromRef(dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x));
    }
    catch (std::exception &e) {
      return dynamic_cast<const PDE_OptVector<Real>&>(x).getField();
    }
  }

public:
  TpetraStefanBoltzmannBranchHelper(const Real tol = 1e-6, const int method = 0)
    : ROL::PEBBL::TpetraBranchHelper<Real>(tol, method) {}

  TpetraStefanBoltzmannBranchHelper(const TpetraStefanBoltzmannBranchHelper &BH)
    : ROL::PEBBL::TpetraBranchHelper<Real>(BH) {}

  //int getMyIndex(const ROL::Vector<Real> &x) const {
  int getMyIndex(const ROL::Vector<Real> &x, const ROL::Vector<Real> &g) const {
    // Use Std implementation
    return ROL::PEBBL::TpetraBranchHelper<Real>::getMyIndex(*getData(x),*getData(g));
  }

  void getMyNumFrac(int &nfrac, Real &integralityMeasure,
                    const ROL::Vector<Real> &x) const {
    // Use Std implementation
    ROL::PEBBL::TpetraBranchHelper<Real>::getMyNumFrac(nfrac, integralityMeasure, *getData(x));
  }

  ROL::Ptr<ROL::PEBBL::IntegerTransformation<Real>> createTransform(void) const {
    return ROL::makePtr<TpetraStefanBoltzmannIntegerTransformation<Real>>();
  }

}; // class Tpetra_StefanBoltzmann_BranchHelper

#endif
