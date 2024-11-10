// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ADVDIFFTEST_BRANCHHELPER_H
#define ROL_ADVDIFFTEST_BRANCHHELPER_H

#include "ROL_PEBBL_StdBranchHelper.hpp"
#include "ROL_PEBBL_TpetraBranchHelper.hpp"
#include "transform.hpp"

template <class Real>
class StdAdvDiffBranchHelper : public ROL::PEBBL::StdBranchHelper<Real> {
private:
  ROL::Ptr<const ROL::StdVector<Real>> getParameter(const ROL::Vector<Real> &x) const {
    try {
      return ROL::makePtrFromRef(dynamic_cast<const ROL::StdVector<Real>&>(x));
    }
    catch (std::exception &e) {
      return dynamic_cast<const PDE_OptVector<Real>&>(x).getParameter();
    }
  }

public:
  StdAdvDiffBranchHelper(const Real tol = 1e-6, const int method = 0)
    : ROL::PEBBL::StdBranchHelper<Real>(tol, method) {}

  StdAdvDiffBranchHelper(const StdAdvDiffBranchHelper &BH)
    : ROL::PEBBL::StdBranchHelper<Real>(BH) {}

  //int getMyIndex(const ROL::Vector<Real> &x) const {
  int getMyIndex(const ROL::Vector<Real> &x, const ROL::Vector<Real> &g) const {
    // Use Std implementation
    return ROL::PEBBL::StdBranchHelper<Real>::getMyIndex(*getParameter(x),*getParameter(g));
  }

  void getMyNumFrac(int &nfrac, Real &integralityMeasure,
                    const ROL::Vector<Real> &x) const {
    // Use Std implementation
    ROL::PEBBL::StdBranchHelper<Real>::getMyNumFrac(nfrac, integralityMeasure, *getParameter(x));
  }

  ROL::Ptr<ROL::PEBBL::IntegerTransformation<Real>> createTransform(void) const {
    return ROL::makePtr<StdAdvDiffIntegerTransformation<Real>>();
  }

}; // class StdAdvDiffBranchHelper

template <class Real>
class TpetraAdvDiffBranchHelper : public ROL::PEBBL::TpetraBranchHelper<Real> {
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
  TpetraAdvDiffBranchHelper(const Real tol = 1e-6, const int method = 0)
    : ROL::PEBBL::TpetraBranchHelper<Real>(tol, method) {}

  TpetraAdvDiffBranchHelper(const TpetraAdvDiffBranchHelper &BH)
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
    return ROL::makePtr<TpetraAdvDiffIntegerTransformation<Real>>();
  }

}; // class Tpetra_AdvDiff_BranchHelper

#endif
