// @ass ROL::PEBBL::TpetraIntegerTransformation
// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ADVDIFF_BRANCHHELPER_H
#define ROL_ADVDIFF_BRANCHHELPER_H

#include "ROL_PEBBL_StdBranchHelper.hpp"
#include "transform.hpp"

template <class Real>
class AdvDiffBranchHelper : public ROL::PEBBL::StdBranchHelper<Real> {
private:
  ROL::Ptr<const ROL::StdVector<Real>> getParameter(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const PDE_OptVector<Real>&>(x).getParameter();
  }

public:
  AdvDiffBranchHelper(const Real tol = 1e-6, const int method = 0)
    : ROL::PEBBL::StdBranchHelper<Real>(tol, method) {}

  AdvDiffBranchHelper(const AdvDiffBranchHelper &BH)
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
    return ROL::makePtr<AdvDiffIntegerTransformation<Real>>();
  }

}; // class AdvDiffBranchHelper

#endif
