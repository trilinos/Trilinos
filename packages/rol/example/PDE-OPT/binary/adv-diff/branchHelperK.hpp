// @ass ROL::PEBBL::TpetraIntegerTransformation
// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ADVDIFF_BRANCHHELPERK_H
#define ROL_ADVDIFF_BRANCHHELPERK_H

#include "ROL_PEBBL_StdBranchHelper.hpp"
#include "transformK.hpp"

template <class Real>
class AdvDiffBranchHelper : public ROL::PEBBL::StdBranchHelper<Real> {
private:
  const ROL::StdVector<Real> getParameter(const ROL::Vector<Real> &x) const {
    return dynamic_cast<const ROL::StdVector<Real>&>(x);
  }

public:
  AdvDiffBranchHelper(const Real tol = 1e-6, const int method = 0)
    : ROL::PEBBL::StdBranchHelper<Real>(tol, method) {}

  AdvDiffBranchHelper(const AdvDiffBranchHelper &BH)
    : ROL::PEBBL::StdBranchHelper<Real>(BH) {}

  //int getMyIndex(const ROL::Vector<Real> &x) const {
  int getMyIndex(const ROL::Vector<Real> &x, const ROL::Vector<Real> &g) const {
    // Use Std implementation
    return ROL::PEBBL::StdBranchHelper<Real>::getMyIndex(x,g);
  }

  void getMyNumFrac(int &nfrac, Real &integralityMeasure,
                    const ROL::Vector<Real> &x) const {
    // Use Std implementation
    ROL::PEBBL::StdBranchHelper<Real>::getMyNumFrac(nfrac, integralityMeasure,x);
  }

  ROL::Ptr<ROL::PEBBL::IntegerTransformation<Real>> createTransform(void) const {
    return ROL::makePtr<AdvDiffIntegerTransformation<Real>>();
  }

}; // class AdvDiffBranchHelper

#endif
