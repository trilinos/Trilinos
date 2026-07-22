// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PEBBL_STDBRANCHHELPER_H
#define ROL_PEBBL_STDBRANCHHELPER_H

#include "ROL_StdVector.hpp"
#include "ROL_PEBBL_BranchHelper.hpp"
#include "ROL_PEBBL_StdIntegerTransformation.hpp"

/** @ingroup func_group
    \class ROL::PEBBL::StdBranchHelper
    \brief Defines the pebbl branch index interface for StdVectors.

    ---
*/


namespace ROL {
namespace PEBBL {

template <class Real>
class StdBranchHelper : public BranchHelper<Real> {
private:
  const Real tol_;
  const int method_;

  using BranchHelper<Real>::getIntegerVector;

  Ptr<const std::vector<Real>> getConstData(const Vector<Real> &x) const {
    return dynamic_cast<const StdVector<Real>&>(*getIntegerVector(x)).getVector();
  }

  // Branching based on distance to integer.
  int getIndex_D(const Vector<Real> &x, const Vector<Real> &g) const {
    // Get index closest to 0.5
    Ptr<const std::vector<Real>> xval = getConstData(x);
    int index(0), cnt(0);
    Real minD(0), minX(ROL_INF<Real>()), half(0.5), fx(0), cx(0);
    for (const auto & x : *xval) {
      fx   = std::floor(x);
      cx   = std::ceil(x);
      minD = std::min(x-fx,cx-x);
      if (std::abs(minD-half) < minX) {
        minX = std::abs(minD-half);
        index = cnt;
      }
      cnt++;
    }
    return index;
  }

  // Branching based on directional derivatives (similar to pseudo costs).
  int getIndex_PC(const Vector<Real> &x, const Vector<Real> &g) const {
    Ptr<const std::vector<Real>> xval = getConstData(x);
    Ptr<const std::vector<Real>> gval = getConstData(g);
    Real maxD(ROL_NINF<Real>()), Li(0), Ui(0), mini(0);
    int index = 0, size = gval->size();
    for (int i = 0; i < size; ++i) {
      Li   = (*gval)[i] * (std::floor((*xval)[i]) - (*xval)[i]);
      Ui   = (*gval)[i] * (std::ceil((*xval)[i])  - (*xval)[i]);
      mini = std::min(std::abs(Li),std::abs(Ui));
      if (mini > maxD) {
        maxD  = mini;
        index = i;
      }
    }
    return index;
  }

public:
  StdBranchHelper(Real tol = 1e-6, int method = 0)
    : tol_(tol), method_(method) {}

  StdBranchHelper(const StdBranchHelper &BH)
    : tol_(BH.tol_), method_(BH.method_) {}

  int getIndex(const Vector<Real> &x, const Vector<Real> &g) const {
    int index(0);
    if (method_ == 1) index = getIndex_D(x,g);
    else              index = getIndex_PC(x,g);
    return index;
  }


  void getNumFrac(int &nfrac, Real &integralityMeasure, const Vector<Real> &x) const {
    // Return number of fractional variables and the
    // sum of the distance to integer for the input vector
    Ptr<const std::vector<Real>> xval = getConstData(x);
    nfrac = 0;
    integralityMeasure = static_cast<Real>(0);
    Real minD(0), fx(0), cx(0);
    for (const auto & x : *xval) {
      fx   = std::floor(x);
      cx   = std::ceil(x);
      minD = std::min(x-fx,cx-x);
      integralityMeasure += minD;
      if (minD > tol_) nfrac++;
    }
  }

  Ptr<IntegerTransformation<Real>> createTransform(void) const {
    return makePtr<StdIntegerTransformation<Real>>();
  }

}; // class StdBranchHelper

} // namespace PEBBL
} // namespace ROL

#endif
