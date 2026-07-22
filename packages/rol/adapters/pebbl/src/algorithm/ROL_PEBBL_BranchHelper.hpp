// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PEBBL_BRANCHHELPER_H
#define ROL_PEBBL_BRANCHHELPER_H

#include "ROL_PartitionedVector.hpp"
#include "ROL_PEBBL_IntegerTransformation.hpp"
#include "ROL_PEBBL_MixedVector.hpp"

/** @ingroup func_group
    \class ROL::PEBBL::BranchHelper
    \brief Defines the pebbl branch index interface.

    ---
*/

namespace ROL {
namespace PEBBL {

template <class Real>
class BranchHelper {
protected:
  Ptr<const Vector<Real>> getOptVector(const Vector<Real> &xs ) const {
    try {
      return dynamic_cast<const PartitionedVector<Real>&>(xs).get(0);
    }
    catch (std::exception &e) {
      return makePtrFromRef(xs);
    }
  }

  Ptr<const Vector<Real>> getIntegerVector(const Vector<Real> &xs) const {
    try {
      return dynamic_cast<const MixedVector<Real>&>(*getOptVector(xs)).getIntegerVariables();
    }
    catch (std::exception &e) {
      return getOptVector(xs);
    }
  }

public:
  virtual ~BranchHelper(void) {}
  BranchHelper(void) {}
  BranchHelper(const BranchHelper &con) {}

  virtual int getIndex(const Vector<Real> &x, const Vector<Real> &g) const = 0;
  virtual void getNumFrac(int &nfrac, Real &integralityMeasure, const Vector<Real> &x) const = 0;
  virtual Ptr<IntegerTransformation<Real>> createTransform(void) const = 0;

}; // class BranchHelper

} // namespace PEBBL
} // namespace ROL

#endif
