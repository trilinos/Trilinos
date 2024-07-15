// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PEBBL_BUILDTRANSFORMATION_H
#define ROL_PEBBL_BUILDTRANSFORMATION_H

#include "ROL_AffineTransformObjective.hpp"
#include "ROL_AffineTransformConstraint.hpp"
#include "ROL_PEBBL_IntegerTransformation.hpp"

/** @ingroup func_group
    \class ROL::PEBBL::BuildTransformation
    \brief Perform integer transformation for PEBBL.

    ---
*/

namespace ROL {
namespace PEBBL {

template <class Real>
class BuildTransformation {
private:
  const Ptr<IntegerTransformation<Real>> trans_;
  const Ptr<const Vector<Real>>          x_;
  const Ptr<VectorController<Real>>      storage_;

public:
  virtual ~BuildTransformation(void) {}

  BuildTransformation(const Ptr<IntegerTransformation<Real>> &trans,
                      const Ptr<const Vector<Real>>          &x)
    : trans_(trans), x_(x),
      storage_(makePtr<VectorController<Real>>()) {}

  Ptr<Objective<Real>> transform(const Ptr<Objective<Real>> &obj) const {
    return makePtr<AffineTransformObjective<Real>>(obj,trans_,*x_,storage_);
  }

  Ptr<Constraint<Real>> transform(const Ptr<Constraint<Real>> &con) const {
    return makePtr<AffineTransformConstraint<Real>>(con,trans_,*x_,storage_);
  }

}; // class BuildTransformation

} // namespace PEBBL
} // namespace ROL

#endif
