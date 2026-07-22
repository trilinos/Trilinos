// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PEBBL_STDINTEGERTRANSFORMATION_H
#define ROL_PEBBL_STDINTEGERTRANSFORMATION_H

#include "ROL_PEBBL_IntegerTransformation.hpp"

/** @ingroup func_group
    \class ROL::PEBBL:IntegerStdTransformation
    \brief Defines the pebbl integer transformation operator interface for StdVectors.

    ROL's pebbl constraint interface is designed to set individual components
    of a vector to a fixed value.  The range space is the same as the domain.

    ---
*/


namespace ROL {
namespace PEBBL {

template <class Real>
class StdIntegerTransformation : public IntegerTransformation<Real> {
private:
  Ptr<std::vector<Real>> getData(Vector<Real> &x) const {
    return dynamic_cast<StdVector<Real>&>(x).getVector();
  }

 using IntegerTransformation<Real>::map_; 

public:
  StdIntegerTransformation(void)
    : IntegerTransformation<Real>() {}

  StdIntegerTransformation(const StdIntegerTransformation &T)
    : IntegerTransformation<Real>(T) {}

  void fixValues(Vector<Real> &c, bool zero = false) const {
    for (auto it=map_.begin(); it!=map_.end(); ++it) {
      (*getData(c))[it->first] = (zero ? static_cast<Real>(0) : it->second);
    }
  }

}; // class StdIntegerTransformation

} // namespace PEBBL
} // namespace ROL

#endif
