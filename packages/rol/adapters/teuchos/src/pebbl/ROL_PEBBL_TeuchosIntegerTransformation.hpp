// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PEBBL_TEUCHOSINTEGERTRANSFORMATION_H
#define ROL_PEBBL_TEUCHOSINTEGERTRANSFORMATION_H

#include "ROL_PEBBL_IntegerTransformation.hpp"
#include "ROL_TeuchosVector.hpp"

/** @ingroup func_group
    \class ROL::PEBBL::TeuchosIntegerTransformation
    \brief Defines the pebbl transform operator interface for TeuchosVectors.

    ROL's pebbl constraint interface is designed to set individual components
    of a vector to a fixed value.  The range space is the same as the domain.

    ---
*/


namespace ROL {
namespace PEBBL {

template <class Ordinal, class Real>
class TeuchosIntegerTransformation : public IntegerTransformation<Real> {
private:
  Ptr<Teuchos::SerialDenseVector<Ordinal,Real>> getData(Vector<Real> &x) const {
    return dynamic_cast<TeuchosVector<Ordinal,Real>&>(x).getVector();
  }

  Ptr<const Teuchos::SerialDenseVector<Ordinal,Real>> getConstData(const Vector<Real> &x) const {
    return dynamic_cast<const TeuchosVector<Ordinal,Real>&>(x).getVector();
  }

 using IntegerTransformation<Real>::map_; 

public:
  TeuchosIntegerTransformation(void)
    : IntegerTransformation<Real>() {}

  TeuchosIntegerTransformation(const TeuchosIntegerTransformation &T)
    : IntegerTransformation<Real>(T) {}

  void fixValues(Vector<Real> &c, bool zero = false) const {
    for (auto it=map_.begin(); it!=map_.end(); ++it) {
      (*getData(c))(it->first) = (zero ? static_cast<Real>(0) : it->second);
    }
  }

}; // class TeuchosIntegerTransformation

} // namespace PEBBL
} // namespace ROL

#endif
