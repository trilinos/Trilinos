// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALARFUNCTION_H
#define ROL_SCALARFUNCTION_H

/** \class ROL::ScalarFunction
    \brief Provides interface for functions that map scalars to scalars.
*/

#include "ROL_Types.hpp"

namespace ROL { 

template<class Real>
class ScalarFunction {
public:
  virtual ~ScalarFunction() {}
  virtual Real value(const Real alpha) = 0;
  virtual Real deriv(const Real alpha) {
    Real val1 = value(alpha);
    Real eta  = std::sqrt(ROL_EPSILON<Real>());
    Real h    = eta*(std::abs(alpha) > eta ? std::abs(alpha) : 1.0);
    Real val0 = value(alpha+h);
    return (val0-val1)/h;
  }
};

}

#endif
