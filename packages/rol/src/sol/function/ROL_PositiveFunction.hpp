// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_POSITIVEFUNCTION_HPP
#define ROL_POSITIVEFUNCTION_HPP

namespace ROL {

template<class Real>
class PositiveFunction {
public: 
  virtual ~PositiveFunction() {}
  virtual Real evaluate(Real input, int deriv) = 0;
};

}

#endif
