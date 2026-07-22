// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALARMINIMIZATION_H
#define ROL_SCALARMINIMIZATION_H

/** \class ROL::ScalarMinimization
    \brief Provides interface for mimimizing functions
           that map scalars to scalars on a prescribed
           interval.
*/

#include "ROL_ScalarFunction.hpp"
#include "ROL_ScalarMinimizationStatusTest.hpp"

namespace ROL { 

template<class Real>
class ScalarMinimization {
public:
  virtual ~ScalarMinimization() {}
  void run(Real &fx, Real &x, int &nfval, int &ngrad,
           ScalarFunction<Real> &f, const Real A, const Real B) const {
    ScalarMinimizationStatusTest<Real> test;
    run(fx,x,nfval,ngrad,f,A,B,test);
  }
  virtual void run(Real &fx, Real &x, int &nfval, int &ngrad,
                   ScalarFunction<Real> &f, const Real A, const Real B,
                   ScalarMinimizationStatusTest<Real> &test) const = 0;
};

}

#endif
