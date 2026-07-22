// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALARMINIMIZATIONTEST_H
#define ROL_SCALARMINIMIZATIONTEST_H

/** \class ROL::ScalarMinimizationTest
    \brief Tests the minimization of scalar functions.
*/

#include "ROL_BrentsScalarMinimization.hpp"
#include "ROL_BisectionScalarMinimization.hpp"
#include "ROL_GoldenSectionScalarMinimization.hpp"
#include "ROL_ScalarFunction.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"
#include <iostream>

namespace ROL { 

template<class Real>
class ScalarMinimizationTest {
private:
  ROL::Ptr<ScalarMinimization<Real> > algo_;

public:
  virtual ~ScalarMinimizationTest(void) {}

  ScalarMinimizationTest(ROL::ParameterList &parlist) {
    std::string type = parlist.sublist("Scalar Minimization").get("Type","Brent's");
    if ( type == "Brent's" ) {
      algo_ = ROL::makePtr<BrentsScalarMinimization<Real>>(parlist);
    }
    else if ( type == "Bisection" ) {
      algo_ = ROL::makePtr<BisectionScalarMinimization<Real>>(parlist);
    }
    else if ( type == "Golden Section" ) {
      algo_ = ROL::makePtr<GoldenSectionScalarMinimization<Real>>(parlist);
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> (ROL::ScalarMinimizationTest): Undefined ScalarMinimization type!");
    }
  }

  virtual bool test(std::ostream &stream = std::cout) = 0;

protected:
  void run(Real &fx, Real &x, int &nfval, int &ngrad,
            ScalarFunction<Real> &f,
           const Real &A, const Real &B) {
    algo_->run(fx, x, nfval, ngrad, f, A, B);
  }
};

}

#endif
