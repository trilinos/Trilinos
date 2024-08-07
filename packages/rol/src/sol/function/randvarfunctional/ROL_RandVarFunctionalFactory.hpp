// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RANDVARFUNCTIONALFACTORY_HPP
#define ROL_RANDVARFUNCTIONALFACTORY_HPP

#include "ROL_RiskMeasureFactory.hpp"
#include "ROL_DeviationMeasureFactory.hpp"
#include "ROL_ErrorMeasureFactory.hpp"
#include "ROL_RegretMeasureFactory.hpp"
#include "ROL_ProbabilityFactory.hpp"

namespace ROL {

  template<class Real>
  inline Ptr<RandVarFunctional<Real> > RandVarFunctionalFactory(ROL::ParameterList &parlist) {
    std::string type = parlist.sublist("SOL").get("Type","Risk Averse");
    if (type == "Risk Averse") {
      return RiskMeasureFactory<Real>(parlist);
    }
    else if (type == "Deviation") {
      return DeviationMeasureFactory<Real>(parlist);
    }
    else if (type == "Error") {
      return ErrorMeasureFactory<Real>(parlist);
    }
    else if (type == "Regret") {
      return RegretMeasureFactory<Real>(parlist);
    }
    else if (type == "Probability") {
      return ProbabilityFactory<Real>(parlist);
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        ">>> (ROL::RandVarFunctionalFactory): Invalid random variable functional type!");
    }
  }
}
#endif
