// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STATUSFACTORY_H
#define ROL_STATUSFACTORY_H

#include "ROL_Types.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"

#include "ROL_StatusTest.hpp"
#include "ROL_BundleStatusTest.hpp"
#include "ROL_ConstraintStatusTest.hpp"
#include "ROL_FletcherStatusTest.hpp"

namespace ROL {
  template<class Real>
  class StatusTestFactory {
    public:
    ~StatusTestFactory(void){}

    ROL::Ptr<StatusTest<Real> > getStatusTest(const std::string step,
                                                  ROL::ParameterList &parlist) {
      EStep els = StringToEStep(step);
      switch(els) {
        case STEP_BUNDLE:              return ROL::makePtr<BundleStatusTest<Real>>(parlist);
        case STEP_AUGMENTEDLAGRANGIAN: return ROL::makePtr<ConstraintStatusTest<Real>>(parlist);
        case STEP_COMPOSITESTEP:       return ROL::makePtr<ConstraintStatusTest<Real>>(parlist);
        case STEP_MOREAUYOSIDAPENALTY: return ROL::makePtr<ConstraintStatusTest<Real>>(parlist);
        case STEP_INTERIORPOINT:       return ROL::makePtr<ConstraintStatusTest<Real>>(parlist);
        case STEP_LINESEARCH:          return ROL::makePtr<StatusTest<Real>>(parlist);
        case STEP_PRIMALDUALACTIVESET: return ROL::makePtr<StatusTest<Real>>(parlist);
        case STEP_TRUSTREGION:         return ROL::makePtr<StatusTest<Real>>(parlist);
        case STEP_FLETCHER:            return ROL::makePtr<FletcherStatusTest<Real>>(parlist);
        default:                       return ROL::nullPtr;
      } 
    }
  };
}

#endif
