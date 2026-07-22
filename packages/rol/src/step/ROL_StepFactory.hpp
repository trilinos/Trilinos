// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STEPFACTORY_H
#define ROL_STEPFACTORY_H

#include "ROL_Types.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"

#include "ROL_Step.hpp"
#include "ROL_LineSearchStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_PrimalDualActiveSetStep.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_AugmentedLagrangianStep.hpp"
#include "ROL_MoreauYosidaPenaltyStep.hpp"
#include "ROL_BundleStep.hpp"
#include "ROL_InteriorPointStep.hpp"
#include "ROL_FletcherStep.hpp"

namespace ROL {

  template<class Real>
  class StepFactory {
    public:
    ~StepFactory(void){}

    ROL::Ptr<Step<Real>> getStep(const std::string &type,
                                 ROL::ParameterList &parlist) const {
      EStep els = StringToEStep(type);
      switch(els) {
        case STEP_AUGMENTEDLAGRANGIAN: return ROL::makePtr<AugmentedLagrangianStep<Real>>(parlist);
        case STEP_BUNDLE:              return ROL::makePtr<BundleStep<Real>>(parlist);
        case STEP_COMPOSITESTEP:       return ROL::makePtr<CompositeStep<Real>>(parlist);
        case STEP_LINESEARCH:          return ROL::makePtr<LineSearchStep<Real>>(parlist);
        case STEP_MOREAUYOSIDAPENALTY: return ROL::makePtr<MoreauYosidaPenaltyStep<Real>>(parlist);
        case STEP_PRIMALDUALACTIVESET: return ROL::makePtr<PrimalDualActiveSetStep<Real>>(parlist);
        case STEP_TRUSTREGION:         return ROL::makePtr<TrustRegionStep<Real>>(parlist);
        case STEP_INTERIORPOINT:       return ROL::makePtr<InteriorPointStep<Real>>(parlist); 
        case STEP_FLETCHER:            return ROL::makePtr<FletcherStep<Real>>(parlist);
        default:                       return ROL::nullPtr;
      }
    }
  };

}

#endif
