// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TRUSTREGIONFACTORY_H
#define ROL_TRUSTREGIONFACTORY_H

#include "ROL_Types.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"

#include "ROL_TrustRegion.hpp"
#include "ROL_CauchyPoint.hpp"
#include "ROL_DogLeg.hpp"
#include "ROL_DoubleDogLeg.hpp"
#include "ROL_TruncatedCG.hpp"
#include "ROL_LinMore.hpp"

namespace ROL {
template<class Real>
  inline ROL::Ptr<TrustRegion<Real> > TrustRegionFactory(ROL::ParameterList &parlist) {
    ETrustRegion etr = StringToETrustRegion(
      parlist.sublist("Step").sublist("Trust Region").get("Subproblem Solver","Dogleg"));
    switch(etr) {
      case TRUSTREGION_CAUCHYPOINT:  return ROL::makePtr<CauchyPoint<Real>>(parlist);
      case TRUSTREGION_DOGLEG:       return ROL::makePtr<DogLeg<Real>>(parlist);
      case TRUSTREGION_DOUBLEDOGLEG: return ROL::makePtr<DoubleDogLeg<Real>>(parlist);
      case TRUSTREGION_TRUNCATEDCG:  return ROL::makePtr<TruncatedCG<Real>>(parlist);
      case TRUSTREGION_LINMORE:      return ROL::makePtr<LinMore<Real>>(parlist);
      default:                       return ROL::nullPtr;
    }
  }
}

#endif
