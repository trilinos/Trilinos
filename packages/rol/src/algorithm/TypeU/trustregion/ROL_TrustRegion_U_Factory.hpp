// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions of enums for trust region algorithms.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_TRUSTREGION_U_FACTORY_HPP
#define ROL_TRUSTREGION_U_FACTORY_HPP

#include "ROL_TrustRegion_U_Types.hpp"
#include "ROL_CauchyPoint_U.hpp"
#include "ROL_DogLeg_U.hpp"
#include "ROL_DoubleDogLeg_U.hpp"
#include "ROL_TruncatedCG_U.hpp"
#include "ROL_SPGTrustRegion_U.hpp"

namespace ROL {
  template<typename Real>
  inline Ptr<TrustRegion_U<Real>> TrustRegionUFactory(ParameterList &list) {
    ETrustRegionU etr = StringToETrustRegionU(
      list.sublist("Step").sublist("Trust Region").get("Subproblem Solver","Dogleg"));
    switch(etr) {
      case TRUSTREGION_U_CAUCHYPOINT:  return makePtr<CauchyPoint_U<Real>>();
      case TRUSTREGION_U_DOGLEG:       return makePtr<DogLeg_U<Real>>();
      case TRUSTREGION_U_DOUBLEDOGLEG: return makePtr<DoubleDogLeg_U<Real>>();
      case TRUSTREGION_U_TRUNCATEDCG:  return makePtr<TruncatedCG_U<Real>>(list);
      case TRUSTREGION_U_SPG:          return makePtr<SPGTrustRegion_U<Real>>(list);
      default:                         return nullPtr;
    }
  }
} // ROL

#endif
