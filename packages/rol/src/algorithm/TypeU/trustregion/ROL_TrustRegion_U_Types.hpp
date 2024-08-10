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

#ifndef ROL_TRUSTREGION_U_TYPES_HPP
#define ROL_TRUSTREGION_U_TYPES_HPP

namespace ROL {

  /** \enum   ROL::ETrustRegionU
      \brief  Enumeration of trust-region solver types.

      \arg    TRUSTREGION_U_CAUCHYPOINT     describe
      \arg    TRUSTREGION_U_TRUNCATEDCG     describe
      \arg    TRUSTREGION_U_SPG             describe
      \arg    TRUSTREGION_U_DOGLEG          describe
      \arg    TRUSTREGION_U_DOUBLEDOGLEG    describe
   */
  enum ETrustRegionU{
    TRUSTREGION_U_CAUCHYPOINT = 0,
    TRUSTREGION_U_TRUNCATEDCG,
    TRUSTREGION_U_SPG,
    TRUSTREGION_U_DOGLEG,
    TRUSTREGION_U_DOUBLEDOGLEG,
    TRUSTREGION_U_LAST
  };

  inline std::string ETrustRegionUToString(ETrustRegionU tr) {
    std::string retString;
    switch(tr) {
      case TRUSTREGION_U_CAUCHYPOINT:   retString = "Cauchy Point";        break;
      case TRUSTREGION_U_TRUNCATEDCG:   retString = "Truncated CG";        break;
      case TRUSTREGION_U_SPG:           retString = "SPG";                 break;
      case TRUSTREGION_U_DOGLEG:        retString = "Dogleg";              break;
      case TRUSTREGION_U_DOUBLEDOGLEG:  retString = "Double Dogleg";       break;
      case TRUSTREGION_U_LAST:          retString = "Last Type (Dummy)";   break;
      default:                          retString = "INVALID ETrustRegionU";
    }
    return retString;
  }

  /** \brief  Verifies validity of a TrustRegionU enum.
    
      \param  tr  [in]  - enum of the TrustRegionU
      \return 1 if the argument is a valid TrustRegionU; 0 otherwise.
    */
  inline int isValidTrustRegionU(ETrustRegionU ls){
    return( (ls == TRUSTREGION_U_CAUCHYPOINT)  ||
            (ls == TRUSTREGION_U_TRUNCATEDCG)  ||
            (ls == TRUSTREGION_U_SPG)          ||
            (ls == TRUSTREGION_U_DOGLEG)       ||
            (ls == TRUSTREGION_U_DOUBLEDOGLEG)
          );
  }

  inline ETrustRegionU & operator++(ETrustRegionU &type) {
    return type = static_cast<ETrustRegionU>(type+1);
  }

  inline ETrustRegionU operator++(ETrustRegionU &type, int) {
    ETrustRegionU oldval = type;
    ++type;
    return oldval;
  }

  inline ETrustRegionU & operator--(ETrustRegionU &type) {
    return type = static_cast<ETrustRegionU>(type-1);
  }

  inline ETrustRegionU operator--(ETrustRegionU &type, int) {
    ETrustRegionU oldval = type;
    --type;
    return oldval;
  }

  inline ETrustRegionU StringToETrustRegionU(std::string s) {
    s = removeStringFormat(s);
    for ( ETrustRegionU tr = TRUSTREGION_U_CAUCHYPOINT; tr < TRUSTREGION_U_LAST; tr++ ) {
      if ( !s.compare(removeStringFormat(ETrustRegionUToString(tr))) ) {
        return tr;
      }
    }
    return TRUSTREGION_U_CAUCHYPOINT;
  }
} // ROL

#endif
