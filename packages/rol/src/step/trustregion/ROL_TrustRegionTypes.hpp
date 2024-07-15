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

#ifndef ROL_TRUSTREGIONTYPES_HPP
#define ROL_TRUSTREGIONTYPES_HPP

#include "ROL_Types.hpp"

namespace ROL {

  /** \enum   ROL::ETrustRegion
      \brief  Enumeration of trust-region solver types.

      \arg    CAUCHYPOINT     describe
      \arg    TRUNCATEDCG     describe
      \arg    DOGLEG          describe
      \arg    DOUBLEDOGLEG    describe
   */
  enum ETrustRegion{
    TRUSTREGION_CAUCHYPOINT = 0,
    TRUSTREGION_TRUNCATEDCG,
    TRUSTREGION_DOGLEG,
    TRUSTREGION_DOUBLEDOGLEG,
    TRUSTREGION_LINMORE,
    TRUSTREGION_LAST
  };

  inline std::string ETrustRegionToString(ETrustRegion tr) {
    std::string retString;
    switch(tr) {
      case TRUSTREGION_CAUCHYPOINT:   retString = "Cauchy Point";        break;
      case TRUSTREGION_TRUNCATEDCG:   retString = "Truncated CG";        break;
      case TRUSTREGION_DOGLEG:        retString = "Dogleg";              break;
      case TRUSTREGION_DOUBLEDOGLEG:  retString = "Double Dogleg";       break;
      case TRUSTREGION_LINMORE:       retString = "Lin-More";            break;
      case TRUSTREGION_LAST:          retString = "Last Type (Dummy)";   break;
      default:                        retString = "INVALID ETrustRegion";
    }
    return retString;
  }

  /** \brief  Verifies validity of a TrustRegion enum.
    
      \param  tr  [in]  - enum of the TrustRegion
      \return 1 if the argument is a valid TrustRegion; 0 otherwise.
    */
  inline int isValidTrustRegion(ETrustRegion ls){
    return( (ls == TRUSTREGION_CAUCHYPOINT)  ||
            (ls == TRUSTREGION_TRUNCATEDCG)  ||
            (ls == TRUSTREGION_DOGLEG)       ||
            (ls == TRUSTREGION_DOUBLEDOGLEG) ||
            (ls == TRUSTREGION_LINMORE)
          );
  }

  inline ETrustRegion & operator++(ETrustRegion &type) {
    return type = static_cast<ETrustRegion>(type+1);
  }

  inline ETrustRegion operator++(ETrustRegion &type, int) {
    ETrustRegion oldval = type;
    ++type;
    return oldval;
  }

  inline ETrustRegion & operator--(ETrustRegion &type) {
    return type = static_cast<ETrustRegion>(type-1);
  }

  inline ETrustRegion operator--(ETrustRegion &type, int) {
    ETrustRegion oldval = type;
    --type;
    return oldval;
  }

  inline ETrustRegion StringToETrustRegion(std::string s) {
    s = removeStringFormat(s);
    for ( ETrustRegion tr = TRUSTREGION_CAUCHYPOINT; tr < TRUSTREGION_LAST; tr++ ) {
      if ( !s.compare(removeStringFormat(ETrustRegionToString(tr))) ) {
        return tr;
      }
    }
    return TRUSTREGION_CAUCHYPOINT;
  }

  /** \enum   ROL::ETrustRegionModel
      \brief  Enumeration of trust-region model types.

      \arg    COLEMANLI   describe
      \arg    KELLEYSACHS describe
   */
  enum ETrustRegionModel{
    TRUSTREGION_MODEL_COLEMANLI = 0,
    TRUSTREGION_MODEL_KELLEYSACHS,
    TRUSTREGION_MODEL_LINMORE,
    TRUSTREGION_MODEL_LAST
  };

  inline std::string ETrustRegionModelToString(ETrustRegionModel tr) {
    std::string retString;
    switch(tr) {
      case TRUSTREGION_MODEL_COLEMANLI:   retString = "Coleman-Li";        break;
      case TRUSTREGION_MODEL_KELLEYSACHS: retString = "Kelley-Sachs";      break;
      case TRUSTREGION_MODEL_LINMORE:     retString = "Lin-More";          break;
      case TRUSTREGION_MODEL_LAST:        retString = "Last Type (Dummy)"; break;
      default:                            retString = "INVALID ETrustRegionModel";
    }
    return retString;
  }

  /** \brief  Verifies validity of a TrustRegionModel enum.
    
      \param  tr  [in]  - enum of the TrustRegionModel
      \return 1 if the argument is a valid TrustRegionModel; 0 otherwise.
    */
  inline int isValidTrustRegionModel(ETrustRegionModel ls){
    return( (ls == TRUSTREGION_MODEL_COLEMANLI)   ||
            (ls == TRUSTREGION_MODEL_KELLEYSACHS) ||
            (ls == TRUSTREGION_MODEL_LINMORE)
          );
  }

  inline ETrustRegionModel & operator++(ETrustRegionModel &type) {
    return type = static_cast<ETrustRegionModel>(type+1);
  }

  inline ETrustRegionModel operator++(ETrustRegionModel &type, int) {
    ETrustRegionModel oldval = type;
    ++type;
    return oldval;
  }

  inline ETrustRegionModel & operator--(ETrustRegionModel &type) {
    return type = static_cast<ETrustRegionModel>(type-1);
  }

  inline ETrustRegionModel operator--(ETrustRegionModel &type, int) {
    ETrustRegionModel oldval = type;
    --type;
    return oldval;
  }

  inline ETrustRegionModel StringToETrustRegionModel(std::string s) {
    s = removeStringFormat(s);
    for ( ETrustRegionModel tr = TRUSTREGION_MODEL_COLEMANLI; tr < TRUSTREGION_MODEL_LAST; tr++ ) {
      if ( !s.compare(removeStringFormat(ETrustRegionModelToString(tr))) ) {
        return tr;
      }
    }
    return TRUSTREGION_MODEL_COLEMANLI;
  }

  inline bool isValidTrustRegionSubproblem(ETrustRegion etr, ETrustRegionModel etrm, bool isBnd) {
    if (etrm != TRUSTREGION_MODEL_LINMORE || !isBnd) {
      switch(etr) {
        case TRUSTREGION_CAUCHYPOINT:  return true;
        case TRUSTREGION_DOGLEG:       return true;
        case TRUSTREGION_DOUBLEDOGLEG: return true;
        case TRUSTREGION_TRUNCATEDCG:  return true;
        case TRUSTREGION_LINMORE:      return false;
        default:                       return false;
      }
    }
    else {
      switch(etr) {
        case TRUSTREGION_CAUCHYPOINT:  return false;
        case TRUSTREGION_DOGLEG:       return false;
        case TRUSTREGION_DOUBLEDOGLEG: return false;
        case TRUSTREGION_TRUNCATEDCG:  return false;
        case TRUSTREGION_LINMORE:      return true;
        default:                       return false;
      }
    }
  }

  /** \enum  ROL::ETrustRegionFlag 
      \brief Enumation of flags used by trust-region solvers.

      \arg TRUSTREGION_FLAG_SUCCESS        Actual and predicted reductions are positive 
      \arg TRUSTREGION_FLAG_POSPREDNEG     Reduction is positive, predicted negative (impossible)
      \arg TRUSTREGION_FLAG_NPOSPREDPOS    Reduction is nonpositive, predicted positive
      \arg TRUSTREGION_FLAG_NPOSPREDNEG    Reduction is nonpositive, predicted negative (impossible)
      \arg TRUSTREGION_FLAG_QMINSUFDEC     Insufficient decrease of the quadratic model (bound constraint only)
      \arg TRUSTREGION_FLAG_NAN            Actual and/or predicted reduction is NaN

  */
  enum ETrustRegionFlag {
    TRUSTREGION_FLAG_SUCCESS = 0,
    TRUSTREGION_FLAG_POSPREDNEG,
    TRUSTREGION_FLAG_NPOSPREDPOS,
    TRUSTREGION_FLAG_NPOSPREDNEG,
    TRUSTREGION_FLAG_QMINSUFDEC,
    TRUSTREGION_FLAG_NAN,
    TRUSTREGION_FLAG_UNDEFINED 
  };
 

  inline std::string ETrustRegionFlagToString(ETrustRegionFlag trf) {
    std::string retString;
    switch(trf) {
      case TRUSTREGION_FLAG_SUCCESS:  
        retString = "Both actual and predicted reductions are positive (success)";
        break;
      case TRUSTREGION_FLAG_POSPREDNEG: 
        retString = "Actual reduction is positive and predicted reduction is negative (impossible)";
        break;
      case TRUSTREGION_FLAG_NPOSPREDPOS: 
        retString = "Actual reduction is nonpositive and predicted reduction is positive";
        break;
      case TRUSTREGION_FLAG_NPOSPREDNEG:
        retString = "Actual reduction is nonpositive and predicted reduction is negative (impossible)";
        break;
      case TRUSTREGION_FLAG_QMINSUFDEC:
        retString = "Sufficient decrease of the quadratic model not met (bound constraints only)";
        break;
      case TRUSTREGION_FLAG_NAN:
        retString = "Actual and/or predicted reduction is a NaN";
        break;
      default:
        retString = "INVALID ETrustRegionFlag";       
    }
    return retString;
  }





}

#endif
