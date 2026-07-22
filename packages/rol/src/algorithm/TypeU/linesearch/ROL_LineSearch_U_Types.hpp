// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINESEARCH_U_TYPES_H
#define ROL_LINESEARCH_U_TYPES_H

#include "ROL_Types.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"

namespace ROL {

  /** \enum   ROL::EDescentU
      \brief  Enumeration of descent direction types.

      \arg    DESCENT_U_STEEPEST        describe
      \arg    DESCENT_U_NONLINEARCG     describe
      \arg    DESCENT_U_SECANT          describe
      \arg    DESCENT_U_NEWTON          describe 
      \arg    DESCENT_U_NEWTONKRYLOV    describe
      \arg    DESCENT_U_SECANTPRECOND   describe
   */
  enum EDescentU{
    DESCENT_U_STEEPEST = 0,
    DESCENT_U_NONLINEARCG,
    DESCENT_U_SECANT,
    DESCENT_U_NEWTON,
    DESCENT_U_NEWTONKRYLOV,
    DESCENT_U_USERDEFINED,
    DESCENT_U_LAST
  };

  inline std::string EDescentUToString(EDescentU tr) {
    std::string retString;
    switch(tr) {
      case DESCENT_U_STEEPEST:     retString = "Steepest Descent";    break;
      case DESCENT_U_NONLINEARCG:  retString = "Nonlinear CG";        break;
      case DESCENT_U_SECANT:       retString = "Quasi-Newton Method"; break;
      case DESCENT_U_NEWTON:       retString = "Newton's Method";     break;
      case DESCENT_U_NEWTONKRYLOV: retString = "Newton-Krylov";       break;
      case DESCENT_U_USERDEFINED:  retString = "User Defined";        break;
      case DESCENT_U_LAST:         retString = "Last Type (Dummy)";   break;
      default:                     retString = "INVALID EDescentU";
    }
    return retString;
  }

  /** \brief  Verifies validity of a DescentU enum.
    
      \param  tr  [in]  - enum of the DescentU
      \return 1 if the argument is a valid DescentU; 0 otherwise.
    */
  inline int isValidDescentU(EDescentU d){
    return( (d == DESCENT_U_STEEPEST)     ||
            (d == DESCENT_U_NONLINEARCG)  ||
            (d == DESCENT_U_SECANT)       ||
            (d == DESCENT_U_NEWTON)       ||
            (d == DESCENT_U_NEWTONKRYLOV) ||
            (d == DESCENT_U_USERDEFINED)
          );
  }

  inline EDescentU & operator++(EDescentU &type) {
    return type = static_cast<EDescentU>(type+1);
  }

  inline EDescentU operator++(EDescentU &type, int) {
    EDescentU oldval = type;
    ++type;
    return oldval;
  }

  inline EDescentU & operator--(EDescentU &type) {
    return type = static_cast<EDescentU>(type-1);
  }

  inline EDescentU operator--(EDescentU &type, int) {
    EDescentU oldval = type;
    --type;
    return oldval;
  }

  inline EDescentU StringToEDescentU(std::string s) {
    s = removeStringFormat(s);
    for ( EDescentU des = DESCENT_U_STEEPEST; des < DESCENT_U_LAST; des++ ) {
      if ( !s.compare(removeStringFormat(EDescentUToString(des))) ) {
        return des;
      }
    }
    return DESCENT_U_SECANT;
  }
  
  /** \enum   ROL::ELineSearchU
      \brief  Enumeration of line-search types.

      \arg    LINESEARCH_U_BACKTRACKING    describe
      \arg    LINESEARCH_U_BISECTION       describe
      \arg    LINESEARCH_U_GOLDENSECTION   describe
      \arg    LINESEARCH_U_CUBICINTERP     describe
      \arg    LINESEARCH_U_BRENTS          describe
      \arg    LINESEARCH_U_USERDEFINED     describe
   */
  enum ELineSearchU{
    LINESEARCH_U_ITERATIONSCALING = 0,
    LINESEARCH_U_PATHBASEDTARGETLEVEL,
    LINESEARCH_U_BACKTRACKING,
    LINESEARCH_U_BISECTION,
    LINESEARCH_U_GOLDENSECTION,
    LINESEARCH_U_CUBICINTERP,
    LINESEARCH_U_BRENTS,
    LINESEARCH_U_USERDEFINED,
    LINESEARCH_U_LAST
  };

  inline std::string ELineSearchUToString(ELineSearchU ls) {
    std::string retString;
    switch(ls) {
      case LINESEARCH_U_ITERATIONSCALING:     retString = "Iteration Scaling";       break;
      case LINESEARCH_U_PATHBASEDTARGETLEVEL: retString = "Path-Based Target Level"; break;
      case LINESEARCH_U_BACKTRACKING:         retString = "Backtracking";            break;
      case LINESEARCH_U_BISECTION:            retString = "Bisection";               break;
      case LINESEARCH_U_GOLDENSECTION:        retString = "Golden Section";          break;
      case LINESEARCH_U_CUBICINTERP:          retString = "Cubic Interpolation";     break;
      case LINESEARCH_U_BRENTS:               retString = "Brent's";                 break;
      case LINESEARCH_U_USERDEFINED:          retString = "User Defined";            break;
      case LINESEARCH_U_LAST:                 retString = "Last Type (Dummy)";       break;
      default:                                retString = "INVALID ELineSearchU";
    }
    return retString;
  }
  
  /** \brief  Verifies validity of a LineSearchU enum.
    
      \param  ls  [in]  - enum of the linesearch
      \return 1 if the argument is a valid linesearch; 0 otherwise.
    */
  inline int isValidLineSearchU(ELineSearchU ls){
    return( (ls == LINESEARCH_U_BACKTRACKING)         ||
            (ls == LINESEARCH_U_ITERATIONSCALING)     ||
            (ls == LINESEARCH_U_PATHBASEDTARGETLEVEL) ||
            (ls == LINESEARCH_U_BISECTION)            ||
            (ls == LINESEARCH_U_GOLDENSECTION)        ||
            (ls == LINESEARCH_U_CUBICINTERP)          ||
            (ls == LINESEARCH_U_BRENTS)               ||
            (ls == LINESEARCH_U_USERDEFINED)
          );
  }

  inline ELineSearchU & operator++(ELineSearchU &type) {
    return type = static_cast<ELineSearchU>(type+1);
  }

  inline ELineSearchU operator++(ELineSearchU &type, int) {
    ELineSearchU oldval = type;
    ++type;
    return oldval;
  }

  inline ELineSearchU & operator--(ELineSearchU &type) {
    return type = static_cast<ELineSearchU>(type-1);
  }

  inline ELineSearchU operator--(ELineSearchU &type, int) {
    ELineSearchU oldval = type;
    --type;
    return oldval;
  }

  inline ELineSearchU StringToELineSearchU(std::string s) {
    s = removeStringFormat(s);
    for ( ELineSearchU ls = LINESEARCH_U_ITERATIONSCALING; ls < LINESEARCH_U_LAST; ls++ ) {
      if ( !s.compare(removeStringFormat(ELineSearchUToString(ls))) ) {
        return ls;
      }
    }
    return LINESEARCH_U_ITERATIONSCALING;
  }

  /** \enum   ROL::ECurvatureConditionU
      \brief  Enumeration of line-search curvature conditions.

      \arg    CURVATURECONDITION_U_WOLFE           describe
      \arg    CURVATURECONDITION_U_STRONGWOLFE     describe
      \arg    CURVATURECONDITION_U_GOLDSTEIN       describe
   */
  enum ECurvatureConditionU{
    CURVATURECONDITION_U_WOLFE = 0,
    CURVATURECONDITION_U_STRONGWOLFE,
    CURVATURECONDITION_U_GENERALIZEDWOLFE,
    CURVATURECONDITION_U_APPROXIMATEWOLFE,
    CURVATURECONDITION_U_GOLDSTEIN,
    CURVATURECONDITION_U_NULL,
    CURVATURECONDITION_U_LAST
  };

  inline std::string ECurvatureConditionUToString(ECurvatureConditionU ls) {
    std::string retString;
    switch(ls) {
      case CURVATURECONDITION_U_WOLFE:            retString = "Wolfe Conditions";             break;
      case CURVATURECONDITION_U_STRONGWOLFE:      retString = "Strong Wolfe Conditions";      break;
      case CURVATURECONDITION_U_GENERALIZEDWOLFE: retString = "Generalized Wolfe Conditions"; break;
      case CURVATURECONDITION_U_APPROXIMATEWOLFE: retString = "Approximate Wolfe Conditions"; break;
      case CURVATURECONDITION_U_GOLDSTEIN:        retString = "Goldstein Conditions";         break;
      case CURVATURECONDITION_U_NULL:             retString = "Null Curvature Condition";     break;
      case CURVATURECONDITION_U_LAST:             retString = "Last Type (Dummy)";            break;
      default:                                    retString = "INVALID ECurvatureConditionU";
    }
    return retString;
  }
  
  /** \brief  Verifies validity of a CurvatureConditionU enum.
    
      \param  ls  [in]  - enum of the Curvature Conditions
      \return 1 if the argument is a valid curvature condition; 0 otherwise.
    */
  inline int isValidCurvatureConditionU(ECurvatureConditionU ls){
    return( (ls == CURVATURECONDITION_U_WOLFE)            ||
            (ls == CURVATURECONDITION_U_STRONGWOLFE)      ||
            (ls == CURVATURECONDITION_U_GENERALIZEDWOLFE) ||
            (ls == CURVATURECONDITION_U_APPROXIMATEWOLFE) ||
            (ls == CURVATURECONDITION_U_GOLDSTEIN)        ||
            (ls == CURVATURECONDITION_U_NULL)
          );
  }

  inline ECurvatureConditionU & operator++(ECurvatureConditionU &type) {
    return type = static_cast<ECurvatureConditionU>(type+1);
  }

  inline ECurvatureConditionU operator++(ECurvatureConditionU &type, int) {
    ECurvatureConditionU oldval = type;
    ++type;
    return oldval;
  }

  inline ECurvatureConditionU & operator--(ECurvatureConditionU &type) {
    return type = static_cast<ECurvatureConditionU>(type-1);
  }

  inline ECurvatureConditionU operator--(ECurvatureConditionU &type, int) {
    ECurvatureConditionU oldval = type;
    --type;
    return oldval;
  }

  inline ECurvatureConditionU StringToECurvatureConditionU(std::string s) {
    s = removeStringFormat(s);
    for ( ECurvatureConditionU cc = CURVATURECONDITION_U_WOLFE; cc < CURVATURECONDITION_U_LAST; cc++ ) {
      if ( !s.compare(removeStringFormat(ECurvatureConditionUToString(cc))) ) {
        return cc;
      }
    }
    return CURVATURECONDITION_U_WOLFE;
  }
} // namespace ROL

#endif
