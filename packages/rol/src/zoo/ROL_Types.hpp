// @HEADER
// ************************************************************************
//
// Questions? Contact Denis Ridzal (dridzal@sandia.gov), or
//                    Drew Kouri   (dpkouri@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions of custom data types in ROL.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_TYPES_HPP
#define ROL_TYPES_HPP

#ifdef  HAVE_ROL_DEBUG
#define ROL_VALIDATE( A )  A
#else
#define ROL_VALIDATE( A ) /* empty */
#endif

#include <Teuchos_ScalarTraits.hpp>

/** \def    ROL_NUM_CHECKDERIV_STEPS
    \brief  Number of steps for derivative checks.
 */
#define ROL_NUM_CHECKDERIV_STEPS 13

namespace ROL {
  
  /** \brief  Platform-dependent machine epsilon. 
   */
  static const double ROL_EPSILON   = std::abs(Teuchos::ScalarTraits<double>::eps());
    
  /** \brief  Tolerance for various equality tests.
   */
  static const double ROL_THRESHOLD = 10.0 * ROL_EPSILON;
  

  /** \enum   ROL::ELineSearch
      \brief  Enumeration of line-search types.

      \arg    BACKTRACKING    describe
      \arg    BISECTION       describe
      \arg    GOLDENSECTION   describe
      \arg    CUBICINTERP     describe
      \arg    BRENTS          describe
   */
  enum ELineSearch{
    LINESEARCH_BACKTRACKING = 0,
    LINESEARCH_BISECTION,
    LINESEARCH_GOLDENSECTION,
    LINESEARCH_CUBICINTERP,
    LINESEARCH_BRENTS,
    LINESEARCH_LAST
  };

  inline std::string ELineSearchToString(ELineSearch ls) {
    std::string retString;
    switch(ls) {
      case LINESEARCH_BACKTRACKING:   retString = "Backtracking";        break;
      case LINESEARCH_BISECTION:      retString = "Bisection";           break;
      case LINESEARCH_GOLDENSECTION:  retString = "Golden Section";      break;
      case LINESEARCH_CUBICINTERP:    retString = "Cubic Interpolation"; break;
      case LINESEARCH_BRENTS:         retString = "Brents";              break;
      case LINESEARCH_LAST:           retString = "Last Type (Dummy)";   break;
      default:                        retString = "INVALID ELineSearch";
    }
    return retString;
  }
  
  /** \brief  Verifies validity of a LineSearch enum.
    
      \param  ls  [in]  - enum of the linesearch
      \return 1 if the argument is a valid linesearch; 0 otherwise.
    */
  inline int isValidLineSearch(ELineSearch ls){
    return( (ls == LINESEARCH_BACKTRACKING)  ||
            (ls == LINESEARCH_BISECTION)     ||
            (ls == LINESEARCH_GOLDENSECTION) ||
            (ls == LINESEARCH_CUBICINTERP)   ||
            (ls == LINESEARCH_BRENTS)
          );
  }

  inline ELineSearch & operator++(ELineSearch &type) {
    return type = static_cast<ELineSearch>(type+1);
  }

  inline ELineSearch operator++(ELineSearch &type, int) {
    ELineSearch oldval = type;
    ++type;
    return oldval;
  }

  inline ELineSearch & operator--(ELineSearch &type) {
    return type = static_cast<ELineSearch>(type-1);
  }

  inline ELineSearch operator--(ELineSearch &type, int) {
    ELineSearch oldval = type;
    --type;
    return oldval;
  }

} // namespace ROL


/*! \mainpage ROL Documentation (Development Version)
 
  \image html rol.png
  \image latex rol.jpg "Rapid Optimization Library" width=1in

  \section intro_sec Introduction

  %ROL, the Rapid Optimization Library, is a Trilinos package for matrix-free
  optimization.
 
  \section overview_sec Overview

  Current release of %ROL includes the following features:
  \li Unconstrained.

  \section quickstart_sec Quick Start

  The following example demonstrates, in 2 steps, the use of %ROL.

  \subsection vector_qs_sec Step 1: Define linear algebra / vector interface.
  
  \code
      ROL::Vector
  \endcode

  We additionally set the number of computational cells \c numCells.


  \subsection objective_qs_sec Step 2: Define objective function.

  \code
      ROL::Objective
  \endcode


  \subsection done_qs_sec Done!


  \section devplans_sec Development Plans

  The next release of %ROL is expected to ...
*/
  

#endif
