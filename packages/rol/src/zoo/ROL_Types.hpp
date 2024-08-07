// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

#include <algorithm>
#include <complex>
#include <exception>
#include <string>
#include <sstream>
#include <limits>
#include <type_traits>
#include <ROL_stacktrace.hpp>
#include "ROL_ScalarTraits.hpp"
#include <ROL_Ptr.hpp>
#include <ROL_Vector.hpp>
#include <ROL_config.h>

/** \def    ROL_NUM_CHECKDERIV_STEPS
    \brief  Number of steps for derivative checks.
 */
#define ROL_NUM_CHECKDERIV_STEPS 13



namespace ROL {

  template<class T>
  std::string NumberToString( T Number )
  {
    std::ostringstream ss;
    ss << Number;
    return ss.str();
  }

  /** \brief  Platform-dependent machine epsilon.
   */
  template<class Real>
  inline Real ROL_EPSILON(void) { return std::abs(ROL::ScalarTraits<Real>::eps()); }
  //static const Real ROL_EPSILON<Real>() = std::abs(ROL::ScalarTraits<Real>::eps());

  /** \brief  Tolerance for various equality tests.
   */
  template<class Real>
  inline Real ROL_THRESHOLD(void) { return 10.0 * ROL_EPSILON<Real>(); }

  /** \brief  Platform-dependent maximum double.
   */
  template<class Real>
  inline Real ROL_OVERFLOW(void) { return std::abs(ROL::ScalarTraits<Real>::rmax()); }

  template<class Real>
  inline Real ROL_INF(void) { return 0.1*ROL_OVERFLOW<Real>(); }

  template<class Real>
  inline Real ROL_NINF(void) { return -ROL_INF<Real>(); }

  /** \brief  Platform-dependent minimum double.
   */
  template<class Real>
  inline Real ROL_UNDERFLOW(void) { return std::abs(ROL::ScalarTraits<Real>::rmin()); }

  /** \brief Enum for algorithm termination.
   */
  enum EExitStatus {
    EXITSTATUS_CONVERGED = 0,
    EXITSTATUS_MAXITER,
    EXITSTATUS_STEPTOL,
    EXITSTATUS_NAN,
    EXITSTATUS_USERDEFINED,
    EXITSTATUS_LAST
  };

  inline std::string EExitStatusToString(EExitStatus tr) {
    std::string retString;
    switch(tr) {
      case EXITSTATUS_CONVERGED:   retString = "Converged";                          break;
      case EXITSTATUS_MAXITER:     retString = "Iteration Limit Exceeded";           break;
      case EXITSTATUS_STEPTOL:     retString = "Step Tolerance Met";                 break;
      case EXITSTATUS_NAN:         retString = "Step and/or Gradient Returned NaN";  break;
      case EXITSTATUS_USERDEFINED: retString = "User Defined";                       break;
      case EXITSTATUS_LAST:        retString = "Last Type (Dummy)";                  break;
      default:                     retString = "INVALID EExitStatus";
    }
    return retString;
  }

  /** \brief  State for algorithm class.  Will be used for restarts.
   */
  template<class Real>
  struct AlgorithmState {
    int  iter;
    int  minIter;
    int  nfval;
    int  ncval;
    int  ngrad;
    Real value;
    Real minValue;
    Real gnorm;
    Real cnorm;
    Real snorm;
    Real aggregateGradientNorm;
    Real aggregateModelError;
    bool flag;
    ROL::Ptr<Vector<Real> > iterateVec;
    ROL::Ptr<Vector<Real> > lagmultVec;
    ROL::Ptr<Vector<Real> > minIterVec;
    EExitStatus statusFlag;

    AlgorithmState(void) : iter(0), minIter(0), nfval(0), ngrad(0), value(0), minValue(0), 
      gnorm(std::numeric_limits<Real>::max()),
      cnorm(std::numeric_limits<Real>::max()),
      snorm(std::numeric_limits<Real>::max()), 
      aggregateGradientNorm(std::numeric_limits<Real>::max()),
      aggregateModelError(std::numeric_limits<Real>::max()),
      flag(false),
      iterateVec(ROL::nullPtr), lagmultVec(ROL::nullPtr), minIterVec(ROL::nullPtr),
      statusFlag(EXITSTATUS_LAST) {}
    
    virtual ~AlgorithmState() {}

    void reset(void) {
      iter                  = 0;
      minIter               = 0;
      nfval                 = 0;
      ncval                 = 0;
      ngrad                 = 0;
      value                 = ROL_INF<Real>();
      minValue              = ROL_INF<Real>();
      gnorm                 = ROL_INF<Real>();
      cnorm                 = ROL_INF<Real>();
      snorm                 = ROL_INF<Real>();
      aggregateGradientNorm = ROL_INF<Real>();
      aggregateModelError   = ROL_INF<Real>();
      flag                  = false;
      if (iterateVec != ROL::nullPtr) {
        iterateVec->zero();
      }
      if (lagmultVec != ROL::nullPtr) {
        lagmultVec->zero();
      }
      if (minIterVec != ROL::nullPtr) {
        minIterVec->zero();
      }
    }
  };

  /** \brief  State for step class.  Will be used for restarts.
   */
  template<class Real>
  struct StepState {
    ROL::Ptr<Vector<Real> > gradientVec;
    ROL::Ptr<Vector<Real> > descentVec;
    ROL::Ptr<Vector<Real> > constraintVec;
    int nfval;
    int ngrad;
    Real searchSize; // line search parameter (alpha) or trust-region radius (delta)
    int flag; // Was step successful?
    int SPiter; // Subproblem iteration count
    int SPflag; // Subproblem termination flag

    StepState(void) : gradientVec(ROL::nullPtr),
                      descentVec(ROL::nullPtr),
                      constraintVec(ROL::nullPtr),
                      nfval(0),
                      ngrad(0),
                      searchSize(0),
                      flag(0),
                      SPiter(0),
                      SPflag(0) {}

    void reset(const Real searchSizeInput = 1.0) {
      if (gradientVec != ROL::nullPtr) {
        gradientVec->zero();
      }
      if (descentVec != ROL::nullPtr) {
        descentVec->zero();
      }
      if (constraintVec != ROL::nullPtr) {
        constraintVec->zero();
      }
      nfval = 0;
      ngrad = 0;
      searchSize = searchSizeInput;
      flag = 0;
      SPiter = 0;
      SPflag = 0;
    }
  };

  struct removeSpecialCharacters {
    bool operator()(char c) {
      return (c ==' ' || c =='-' || c == '(' || c == ')' || c=='\'' || c=='\r' || c=='\n' || c=='\t');
    }
  };

  inline std::string removeStringFormat( std::string s ) {
    std::string output = s;
    output.erase( std::remove_if( output.begin(), output.end(), removeSpecialCharacters()), output.end() );
    std::transform( output.begin(), output.end(), output.begin(), ::tolower );
    return output;
  }

  // Types of optimization problem
  enum EProblem {
    TYPE_U = 0,
    TYPE_B,
    TYPE_E,
    TYPE_EB,
    TYPE_LAST
  };

  /** \enum   ROL::EStep
      \brief  Enumeration of step types.

      \arg    AUGMENTEDLAGRANGIAN     describe
      \arg    BUNDLE                  describe
      \arg    COMPOSITESTEP           describe
      \arg    LINESEARCH              describe
      \arg    MOREAUYOSIDAPENALTY     describe
      \arg    PRIMALDUALACTIVESET     describe
      \arg    TRUSTREGION             describe
   */
  enum EStep{
    STEP_AUGMENTEDLAGRANGIAN = 0,
    STEP_BUNDLE,
    STEP_COMPOSITESTEP,
    STEP_LINESEARCH,
    STEP_MOREAUYOSIDAPENALTY,
    STEP_PRIMALDUALACTIVESET,
    STEP_TRUSTREGION,
    STEP_INTERIORPOINT,
    STEP_FLETCHER,
    STEP_LAST
  };

  inline std::string EStepToString(EStep tr) {
    std::string retString;
    switch(tr) {
      case STEP_AUGMENTEDLAGRANGIAN: retString = "Augmented Lagrangian";   break;
      case STEP_BUNDLE:              retString = "Bundle";                 break;
      case STEP_COMPOSITESTEP:       retString = "Composite Step";         break;
      case STEP_LINESEARCH:          retString = "Line Search";            break;
      case STEP_MOREAUYOSIDAPENALTY: retString = "Moreau-Yosida Penalty";  break;
      case STEP_PRIMALDUALACTIVESET: retString = "Primal Dual Active Set"; break;
      case STEP_TRUSTREGION:         retString = "Trust Region";           break;
      case STEP_INTERIORPOINT:       retString = "Interior Point";         break;
      case STEP_FLETCHER:            retString = "Fletcher";               break;
      case STEP_LAST:                retString = "Last Type (Dummy)";      break;
      default:                       retString = "INVALID EStep";
    }
    return retString;
  }

  inline bool isCompatibleStep( EProblem p, EStep s ) {
    bool comp = false;
    switch(p) {

      case TYPE_U:    comp = ( (s == STEP_LINESEARCH) ||
                               (s == STEP_TRUSTREGION) ||
                               (s == STEP_BUNDLE) );
        break;

      case TYPE_B:    comp = ( (s == STEP_LINESEARCH)  ||
                               (s == STEP_TRUSTREGION) ||
                               (s == STEP_MOREAUYOSIDAPENALTY) ||
                               (s == STEP_PRIMALDUALACTIVESET) ||
                               (s == STEP_INTERIORPOINT) );
        break;

      case TYPE_E:    comp = ( (s == STEP_COMPOSITESTEP) ||
                               (s == STEP_AUGMENTEDLAGRANGIAN) ||
			       (s == STEP_FLETCHER) );
        break;

      case TYPE_EB:   comp = ( (s == STEP_AUGMENTEDLAGRANGIAN) ||
                               (s == STEP_MOREAUYOSIDAPENALTY) ||
                               (s == STEP_INTERIORPOINT) ||
			       (s == STEP_FLETCHER) );
        break;

      case TYPE_LAST: comp = false; break;
      default:        comp = false;
    }
    return comp;
  }

  inline std::string EProblemToString( EProblem p ) {
    std::string retString;
    switch(p) {
      case TYPE_U:     retString = "Type-U";             break;
      case TYPE_E:     retString = "Type-E";             break;
      case TYPE_B:     retString = "Type-B";             break;
      case TYPE_EB:    retString = "Type-EB";            break;
      case TYPE_LAST:  retString = "Type-Last (Dummy)";  break;
      default:         retString = "Invalid EProblem";
    }
    return retString;
  }
  
 
  /** \brief  Verifies validity of a TrustRegion enum.
    
      \param  tr  [in]  - enum of the TrustRegion
      \return 1 if the argument is a valid TrustRegion; 0 otherwise.
    */
  inline int isValidStep(EStep ls) {
    return( (ls == STEP_AUGMENTEDLAGRANGIAN) ||
            (ls == STEP_BUNDLE) ||
            (ls == STEP_COMPOSITESTEP) ||
            (ls == STEP_LINESEARCH) ||
            (ls == STEP_MOREAUYOSIDAPENALTY) ||
            (ls == STEP_PRIMALDUALACTIVESET) ||
            (ls == STEP_TRUSTREGION) || 
            (ls == STEP_INTERIORPOINT) ||
	    (ls == STEP_FLETCHER) ) ;
  }

  inline EStep & operator++(EStep &type) {
    return type = static_cast<EStep>(type+1);
  }

  inline EStep operator++(EStep &type, int) {
    EStep oldval = type;
    ++type;
    return oldval;
  }

  inline EStep & operator--(EStep &type) {
    return type = static_cast<EStep>(type-1);
  }

  inline EStep operator--(EStep &type, int) {
    EStep oldval = type;
    --type;
    return oldval;
  }

  inline EStep StringToEStep(std::string s) {
    s = removeStringFormat(s);
    for ( EStep st = STEP_AUGMENTEDLAGRANGIAN; st < STEP_LAST; ++st ) {
      if ( !s.compare(removeStringFormat(EStepToString(st))) ) {
        return st;
      }
    }
    return STEP_LAST;
  }

  /** \enum   ROL::EDescent
      \brief  Enumeration of descent direction types.

      \arg    STEEPEST        describe
      \arg    NONLINEARCG     describe
      \arg    SECANT          describe
      \arg    NEWTON          describe 
      \arg    NEWTONKRYLOV    describe
      \arg    SECANTPRECOND   describe
   */
  enum EDescent{
    DESCENT_STEEPEST = 0,
    DESCENT_NONLINEARCG,
    DESCENT_SECANT,
    DESCENT_NEWTON,
    DESCENT_NEWTONKRYLOV,
    DESCENT_LAST
  };

  inline std::string EDescentToString(EDescent tr) {
    std::string retString;
    switch(tr) {
      case DESCENT_STEEPEST:             retString = "Steepest Descent";                          break;
      case DESCENT_NONLINEARCG:          retString = "Nonlinear CG";                              break;
      case DESCENT_SECANT:               retString = "Quasi-Newton Method";                       break;
      case DESCENT_NEWTON:               retString = "Newton's Method";                           break;
      case DESCENT_NEWTONKRYLOV:         retString = "Newton-Krylov";                             break;
      case DESCENT_LAST:                 retString = "Last Type (Dummy)";                         break;
      default:                           retString = "INVALID ESecant";
    }
    return retString;
  }

  /** \brief  Verifies validity of a Secant enum.
    
      \param  tr  [in]  - enum of the Secant
      \return 1 if the argument is a valid Secant; 0 otherwise.
    */
  inline int isValidDescent(EDescent d){
    return( (d == DESCENT_STEEPEST)      ||
            (d == DESCENT_NONLINEARCG)   ||
            (d == DESCENT_SECANT)        ||
            (d == DESCENT_NEWTON)        ||
            (d == DESCENT_NEWTONKRYLOV)
          );
  }

  inline EDescent & operator++(EDescent &type) {
    return type = static_cast<EDescent>(type+1);
  }

  inline EDescent operator++(EDescent &type, int) {
    EDescent oldval = type;
    ++type;
    return oldval;
  }

  inline EDescent & operator--(EDescent &type) {
    return type = static_cast<EDescent>(type-1);
  }

  inline EDescent operator--(EDescent &type, int) {
    EDescent oldval = type;
    --type;
    return oldval;
  }

  inline EDescent StringToEDescent(std::string s) {
    s = removeStringFormat(s);
    for ( EDescent des = DESCENT_STEEPEST; des < DESCENT_LAST; des++ ) {
      if ( !s.compare(removeStringFormat(EDescentToString(des))) ) {
        return des;
      }
    }
    return DESCENT_SECANT;
  }
  
  /** \enum   ROL::ESecant
      \brief  Enumeration of secant update algorithms.

      \arg    LBFGS           describe
      \arg    LDFP            describe
      \arg    LSR1            describe 
      \arg    BARZILAIBORWEIN describe
   */
  enum ESecant{
    SECANT_LBFGS = 0,
    SECANT_LDFP,
    SECANT_LSR1,
    SECANT_BARZILAIBORWEIN,
    SECANT_USERDEFINED,
    SECANT_LAST
  };

  inline std::string ESecantToString(ESecant tr) {
    std::string retString;
    switch(tr) {
      case SECANT_LBFGS:           retString = "Limited-Memory BFGS"; break;
      case SECANT_LDFP:            retString = "Limited-Memory DFP";  break;
      case SECANT_LSR1:            retString = "Limited-Memory SR1";  break;
      case SECANT_BARZILAIBORWEIN: retString = "Barzilai-Borwein";    break;
      case SECANT_USERDEFINED:     retString = "User-Defined";        break;
      case SECANT_LAST:            retString = "Last Type (Dummy)";   break;
      default:                     retString = "INVALID ESecant";
    }
    return retString;
  }
  
  /** \brief  Verifies validity of a Secant enum.
    
      \param  tr  [in]  - enum of the Secant
      \return 1 if the argument is a valid Secant; 0 otherwise.
    */
  inline int isValidSecant(ESecant s) {
    return( (s == SECANT_LBFGS)           ||
            (s == SECANT_LDFP)            ||
            (s == SECANT_LSR1)            ||
            (s == SECANT_BARZILAIBORWEIN) ||
            (s == SECANT_USERDEFINED)
          );
  }

  inline ESecant & operator++(ESecant &type) {
    return type = static_cast<ESecant>(type+1);
  }

  inline ESecant operator++(ESecant &type, int) {
    ESecant oldval = type;
    ++type;
    return oldval;
  }

  inline ESecant & operator--(ESecant &type) {
    return type = static_cast<ESecant>(type-1);
  }

  inline ESecant operator--(ESecant &type, int) {
    ESecant oldval = type;
    --type;
    return oldval;
  }

  inline ESecant StringToESecant(std::string s) {
    s = removeStringFormat(s);
    for ( ESecant sec = SECANT_LBFGS; sec < SECANT_LAST; sec++ ) {
      if ( !s.compare(removeStringFormat(ESecantToString(sec))) ) {
        return sec;
      }
    }
    return SECANT_LBFGS;
  }
  
  /** \enum   ROL::ENonlinearCG
      \brief  Enumeration of nonlinear CG algorithms.

      \arg    HESTENES_STIEFEL   \f$ \frac{g_{k+1}^\top y_k}{d_k^\top y_k } \f$
      \arg    FLETCHER_REEVES    \f$ \frac{\|g_{k+1}\|^2}{\|g_k\|^2} \f$
      \arg    DANIEL             \f$ \frac{g_{k+1}^\top \nabla^2 f(x_k) d_k}{d_k^\top \nabla^2 f(x_k) d_k} \f$
      \arg    POLAK_RIBIERE      \f$ \frac{g_{k+1}^\top y_k}{\|g_k\|^2} \f$
      \arg    FLETCHER_CONJDESC  \f$ -\frac{\|g_{k+1}\|^2}{d_k^\top g_k} \f$
      \arg    LIU_STOREY         \f$ -\frac{g_k^\top y_{k-1} }{d_{k-1}^\top g_{k-1} \f$
      \arg    DAI_YUAN           \f$ \frac{\|g_{k+1}\|^2}{d_k^\top y_k} \f$
      \arg    HAGER_ZHANG        \f$ \frac{g_{k+1}^\top y_k}{d_k^\top y_k} - 2 \frac{\|y_k\|^2}{d_k^\top y_k} \frac{g_{k+1}^\top d_k}{d_k^\top y_k} \f$
      \arg    OREN_LUENBERGER    \f$ \frac{g_{k+1}^\top y_k}{d_k^\top y_k} - \frac{\|y_k\|^2}{d_k^\top y_k} \frac{g_{k+1}^\top d_k}{d_k^\top y_k} \f$ 
   */
  enum ENonlinearCG{
    NONLINEARCG_HESTENES_STIEFEL = 0,
    NONLINEARCG_FLETCHER_REEVES,
    NONLINEARCG_DANIEL,
    NONLINEARCG_POLAK_RIBIERE,
    NONLINEARCG_FLETCHER_CONJDESC,
    NONLINEARCG_LIU_STOREY,
    NONLINEARCG_DAI_YUAN,
    NONLINEARCG_HAGER_ZHANG,
    NONLINEARCG_OREN_LUENBERGER,
    NONLINEARCG_USERDEFINED,
    NONLINEARCG_LAST
  };

  inline std::string ENonlinearCGToString(ENonlinearCG tr) {
    std::string retString;
    switch(tr) {
      case NONLINEARCG_HESTENES_STIEFEL:      retString = "Hestenes-Stiefel";            break;
      case NONLINEARCG_FLETCHER_REEVES:       retString = "Fletcher-Reeves";             break;
      case NONLINEARCG_DANIEL:                retString = "Daniel (uses Hessian)";       break;
      case NONLINEARCG_POLAK_RIBIERE:         retString = "Polak-Ribiere";               break;
      case NONLINEARCG_FLETCHER_CONJDESC:     retString = "Fletcher Conjugate Descent";  break;
      case NONLINEARCG_LIU_STOREY:            retString = "Liu-Storey";                  break;
      case NONLINEARCG_DAI_YUAN:              retString = "Dai-Yuan";                    break;
      case NONLINEARCG_HAGER_ZHANG:           retString = "Hager-Zhang";                 break;
      case NONLINEARCG_OREN_LUENBERGER:       retString = "Oren-Luenberger";             break;
      case NONLINEARCG_USERDEFINED:           retString = "User Defined";                break;
      case NONLINEARCG_LAST:                  retString = "Last Type (Dummy)";           break;
      default:                                retString = "INVALID ENonlinearCG";
    }
    return retString;
  }
  
  /** \brief  Verifies validity of a NonlinearCG enum.
    
      \param  tr  [in]  - enum of the NonlinearCG
      \return 1 if the argument is a valid NonlinearCG; 0 otherwise.
    */
  inline int isValidNonlinearCG(ENonlinearCG s) {
    return( (s == NONLINEARCG_HESTENES_STIEFEL)  ||
            (s == NONLINEARCG_FLETCHER_REEVES)   ||
            (s == NONLINEARCG_DANIEL)            ||
            (s == NONLINEARCG_POLAK_RIBIERE)     ||
            (s == NONLINEARCG_FLETCHER_CONJDESC) ||
            (s == NONLINEARCG_LIU_STOREY)        ||
            (s == NONLINEARCG_DAI_YUAN)          ||
            (s == NONLINEARCG_HAGER_ZHANG)       ||
            (s == NONLINEARCG_OREN_LUENBERGER)   ||
            (s == NONLINEARCG_USERDEFINED)
          );
  }

  inline ENonlinearCG & operator++(ENonlinearCG &type) {
    return type = static_cast<ENonlinearCG>(type+1);
  }

  inline ENonlinearCG operator++(ENonlinearCG &type, int) {
    ENonlinearCG oldval = type;
    ++type;
    return oldval;
  }

  inline ENonlinearCG & operator--(ENonlinearCG &type) {
    return type = static_cast<ENonlinearCG>(type-1);
  }

  inline ENonlinearCG operator--(ENonlinearCG &type, int) {
    ENonlinearCG oldval = type;
    --type;
    return oldval;
  }

  inline ENonlinearCG StringToENonlinearCG(std::string s) {
    s = removeStringFormat(s);
    for ( ENonlinearCG nlcg = NONLINEARCG_HESTENES_STIEFEL; nlcg < NONLINEARCG_LAST; nlcg++ ) {
      if ( !s.compare(removeStringFormat(ENonlinearCGToString(nlcg))) ) {
        return nlcg;
      }
    }
    return NONLINEARCG_HESTENES_STIEFEL;
  }
  
  /** \enum   ROL::ELineSearch
      \brief  Enumeration of line-search types.

      \arg    BACKTRACKING    describe
      \arg    BISECTION       describe
      \arg    GOLDENSECTION   describe
      \arg    CUBICINTERP     describe
      \arg    BRENTS          describe
      \arg    USERDEFINED     describe
   */
  enum ELineSearch{
    LINESEARCH_ITERATIONSCALING = 0,
    LINESEARCH_PATHBASEDTARGETLEVEL,
    LINESEARCH_BACKTRACKING,
    LINESEARCH_BISECTION,
    LINESEARCH_GOLDENSECTION,
    LINESEARCH_CUBICINTERP,
    LINESEARCH_BRENTS,
    LINESEARCH_USERDEFINED,
    LINESEARCH_LAST
  };

  inline std::string ELineSearchToString(ELineSearch ls) {
    std::string retString;
    switch(ls) {
      case LINESEARCH_ITERATIONSCALING:     retString = "Iteration Scaling";       break;
      case LINESEARCH_PATHBASEDTARGETLEVEL: retString = "Path-Based Target Level"; break;
      case LINESEARCH_BACKTRACKING:         retString = "Backtracking";            break;
      case LINESEARCH_BISECTION:            retString = "Bisection";               break;
      case LINESEARCH_GOLDENSECTION:        retString = "Golden Section";          break;
      case LINESEARCH_CUBICINTERP:          retString = "Cubic Interpolation";     break;
      case LINESEARCH_BRENTS:               retString = "Brent's";                 break;
      case LINESEARCH_USERDEFINED:          retString = "User Defined";            break;
      case LINESEARCH_LAST:                 retString = "Last Type (Dummy)";       break;
      default:                              retString = "INVALID ELineSearch";
    }
    return retString;
  }
  
  /** \brief  Verifies validity of a LineSearch enum.
    
      \param  ls  [in]  - enum of the linesearch
      \return 1 if the argument is a valid linesearch; 0 otherwise.
    */
  inline int isValidLineSearch(ELineSearch ls){
    return( (ls == LINESEARCH_BACKTRACKING)         ||
            (ls == LINESEARCH_ITERATIONSCALING)     ||
            (ls == LINESEARCH_PATHBASEDTARGETLEVEL) ||
            (ls == LINESEARCH_BISECTION)            ||
            (ls == LINESEARCH_GOLDENSECTION)        ||
            (ls == LINESEARCH_CUBICINTERP)          ||
            (ls == LINESEARCH_BRENTS)               ||
            (ls == LINESEARCH_USERDEFINED)
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

  inline ELineSearch StringToELineSearch(std::string s) {
    s = removeStringFormat(s);
    for ( ELineSearch ls = LINESEARCH_ITERATIONSCALING; ls < LINESEARCH_LAST; ls++ ) {
      if ( !s.compare(removeStringFormat(ELineSearchToString(ls))) ) {
        return ls;
      }
    }
    return LINESEARCH_ITERATIONSCALING;
  }

  /** \enum   ROL::ECurvatureCondition
      \brief  Enumeration of line-search curvature conditions.

      \arg    WOLFE           describe
      \arg    STRONGWOLFE     describe
      \arg    GOLDSTEIN       describe
   */
  enum ECurvatureCondition{
    CURVATURECONDITION_WOLFE = 0,
    CURVATURECONDITION_STRONGWOLFE,
    CURVATURECONDITION_GENERALIZEDWOLFE,
    CURVATURECONDITION_APPROXIMATEWOLFE,
    CURVATURECONDITION_GOLDSTEIN,
    CURVATURECONDITION_NULL,
    CURVATURECONDITION_LAST
  };

  inline std::string ECurvatureConditionToString(ECurvatureCondition ls) {
    std::string retString;
    switch(ls) {
      case CURVATURECONDITION_WOLFE:            retString = "Wolfe Conditions";             break;
      case CURVATURECONDITION_STRONGWOLFE:      retString = "Strong Wolfe Conditions";      break;
      case CURVATURECONDITION_GENERALIZEDWOLFE: retString = "Generalized Wolfe Conditions"; break;
      case CURVATURECONDITION_APPROXIMATEWOLFE: retString = "Approximate Wolfe Conditions"; break;
      case CURVATURECONDITION_GOLDSTEIN:        retString = "Goldstein Conditions";         break;
      case CURVATURECONDITION_NULL:             retString = "Null Curvature Condition";     break;
      case CURVATURECONDITION_LAST:             retString = "Last Type (Dummy)";            break;
      default:                                  retString = "INVALID ECurvatureCondition";
    }
    return retString;
  }
  
  /** \brief  Verifies validity of a CurvatureCondition enum.
    
      \param  ls  [in]  - enum of the Curvature Conditions
      \return 1 if the argument is a valid curvature condition; 0 otherwise.
    */
  inline int isValidCurvatureCondition(ECurvatureCondition ls){
    return( (ls == CURVATURECONDITION_WOLFE)            ||
            (ls == CURVATURECONDITION_STRONGWOLFE)      ||
            (ls == CURVATURECONDITION_GENERALIZEDWOLFE) ||
            (ls == CURVATURECONDITION_APPROXIMATEWOLFE) ||
            (ls == CURVATURECONDITION_GOLDSTEIN)        ||
            (ls == CURVATURECONDITION_NULL)
          );
  }

  inline ECurvatureCondition & operator++(ECurvatureCondition &type) {
    return type = static_cast<ECurvatureCondition>(type+1);
  }

  inline ECurvatureCondition operator++(ECurvatureCondition &type, int) {
    ECurvatureCondition oldval = type;
    ++type;
    return oldval;
  }

  inline ECurvatureCondition & operator--(ECurvatureCondition &type) {
    return type = static_cast<ECurvatureCondition>(type-1);
  }

  inline ECurvatureCondition operator--(ECurvatureCondition &type, int) {
    ECurvatureCondition oldval = type;
    --type;
    return oldval;
  }

  inline ECurvatureCondition StringToECurvatureCondition(std::string s) {
    s = removeStringFormat(s);
    for ( ECurvatureCondition cc = CURVATURECONDITION_WOLFE; cc < CURVATURECONDITION_LAST; cc++ ) {
      if ( !s.compare(removeStringFormat(ECurvatureConditionToString(cc))) ) {
        return cc;
      }
    }
    return CURVATURECONDITION_WOLFE;
  }

  /** \enum  ROL::ECGFlag 
      \brief Enumation of flags used by conjugate gradient methods.

    \arg CG_FLAG_SUCCESS     Residual Tolerance Met
    \arg CG_FLAG_ITEREXCEED  Iteration Limit Exceeded
    \arg CG_FLAG_NEGCURVE    Negative Curvature Detected
    \arh CG_FLAG_TRRADEX     Trust-Region Radius Exceeded
    \arh CG_FLAG_ZERORHS     Initiali Right Hand Side is Zero

  */
  enum ECGFlag {
    CG_FLAG_SUCCESS = 0,
    CG_FLAG_ITEREXCEED,
    CG_FLAG_NEGCURVE,
    CG_FLAG_TRRADEX,
    CG_FLAG_ZERORHS,
    CG_FLAG_UNDEFINED 
  };


  inline std::string ECGFlagToString(ECGFlag cgf) {
    std::string retString;
    switch(cgf) {
      case CG_FLAG_SUCCESS:
        retString = "Residual tolerance met";
        break;
      case CG_FLAG_ITEREXCEED:
        retString = "Iteration limit exceeded";
        break;
      case CG_FLAG_NEGCURVE:
        retString = "Negative curvature detected";
        break;
      case CG_FLAG_TRRADEX:   
        retString = "Trust-Region radius exceeded";
        break;
      case CG_FLAG_ZERORHS:
        retString = "Initial right hand side is zero";
        break;
      default:
        retString = "INVALID ECGFlag";  
    }
    return retString;
  }



  // For use in gradient and Hessian checks
  namespace Finite_Difference_Arrays {

    // Finite difference steps in axpy form    
    const int shifts[4][4] = { {  1,  0,  0, 0 },  // First order
                               { -1,  2,  0, 0 },  // Second order
                               { -1,  2,  1, 0 },  // Third order
                               { -1, -1,  3, 1 }   // Fourth order
                             };

      // Finite difference weights     
     const double weights[4][5] = { { -1.0,          1.0, 0.0,      0.0,      0.0      },  // First order
                                    {  0.0,     -1.0/2.0, 1.0/2.0,  0.0,      0.0      },  // Second order
                                    { -1.0/2.0, -1.0/3.0, 1.0,     -1.0/6.0,  0.0      },  // Third order
                                    {  0.0,     -2.0/3.0, 1.0/12.0, 2.0/3.0, -1.0/12.0 }   // Fourth order
                                  };

  }


// Generic conversion from Element type to Real type
template<class Real, class Element>
struct TypeCaster {
  static Real ElementToReal( const Element &val ) {
    return Real(0);
  }
};

// Partially specialize for complex<Real>
template<class Real>
struct TypeCaster<Real, std::complex<Real> > {
  static Real ElementToReal( const std::complex<Real> &val ) {
    return val.real();
  } 
};

// Fully specialize for double,float
template<>
struct TypeCaster<double,float> {
  static double ElementToReal( const float &val ) {
    return static_cast<double>(val);
  }
};

// Cast from Element type to Real type
template<class Element, class Real>
Real rol_cast(const Element &val) {
  return TypeCaster<Real,Element>::ElementToReal(val);
}






namespace Exception {

class NotImplemented : public std::logic_error {
public:
  NotImplemented( const std::string& what_arg ) :
    std::logic_error(what_arg) {}


}; // class NotImplemented
 

#if __cplusplus >= 201402L // using C++14

using std::enable_if_t;

#else // No C++14

template<bool B, class T=void>
using enable_if_t = typename std::enable_if<B,T>::type;

#endif





} // namespace Exception


} // namespace ROL


/*! \mainpage %ROL Documentation (Development Version)
 
  \image html rol.png "Rapid Optimization Library" width=1in
  \image latex rol.pdf "Rapid Optimization Library" width=1in

  \section intro_sec Introduction

  Rapid Optimization Library (%ROL) is a C++ package for large-scale
  optimization.
  It is used for the solution of optimal design, optimal control and
  inverse problems in large-scale engineering applications.
  Other uses include mesh optimization and image processing. 

 
  \section overview_sec Overview

  %ROL aims to combine flexibility, efficiency and robustness.  Key features:

  \li  Matrix-free application programming interfaces (APIs) ---enable direct
       use of application data structures and memory spaces, linear solvers,
       nonlinear solvers and preconditioners.
  \li  State-of-the-art algorithms for unconstrained optimization,
       constrained optimization and optimization under uncertainty ---enable
       inexact and adaptive function evaluations and iterative linear
       system solves.
  \li  Special APIs for simulation-based optimization ---enable a
       streamlined embedding into engineering applications, rigorous
       implementation verification and efficient use.
  \li  Modular interfaces throughout the optimization process ---enable custom
       and user-defined algorithms, stopping criteria, hierarchies of
       algorithms, and selective use of a variety of tools and components.

  For a detailed description of user interfaces and algorithms, see the
  presentations ROL-Trilinos-xx.x.pptx (or .pdf) in the doc/presentations
  directory.

  To start using %ROL, including all its advanced algorithms and features,
  jump to the <a href="modules.html">Modules</a> page.

  For a basic example, see below.

  \section quickstart_sec Quick Start

  The Rosenbrock example (rol/example/rosenbrock/example_01.cpp) demonstrates
  the basic use of %ROL.
  It amounts to six steps:

  \subsection vector_qs_sec Step 1: Implement linear algebra / vector interface.
  --- or try one of the provided implementations, such as ROL::StdVector in rol/vector.
  
  ~~~{.hpp}
      ROL::Vector
  ~~~

  \subsection objective_qs_sec Step 2: Implement objective function interface.
  --- or try one of the provided functions, such as @b ROL::ZOO::Objective_Rosenbrock in rol/zoo.

  \code
      ROL::Objective
  \endcode

  \subsection step_qs_sec Step 3: Choose optimization algorithm.
  ---  with @b ROL::ParameterList settings in the variable @b parlist.

  \code
      ROL::Algorithm<RealT> algo("Line Search",parlist);
  \endcode

  \subsection run_qs_sec Step 4: Run algorithm.
  ---  starting from the initial iterate @b x, applied to objective function @b obj.

  \code
      algo.run(x, obj);
  \endcode
*/

/** @defgroup interface_group User Interface
 *  \brief ROL's interface, to be implemented by the user.
 */

/** @defgroup la_group Linear Algebra Interface
 *  @ingroup interface_group
 *  \brief ROL's linear algebra or vector space interface.
 */

/** @defgroup func_group Functional Interface
 *  @ingroup interface_group
 *  \brief ROL's functional interface.

    ROL is used for the numerical solution of smooth optimization problems
    \f[
      \begin{array}{rl}
        \displaystyle \min_{x} & f(x) \\
        \mbox{subject to} & c(x) = 0 \,, \\
                          & a \le x \le b \,,
      \end{array}
    \f]
    where:
      \li \f$f : \mathcal{X} \rightarrow \mathbb{R}\f$ is a Fr&eacute;chet differentiable functional,
      \li \f$c : \mathcal{X} \rightarrow \mathcal{C}\f$ is a Fr&eacute;chet differentiable operator,
      \li \f$\mathcal{X}\f$ and \f$\mathcal{C}\f$ are Banach spaces of functions, and
      \li \f$a \le x \le b\f$ defines pointwise (componentwise) bounds on \f$x\f$.

    This formulation encompasses a variety of useful problem scenarios.

    First, the vector spaces \f$\mathcal{X}\f$ and \f$\mathcal{C}\f$, to be defined by the user
    through the \ref la_group, can represent real spaces, such as \f$\mathcal{X} = \mathbb{R}^n\f$ and
    \f$\mathcal{C} = \mathbb{R}^m\f$, and function spaces, such as Hilbert and Banach function spaces.
    ROL's vector space abstractions enable efficent implementations of general optimization algorithms,
    spanning traditional nonlinear programming (NLP), optimization with partial differential
    equation (PDE) or differential algebraic equation (DAE) constraints, and
    stochastic optimization.

    Second, ROL's core methods can solve four types of smooth optimization problems, depending on the
    structure of the constraints.
      \li @b Type-U. No constraints (where \f$c(x) = 0\f$ and \f$a \le x \le b\f$ are absent):
          \f[
            \begin{array}{rl}
              \displaystyle \min_{x} & f(x)
            \end{array}
          \f]
          These problems are known as unconstrained optimization problems.
          The user implements the methods of the #ROL::Objective interface.
      \li @b Type-B. Bound constraints (where \f$c(x) = 0\f$ is absent):
          \f[
            \begin{array}{rl}
              \displaystyle \min_{x} & f(x) \\
              \mbox{subject to} & a \le x \le b \,.
            \end{array}
          \f]
          This problem is typically handled using projections on active sets or primal-dual active-set methods.
          ROL provides example implementations of the projections for simple box constraints.
          Other projections are defined by the user.
          The user implements the methods of the #ROL::BoundConstraint interface.
      \li @b Type-E. Equality constraints, generally nonlinear and nonconvex (where \f$a \le x \le b\f$ is absent):
          \f[
            \begin{array}{rl}
              \displaystyle \min_{x} & f(x) \\
              \mbox{subject to} & c(x) = 0 \,.
            \end{array}
          \f]
          Equality constraints are handled in ROL using, for example, matrix-free composite-step methods, including sequential quadratic programming (SQP).
          The user implements the methods of the #ROL::EqualityConstraint interface.
      \li @b Type-EB. Equality and bound constraints:
          \f[
            \begin{array}{rl}
              \displaystyle \min_{x} & f(x) \\
              \mbox{subject to} & c(x) = 0 \\
                                & a \le x \le b \,.
            \end{array}
          \f]
          This formulation includes general inequality constraints.
          For example, we can consider the reformulation:
          \f[
            \begin{array}{rlcccrl}
              \displaystyle \min_{x} & f(x) &&&& \displaystyle \min_{x,s} & f(x) \\
              \mbox{subject to} & c(x) \le 0 & & \quad \longleftrightarrow \quad & & \mbox{subject to} & c(x) + s = 0 \,, \\
              &&&&&& s \ge 0 \,.
            \end{array}
          \f]
          ROL uses a combination of matrix-free composite-step methods, projection methods and primal-dual active set methods to solve these problems.
          The user implements the methods of the #ROL::EqualityConstraint and the #ROL::BoundConstraint interfaces.

    Third, ROL's design enables streamlined algorithmic extensions, such as the \ref stochastic_group capability.

    ---
 */

/** @defgroup step_group Optimization Steps
 *  \brief ROL's optimization steps.
 */

/** @defgroup extensions_group Algorithmic Extensions
 *  \brief ROL's algorithmic extensions.
 */
  
/** @defgroup stochastic_group Stochastic Optimization
 *  @ingroup extensions_group
 *  \brief ROL's stochastic optimization capability.
 */

/** @defgroup risk_group Risk Measures
 *  @ingroup stochastic_group
 * \brief ROL's risk measure implementations.
*/ 

/** @defgroup dynamic_group Dynamic functions
 *  @ingroup interface_group
 *  \brief ROL's interfaces for time-dependent constraints and objectives
 */

/** @defgroup examples_group Examples
 *  \brief ROL's examples
 *  <ul>
 *  <li><b>Unconstrained Examples</b>   
 *  <ol>
 *  <li>\link rol/example/rosenbrock/example_01.cpp Minimizing the Rosenbrock function\endlink</li>
 *  <li>\link rol/example/zakharov/example_01.cpp Minimizing the Zakharov function\endlink</li>
 *  <li>\link rol/example/sacado/example_01.hpp Using Sacado with ROL\endlink</li>
 *  <li>\link rol/example/dual-spaces/rosenbrock-1/example_01.cpp Using Dual Spaces\endlink</li>
 *  </ol>
 *  <li><b>Constrained Examples</b></li>
 *  <ol>
 *  <li>\link rol/example/dual-spaces/simple-eq-constr-1/example_01.cpp Using Dual Spaces\endlink</li>
 *  <li>\link rol/example/sacado/example_02.hpp Using Sacado with ROL\endlink</li>
 *  <li>\link rol/example/poisson-control/example_01.cpp Poisson control\endlink</li>
 *  <li>\link rol/example/poisson-inversion/example_01.cpp Poisson inversion\endlink</li>
 *  <li>\link rol/example/burgers-control/example_01.cpp Burgers control\endlink</li>
 *  <li>\link rol/example/gross-pitaevskii/example_01.hpp Minimizing the Gross-Pitaevskii functional \endlink</li>
 *  <li>\link rol/example/gross-pitaevskii/example_02.hpp Gross-Pitaevskii functional with \f$H^1\f$ gradient \endlink</li>
 *  </ol>
 *  <li><b>Examples using Third Party Libraries</b></li> 
 *  <ol>
 *  <li>\link rol/example/json/example_01.cpp Using a JSON file to provide ROL with parameters\endlink</li> 
 *  </ol>
 *  </ul> 
*/  

#endif
