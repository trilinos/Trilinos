// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_KRYLOVFACTORY_H
#define ROL_KRYLOVFACTORY_H

#include "ROL_Ptr.hpp"
#include "ROL_Types.hpp"

#include "ROL_ConjugateGradients.hpp"
#include "ROL_ConjugateResiduals.hpp"
#include "ROL_GMRES.hpp"
#include "ROL_MINRES.hpp"

namespace ROL {
  /** \enum   ROL::EKrylov
      \brief  Enumeration of Krylov methods.

      \arg    CG          Conjugate Gradient Method
      \arg    CR          Conjugate Residual Method
      \arg    GMRES       Generalized Minimum Residual Method
      \arg    MINRES      Minimum Residual Method
      \arg    USERDEFINED User defined Krylov method
      \arg    LAST        Dummy type
   */
  enum EKrylov{
    KRYLOV_CG = 0,
    KRYLOV_CR,
    KRYLOV_GMRES,
    KRYLOV_MINRES,
    KRYLOV_USERDEFINED,
    KRYLOV_LAST
  };

  inline std::string EKrylovToString(EKrylov type) {
    std::string retString;
    switch(type) {
      case KRYLOV_CG:          retString = "Conjugate Gradients"; break;
      case KRYLOV_CR:          retString = "Conjugate Residuals"; break;
      case KRYLOV_GMRES:       retString = "GMRES";               break;
      case KRYLOV_MINRES:      retString = "MINRES";              break;
      case KRYLOV_USERDEFINED: retString = "User Defined";        break;
      case KRYLOV_LAST:        retString = "Last Type (Dummy)";   break;
      default:                 retString = "INVALID EKrylov";
    }
    return retString;
  }

  /** \brief  Verifies validity of a Krylov enum.
    
      \param  type  [in]  - enum of the Krylov
      \return 1 if the argument is a valid Secant; 0 otherwise.
    */
  inline int isValidKrylov(EKrylov type){
    return( (type == KRYLOV_CG)      ||
            (type == KRYLOV_CR)      ||
            (type == KRYLOV_GMRES)   ||
            (type == KRYLOV_USERDEFINED) );
  }

  inline EKrylov & operator++(EKrylov &type) {
    return type = static_cast<EKrylov>(type+1);
  }

  inline EKrylov operator++(EKrylov &type, int) {
    EKrylov oldval = type;
    ++type;
    return oldval;
  }

  inline EKrylov & operator--(EKrylov &type) {
    return type = static_cast<EKrylov>(type-1);
  }

  inline EKrylov operator--(EKrylov &type, int) {
    EKrylov oldval = type;
    --type;
    return oldval;
  }

  inline EKrylov StringToEKrylov(std::string s) {
    s = removeStringFormat(s);
    for ( EKrylov type = KRYLOV_CG; type < KRYLOV_LAST; type++ ) {
      if ( !s.compare(removeStringFormat(EKrylovToString(type))) ) {
        return type;
      }
    }
    return KRYLOV_CG;
  }

  template<class Real>
  inline Ptr<Krylov<Real>> KrylovFactory( ParameterList &parlist ) {
    Real em4(1e-4), em2(1e-2);
    EKrylov ekv = StringToEKrylov(
                   parlist.sublist("General").sublist("Krylov").get("Type","GMRES"));
    Real absTol  = parlist.sublist("General").sublist("Krylov").get("Absolute Tolerance", em4);
    Real relTol  = parlist.sublist("General").sublist("Krylov").get("Relative Tolerance", em2);
    int maxit    = parlist.sublist("General").sublist("Krylov").get("Iteration Limit", 20);
    bool inexact = parlist.sublist("General").get("Inexact Hessian-Times-A-Vector",false);
    switch(ekv) {
      case KRYLOV_CR: 
        return makePtr<ConjugateResiduals<Real>>(absTol,relTol,maxit,inexact);
      case KRYLOV_CG: 
        return makePtr<ConjugateGradients<Real>>(absTol,relTol,maxit,inexact);
      case KRYLOV_MINRES:
        return makePtr<MINRES<Real>>(absTol,relTol,maxit,inexact);
      case KRYLOV_GMRES:
        return makePtr<GMRES<Real>>(parlist);
      default:
        return nullPtr;
    }
  }

}

#endif
