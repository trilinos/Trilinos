// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
