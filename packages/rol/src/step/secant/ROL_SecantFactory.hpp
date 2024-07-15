// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SECANTFACTORY_H
#define ROL_SECANTFACTORY_H

#include "ROL_Types.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_Ptr.hpp"

#include "ROL_Secant.hpp"
#include "ROL_lBFGS.hpp"
#include "ROL_lDFP.hpp"
#include "ROL_lSR1.hpp"
#include "ROL_BarzilaiBorwein.hpp"

namespace ROL {
  template<class Real>
  inline ROL::Ptr<Secant<Real> > getSecant( ESecant esec = SECANT_LBFGS, int L = 10, int BBtype = 1 ) {
    switch (esec) {
      case SECANT_LBFGS:           return ROL::makePtr<lBFGS<Real>>(L);
      case SECANT_LDFP:            return ROL::makePtr<lDFP<Real>>(L);
      case SECANT_LSR1:            return ROL::makePtr<lSR1<Real>>(L);
      case SECANT_BARZILAIBORWEIN: return ROL::makePtr<BarzilaiBorwein<Real>>(BBtype);
      default:                     return ROL::nullPtr; 
    }
  }

  template<class Real>
  inline ROL::Ptr<Secant<Real> > SecantFactory( ROL::ParameterList &parlist, ESecantMode mode = SECANTMODE_BOTH ) {
    ESecant esec = StringToESecant(
             parlist.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS") );
    int L    = parlist.sublist("General").sublist("Secant").get("Maximum Storage",10);
    int BB   = parlist.sublist("General").sublist("Secant").get("Barzilai-Borwein",1);
    bool uds = parlist.sublist("General").sublist("Secant").get("Use Default Scaling",true);
    Real s   = parlist.sublist("General").sublist("Secant").get("Initial Hessian Scale",1.0);
    switch (esec) {
      case SECANT_LBFGS:           return ROL::makePtr<lBFGS<Real>>(L,uds,s);
      case SECANT_LDFP:            return ROL::makePtr<lDFP<Real>>(L,uds,s);
      case SECANT_LSR1:            return ROL::makePtr<lSR1<Real>>(L,uds,s,mode);
      case SECANT_BARZILAIBORWEIN: return ROL::makePtr<BarzilaiBorwein<Real>>(BB);
      default:                     return ROL::nullPtr;
    }
  }
}

#endif
