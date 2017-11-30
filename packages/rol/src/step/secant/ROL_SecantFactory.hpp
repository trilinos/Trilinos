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

#ifndef ROL_SECANTFACTORY_H
#define ROL_SECANTFACTORY_H

#include "ROL_Types.hpp"

#include "Teuchos_ParameterList.hpp"
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
  inline ROL::Ptr<Secant<Real> > SecantFactory( Teuchos::ParameterList &parlist ) {
    ESecant esec = StringToESecant(
             parlist.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS") );
    int L  = parlist.sublist("General").sublist("Secant").get("Maximum Storage",10);
    int BB = parlist.sublist("General").sublist("Secant").get("Barzilai-Borwein",1);
    switch (esec) {
      case SECANT_LBFGS:           return ROL::makePtr<lBFGS<Real>>(L);
      case SECANT_LDFP:            return ROL::makePtr<lDFP<Real>>(L);
      case SECANT_LSR1:            return ROL::makePtr<lSR1<Real>>(L);
      case SECANT_BARZILAIBORWEIN: return ROL::makePtr<BarzilaiBorwein<Real>>(BB);
      default:                     return ROL::nullPtr;
    }
  }
}

#endif
