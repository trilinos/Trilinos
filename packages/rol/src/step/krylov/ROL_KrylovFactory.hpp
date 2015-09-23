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

#include "ROL_Types.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "ROL_Krylov.hpp"
#include "ROL_ConjugateGradients.hpp"
#include "ROL_ConjugateResiduals.hpp"

namespace ROL {
  template<class Real>
  inline Teuchos::RCP<Krylov<Real> > KrylovFactory( Teuchos::ParameterList &parlist ) {
    EKrylov ekv = StringToEKrylov(
                   parlist.sublist("General").sublist("Krylov").get("Type","Conjugate Gradients"));
    Real absTol  = parlist.sublist("General").sublist("Krylov").get("Absolute Tolerance", 1.e-4);
    Real relTol  = parlist.sublist("General").sublist("Krylov").get("Relative Tolerance", 1.e-2);
    int maxit    = parlist.sublist("General").sublist("Krylov").get("Iteration Limit", 20);
    bool inexact = parlist.sublist("General").get("Inexact Hessian-Times-A-Vector",false);
    switch(ekv) {
      case KRYLOV_CR: return Teuchos::rcp( new ConjugateResiduals<Real>(absTol,relTol,maxit,inexact) );
      case KRYLOV_CG: return Teuchos::rcp( new ConjugateGradients<Real>(absTol,relTol,maxit,inexact) );
      default:        return Teuchos::null;
    }
  }
}

#endif
