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

#ifndef ROL_KRYLOV_H
#define ROL_KRYLOV_H

/** \class ROL::Krylov
    \brief Provides definitions for Krylov solvers.
*/

#include "ROL_Vector.hpp"
#include "ROL_LinearOperator.hpp"

namespace ROL {

template<class Real>
class Krylov {

  Real absTol_;      // Absolute residual tolerance
  Real relTol_;      // Relative residual tolerance
  unsigned  maxit_;  // Maximum number of iterations

public:
  virtual ~Krylov(void) {}

  Krylov( Real absTol = 1.e-4, Real relTol = 1.e-2, unsigned maxit = 100 ) 
    : absTol_(absTol), relTol_(relTol), maxit_(maxit) {}

  Krylov( Teuchos::ParameterList &parlist ) {
    Teuchos::ParameterList &krylovList = parlist.sublist("General").sublist("Krylov");
    absTol_ = krylovList.get("Absolute Tolerance", 1.e-4);
    relTol_ = krylovList.get("Relative Tolerance", 1.e-2);
    maxit_  = krylovList.get("Iteration Limit", 100);
  }

  // Run Krylov Method
  virtual Real run( Vector<Real> &x, LinearOperator<Real> &A, const Vector<Real> &b, LinearOperator<Real> &M, 
                    int &iter, int &flag ) = 0;

  void resetAbsoluteTolerance(const Real absTol) {
    absTol_ = absTol;
  }
  void resetRelativeTolerance(const Real relTol) {
    relTol_ = relTol;
  }
  void resetMaximumIteration(const unsigned maxit) {
    maxit_ = maxit;
  }
  Real getAbsoluteTolerance(void) const {
    return absTol_;
  }
  Real getRelativeTolerance(void) const {
    return relTol_;
  }
  unsigned getMaximumIteration(void) const {
    return maxit_;
  }
};

}

#include "ROL_KrylovFactory.hpp"

#endif
