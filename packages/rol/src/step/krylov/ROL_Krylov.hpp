// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_KRYLOV_H
#define ROL_KRYLOV_H

/** \class ROL::Krylov
    \brief Provides definitions for Krylov solvers.
*/

#include "ROL_Vector.hpp"
#include "ROL_LinearOperator.hpp"
#include "ROL_ParameterList.hpp"

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

  Krylov( ROL::ParameterList &parlist ) {
    ROL::ParameterList &krylovList = parlist.sublist("General").sublist("Krylov");
    absTol_ = krylovList.get("Absolute Tolerance", 1.e-4);
    relTol_ = krylovList.get("Relative Tolerance", 1.e-2);
    maxit_  = krylovList.get("Iteration Limit", 100);
  }

  // Run Krylov Method
  virtual Real run( Vector<Real> &x, LinearOperator<Real> &A,
              const Vector<Real> &b, LinearOperator<Real> &M, 
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

#endif
