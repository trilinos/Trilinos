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
G// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#pragma once
#ifndef ROL2_KRYLOV_DECL_HPP
#define ROL2_KRYLOV_DECL_HPP

/** \class ROL2::Krylov
    \brief Provides definitions for Krylov solvers.
*/

namespace ROL2 {

template<class Real>
class Krylov {
public:

  enum class Type : std::int16_t { 
    CG = 0,
    CR,
    GMRES,
    MINRES,
    UserDefined,
    Last
  };

  virtual ~Krylov() = default;

  Krylov( Real          absTol = 1.e-4, 
          Real          relTol = 1.e-2,
          std::uint16_t maxit = 100 )
    : absTol_(absTol), relTol_(relTol), maxit_(maxit) {}

  Krylov( ParameterList& parlist ) {
    ROL::ParameterList& krylovList = parlist.sublist("General").sublist("Krylov");
    absTol_ = krylovList.get("Absolute Tolerance", 1.e-4);
    relTol_ = krylovList.get("Relative Tolerance", 1.e-2);
    maxit_  = krylovList.get("Iteration Limit", 100);
  }

  // Run Krylov Method
  virtual Real run(       Vector<Real>&         x, 
                          LinearOperator<Real>& A,
                    const Vector<Real>&         b, 
                          LinearOperator<Real>& M, 
                          int&                  iter, 
                          int&                  flag ) = 0;

  void resetAbsoluteTolerance(Real absTol) {
    absTol_ = absTol;
  }

  void resetRelativeTolerance(Real relTol) {
    relTol_ = relTol;
  }

  void resetMaximumIteration(std::uint16_t maxit) {
    maxit_ = maxit;
  }

  Real getAbsoluteTolerance() const {
    return absTol_;
  }

  Real getRelativeTolerance() const {
    return relTol_;
  }

  std::uint16_t getMaximumIteration() const {
    return maxit_;
  }

  static EnumMap<Type> type_dict;

  static Ptr<Krylov> create( ParameterList& parlist );

private:
  Real           absTol_;      // Absolute residual tolerance
  Real           relTol_;      // Relative residual tolerance
  std::uint16_t  maxit_;  // Maximum number of iterations
}; // class Krylov

template<class Real>
EnumMap<Krylov<Real>::Type>
Krylov<Real>::type_dict = { "Conjugate Gradient",
                            "Conjugate Residuals",
                            "GMRES",
                            "MINRES",
                            "User Defined" };

} // namespace ROL2

#endif // ROL2_KRYLOV_DECL_HPP
