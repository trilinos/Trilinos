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

#ifndef ROL_QUASINEWTONALGORITHM_B_H
#define ROL_QUASINEWTONALGORITHM_B_H

#include "ROL_Algorithm_B.hpp"
#include "ROL_SecantFactory.hpp"

/** \class ROL::QuasiNewtonAlgorithm_B
    \brief Provides an interface to run the projected secant algorithm.
*/

namespace ROL {

template<typename Real>
class QuasiNewtonAlgorithm_B : public Algorithm_B<Real> {
private:
  Ptr<Secant<Real>> secant_; ///< Secant object (used for quasi-Newton)
  ESecant esec_;             ///< Secant type
  std::string secantName_;   ///< Secant name

  int maxit_;         ///< Maximum number of line search steps (default: 20)
  Real rhodec_;       ///< Backtracking rate (default: 0.5)
  Real c1_;           ///< Sufficient Decrease Parameter (default: 1e-4)
  Real sigma1_;       ///< Lower safeguard for quadratic line search (default: 0.1)
  Real sigma2_;       ///< Upper safeguard for quadratic line search (default: 0.9)
  std::string algoName_;

  ParameterList list_;

  bool hasLEC_;
  int ls_nfval_, spgIter_;
  int verbosity_;
  bool printHeader_;

  using Algorithm_B<Real>::status_;
  using Algorithm_B<Real>::state_;
  using Algorithm_B<Real>::proj_;

  void initialize(Vector<Real>          &x,
                  const Vector<Real>    &g,
                  Objective<Real>       &obj,
                  BoundConstraint<Real> &bnd,
                  std::ostream &outStream = std::cout); 

public:

  QuasiNewtonAlgorithm_B(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

  using Algorithm_B<Real>::run;
  std::vector<std::string> run( Vector<Real>          &x,
                                const Vector<Real>    &g, 
                                Objective<Real>       &obj,
                                BoundConstraint<Real> &bnd,
                                std::ostream          &outStream = std::cout);

  std::string printHeader( void ) const override;

  std::string printName( void ) const override;

  std::string print( const bool print_header = false ) const override;

}; // class ROL::QuasiNewtonAlgorithm_B

} // namespace ROL

#include "ROL_QuasiNewtonAlgorithm_B_Def.hpp"

#endif
