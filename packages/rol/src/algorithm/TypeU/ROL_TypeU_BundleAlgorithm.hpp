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

#ifndef ROL_TYPEU_BUNDLEALGORITHM_H
#define ROL_TYPEU_BUNDLEALGORITHM_H

#include "ROL_TypeU_Algorithm.hpp"
#include "ROL_Bundle_U.hpp"
#include "ROL_LineSearch_U.hpp"

/** \class ROL::TypeU::BundleAlgorithm
    \brief Provides an interface to run trust-bundle methods for unconstrained
           optimization algorithms.
*/

namespace ROL {
namespace TypeU {

template<typename Real>
class BundleAlgorithm : public Algorithm<Real> {
private:
  // Bundle
  Ptr<Bundle_U<Real>>     bundle_;     // Bundle of subgradients and linearization errors
  Ptr<LineSearch_U<Real>> lineSearch_; // Line-search object for nonconvex problems

  // Dual cutting plane solution
  unsigned QPiter_;  // Number of QP solver iterations
  unsigned QPmaxit_; // Maximum number of QP iterations
  Real QPtol_;       // QP subproblem tolerance

  // Step flag
  int step_flag_; // Whether serious or null step

  // Aggregate subgradients, linearizations, and distance measures

  // Algorithmic parameters
  Real T_;
  Real tol_;
  Real m1_;
  Real m2_;
  Real m3_;
  Real nu_;

  // Line-search parameters
  int ls_maxit_;

  bool first_print_;
  bool isConvex_;

  int verbosity_;
  bool printHeader_;

  using Algorithm<Real>::state_;
  using Algorithm<Real>::status_;
  using Algorithm<Real>::initialize;

public:

  BundleAlgorithm( ParameterList &parlist,
                   const Ptr<LineSearch_U<Real>> &lineSearch = nullPtr );

  void run( Vector<Real>       &x,
            const Vector<Real> &g, 
            Objective<Real>    &obj,
            std::ostream       &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os) const override;

  void writeOutput( std::ostream& os, bool print_header = false ) const override;

private:

  void initialize(const Vector<Real> &x, const Vector<Real> &g,
                  Objective<Real> &obj, std::ostream &outStream = std::cout);

}; // class ROL::BundleAlgorithm
} // namespace TypeU
} // namespace ROL

#include "ROL_TypeU_BundleAlgorithm_Def.hpp"

#endif
