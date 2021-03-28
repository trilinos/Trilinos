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

#pragma once
#ifndef ROL2_TYPEU_TRUSTREGION_DECL_H
#define ROL2_TYPEU_TRUSTREGION_DECL_H

/** \class ROL2::TypeU::TrustRegion
    \brief Provides interface for and implements trust-region subproblem solvers.
*/

namespace ROL2 {
namespace TypeU {

template<class Real>
class TrustRegion {
public:

  enum class Type : std::int16_t {
    CauchyPoint = 0,
    TruncatedCG,
    SPG,
    DogLeg,
    DoubleDogLeg,
    Last
  };

  enum class Flag : std::int16_t {
    Success = 0,
    PosPredNeg,
    NPosPredPos,
    NPosPredNeg,
    TRNaN,
    QMinSufDec,
    Undefined,
    Last
  };

  static EnumMap<Type> type_dict;
  static EnumMap<Flag> flag_dict;

  virtual ~TrustRegion() = default;

  virtual void initialize( const Vector<Real>& x, 
                           const Vector<Real>& g ) {}

  virtual void solve( Vector<Real>&           s,          // Step (to be computed)
                      Real&                   snorm,      // Step norm (to be computed)
                      Real&                   pRed,       // Predicted reduction (to be computed)
                      int&                    iflag,      // Exit flag (to be computed)
                      int&                    iter,       // Iteration count (to be computed)
                      Real                    del,        // Trust-region radius
                      TrustRegionModel<Real>& model) = 0; // Trust-region model

  static Real initialRadius(       int&                    nfval,
                             const Vector<Real>&           x,
                             const Vector<Real>&           g,
                                   Vector<Real>&           Bg,
                                   Real                    fx,
                                   Real                    gnorm,
                                   Objective<Real>&        obj,
                                   TrustRegionModel<Real>& model,
                                   Real                    delMax,
                                   std::ostream&           os,
                                   bool                    print = false);

  static void analyzeRatio( Real&         rho,
                            Flag&         flag,
                            Real          fold,
                            Real          ftrial,
                            Real          pRed,
                            Real          epsi,
                            std::ostream& os = std::cout,
                            bool          print = false);

  static Real interpolateRadius( const Vector<Real>& g,
                                 const Vector<Real>& s,
                                       Real          snorm,
                                       Real          pRed,
                                       Real          fold,
                                       Real          ftrial,
                                       Real          del,
                                       Real          gamma0,
                                       Real          gamma1,
                                       Real          eta2,
                                       std::ostream& os = std::cout,
                                       bool          print = false);

  static Ptr<TrustRegion<Real>> create( ParameterList& parlist );

}; // TrustRegion

template<class Real>
EnumMap<typename TrustRegion<Real>::Type> 
TrustRegion<Real>::type_dict = { "Cauchy Point",
                                 "Truncated CG",
                                 "SPG",
                                 "DogLeg",
                                 "Double DogLeg" };
template<class Real>
EnumMap<typename TrustRegion<Real>::Flag> 
TrustRegion<Real>::flag_dict = { 
  "Both actual and predicted reductions are positive (success)",
  "Actual reduction is positive and predicted reduction is negative (impossible)",
  "Actual reduction is nonpositive and predicted reduction is positive",
  "Actual reduction is nonpositive and predicted reduction is negative (impossible)",
  "Actual and/or predicted reduction is a NaN",
  "Subproblem solution did not produce sufficient decrease",
  "Undefined"
};

} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_TRUSTREGION_DECL_H
