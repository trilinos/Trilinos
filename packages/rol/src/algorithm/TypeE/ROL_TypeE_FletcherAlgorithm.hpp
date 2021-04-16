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

#ifndef ROL_TYPEE_FLETCHERALGORITHM_H
#define ROL_TYPEE_FLETCHERALGORITHM_H

#include "ROL_TypeE_Algorithm.hpp"
#include "ROL_FletcherObjectiveE.hpp"

/** \class ROL::TypeE::FletcherAlgorithm
    \brief Provides an interface to run equality constrained optimization algorithms
           using Fletcher's exact penalty.
*/

namespace ROL {
namespace TypeE {

template<typename Real>
class FletcherAlgorithm : public TypeE::Algorithm<Real> {
private:

  ParameterList list_;
  std::string subStep_;
  // Penalty function data
  Real merit_, gpnorm_;
  Real sigma_, delta_;
  Real minSigma_, maxSigma_, sigmaUpdate_;
  Real minDelta_, deltaUpdate_;
  bool modifySigma_;
  int subproblemIter_;
  // Verbosity flag
  int verbosity_;
  bool printHeader_;

  using TypeE::Algorithm<Real>::state_;
  using TypeE::Algorithm<Real>::status_;

public:

  FletcherAlgorithm(ParameterList &list);

  void initialize( Vector<Real>             &x,
                   const Vector<Real>       &g,
                   const Vector<Real>       &l,
                   const Vector<Real>       &c,
                   FletcherObjectiveE<Real> &fobj,
                   Constraint<Real>         &con,
                   std::ostream             &outStream);

  using TypeE::Algorithm<Real>::run;
  virtual void run( Vector<Real>       &x,
                    const Vector<Real> &g, 
                    Objective<Real>    &obj,
                    Constraint<Real>   &econ,
                    Vector<Real>       &emul,
                    const Vector<Real> &eres,
                    std::ostream       &outStream = std::cout) override;

  virtual void writeHeader( std::ostream& os ) const override;

  virtual void writeName( std::ostream& os ) const override;

  virtual void writeOutput( std::ostream& os, const bool print_header = false ) const override;

}; // class ROL::TypeE::FletcherAlgorithm

} // namespace TypeE
} // namespace ROL

#include "ROL_TypeE_FletcherAlgorithm_Def.hpp"

#endif
