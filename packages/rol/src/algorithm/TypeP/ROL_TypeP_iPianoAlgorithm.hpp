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

#ifndef ROL_TYPEP_IPIANOALGORITHM_HPP
#define ROL_TYPEP_IPIANOALGORITHM_HPP

#include "ROL_TypeP_Algorithm.hpp"

/** \class ROL::TypeP::iPianoAlgorithm
    \brief Provides an interface to run the proximal gradient algorithm.
*/

namespace ROL {
namespace TypeP {

template<typename Real>
class iPianoAlgorithm : public TypeP::Algorithm<Real> {
private:
  int maxit_;
  Real t0_, alpha_, beta_, rhodec_, rhoinc_, c1_, c2_, L_;
  bool useConstBeta_, initProx_;
  int verbosity_;
  bool writeHeader_;

  using TypeP::Algorithm<Real>::status_;
  using TypeP::Algorithm<Real>::state_;
  using TypeP::Algorithm<Real>::pgstep;

  void initialize(Vector<Real>          &x,
                  const Vector<Real>    &g,
                  Objective<Real>       &sobj,
                  Objective<Real>       &nobj,
                  Vector<Real>          &px,
                  Vector<Real>          &dg,
                  std::ostream &outStream = std::cout); 
public:

  iPianoAlgorithm(ParameterList &list);

  using TypeP::Algorithm<Real>::run;
  void run( Vector<Real>          &x,
            const Vector<Real>    &g, 
            Objective<Real>       &sobj,
            Objective<Real>       &nobj,
            std::ostream          &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream& os, bool write_header = false ) const override;

}; // class ROL::TypeP::iPianoAlgorithm

} // namespace TypeP
} // namespace ROL

#include "ROL_TypeP_iPianoAlgorithm_Def.hpp"

#endif
