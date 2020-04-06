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

#ifndef ROL_MOREAUYOSIDAALGORITHM_G_H
#define ROL_MOREAUYOSIDAALGORITHM_G_H

#include "ROL_Algorithm_G.hpp"
#include "ROL_MoreauYosidaObjective.hpp"

/** \class ROL::MoreauYosidaAlgorithm_G
    \brief Provides an interface to run the Moreau-Yosida algorithm.
*/

namespace ROL {

template<typename Real>
class MoreauYosidaAlgorithm_G : public Algorithm_G<Real> {
private:
  Real compViolation_;
  Real gnorm_;
  Real maxPenalty_;
  Real tau_;
  bool print_;
  bool updatePenalty_;
  bool updateMultiplier_;

  ROL::ParameterList list_;
  int subproblemIter_;

  std::string stepname_;

  int verbosity_;
  bool printHeader_;

  bool hasPolyProj_;

  using Algorithm_G<Real>::status_;
  using Algorithm_G<Real>::state_;
  using Algorithm_G<Real>::proj_;

  void initialize(Vector<Real>                &x,
                  const Vector<Real>          &g,
                  const Vector<Real>          &l,
                  const Vector<Real>          &c,
                  MoreauYosidaObjective<Real> &myobj,
                  BoundConstraint<Real>       &bnd,
                  Constraint<Real>            &con,
                  Vector<Real>                &pwa,
                  Vector<Real>                &dwa,
                  std::ostream &outStream = std::cout); 

  void updateState(const Vector<Real>          &x,
                   const Vector<Real>          &l,
                   MoreauYosidaObjective<Real> &myobj,
                   BoundConstraint<Real>       &bnd,
                   Constraint<Real>            &con,
                   Vector<Real>                &pwa,
                   Vector<Real>                &dwa,
                   std::ostream &outStream = std::cout);
public:

  MoreauYosidaAlgorithm_G(ParameterList &list);

  using Algorithm_G<Real>::run;
  std::vector<std::string> run( Vector<Real>          &x,
                                const Vector<Real>    &g, 
                                Objective<Real>       &obj,
                                BoundConstraint<Real> &bnd,
                                Constraint<Real>      &econ,
                                Vector<Real>          &emul,
                                const Vector<Real>    &eres,
                                std::ostream          &outStream = std::cout);

  std::string printHeader( void ) const override;

  std::string printName( void ) const override;

  std::string print( const bool print_header = false ) const override;

}; // class ROL::MoreauYosidaAlgorithm_G

} // namespace ROL

#include "ROL_MoreauYosidaAlgorithm_G_Def.hpp"

#endif
