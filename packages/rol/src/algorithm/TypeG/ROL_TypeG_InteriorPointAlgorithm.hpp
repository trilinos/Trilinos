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

#ifndef ROL_TYPEG_INTERIORPOINTALGORITHM_H
#define ROL_TYPEG_INTERIORPOINTALGORITHM_H

#include "ROL_TypeG_Algorithm.hpp"
#include "ROL_InteriorPointObjective.hpp"

/** \class ROL::TypeG::InteriorPointAlgorithm
    \brief Provides an interface to run the interior point algorithm.
*/

namespace ROL {
namespace TypeG {

template<typename Real>
class InteriorPointAlgorithm : public TypeG::Algorithm<Real> {
private:
  Real mumin_;
  Real mumax_;
  Real rho_;
  bool useLinearDamping_;
  Real kappaD_;
  Real gtol_;
  Real ctol_;
  Real stol_;
  Real gtolrate_;
  Real ctolrate_;
  Real mingtol_;
  Real minctol_;

  ROL::ParameterList list_;
  int subproblemIter_;

  std::string stepname_;

  bool print_;
  int verbosity_;
  bool printHeader_;

  bool hasPolyProj_;

  using TypeG::Algorithm<Real>::status_;
  using TypeG::Algorithm<Real>::state_;
  using TypeG::Algorithm<Real>::proj_;

  void initialize(Vector<Real>                 &x,
                  const Vector<Real>           &g,
                  const Vector<Real>           &l,
                  const Vector<Real>           &c,
                  InteriorPointObjective<Real> &ipobj,
                  BoundConstraint<Real>        &bnd,
                  Constraint<Real>             &con,
                  Vector<Real>                 &pwa,
                  Vector<Real>                 &dwa,
                  std::ostream &outStream = std::cout); 

  void updateState(const Vector<Real>           &x,
                   const Vector<Real>           &l,
                   InteriorPointObjective<Real> &ipobj,
                   BoundConstraint<Real>        &bnd,
                   Constraint<Real>             &con,
                   Vector<Real>                 &pwa,
                   Vector<Real>                 &dwa,
                   std::ostream                 &outStream = std::cout);
public:

  InteriorPointAlgorithm(ParameterList &list);

  using TypeG::Algorithm<Real>::run;
  void run( Vector<Real>          &x,
            const Vector<Real>    &g, 
            Objective<Real>       &obj,
            BoundConstraint<Real> &bnd,
            Constraint<Real>      &econ,
            Vector<Real>          &emul,
            const Vector<Real>    &eres,
            std::ostream          &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream& os, const bool print_header = false ) const override;

}; // class ROL::TypeG::InteriorPointAlgorithm

} // namespace TypeG
} // namespace ROL

#include "ROL_TypeG_InteriorPointAlgorithm_Def.hpp"

#endif
