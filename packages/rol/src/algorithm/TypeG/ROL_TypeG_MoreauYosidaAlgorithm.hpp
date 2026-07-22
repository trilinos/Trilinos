// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEG_MOREAUYOSIDAALGORITHM_H
#define ROL_TYPEG_MOREAUYOSIDAALGORITHM_H

#include "ROL_TypeG_Algorithm.hpp"
#include "ROL_MoreauYosidaObjective.hpp"
#include "ROL_Secant.hpp"

/** \class ROL::TypeG::MoreauYosidaAlgorithm
    \brief Provides an interface to run the Moreau-Yosida algorithm.
*/

namespace ROL {
namespace TypeG {

template<typename Real>
class MoreauYosidaAlgorithm : public TypeG::Algorithm<Real> {
private:
  const Ptr<Secant<Real>> secant_;

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

  using TypeG::Algorithm<Real>::status_;
  using TypeG::Algorithm<Real>::state_;
  using TypeG::Algorithm<Real>::proj_;

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

  MoreauYosidaAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

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

}; // class ROL::TypeG::MoreauYosidaAlgorithm

} // namespace TypeG
} // namespace ROL

#include "ROL_TypeG_MoreauYosidaAlgorithm_Def.hpp"

#endif
