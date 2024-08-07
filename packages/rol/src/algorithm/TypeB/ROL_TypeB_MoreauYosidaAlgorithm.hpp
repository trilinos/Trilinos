// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_MOREAUYOSIDAALGORITHM_HPP
#define ROL_TYPEB_MOREAUYOSIDAALGORITHM_HPP

#include "ROL_TypeB_Algorithm.hpp"
#include "ROL_MoreauYosidaObjective.hpp"
#include "ROL_Secant.hpp"

/** \class ROL::TypeB::MoreauYosidaAlgorithm
    \brief Provides an interface to run the Moreau-Yosida algorithm.
*/

namespace ROL {
namespace TypeB {

template<typename Real>
class MoreauYosidaAlgorithm : public TypeB::Algorithm<Real> {
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
  bool writeHeader_;

  bool hasEcon_;

  using TypeB::Algorithm<Real>::status_;
  using TypeB::Algorithm<Real>::state_;
  using TypeB::Algorithm<Real>::proj_;

  void initialize(Vector<Real>                &x,
                  const Vector<Real>          &g,
                  MoreauYosidaObjective<Real> &myobj,
                  BoundConstraint<Real>       &bnd,
                  Vector<Real>                &pwa,
                  std::ostream &outStream = std::cout); 

  void updateState(const Vector<Real>          &x,
                   MoreauYosidaObjective<Real> &myobj,
                   BoundConstraint<Real>       &bnd,
                   Vector<Real>                &pwa,
                   std::ostream &outStream = std::cout);
public:

  MoreauYosidaAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

  using TypeB::Algorithm<Real>::run;
  void run( Vector<Real>          &x,
            const Vector<Real>    &g, 
            Objective<Real>       &obj,
            BoundConstraint<Real> &bnd,
            std::ostream          &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream& os, const bool write_header = false ) const override;

}; // class ROL::TypeB::MoreauYosidaAlgorithm

} // namespace TypeB
} // namespace ROL

#include "ROL_TypeB_MoreauYosidaAlgorithm_Def.hpp"

#endif
