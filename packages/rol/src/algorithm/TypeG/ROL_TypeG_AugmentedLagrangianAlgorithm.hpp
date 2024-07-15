// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEG_AUGMENTEDLAGRANGIANALGORITHM_H
#define ROL_TYPEG_AUGMENTEDLAGRANGIANALGORITHM_H

#include "ROL_TypeG_Algorithm.hpp"
#include "ROL_AugmentedLagrangianObjective.hpp"
#include "ROL_Secant.hpp"

/** \class ROL::TypeG::AugmentedLagrangianAlgorithm
    \brief Provides an interface to run general constrained optimization algorithms
           using Augmented Lagrangians.
*/

namespace ROL {
namespace TypeG {


template<typename Real>
class AugmentedLagrangianAlgorithm : public TypeG::Algorithm<Real> {
private:
  const Ptr<Secant<Real>> secant_;
  ParameterList list_;
  bool useRelTol_;
  // Lagrange multiplier update
  bool useDefaultInitPen_;
  bool scaleLagrangian_;
  Real minPenaltyReciprocal_;
  Real minPenaltyLowerBound_;
  Real penaltyUpdate_;
  Real maxPenaltyParam_;
  // Optimality tolerance update
  Real optIncreaseExponent_;
  Real optDecreaseExponent_;
  Real optToleranceInitial_;
  Real optTolerance_;
  // Feasibility tolerance update
  Real feasIncreaseExponent_;
  Real feasDecreaseExponent_;
  Real feasToleranceInitial_;
  Real feasTolerance_;
  // Subproblem information
  bool print_;
  int maxit_;
  int subproblemIter_;
  std::string subStep_;
  int HessianApprox_;
  Real outerOptTolerance_;
  Real outerFeasTolerance_;
  Real outerStepTolerance_;
  // Scaling information
  bool useDefaultScaling_;
  Real fscale_;
  Real cscale_;
  // Verbosity flag
  int verbosity_;
  bool printHeader_;
  bool hasPolyProj_;

  using TypeG::Algorithm<Real>::state_;
  using TypeG::Algorithm<Real>::status_;
  using TypeG::Algorithm<Real>::proj_;

  void initialize(Vector<Real>                       &x,
                  const Vector<Real>                 &g,
                  const Vector<Real>                 &l,
                  const Vector<Real>                 &c,
                  AugmentedLagrangianObjective<Real> &alobj,
                  BoundConstraint<Real>              &bnd,
                  Constraint<Real>                   &con,
                  std::ostream                       &outStream = std::cout);

public:

  AugmentedLagrangianAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

  using TypeG::Algorithm<Real>::run;
  void run( Vector<Real>          &x,
            const Vector<Real>    &g,
            Objective<Real>       &obj,
            BoundConstraint<Real> &bnd,
            Constraint<Real>      &econ,
            Vector<Real>          &emul,
            const Vector<Real>    &eres,
            std::ostream          &outStream = std::cout) override ;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream& os, const bool print_header = false ) const override;

}; // class ROL::TypeG::AugmentedLagrangianAlgorithm

} // namespace TypeG
} // namespace ROL

#include "ROL_TypeG_AugmentedLagrangianAlgorithm_Def.hpp"

#endif
