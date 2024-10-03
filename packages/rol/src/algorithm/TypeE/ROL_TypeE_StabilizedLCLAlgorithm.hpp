// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEE_STABILIZEDLCLALGORITHM_H
#define ROL_TYPEE_STABILIZEDLCLALGORITHM_H

#include "ROL_TypeE_Algorithm.hpp"
#include "ROL_AugmentedLagrangianObjective.hpp"
#include "ROL_ElasticLinearConstraint.hpp"
#include "ROL_ElasticObjective.hpp"
#include "ROL_Secant.hpp"

/** \class ROL::TypeE::StabilizedLCLAlgorithm
    \brief Provides an interface to run equality constrained optimization algorithms
           using Stabilized LCL.
*/

namespace ROL {
namespace TypeE {


template<typename Real>
class StabilizedLCLAlgorithm : public TypeE::Algorithm<Real> {
private:
  const Ptr<Secant<Real>> secant_;
  ParameterList list_;
  // Lagrange multiplier update
  bool useDefaultInitPen_;
  bool scaleLagrangian_;
  Real penaltyUpdate_;
  Real maxPenaltyParam_;
  Real sigma_;
  Real sigmaMax_;
  Real sigmaUpdate_;
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

  using TypeE::Algorithm<Real>::state_;
  using TypeE::Algorithm<Real>::status_;

  void initialize(Vector<Real>           &x,
                  const Vector<Real>     &g,
                  const Vector<Real>     &l,
                  const Vector<Real>     &c,
                  ElasticObjective<Real> &alobj,
                  Constraint<Real>       &con,
                  std::ostream           &outStream = std::cout);

public:

  StabilizedLCLAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

  using TypeE::Algorithm<Real>::run;
  virtual void run( Problem<Real> &problem,
                    std::ostream  &outStream = std::cout) override;

  virtual void run( Vector<Real>          &x,
                    const Vector<Real>    &g, 
                    Objective<Real>       &obj,
                    Constraint<Real>      &econ,
                    Vector<Real>          &emul,
                    const Vector<Real>    &eres,
                    std::ostream          &outStream = std::cout) override;

  virtual void writeHeader( std::ostream& os ) const override;

  virtual void writeName( std::ostream& os ) const override;

  virtual void writeOutput( std::ostream& os, const bool print_header = false ) const override;

}; // class ROL::TypeE::StabilizedLCLAlgorithm

} // namespace TypeE
} // namespace ROL

#include "ROL_TypeE_StabilizedLCLAlgorithm_Def.hpp"

#endif
