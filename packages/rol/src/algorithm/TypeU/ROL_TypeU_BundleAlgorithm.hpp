// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

  void writeOutput( std::ostream& os, const bool print_header = false ) const override;

private:

  void initialize(const Vector<Real> &x, const Vector<Real> &g,
                  Objective<Real> &obj, std::ostream &outStream = std::cout);

}; // class ROL::BundleAlgorithm
} // namespace TypeU
} // namespace ROL

#include "ROL_TypeU_BundleAlgorithm_Def.hpp"

#endif
