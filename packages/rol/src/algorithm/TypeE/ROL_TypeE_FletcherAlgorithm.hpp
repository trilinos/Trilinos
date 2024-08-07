// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEE_FLETCHERALGORITHM_H
#define ROL_TYPEE_FLETCHERALGORITHM_H

#include "ROL_TypeE_Algorithm.hpp"
#include "ROL_FletcherObjectiveE.hpp"
#include "ROL_Secant.hpp"

/** \class ROL::TypeE::FletcherAlgorithm
    \brief Provides an interface to run equality constrained optimization algorithms
           using Fletcher's exact penalty.
*/

namespace ROL {
namespace TypeE {

template<typename Real>
class FletcherAlgorithm : public TypeE::Algorithm<Real> {
private:

  const Ptr<Secant<Real>> secant_;
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

  FletcherAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

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
