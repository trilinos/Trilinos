// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_GRADIENTALGORITHM_HPP
#define ROL_TYPEB_GRADIENTALGORITHM_HPP

#include "ROL_TypeB_Algorithm.hpp"

/** \class ROL::TypeB::GradientAlgorithm
    \brief Provides an interface to run the projected gradient algorithm.
*/

namespace ROL {
namespace TypeB {

template<typename Real>
class GradientAlgorithm : public TypeB::Algorithm<Real> {
private:
  int maxit_;
  Real alpha0_, alpha0bnd_, rhodec_, rhoinc_, c1_, maxAlpha_;
  bool useralpha_, usePrevAlpha_, useAdapt_, normAlpha_;
  int verbosity_;
  bool writeHeader_;

  using TypeB::Algorithm<Real>::status_;
  using TypeB::Algorithm<Real>::state_;
  using TypeB::Algorithm<Real>::proj_;

  void initialize(Vector<Real>          &x,
                  const Vector<Real>    &g,
                  Objective<Real>       &obj,
                  BoundConstraint<Real> &bnd,
                  std::ostream &outStream = std::cout); 

public:

  GradientAlgorithm(ParameterList &list);

  using TypeB::Algorithm<Real>::run;
  void run( Vector<Real>          &x,
            const Vector<Real>    &g, 
            Objective<Real>       &obj,
            BoundConstraint<Real> &bnd,
            std::ostream          &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream& os, const bool write_header = false ) const override;

}; // class ROL::TypeB::GradientAlgorithm

} // namespace TypeB
} // namespace ROL

#include "ROL_TypeB_GradientAlgorithm_Def.hpp"

#endif
