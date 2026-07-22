// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
