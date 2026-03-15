// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_UNIFORMSIMPLEXGENERATOR_HPP
#define ROL_UNIFORMSIMPLEXGENERATOR_HPP

#include "ROL_MonteCarloGenerator.hpp"

namespace ROL {

template<typename Real>
class UniformSimplexGenerator : public MonteCarloGenerator<Real> {
private: 
  const int dim_;

protected:
  std::vector<std::vector<Real>> sample(int nSamp, bool store = true, bool refine = false);

public:
  UniformSimplexGenerator(int nSamp,
                          int dim,
                          const Ptr<BatchManager<Real>> &bman, 
                          bool use_SA = false,
                          bool adaptive = false,
                          int numNewSamps = 0);

};

}

#include "ROL_UniformSimplexGeneratorDef.hpp"

#endif
