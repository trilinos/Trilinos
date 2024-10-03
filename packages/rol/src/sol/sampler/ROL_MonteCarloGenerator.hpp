// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MONTECARLOGENERATOR_HPP
#define ROL_MONTECARLOGENERATOR_HPP

#include "ROL_SampleGenerator.hpp"
#include "ROL_Distribution.hpp"

namespace ROL {

template<typename Real>
class MonteCarloGenerator : public SampleGenerator<Real> {
private: 
  std::vector<std::vector<Real>> data_; 
  const std::vector<Ptr<Distribution<Real>>> dist_;

  int  nSamp_;
  const bool use_normal_;
  const bool use_SA_;
  const bool adaptive_;
  const int  numNewSamps_;
  const bool useDist_;
  const int seed_;

  Real sum_val_;
  Real sum_val2_;
  Real sum_ng_;
  Real sum_ng2_;

  Real ierf(Real input) const;
  Real random(void) const;
  std::vector<std::vector<Real>> sample(int nSamp, bool store = true, bool refine = false);

public:
  MonteCarloGenerator(int nSamp,
                      const std::vector<Ptr<Distribution<Real>>> &dist, 
                      const Ptr<BatchManager<Real>> &bman, 
                      bool use_SA = false,
                      bool adaptive = false,
                      int numNewSamps = 0,
                      int seed = 123454321);

  MonteCarloGenerator(int nSamp,
                      std::vector<std::vector<Real>> &bounds, 
                      const Ptr<BatchManager<Real>> &bman,  
                      bool use_SA = false,
                      bool adaptive = false,
                      int numNewSamps = 0,
                      int seed = 123454321);

  MonteCarloGenerator(int nSamp,
                      const std::vector<Real> &mean,
                      const std::vector<Real> &std, 
                      const Ptr<BatchManager<Real>> &bman,
                      bool use_SA = false,
                      bool adaptive = false,
                      int numNewSamps = 0,
                      int seed = 123454321);

  void update( const Vector<Real> &x ) override;
  Real computeError( std::vector<Real> &vals ) override;
  Real computeError( std::vector<Ptr<Vector<Real>>> &vals, const Vector<Real> &x ) override;
  void refine(void) override;
  int numGlobalSamples(void) const override;

};

}

#include "ROL_MonteCarloGeneratorDef.hpp"

#endif
