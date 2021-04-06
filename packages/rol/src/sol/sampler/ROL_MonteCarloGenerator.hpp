// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
