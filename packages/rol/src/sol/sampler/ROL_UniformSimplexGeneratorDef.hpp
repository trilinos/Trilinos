// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_UNIFORMSIMPLEXGENERATORDEF_HPP
#define ROL_UNIFORMSIMPLEXGENERATORDEF_HPP

namespace ROL {

template<typename Real>
std::vector<std::vector<Real>> UniformSimplexGenerator<Real>::sample(int nSamp, bool store, bool refine) {
  if (!refine) srand(MonteCarloGenerator<Real>::seed_);
  const Real zero(0), one(1);
  const int rank = SampleGenerator<Real>::batchID();
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<Real> real_dist(zero, one);

  std::vector<Real> pts(nSamp*dim_, zero), tmp(dim_+1, zero);
  tmp[dim_] = one;
  if (rank == 0) {
    // Generate samples
    for (int i = 0; i < nSamp; ++i) {
      // Generate simplex sample
      std::generate(tmp.begin()+1,tmp.end()-1,[&]() { return real_dist(gen); });
      std::sort(tmp.begin(),tmp.end());
      std::adjacent_difference(tmp.begin(),tmp.end(),tmp.begin());
      for (int j = 0; j < dim_; ++j) pts[j + i*dim_] = tmp[j+1];
    }
  }
  SampleGenerator<Real>::broadcast(&pts[0],nSamp*dim_,0);
  // Separate samples across processes
  int nProc  = SampleGenerator<Real>::numBatches();
  int frac   = nSamp / nProc;
  int rem    = nSamp % nProc;
  int N      = frac + ((rank < rem) ? 1 : 0);
  int offset = 0;
  for (int i = 0; i < rank; ++i) offset += frac + ((i < rem) ? 1 : 0);
  std::vector<std::vector<Real>> mypts;
  std::vector<Real> pt(dim_);
  for (int i = 0; i < N; ++i) {
    int I = offset+i;
    for (int j = 0; j < dim_; ++j) pt[j] = pts[j + I*dim_];
    mypts.push_back(pt);
  }
  if ( store ) {
    std::vector<Real> mywts(N, one/static_cast<Real>(nSamp));
    SampleGenerator<Real>::setPoints(mypts);
    SampleGenerator<Real>::setWeights(mywts);
  }
  return mypts;
}

template<typename Real>
UniformSimplexGenerator<Real>::UniformSimplexGenerator(int nSamp, int dim,
                                                       const Ptr<BatchManager<Real>> &bman, 
                                                       bool use_SA,
                                                       bool adaptive,
                                                       int numNewSamps)
  : MonteCarloGenerator<Real>(bman,use_SA,adaptive,numNewSamps),
    dim_(dim) {
  int nProc = SampleGenerator<Real>::numBatches();
  ROL_TEST_FOR_EXCEPTION( nSamp < nProc, std::invalid_argument,
    ">>> ERROR (ROL::UniformSimplexGenerator): Total number of samples is less than the number of batches!");
  MonteCarloGenerator<Real>::nSamp_ = nSamp;
  sample(nSamp,true,false);
}

}
#endif
