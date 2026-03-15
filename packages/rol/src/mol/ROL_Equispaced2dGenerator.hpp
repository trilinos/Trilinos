// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_EQUISPACED2DGENERATOR_HPP
#define ROL_EQUISPACED2DGENERATOR_HPP

#include "ROL_SampleGenerator.hpp"

namespace ROL {

template<class Real>
class Equispaced2dGenerator : public SampleGenerator<Real> {
private:
  int nsamp_;

public:
  Equispaced2dGenerator(int nsamp, const ROL::Ptr<BatchManager<Real>> &bman, bool useTrig = false)
    : SampleGenerator<Real>(bman), nsamp_(nsamp) {
    const Real zero(0), one(1), pi2(M_PI*0.5);
    int rank = SampleGenerator<Real>::batchID();
    std::vector<Real> pts(2*nsamp, zero);
    if (rank == 0) {
      // Generate samples
      Real lam(0);
      const Real dlam = static_cast<Real>(1)/static_cast<Real>(nsamp+1);
      for (int i = 0; i < nsamp; ++i) {
        lam = static_cast<Real>(i+1)*dlam; 
	if (useTrig) pts[1 + i*2] = std::sin(pi2*lam)/(std::cos(pi2*lam)+std::sin(pi2*lam));
	else         pts[1 + i*2] = lam;
        pts[0 + i*2] = one - pts[1 + i*2];
      }
    }
    SampleGenerator<Real>::broadcast(&pts[0],nsamp*2,0);
    // Separate samples across processes
    int nproc  = SampleGenerator<Real>::numBatches();
    int frac   = nsamp / nproc;
    int rem    = nsamp % nproc;
    int N      = frac + ((rank < rem) ? 1 : 0);
    int offset = 0;
    for (int i = 0; i < rank; ++i) offset += frac + ((i < rem) ? 1 : 0);
    std::vector<std::vector<Real>> mypts;
    std::vector<Real> pt(2);
    for (int i = 0; i < N; ++i) {
      int I = offset+i;
      for (int j = 0; j < 2; ++j) pt[j] = pts[j + 2*I];
      mypts.push_back(pt);
    }
    std::vector<Real> mywts(N, one/static_cast<Real>(nsamp));
    SampleGenerator<Real>::setPoints(mypts);
    SampleGenerator<Real>::setWeights(mywts);
  }
};

}

#endif
