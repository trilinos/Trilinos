// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_PRINTDESIGN_HPP
#define ROL_OED_PRINTDESIGN_HPP

namespace ROL::OED {
template<typename Real>
void PrintDesign(const Ptr<const Vector<Real>>& input,
                 const Ptr<const SampleGenerator<Real>>& sampler,
	         const std::string& name, const std::string& ext=".txt") {
  int nBatch = sampler->numBatches();
  int myRank = sampler->batchID();
  int nsamp  = sampler->numMySamples();
  int dim    = sampler->getMyPoint(0).size();
  std::stringstream filename;
  filename << name;
  if (nBatch > 1) filename << "_" << myRank;
  filename << ext;
  std::ofstream file;
  file.open(filename.str());
  file << std::scientific << std::setprecision(15);
  for (int i = 0; i < nsamp; ++i) {
    for (int j = 0; j < dim; ++j) {
      file << std::right << std::setw(25)
           << sampler->getMyPoint(i)[j];
    }
    file << std::right << std::setw(25)
         << dynamicPtrCast<const ProbabilityVector<Real>>(input)->getProbability(i)
         << std::endl;
  }
  file.close();
}
}

#endif
