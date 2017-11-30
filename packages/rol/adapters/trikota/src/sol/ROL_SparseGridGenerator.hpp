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

#ifndef ROL_SPARSEGRIDGENERATOR_HPP
#define ROL_SPARSEGRIDGENERATOR_HPP

#include "ROL_SampleGenerator.hpp"
#include "ROL_Quadrature.hpp"

namespace ROL {

template<class Real>
class SparseGridGenerator : public SampleGenerator<Real> {
private:
  ROL::Ptr<Quadrature<Real> > grid_;
  ROL::Ptr<SparseGridIndexSet<Real> > indices_;
  bool adaptive_;
  QuadratureInfo info_;
  Real error_;
  int npts_;
  std::vector<int> index_;
  std::vector<int> search_index_;
  int direction_;

  ROL::Ptr<Vector<Real> > mydiff_, diff_;
  bool isVectorInit_;

  void buildDiffRule(Quadrature<Real> &outRule, const std::vector<int> &index) const;
  void splitSamples(std::vector<std::vector<Real> > &mypts, std::vector<Real> &mywts);
  void updateSamples(Quadrature<Real> &grid);

public:
  SparseGridGenerator(const ROL::Ptr<BatchManager<Real> > &bman,
                      const QuadratureInfo &info,
                      const bool adaptive = false);

  SparseGridGenerator(const ROL::Ptr<BatchManager<Real> > &bman,
                      const char* SGinfo,
                      const char* SGdata,
                      const bool isNormalized = true);

  void update(const Vector<Real> &x);
  Real computeError(std::vector<Real> &vals);
  Real computeError(std::vector<ROL::Ptr<Vector<Real> > > &vals, const Vector<Real> &x);
  void refine(void);
  void setSamples(bool inConstructor = false);
  void printIndexSet(void) const;
}; // class SparseGridGenerator

} // namespace ROL

#include <ROL_SparseGridGeneratorDef.hpp>

#endif
