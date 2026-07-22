// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  int numGlobalSamples(void) const {
    return npts_;
  }
}; // class SparseGridGenerator

} // namespace ROL

#include <ROL_SparseGridGeneratorDef.hpp>

#endif
