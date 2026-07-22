// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PROJECTION_PARTITIONED_H
#define ROL_PROJECTION_PARTITIONED_H

#include "ROL_PartitionedVector.hpp"
#include <iostream>

namespace ROL {

template<typename Real>
class Projection_Partitioned : public Projection<Real> {
private:
  std::vector<Ptr<Projection<Real>>> pvec_;

public:
  Projection_Partitioned(const std::vector<Ptr<Projection<Real>>> &pvec)
  : pvec_(pvec)
  {}


  virtual void project(Vector<Real> &x, std::ostream &stream = std::cout) {
    PartitionedVector<Real> &xpv = static_cast<PartitionedVector<Real>&>(x);
    for (unsigned i = 0; i < pvec_.size(); ++i) pvec_[i]->project(*xpv.get(i),stream);
  }

  virtual void applyJacobian(Vector<Real> &v, const Vector<Real> &x) {
    PartitionedVector<Real> &vpv = dynamic_cast<PartitionedVector<Real>&>(v);
    const PartitionedVector<Real> &xpv = dynamic_cast<const PartitionedVector<Real>&>(x);
    for (unsigned i = 0; i < pvec_.size(); ++i)
      pvec_[i]->applyJacobian(*vpv.get(i),*xpv.get(i));
  }

}; // class Projection_Partitioned

} // namespace ROL

#endif
