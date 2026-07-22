// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BATCHSTDVECTOR_H
#define ROL_BATCHSTDVECTOR_H

#include "ROL_StdVector.hpp"
#include "ROL_BatchManager.hpp"

/** \class ROL::BatchStdVector
    \brief Provides the std::vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template <class Real>
class BatchStdVector : public StdVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  const Ptr<BatchManager<Real>> bman_;

public:
  BatchStdVector(const Ptr<std::vector<Real>> &vec,
                 const Ptr<BatchManager<Real>> &bman)
   : StdVector<Real>(vec), bman_(bman) {}
   
  virtual Real dot(const Vector<Real> &x) const {
    const std::vector<Real> &xval = *(dynamic_cast<const StdVector<Real>&>(x).getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint numMySamples = yval.size();
    ROL_TEST_FOR_EXCEPTION( xval.size() != numMySamples, std::invalid_argument,
      "Error: Vectors must have the same dimension." );
    Real val(0), sum_val(0);
    for (uint i = 0; i < numMySamples; ++i) {
      val += xval[i] * yval[i];
    }
    // Global sum
    bman_->sumAll(&val,&sum_val,1);
    return sum_val;
  }

  virtual Ptr<Vector<Real>> clone(void) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint numMySamples = yval.size();
    return makePtr<BatchStdVector>(
           makePtr<std::vector<Real>>(numMySamples),bman_);
  }

  int dimension(void) const {
    Real dim = (Real)StdVector<Real>::dimension();
    Real sum = 0.;
    bman_->sumAll(&dim,&sum,1);
    return (int)sum;
  }

  Real reduce(const Elementwise::ReductionOp<Real> &r) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint numMySamples = yval.size();
    Real result = r.initialValue();
    for (uint i = 0; i < numMySamples; i++) {
      r.reduce(yval[i],result);
    }
    // Global sum
    Real sum = 0.;
    bman_->reduceAll(&result,&sum,1,r);
    return sum;
  }

  const Ptr<BatchManager<Real>> getBatchManager(void) const {
    return bman_;
  }
};

} // namespace ROL

#endif
