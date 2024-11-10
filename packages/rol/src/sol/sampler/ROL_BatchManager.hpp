// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_BATCHMANAGER_HPP
#define ROL_BATCHMANAGER_HPP

#include "ROL_Vector.hpp"
#include "ROL_Elementwise_Function.hpp"

namespace ROL {

template<class Real> 
class BatchManager {
public:
  virtual ~BatchManager() {}

  virtual int batchID() {
    return 0;
  }

  virtual int numBatches() {
    return 1;
  }

  virtual void sumAll(Real* input, Real* output, int dim) {
    for (int i=0; i<dim; ++i) output[i] = input[i];
  }

  virtual void sumAll(Vector<Real> &input, Vector<Real> &output) {
    output.set(input);
  }

  virtual void reduceAll(Real *input, Real* output, int dim,
                         const Elementwise::ReductionOp<Real> &r) {
    for (int i = 0; i < dim; ++i) output[i] = input[i];
  }

  virtual void gatherAll(const Real *send, const int ssize, Real *receive, int const rsize) const {
    for (int i = 0; i < rsize; ++i) receive[i] = send[i];
  }

  virtual void broadcast(Real *input, int cnt, int root) {}

  virtual void barrier() {}
};

}

#endif
