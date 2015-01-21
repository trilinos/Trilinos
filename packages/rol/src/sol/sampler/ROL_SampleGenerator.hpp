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

#ifndef ROL_SAMPLEGENERATOR_HPP
#define ROL_SAMPLEGENERATOR_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Comm.h"
#include "ROL_BatchManager.hpp"
#include "ROL_Vector.hpp"

namespace ROL {

template<class Real> 
class SampleGenerator {
private:
  int begin_;
  Teuchos::RCP<BatchManager<Real> > bman_;
  std::vector<std::vector<Real> > points_;
  std::vector<Real> weights_;

protected:
  void setPoints(std::vector<std::vector<Real> > &p) {
    points_.clear();
    points_.insert(points_.begin(),p.begin(),p.end());
  }
  void setWeights(std::vector<Real> &w) {
    weights_.clear();
    weights_.insert(weights_.begin(),w.begin(),w.end());
  }

public:
  virtual ~SampleGenerator() {}
  SampleGenerator(Teuchos::RCP<BatchManager<Real> > &bman) : begin_(0), bman_(bman) {}

  virtual void update(const Vector<Real> &x) {
    begin_ = 0;
  } 

  virtual int start(void) {
    return begin_;
  }

  virtual Real computeError(std::vector<Real> &vals) {
    return 0.0;
  }

  virtual Real computeError(std::vector<Teuchos::RCP<Vector<Real> > > &vals, const Vector<Real> &x) {
    return 0.0;
  }

  virtual void refine(void) {
    begin_ = numMySamples();
  }

  virtual void setSamples(bool inConstructor = false) {}

  virtual int numMySamples(void) {
    return weights_.size();
  }

  virtual std::vector<Real> getMyPoint(const int i) {
    return points_[i];
  }  

  virtual Real getMyWeight(const int i) {
    return weights_[i];
  }

  int batchID(void) {
    return bman_->batchID();
  }

  int numBatches(void) {
    return bman_->numBatches();
  }

  void sumAll(Real *input, Real *output, int dim) {
    bman_->sumAll(input, output, dim);
  }

  void sumAll(Vector<Real> &input, Vector<Real> &output) {
    bman_->sumAll(input,output);
  }

  void barrier(void) {
    bman_->barrier();
  }
};

}

#endif
