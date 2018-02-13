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
  const ROL::Ptr<BatchManager<Real> > bman_;

protected:
  const ROL::Ptr<BatchManager<Real> > getBatchManager(void) const {
    return bman_;
  }

public:
  BatchStdVector(const ROL::Ptr<std::vector<Real> > &vec,
                 const ROL::Ptr<BatchManager<Real> > &bman)
   : StdVector<Real>(vec), bman_(bman) {}
   
  virtual Real dot(const Vector<Real> &x) const {
    const std::vector<Real> &xval = *(dynamic_cast<const StdVector<Real>&>(x).getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint numMySamples = yval.size();
    TEUCHOS_TEST_FOR_EXCEPTION( xval.size() != numMySamples, std::invalid_argument,
      "Error: Vectors must have the same dimension." );
    Real val(0), sum_val(0);
    for (uint i = 0; i < numMySamples; ++i) {
      val += xval[i] * yval[i];
    }
    // Global sum
    bman_->sumAll(&val,&sum_val,1);
    return sum_val;
  }

  virtual ROL::Ptr<Vector<Real> > clone(void) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint numMySamples = yval.size();
    return ROL::makePtr<BatchStdVector>(
           ROL::makePtr<std::vector<Real>>(numMySamples),bman_);
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
};

} // namespace ROL

#endif
