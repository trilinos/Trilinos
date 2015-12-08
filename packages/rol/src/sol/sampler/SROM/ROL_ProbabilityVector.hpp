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

#ifndef ROL_PROBABILITYVECTOR_H
#define ROL_PROBABILITYVECTOR_H

#include "ROL_StdVector.hpp"
#include "ROL_BatchManager.hpp"

/** \class ROL::ProbabilityVector
    \brief Provides the std::vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template <class Real>
class PrimalProbabilityVector;

template <class Real>
class DualProbabilityVector;

template <class Real>
class ProbabilityVector : public StdVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  uint numMySamples_;

public:
  ProbabilityVector(const Teuchos::RCP<std::vector<Real> > &vec)
    : StdVector<Real>(vec), numMySamples_(vec->size()) {}

  const Real getProbability(const int i) const {
    if ( i >= 0 && i < (int)numMySamples_ ) {
      const std::vector<Real> &yval = *(StdVector<Real>::getVector());
      return yval[i];
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> ERROR (ROL::ProbabilityVector): index out of bounds in getProbability!");
    }
    return 0;
  }

  void setProbability(const int i, const Real wt) {
    if ( i >= 0 && i < (int)numMySamples_ ) {
      std::vector<Real> &yval = *(StdVector<Real>::getVector());
      yval[i] = wt;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> ERROR (ROL::ProbabilityVector): index out of bounds in setProbability!");
    }
  }

  int getNumMyAtoms(void) const {
    return StdVector<Real>::dimension();
  }
};

template<class Real>
class PrimalProbabilityVector : public ProbabilityVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  uint numMySamples_;
  Teuchos::RCP<std::vector<Real> > scale_;
  Teuchos::RCP<BatchManager<Real> > bman_;

  mutable Teuchos::RCP<DualProbabilityVector<Real> > dual_vec_;

public:
  PrimalProbabilityVector(const Teuchos::RCP<std::vector<Real> > &vec,
                          const Teuchos::RCP<std::vector<Real> > &scale,
                          const Teuchos::RCP<BatchManager<Real> > &bman)
    : ProbabilityVector<Real>(vec), numMySamples_(vec->size()),
      scale_(scale), bman_(bman) {}

  Real dot(const Vector<Real> &x) const {
    const PrimalProbabilityVector<Real> &ex = Teuchos::dyn_cast<const PrimalProbabilityVector>(x);
    const std::vector<Real> &xval = *(ex.getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    Real val = 0;
    for (uint i = 0; i < numMySamples_; i++) {
      val += xval[i] * (*scale_)[i] * yval[i];
    }
    // Global sum
    Real sum_val = 0;
    bman_->sumAll(&val,&sum_val,1);
    return sum_val;
  }

  Teuchos::RCP<Vector<Real> > clone(void) const {
    return Teuchos::rcp(new PrimalProbabilityVector(
           Teuchos::rcp(new std::vector<Real>(numMySamples_)),scale_,bman_));
  }

  const Vector<Real> & dual(void) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    std::vector<Real> tmp(yval);
    for (uint i = 0; i < numMySamples_; i++) {
      tmp[i] *= (*scale_)[i];
    }
    dual_vec_ = Teuchos::rcp(new DualProbabilityVector<Real>(
                Teuchos::rcp(new std::vector<Real>(tmp)),scale_,bman_));
    return *dual_vec_;
  }

  int dimension(void) const {
    Real dim = (Real)StdVector<Real>::dimension();
    Real sum = 0.;
    bman_->sumAll(&dim,&sum,1);
    return (int)sum;
  }

  Real reduce(const Elementwise::ReductionOp<Real> &r) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    Real result = r.initialValue();
    for (uint i = 0; i < numMySamples_; i++) {
      r.reduce(yval[i],result);
    }
    // Global sum
    Real sum = 0.;
    bman_->reduceAll(&result,&sum,r);
    return sum;
  }
};

template<class Real>
class DualProbabilityVector : public ProbabilityVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  uint numMySamples_;
  Teuchos::RCP<std::vector<Real> > scale_;
  Teuchos::RCP<BatchManager<Real> > bman_;

  mutable Teuchos::RCP<PrimalProbabilityVector<Real> > dual_vec_;

public:
  DualProbabilityVector(const Teuchos::RCP<std::vector<Real> > &vec,
                        const Teuchos::RCP<std::vector<Real> > &scale,
                        const Teuchos::RCP<BatchManager<Real> > &bman)
    : ProbabilityVector<Real>(vec), numMySamples_(vec->size()),
      scale_(scale), bman_(bman) {}

  Real dot(const Vector<Real> &x) const {
    const DualProbabilityVector<Real> &ex = Teuchos::dyn_cast<const DualProbabilityVector>(x);
    const std::vector<Real> &xval = *(ex.getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    Real val = 0;
    for (uint i = 0; i < numMySamples_; i++) {
      val += xval[i] * yval[i] / (*scale_)[i];
    }
    // Global sum
    Real sum_val = 0;
    bman_->sumAll(&val,&sum_val,1);
    return sum_val;
  }

  Teuchos::RCP<Vector<Real> > clone(void) const {
    return Teuchos::rcp(new DualProbabilityVector(
           Teuchos::rcp(new std::vector<Real>(numMySamples_)),scale_,bman_));
  }

  const Vector<Real> & dual(void) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    std::vector<Real> tmp(yval);
    for (uint i = 0; i < numMySamples_; i++) {
      tmp[i] /= (*scale_)[i];
    }
    dual_vec_ = Teuchos::rcp(new PrimalProbabilityVector<Real>(
                Teuchos::rcp(new std::vector<Real>(tmp)),scale_,bman_));
    return *dual_vec_;
  }

  int dimension(void) const {
    Real dim = (Real)StdVector<Real>::dimension();
    Real sum = 0.;
    bman_->sumAll(&dim,&sum,1);
    return (int)sum;
  }

  Real reduce(const Elementwise::ReductionOp<Real> &r) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    Real result = r.initialValue();
    for (uint i = 0; i < numMySamples_; i++) {
      r.reduce(yval[i],result);
    }
    // Global sum
    Real sum = 0.;
    bman_->reduceAll(&result,&sum,r);
    return sum;
  }
};

} // namespace ROL

#endif
