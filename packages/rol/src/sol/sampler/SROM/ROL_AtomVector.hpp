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

#ifndef ROL_ATOMVECTOR_H
#define ROL_ATOMVECTOR_H

#include "ROL_StdVector.hpp"
#include "ROL_BatchManager.hpp"

/** \class ROL::AtomVector
    \brief Provides the std::vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template <class Real>
class PrimalAtomVector;

template <class Real>
class DualAtomVector;

template <class Real>
class AtomVector : public StdVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  const uint numMySamples_;
  const uint dimension_;

public:
  AtomVector(const Teuchos::RCP<std::vector<Real> > &vec,
             const int numMySamples, const int dimension)
    : StdVector<Real>(vec), numMySamples_(numMySamples), dimension_(dimension) {}

  Teuchos::RCP<const std::vector<Real> > getAtom(const int i) const {
    std::vector<Real> pt(dimension_,0.);
    if ( i >= 0 && i < (int)numMySamples_ ) {
      const std::vector<Real> &yval = *(StdVector<Real>::getVector());
      for (uint j = 0; j < dimension_; j++) {
        pt[j] = yval[i*dimension_ + j];
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> ERROR (ROL::AtomVector): index out of bounds in getAtom!");
    }
    return Teuchos::rcp(new std::vector<Real>(pt));
  }

  void setAtom(const int i, const std::vector<Real> &pt) {
    if ( i >= 0 && i < (int)numMySamples_ ) {
      std::vector<Real> &yval = *(StdVector<Real>::getVector());
      for (uint j = 0; j < dimension_; j++) {
        yval[i*dimension_ + j] = pt[j];
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> ERROR (ROL::AtomVector): index out of bounds in setAtom!");
    }
  }
};

template<class Real>
class PrimalAtomVector : public AtomVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  const uint numMySamples_;
  const uint dimension_;
  const Teuchos::RCP<std::vector<Real> > scale_;
  const Teuchos::RCP<BatchManager<Real> > bman_;

  mutable Teuchos::RCP<DualAtomVector<Real> > dual_vec_;

public:
  PrimalAtomVector(const Teuchos::RCP<std::vector<Real> > &vec,
                   const int numMySamples,
                   const int dimension,
                   const Teuchos::RCP<std::vector<Real> > &scale,
                   const Teuchos::RCP<BatchManager<Real> > &bman)
    : AtomVector<Real>(vec,numMySamples,dimension),
      numMySamples_((uint)numMySamples), dimension_((uint)dimension),
      scale_(scale), bman_(bman) {}

  Real dot(const Vector<Real> &x) const {
    const PrimalAtomVector<Real> &ex = Teuchos::dyn_cast<const PrimalAtomVector>(x);
    const std::vector<Real> &xval = *(ex.getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint index = 0;
    Real val = 0;
    for (uint i = 0; i < numMySamples_; i++) {
      for (uint j = 0; j < dimension_; j++) {
        index = i*dimension_ + j;
        val += xval[index] * (*scale_)[index] * yval[index];
      }
    }
    // Global sum
    Real sum_val = 0;
    bman_->sumAll(&val,&sum_val,1);
    return sum_val;
  }

  Teuchos::RCP<Vector<Real> > clone(void) const {
    return Teuchos::rcp(new PrimalAtomVector(
           Teuchos::rcp(new std::vector<Real>(numMySamples_*dimension_)),
                        numMySamples_,dimension_,scale_,bman_));
  }

  const Vector<Real> & dual(void) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint index = 0;
    std::vector<Real> tmp(yval);
    for (uint i = 0; i < numMySamples_; i++) {
      for (uint j = 0; j < dimension_; j++) {
        index = i*dimension_ + j;
        tmp[index] *= (*scale_)[index];
      }
    }
    dual_vec_ = Teuchos::rcp(new DualAtomVector<Real>(
                Teuchos::rcp(new std::vector<Real>(tmp)),
                             numMySamples_,dimension_,scale_,bman_));
    return *dual_vec_;
  }

  Real reduce(const Elementwise::ReductionOp<Real> &r) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint index = 0;
    Real result = r.initialValue();
    for (uint i = 0; i < numMySamples_; i++) {
      for (uint j = 0; j < dimension_; j++) {
        index = i*dimension_ + j;
        r.reduce(yval[index],result);
      }
    }
    // Global sum
    Real sum = 0.;
    bman_->reduceAll(&result,&sum,r);
    return sum;
  }
};

template<class Real>
class DualAtomVector : public AtomVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  const uint numMySamples_;
  const uint dimension_;
  const Teuchos::RCP<std::vector<Real> > scale_;
  const Teuchos::RCP<BatchManager<Real> > bman_;

  mutable Teuchos::RCP<PrimalAtomVector<Real> > dual_vec_;

public:
  DualAtomVector(const Teuchos::RCP<std::vector<Real> > &vec,
                 const int numMySamples,
                 const int dimension,
                 const Teuchos::RCP<std::vector<Real> > &scale,
                 const Teuchos::RCP<BatchManager<Real> > &bman)
    : AtomVector<Real>(vec,numMySamples,dimension),
      numMySamples_((uint)numMySamples), dimension_((uint)dimension),
      scale_(scale), bman_(bman) {}

  Real dot(const Vector<Real> &x) const {
    const DualAtomVector<Real> &ex = Teuchos::dyn_cast<const DualAtomVector>(x);
    const std::vector<Real> &xval = *(ex.getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint index = 0;
    Real val = 0;
    for (uint i = 0; i < numMySamples_; i++) {
      for (uint j = 0; j < dimension_; j++) {
        index = i*dimension_ + j;
        val += xval[index] * yval[index] / (*scale_)[index];
      }
    }
    // Global sum
    Real sum_val = 0;
    bman_->sumAll(&val,&sum_val,1);
    return sum_val;
  }

  Teuchos::RCP<Vector<Real> > clone(void) const {
    return Teuchos::rcp(new DualAtomVector(
           Teuchos::rcp(new std::vector<Real>(numMySamples_*dimension_)),
                        numMySamples_,dimension_,scale_,bman_));
  }

  const Vector<Real> & dual(void) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint index = 0;
    std::vector<Real> tmp(yval);
    for (uint i = 0; i < numMySamples_; i++) {
      for (uint j = 0; j < dimension_; j++) {
        index = i*dimension_ + j;
        tmp[index] /= (*scale_)[index];
      }
    }
    dual_vec_ = Teuchos::rcp(new PrimalAtomVector<Real>(
                Teuchos::rcp(new std::vector<Real>(tmp)),
                             numMySamples_,dimension_,scale_,bman_));
    return *dual_vec_;
  }

  Real reduce(const Elementwise::ReductionOp<Real> &r) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint index = 0;
    Real result = r.initialValue();
    for (uint i = 0; i < numMySamples_; i++) {
      for (uint j = 0; j < dimension_; j++) {
        index = i*dimension_ + j;
        r.reduce(yval[index],result);
      }
    }
    // Global sum
    Real sum = 0.;
    bman_->reduceAll(&result,&sum,r);
    return sum;
  }
};

} // namespace ROL

#endif
