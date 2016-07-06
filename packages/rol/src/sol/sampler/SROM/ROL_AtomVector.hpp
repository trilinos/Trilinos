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

#include "ROL_BatchStdVector.hpp"

/** \class ROL::AtomVector
    \brief Provides the std::vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template <class Real>
class PrimalAtomVector;

template <class Real>
class DualAtomVector;

template <class Real>
class AtomVector : public BatchStdVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  const int numMySamples_;
  const int dimension_;

public:
  AtomVector(const Teuchos::RCP<std::vector<Real> > &vec,
             const Teuchos::RCP<BatchManager<Real> > &bman,
             const int numMySamples, const int dimension)
    : BatchStdVector<Real>(vec,bman),
      numMySamples_(numMySamples), dimension_(dimension) {}

  Teuchos::RCP<const std::vector<Real> > getAtom(const int i) const {
    TEUCHOS_TEST_FOR_EXCEPTION((i < 0 || i > numMySamples_), std::invalid_argument,
      ">>> ERROR (ROL::AtomVector): index out of bounds in getAtom!");
    uint dim = static_cast<uint>(dimension_), I = static_cast<uint>(i);
    std::vector<Real> pt(dim,0);
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    for (uint j = 0; j < dim; ++j) {
      pt[j] = yval[I*dim + j];
    }
    return Teuchos::rcp(new std::vector<Real>(pt));
  }

  void setAtom(const int i, const std::vector<Real> &pt) {
    TEUCHOS_TEST_FOR_EXCEPTION((i < 0 || i > numMySamples_), std::invalid_argument,
      ">>> ERROR (ROL::AtomVector): index out of bounds in setAtom!");
    uint dim = static_cast<uint>(dimension_), I = static_cast<uint>(i);
    std::vector<Real> &yval = *(StdVector<Real>::getVector());
    for (uint j = 0; j < dim; ++j) {
      yval[I*dim + j] = pt[j];
    }
  }

  int getNumMyAtoms(void) const {
    return numMySamples_;
  }

  int getDimension(void) const {
    return dimension_;
  }
};

template<class Real>
class PrimalAtomVector : public AtomVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  const Teuchos::RCP<std::vector<Real> > scale_;
  mutable Teuchos::RCP<DualAtomVector<Real> > dual_vec_;
  mutable bool isDualInitialized_;

public:
  PrimalAtomVector(const Teuchos::RCP<std::vector<Real> > &vec,
                   const Teuchos::RCP<BatchManager<Real> > &bman,
                   const int numMySamples, const int dimension,
                   const Teuchos::RCP<std::vector<Real> > &scale)
    : AtomVector<Real>(vec,bman,numMySamples,dimension),
      scale_(scale), isDualInitialized_(false) {}

  Real dot(const Vector<Real> &x) const {
    const std::vector<Real> &xval = *(Teuchos::dyn_cast<const StdVector<Real> >(x).getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint ysize = yval.size();
    TEUCHOS_TEST_FOR_EXCEPTION( xval.size() != ysize, std::invalid_argument,
      "Error: Vectors must have the same dimension." );
    uint index        = 0;
    uint numMySamples = static_cast<uint>(AtomVector<Real>::getNumMyAtoms());
    uint dimension    = static_cast<uint>(AtomVector<Real>::getDimension());
    Real val(0), sum_val(0);
    for (uint i = 0; i < numMySamples; i++) {
      for (uint j = 0; j < dimension; j++) {
        index = i*dimension + j;
        val += xval[index] * (*scale_)[index] * yval[index];
      }
    }
    // Global sum
    BatchStdVector<Real>::getBatchManager()->sumAll(&val,&sum_val,1);
    return sum_val;
  }

  Teuchos::RCP<Vector<Real> > clone(void) const {
    uint numMySamples = static_cast<uint>(AtomVector<Real>::getNumMyAtoms());
    uint dimension    = static_cast<uint>(AtomVector<Real>::getDimension());
    return Teuchos::rcp(new PrimalAtomVector(
           Teuchos::rcp(new std::vector<Real>(numMySamples*dimension)),
                        BatchStdVector<Real>::getBatchManager(),
                        numMySamples,dimension,scale_));
  }

  const Vector<Real> & dual(void) const {
    uint numMySamples = static_cast<uint>(AtomVector<Real>::getNumMyAtoms());
    uint dimension    = static_cast<uint>(AtomVector<Real>::getDimension());
    if ( !isDualInitialized_ ) {
      dual_vec_ = Teuchos::rcp(new DualAtomVector<Real>(
                  Teuchos::rcp(new std::vector<Real>(numMySamples*dimension)),
                               BatchStdVector<Real>::getBatchManager(),
                               numMySamples,dimension,scale_));
      isDualInitialized_ = true;
    }
    uint index = 0;
    for (uint i = 0; i < numMySamples; i++) {
      for (uint j = 0; j < dimension; j++) {
        index = i*dimension + j;
        (*(dual_vec_->getVector()))[index]
          = (*scale_)[index] * (*(StdVector<Real>::getVector()))[index];
      }
    }
    return *dual_vec_;
  }
};

template<class Real>
class DualAtomVector : public AtomVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  const Teuchos::RCP<std::vector<Real> > scale_;
  mutable Teuchos::RCP<PrimalAtomVector<Real> > primal_vec_;
  mutable bool isDualInitialized_;

public:
  DualAtomVector(const Teuchos::RCP<std::vector<Real> > &vec,
                 const Teuchos::RCP<BatchManager<Real> > &bman,
                 const int numMySamples, const int dimension,
                 const Teuchos::RCP<std::vector<Real> > &scale)
    : AtomVector<Real>(vec,bman,numMySamples,dimension),
      scale_(scale), isDualInitialized_(false) {}

  Real dot(const Vector<Real> &x) const {
    const std::vector<Real> &xval = *(Teuchos::dyn_cast<const StdVector<Real> >(x).getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint ysize = yval.size();
    TEUCHOS_TEST_FOR_EXCEPTION( xval.size() != ysize, std::invalid_argument,
      "Error: Vectors must have the same dimension." );
    uint index        = 0;
    uint numMySamples = static_cast<uint>(AtomVector<Real>::getNumMyAtoms());
    uint dimension    = static_cast<uint>(AtomVector<Real>::getDimension());
    Real val(0), sum_val(0);
    for (uint i = 0; i < numMySamples; i++) {
      for (uint j = 0; j < dimension; j++) {
        index = i*dimension + j;
        val += xval[index] * yval[index] / (*scale_)[index];
      }
    }
    // Global sum
    BatchStdVector<Real>::getBatchManager()->sumAll(&val,&sum_val,1);
    return sum_val;
  }

  Teuchos::RCP<Vector<Real> > clone(void) const {
    uint numMySamples = static_cast<uint>(AtomVector<Real>::getNumMyAtoms());
    uint dimension    = static_cast<uint>(AtomVector<Real>::getDimension());
    return Teuchos::rcp(new DualAtomVector(
           Teuchos::rcp(new std::vector<Real>(numMySamples*dimension)),
                        BatchStdVector<Real>::getBatchManager(),
                        numMySamples,dimension,scale_));
  }

  const Vector<Real> & dual(void) const {
    uint numMySamples = static_cast<uint>(AtomVector<Real>::getNumMyAtoms());
    uint dimension    = static_cast<uint>(AtomVector<Real>::getDimension());
    if ( !isDualInitialized_ ) {
      primal_vec_ = Teuchos::rcp(new PrimalAtomVector<Real>(
                    Teuchos::rcp(new std::vector<Real>(numMySamples*dimension)),
                               BatchStdVector<Real>::getBatchManager(),
                               numMySamples,dimension,scale_));
      isDualInitialized_ = true;
    }
    uint index = 0;
    for (uint i = 0; i < numMySamples; i++) {
      for (uint j = 0; j < dimension; j++) {
        index = i*dimension + j;
        (*(primal_vec_->getVector()))[index]
          = (*(StdVector<Real>::getVector()))[index] / (*scale_)[index];
      }
    }
    return *primal_vec_;
  }
};

} // namespace ROL

#endif
