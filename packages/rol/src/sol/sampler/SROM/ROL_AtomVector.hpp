// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  AtomVector(const ROL::Ptr<std::vector<Real> > &vec,
             const ROL::Ptr<BatchManager<Real> > &bman,
             const int numMySamples, const int dimension)
    : BatchStdVector<Real>(vec,bman),
      numMySamples_(numMySamples), dimension_(dimension) {}

  ROL::Ptr<const std::vector<Real> > getAtom(const int i) const {
    ROL_TEST_FOR_EXCEPTION((i < 0 || i > numMySamples_), std::invalid_argument,
      ">>> ERROR (ROL::AtomVector): index out of bounds in getAtom!");
    uint dim = static_cast<uint>(dimension_), I = static_cast<uint>(i);
    std::vector<Real> pt(dim,0);
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    for (uint j = 0; j < dim; ++j) {
      pt[j] = yval[I*dim + j];
    }
    return ROL::makePtr<std::vector<Real>>(pt);
  }

  void setAtom(const int i, const std::vector<Real> &pt) {
    ROL_TEST_FOR_EXCEPTION((i < 0 || i > numMySamples_), std::invalid_argument,
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
  const ROL::Ptr<std::vector<Real> > scale_;
  mutable ROL::Ptr<DualAtomVector<Real> > dual_vec_;
  mutable bool isDualInitialized_;

public:
  PrimalAtomVector(const ROL::Ptr<std::vector<Real> > &vec,
                   const ROL::Ptr<BatchManager<Real> > &bman,
                   const int numMySamples, const int dimension,
                   const ROL::Ptr<std::vector<Real> > &scale)
    : AtomVector<Real>(vec,bman,numMySamples,dimension),
      scale_(scale), isDualInitialized_(false) {}

  Real dot(const Vector<Real> &x) const {
    const std::vector<Real> &xval = *(dynamic_cast<const StdVector<Real>&>(x).getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint ysize = yval.size();
    ROL_TEST_FOR_EXCEPTION( xval.size() != ysize, std::invalid_argument,
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

  ROL::Ptr<Vector<Real> > clone(void) const {
    uint numMySamples = static_cast<uint>(AtomVector<Real>::getNumMyAtoms());
    uint dimension    = static_cast<uint>(AtomVector<Real>::getDimension());
    return ROL::makePtr<PrimalAtomVector>(
           ROL::makePtr<std::vector<Real>>(numMySamples*dimension),
                        BatchStdVector<Real>::getBatchManager(),
                        numMySamples,dimension,scale_);
  }

  const Vector<Real> & dual(void) const {
    uint numMySamples = static_cast<uint>(AtomVector<Real>::getNumMyAtoms());
    uint dimension    = static_cast<uint>(AtomVector<Real>::getDimension());
    if ( !isDualInitialized_ ) {
      dual_vec_ = ROL::makePtr<DualAtomVector<Real>>(
                  ROL::makePtr<std::vector<Real>>(numMySamples*dimension),
                               BatchStdVector<Real>::getBatchManager(),
                               numMySamples,dimension,scale_);
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
  const ROL::Ptr<std::vector<Real> > scale_;
  mutable ROL::Ptr<PrimalAtomVector<Real> > primal_vec_;
  mutable bool isDualInitialized_;

public:
  DualAtomVector(const ROL::Ptr<std::vector<Real> > &vec,
                 const ROL::Ptr<BatchManager<Real> > &bman,
                 const int numMySamples, const int dimension,
                 const ROL::Ptr<std::vector<Real> > &scale)
    : AtomVector<Real>(vec,bman,numMySamples,dimension),
      scale_(scale), isDualInitialized_(false) {}

  Real dot(const Vector<Real> &x) const {
    const std::vector<Real> &xval = *(dynamic_cast<const StdVector<Real>&>(x).getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint ysize = yval.size();
    ROL_TEST_FOR_EXCEPTION( xval.size() != ysize, std::invalid_argument,
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

  ROL::Ptr<Vector<Real> > clone(void) const {
    uint numMySamples = static_cast<uint>(AtomVector<Real>::getNumMyAtoms());
    uint dimension    = static_cast<uint>(AtomVector<Real>::getDimension());
    return ROL::makePtr<DualAtomVector>(
           ROL::makePtr<std::vector<Real>>(numMySamples*dimension),
                        BatchStdVector<Real>::getBatchManager(),
                        numMySamples,dimension,scale_);
  }

  const Vector<Real> & dual(void) const {
    uint numMySamples = static_cast<uint>(AtomVector<Real>::getNumMyAtoms());
    uint dimension    = static_cast<uint>(AtomVector<Real>::getDimension());
    if ( !isDualInitialized_ ) {
      primal_vec_ = ROL::makePtr<PrimalAtomVector<Real>>(
                    ROL::makePtr<std::vector<Real>>(numMySamples*dimension),
                               BatchStdVector<Real>::getBatchManager(),
                               numMySamples,dimension,scale_);
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
