// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PROBABILITYVECTOR_H
#define ROL_PROBABILITYVECTOR_H

#include "ROL_BatchStdVector.hpp"

/** \class ROL::ProbabilityVector
    \brief Provides the std::vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template <class Real>
class PrimalProbabilityVector;

template <class Real>
class DualProbabilityVector;

template <class Real>
class ProbabilityVector : public BatchStdVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
public:
  ProbabilityVector(const ROL::Ptr<std::vector<Real>> &vec,
                    const ROL::Ptr<BatchManager<Real>> &bman)
   : BatchStdVector<Real>(vec,bman) {}

  const Real getProbability(const int i) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    int numMySamples = static_cast<int>(yval.size());
    ROL_TEST_FOR_EXCEPTION((i < 0 || i >= numMySamples), std::invalid_argument,
      ">>> ERROR (ROL::ProbabilityVector): index out of bounds in getProbability!");
    return yval[i];
  }

  void setProbability(const int i, const Real wt) {
    std::vector<Real> &yval = *(StdVector<Real>::getVector());
    int numMySamples = static_cast<int>(yval.size());
    ROL_TEST_FOR_EXCEPTION((i < 0 || i >= numMySamples), std::invalid_argument,
      ">>> ERROR (ROL::ProbabilityVector): index out of bounds in setProbability!");
    yval[i] = wt;
  }

  int getNumMyAtoms(void) const {
    int numMySamples = static_cast<int>(StdVector<Real>::getVector()->size());
    return numMySamples;
  }

  ROL::Ptr<Vector<Real> > clone(void) const {
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint numMySamples = yval.size();
    return ROL::makePtr<ProbabilityVector>(
           ROL::makePtr<std::vector<Real>>(numMySamples),BatchStdVector<Real>::getBatchManager());
  }
};

template<class Real>
class PrimalProbabilityVector : public ProbabilityVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  ROL::Ptr<std::vector<Real> > scale_;
  mutable ROL::Ptr<DualProbabilityVector<Real> > dual_vec_;
  mutable bool isDualInitialized_;

public:
  PrimalProbabilityVector(const ROL::Ptr<std::vector<Real> > &vec,
                          const ROL::Ptr<BatchManager<Real> > &bman,
                          const ROL::Ptr<std::vector<Real> > &scale)
    : ProbabilityVector<Real>(vec,bman), scale_(scale),
      isDualInitialized_(false) {}

  Real dot(const Vector<Real> &x) const {
    const std::vector<Real> &xval = *(dynamic_cast<const StdVector<Real>&>(x).getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint numMySamples = static_cast<uint>(yval.size());
    ROL_TEST_FOR_EXCEPTION( xval.size() != numMySamples, std::invalid_argument,
      "Error: Vectors must have the same dimension." );
    Real val(0), sum_val(0);
    for (uint i = 0; i < numMySamples; i++) {
      val += xval[i] * (*scale_)[i] * yval[i];
    }
    // Global sum
    BatchStdVector<Real>::getBatchManager()->sumAll(&val,&sum_val,1);
    return sum_val;
  }

  ROL::Ptr<Vector<Real> > clone(void) const {
    uint numMySamples = static_cast<uint>(StdVector<Real>::getVector()->size());
    return ROL::makePtr<PrimalProbabilityVector>(
           ROL::makePtr<std::vector<Real>>(numMySamples),
           BatchStdVector<Real>::getBatchManager(),scale_);
  }

  const Vector<Real> & dual(void) const {
    uint numMySamples = static_cast<uint>(StdVector<Real>::getVector()->size());
    if ( !isDualInitialized_ ) {
      dual_vec_ = ROL::makePtr<DualProbabilityVector<Real>>(
                  ROL::makePtr<std::vector<Real>>(numMySamples),
                  BatchStdVector<Real>::getBatchManager(),scale_);
      isDualInitialized_ = true;
    }
    for (uint i = 0; i < numMySamples; ++i) {
      (*(dual_vec_->getVector()))[i]
        = (*scale_)[i]*(*StdVector<Real>::getVector())[i];
    }
    return *dual_vec_;
  }
};

template<class Real>
class DualProbabilityVector : public ProbabilityVector<Real> {
  typedef typename std::vector<Real>::size_type uint;
private:
  ROL::Ptr<std::vector<Real> > scale_;
  mutable ROL::Ptr<PrimalProbabilityVector<Real> > primal_vec_;
  mutable bool isDualInitialized_;

public:
  DualProbabilityVector(const ROL::Ptr<std::vector<Real> > &vec,
                        const ROL::Ptr<BatchManager<Real> > &bman,
                        const ROL::Ptr<std::vector<Real> > &scale)
    : ProbabilityVector<Real>(vec,bman), scale_(scale),
      isDualInitialized_(false) {}

  Real dot(const Vector<Real> &x) const {
    const std::vector<Real> &xval = *(dynamic_cast<const StdVector<Real>&>(x).getVector());
    const std::vector<Real> &yval = *(StdVector<Real>::getVector());
    uint numMySamples = static_cast<uint>(yval.size());
    ROL_TEST_FOR_EXCEPTION( xval.size() != numMySamples, std::invalid_argument,
      "Error: Vectors must have the same dimension." );
    Real val(0), sum_val(0);
    for (uint i = 0; i < numMySamples; i++) {
      val += xval[i] * yval[i] / (*scale_)[i];
    }
    // Global sum
    BatchStdVector<Real>::getBatchManager()->sumAll(&val,&sum_val,1);
    return sum_val;
  }

  ROL::Ptr<Vector<Real> > clone(void) const {
    uint numMySamples = static_cast<uint>(StdVector<Real>::getVector()->size());
    return ROL::makePtr<DualProbabilityVector>(
           ROL::makePtr<std::vector<Real>>(numMySamples),
           BatchStdVector<Real>::getBatchManager(),scale_);
  }

  const Vector<Real> & dual(void) const {
    uint numMySamples = static_cast<uint>(StdVector<Real>::getVector()->size());
    if ( !isDualInitialized_ ) {
      primal_vec_ = ROL::makePtr<PrimalProbabilityVector<Real>>(
                    ROL::makePtr<std::vector<Real>>(numMySamples),
                    BatchStdVector<Real>::getBatchManager(),scale_);
      isDualInitialized_ = true;
    }
    for (uint i = 0; i < numMySamples; i++) {
      (*(primal_vec_->getVector()))[i]
        = (*StdVector<Real>::getVector())[i]/(*scale_)[i];
    }
    return *primal_vec_;
  }
};

} // namespace ROL

#endif
