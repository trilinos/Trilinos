// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SCALEDSTDVECTOR_H
#define ROL_SCALEDSTDVECTOR_H

#include "ROL_StdVector.hpp"

/** \class ROL::PrimalScaledStdVector
    \brief Provides the std::vector implementation of the ROL::Vector interface
           that handles scalings in the inner product.  Also see ROL::DualScaledStdVector.
*/

/** \class ROL::DualScaledStdVector
    \brief Provides the std::vector implementation of the ROL::Vector interface
           that handles scalings in the inner product.  Also see ROL::PrimalScaledStdVector.
*/

namespace ROL {

template <class Real, class Element=Real>
class PrimalScaledStdVector;

template <class Real, class Element=Real>
class DualScaledStdVector;

template <class Real, class Element>
class PrimalScaledStdVector : public StdVector<Real> {

  typedef typename std::vector<Element>::size_type uint;

private:

  Ptr<std::vector<Element> >               scaling_vec_;
  mutable Ptr<DualScaledStdVector<Real> >  dual_vec_;
  mutable bool isDualInitialized_;

public:

  PrimalScaledStdVector(const Ptr<std::vector<Element> > & std_vec,
                        const Ptr<std::vector<Element> > & scaling_vec) :
    StdVector<Real>(std_vec), scaling_vec_(scaling_vec),
    isDualInitialized_(false) {}

  Real dot( const Vector<Real> &x ) const {
    const PrimalScaledStdVector & ex = dynamic_cast<const PrimalScaledStdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    const std::vector<Element>& yval = *(StdVector<Real>::getVector());
    uint dimension = yval.size();
    Real val = 0;
    for (uint i=0; i<dimension; i++) {
      val += yval[i]*xval[i]*(*scaling_vec_)[i];
    }
    return val;
  }

  Ptr<Vector<Real> > clone() const {
    uint dimension = (StdVector<Real>::getVector())->size();
    return makePtr<PrimalScaledStdVector>(
           makePtr<std::vector<Element>>(dimension), scaling_vec_ );
  }

  const Vector<Real> & dual() const {
    uint n = StdVector<Real>::getVector()->size();
    if ( !isDualInitialized_ ) {
      dual_vec_ = makePtr<DualScaledStdVector<Real>>(
                  makePtr<std::vector<Element>>(n),
                  scaling_vec_);
      isDualInitialized_ = true;
    }
    for (uint i = 0; i < n; i++) {
      (*(dual_vec_->getVector()))[i]
        = (*scaling_vec_)[i]*(*StdVector<Real>::getVector())[i];
    }
    return *dual_vec_;
  }

  Real apply( const Vector<Real> &x ) const {
    const DualScaledStdVector<Real> & ex = dynamic_cast<const DualScaledStdVector<Real>&>(x);
    return StdVector<Real>::dot(ex);
    //const std::vector<Element>& xval = *ex.getVector();
    //const std::vector<Element>& yval = *(StdVector<Real>::getVector());
    //uint dimension = yval.size();
    //Real val = 0;
    //for (uint i=0; i<dimension; i++) {
    //  val += yval[i]*xval[i];
    //}
    //return val;
  }

}; // class PrimalScaledStdVector



template <class Real, class Element>
class DualScaledStdVector : public StdVector<Real> {

  typedef typename std::vector<Element>::size_type uint;

private:

  Ptr<std::vector<Element> >                 scaling_vec_;
  mutable Ptr<PrimalScaledStdVector<Real> >  primal_vec_;
  mutable bool isDualInitialized_;

public:

  DualScaledStdVector(const Ptr<std::vector<Element> > & std_vec,
                      const Ptr<std::vector<Element> > & scaling_vec) :
    StdVector<Real>(std_vec), scaling_vec_(scaling_vec),
    isDualInitialized_(false) {}

  Real dot( const Vector<Real> &x ) const {
    const DualScaledStdVector & ex = dynamic_cast<const DualScaledStdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    const std::vector<Element>& yval = *(StdVector<Real>::getVector());
    uint dimension = yval.size();
    Real val = 0;
    for (uint i=0; i<dimension; i++) {
      val += yval[i]*xval[i]/(*scaling_vec_)[i];
    }
    return val;
  }

  Ptr<Vector<Real> > clone() const {
    uint dimension = (StdVector<Real>::getVector())->size();
    return makePtr<DualScaledStdVector>(
           makePtr<std::vector<Element>>(dimension), scaling_vec_ );
  }

  const Vector<Real> & dual() const {
    uint n = StdVector<Real>::getVector()->size();
    if ( !isDualInitialized_ ) {
      primal_vec_ = makePtr<PrimalScaledStdVector<Real>>(
                    makePtr<std::vector<Element>>(n),
                    scaling_vec_);
      isDualInitialized_ = true;
    }
    for (uint i = 0; i < n; i++) {
      (*(primal_vec_->getVector()))[i]
        = (*StdVector<Real>::getVector())[i]/(*scaling_vec_)[i];
    }
    return *primal_vec_;
  }

  Real apply( const Vector<Real> &x ) const {
    const PrimalScaledStdVector<Real> & ex = dynamic_cast<const PrimalScaledStdVector<Real>&>(x);
    return StdVector<Real>::dot(ex);
//    const std::vector<Element>& xval = *ex.getVector();
//    const std::vector<Element>& yval = *(StdVector<Real>::getVector());
//    uint dimension = yval.size();
//    Real val = 0;
//    for (uint i=0; i<dimension; i++) {
//      val += yval[i]*xval[i];
//    }
//    return val;
  }

}; // class DualScaledStdVector

} // namespace ROL

#endif
