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

  ROL::Ptr<std::vector<Element> >               scaling_vec_;
  mutable ROL::Ptr<DualScaledStdVector<Real> >  dual_vec_;
  mutable bool isDualInitialized_;

public:

  PrimalScaledStdVector(const ROL::Ptr<std::vector<Element> > & std_vec,
                        const ROL::Ptr<std::vector<Element> > & scaling_vec) :
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

  ROL::Ptr<Vector<Real> > clone() const {
    uint dimension = (StdVector<Real>::getVector())->size();
    return ROL::makePtr<PrimalScaledStdVector>(
           ROL::makePtr<std::vector<Element>>(dimension), scaling_vec_ );
  }

  const ROL::Vector<Real> & dual() const {
    uint n = StdVector<Real>::getVector()->size();
    if ( !isDualInitialized_ ) {
      dual_vec_ = ROL::makePtr<DualScaledStdVector<Real>>(
                  ROL::makePtr<std::vector<Element>>(n),
                  scaling_vec_);
      isDualInitialized_ = true;
    }
    for (uint i = 0; i < n; i++) {
      (*(dual_vec_->getVector()))[i]
        = (*scaling_vec_)[i]*(*StdVector<Real>::getVector())[i];
    }
    return *dual_vec_;
  }

}; // class PrimalScaledStdVector



template <class Real, class Element>
class DualScaledStdVector : public StdVector<Real> {

  typedef typename std::vector<Element>::size_type uint;

private:

  ROL::Ptr<std::vector<Element> >                 scaling_vec_;
  mutable ROL::Ptr<PrimalScaledStdVector<Real> >  primal_vec_;
  mutable bool isDualInitialized_;

public:

  DualScaledStdVector(const ROL::Ptr<std::vector<Element> > & std_vec,
                      const ROL::Ptr<std::vector<Element> > & scaling_vec) :
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

  ROL::Ptr<Vector<Real> > clone() const {
    uint dimension = (StdVector<Real>::getVector())->size();
    return ROL::makePtr<DualScaledStdVector>(
           ROL::makePtr<std::vector<Element>>(dimension), scaling_vec_ );
  }

  const ROL::Vector<Real> & dual() const {
    uint n = StdVector<Real>::getVector()->size();
    if ( !isDualInitialized_ ) {
      primal_vec_ = ROL::makePtr<PrimalScaledStdVector<Real>>(
                    ROL::makePtr<std::vector<Element>>(n),
                    scaling_vec_);
      isDualInitialized_ = true;
    }
    for (uint i = 0; i < n; i++) {
      (*(primal_vec_->getVector()))[i]
        = (*StdVector<Real>::getVector())[i]/(*scaling_vec_)[i];
    }
    return *primal_vec_;
  }

}; // class DualScaledStdVector

} // namespace ROL

#endif
