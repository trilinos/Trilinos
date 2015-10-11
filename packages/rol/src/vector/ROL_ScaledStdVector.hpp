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

#include <algorithm>
#include <cstdlib>

#include "ROL_Vector.hpp"

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
class PrimalScaledStdVector : public Vector<Real> {

  typedef typename std::vector<Element>::size_type uint;

private:

  Teuchos::RCP<std::vector<Element> >               std_vec_;
  Teuchos::RCP<std::vector<Element> >               scaling_vec_;
  mutable Teuchos::RCP<DualScaledStdVector<Real> >  dual_vec_;

public:

  PrimalScaledStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec,
                        const Teuchos::RCP<std::vector<Element> > & scaling_vec) :
    std_vec_(std_vec), scaling_vec_(scaling_vec) {}

  void set( const Vector<Real> &x ) {
    const PrimalScaledStdVector &ex = Teuchos::dyn_cast<const PrimalScaledStdVector>(x);
    const std::vector<Element>& xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),std_vec_->begin());   
  }

  void plus( const Vector<Real> &x ) {
    const PrimalScaledStdVector &ex = Teuchos::dyn_cast<const PrimalScaledStdVector>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] += xval[i];
    }
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const PrimalScaledStdVector &ex = Teuchos::dyn_cast<const PrimalScaledStdVector>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] += alpha*xval[i];
    }
  }

  void scale( const Real alpha ) {
    uint dimension = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] *= alpha;
    }
  }

  Real dot( const Vector<Real> &x ) const {
    const PrimalScaledStdVector & ex = Teuchos::dyn_cast<const PrimalScaledStdVector>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    Real val = 0;
    for (uint i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*xval[i]*(*scaling_vec_)[i];
    }
    return val;
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  Teuchos::RCP<Vector<Real> > clone() const {
    return Teuchos::rcp( new PrimalScaledStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size())), scaling_vec_ ));
  }

  Teuchos::RCP<const std::vector<Element> > getVector() const {
    return std_vec_;
  }

  Teuchos::RCP<std::vector<Element> > getVector() {
    return std_vec_;
  }

  const ROL::Vector<Real> & dual() const {
    std::vector<Element> tmp_vec(*std_vec_);
    tmp_vec[0] *= (*scaling_vec_)[0];
    tmp_vec[1] *= (*scaling_vec_)[1];
    dual_vec_ = Teuchos::rcp( new DualScaledStdVector<Real>( Teuchos::rcp(new std::vector<Element>(tmp_vec)), scaling_vec_ ) );
    return *dual_vec_;
  }

  Teuchos::RCP<Vector<Real> > basis( const int i ) const {
    Teuchos::RCP<PrimalScaledStdVector> e = Teuchos::rcp( new PrimalScaledStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size(), 0.0)), scaling_vec_ ));
    (*e->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return std_vec_->size();
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    uint dimension  = std_vec_->size();
    for(uint i=0; i<dimension; ++i) {
      (*std_vec_)[i] = f.apply((*std_vec_)[i]);
    }

  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const PrimalScaledStdVector & ex = Teuchos::dyn_cast<const PrimalScaledStdVector>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] = f.apply((*std_vec_)[i],xval[i]);
    }

  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    uint dimension  = std_vec_->size();
    for(uint i=0; i<dimension; ++i) {
      r.reduce((*std_vec_)[i],result);
    }
    return result;
  }

}; // class PrimalScaledStdVector



template <class Real, class Element>
class DualScaledStdVector : public Vector<Real> {

  typedef typename std::vector<Element>::size_type uint;

private:

  Teuchos::RCP<std::vector<Element> >                 std_vec_;
  Teuchos::RCP<std::vector<Element> >                 scaling_vec_;
  mutable Teuchos::RCP<PrimalScaledStdVector<Real> >  dual_vec_;

public:

  DualScaledStdVector(const Teuchos::RCP<std::vector<Element> > & std_vec,
                      const Teuchos::RCP<std::vector<Element> > & scaling_vec) :
    std_vec_(std_vec), scaling_vec_(scaling_vec) {}

  void set( const Vector<Real> &x ) {
    const DualScaledStdVector &ex = Teuchos::dyn_cast<const DualScaledStdVector>(x);
    const std::vector<Element>& xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),std_vec_->begin());   
  }

  void plus( const Vector<Real> &x ) {
    const DualScaledStdVector &ex = Teuchos::dyn_cast<const DualScaledStdVector>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] += xval[i];
    }
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    const DualScaledStdVector &ex = Teuchos::dyn_cast<const DualScaledStdVector>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] += alpha*xval[i];
    }
  }

  void scale( const Real alpha ) {
    uint dimension = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] *= alpha;
    }
  }

  Real dot( const Vector<Real> &x ) const {
    const DualScaledStdVector & ex = Teuchos::dyn_cast<const DualScaledStdVector>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    Real val = 0;
    for (uint i=0; i<dimension; i++) {
      val += (*std_vec_)[i]*xval[i]/(*scaling_vec_)[i];
    }
    return val;
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  Teuchos::RCP<Vector<Real> > clone() const {
    return Teuchos::rcp( new DualScaledStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size())), scaling_vec_ ));
  }

  Teuchos::RCP<const std::vector<Element> > getVector() const {
    return std_vec_;
  }

  Teuchos::RCP<std::vector<Element> > getVector() {
    return std_vec_;
  }

  const ROL::Vector<Real> & dual() const {
    std::vector<Element> tmp_vec(*std_vec_);
    tmp_vec[0] /= (*scaling_vec_)[0];
    tmp_vec[1] /= (*scaling_vec_)[1];
    dual_vec_ = Teuchos::rcp( new PrimalScaledStdVector<Real>( Teuchos::rcp(new std::vector<Element>(tmp_vec)), scaling_vec_ ) );
    return *dual_vec_;
  }

  Teuchos::RCP<Vector<Real> > basis( const int i ) const {
    Teuchos::RCP<DualScaledStdVector> e = Teuchos::rcp( new DualScaledStdVector( Teuchos::rcp(new std::vector<Element>(std_vec_->size(), 0.0)), scaling_vec_ ));
    (*e->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return std_vec_->size();
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    uint dimension  = std_vec_->size();
    for(uint i=0; i<dimension; ++i) {
      (*std_vec_)[i] = f.apply((*std_vec_)[i]);
    }

  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {
    const DualScaledStdVector & ex = Teuchos::dyn_cast<const DualScaledStdVector>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dimension  = std_vec_->size();
    for (uint i=0; i<dimension; i++) {
      (*std_vec_)[i] = f.apply((*std_vec_)[i],xval[i]);
    }

  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    uint dimension  = std_vec_->size();
    for(uint i=0; i<dimension; ++i) {
      r.reduce((*std_vec_)[i],result);
    }
    return result;
  }

}; // class DualScaledStdVector


} // namespace ROL

#endif
