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

#ifndef ROL_STDVECTOR_H
#define ROL_STDVECTOR_H

#include "ROL_Vector.hpp"

/** \class ROL::StdVector
    \brief Provides the std::vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template <class Real, class Element=Real>
class StdVector : public Vector<Real> {

private:

  Teuchos::RCP<std::vector<Element> >  std_vec_;

public:

  StdVector(const Teuchos::RCP<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

  void plus( const Vector<Real> &x ) {
    const StdVector &ex = Teuchos::dyn_cast<const StdVector>(x);
    int dimension  = (int)(this->std_vec_->size());
    for (int i=0; i<dimension; i++) {
      (*this->std_vec_)[i] += (*ex.getVector())[i];
    }
  }

  void scale( const Real alpha ) {
    int dimension = (int)(this->std_vec_->size());
    for (int i=0; i<dimension; i++) {
      (*this->std_vec_)[i] *= alpha;
    }
  }

  Real dot( const Vector<Real> &x ) const {
    const StdVector & ex = Teuchos::dyn_cast<const StdVector>(x);
    int dimension  = (int)(this->std_vec_->size());
    Real val = 0;

    for (int i=0; i<dimension; i++) {
      val += (*this->std_vec_)[i]*(*ex.getVector())[i];
    }
    return val;
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( this->dot(*this) );
    return val;
  }

  Teuchos::RCP<Vector<Real> > clone() const {
    return Teuchos::rcp( new StdVector( Teuchos::rcp(new std::vector<Element>(this->std_vec_->size())) ));
  }

  Teuchos::RCP<const std::vector<Element> > getVector() const {
    return this->std_vec_;
  }

  Teuchos::RCP<std::vector<Element> > getNonConstVector() {
    return this->std_vec_;
  }

  Teuchos::RCP<Vector<Real> > basis( const int i ) const {
    Teuchos::RCP<StdVector> e = Teuchos::rcp( new StdVector( Teuchos::rcp(new std::vector<Element>(this->std_vec_->size(), 0.0)) ));
    (*e->getNonConstVector())[i]= 1.0;
    return e;
  }

  int dimension() const {return this->std_vec_->size();}

}; // class StdVector

namespace StdVector_Helper {

template<class Real>
Teuchos::RCP<const std::vector<Real> > constDownCast( const Vector<Real> &x ) {
  return (Teuchos::dyn_cast<const StdVector<Real> >(x)).getVector();
}

template<class Real>
Teuchos::RCP<std::vector<Real> > downCast( Vector<Real> &x ) {
  return (Teuchos::dyn_cast<StdVector<Real> >(x)).getNonConstVector();
}

} // namespace StdVector_Helper

} // namespace ROL

#endif
