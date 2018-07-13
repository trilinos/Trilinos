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

#include <algorithm>
#include <cstdlib>

#include "ROL_Vector.hpp"

/** \class ROL::StdVector
    \brief Provides the std::vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template <class Real, class Element=Real>
class StdVector : public Vector<Real> {

  typedef typename std::vector<Real>::size_type uint;

private:

  Ptr<std::vector<Element> >  std_vec_;

public:

  StdVector(const Ptr<std::vector<Element> > & std_vec) : std_vec_(std_vec) {}

  void set( const Vector<Real> &x ) {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

    const StdVector &ex = dynamic_cast<const StdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),std_vec_->begin());
  }

  void plus( const Vector<Real> &x ) {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

    const StdVector &ex = dynamic_cast<const StdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dim  = std_vec_->size();
    for (uint i=0; i<dim; i++) {
      (*std_vec_)[i] += xval[i];
    }
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

    const StdVector &ex = dynamic_cast<const StdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dim  = std_vec_->size();
    for (uint i=0; i<dim; i++) {
      (*std_vec_)[i] += alpha*xval[i];
    }
  }

  void scale( const Real alpha ) {
    uint dim = std_vec_->size();
    for (uint i=0; i<dim; i++) {
      (*std_vec_)[i] *= alpha;
    }
  }

  virtual Real dot( const Vector<Real> &x ) const {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

    const StdVector & ex = dynamic_cast<const StdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dim  = std_vec_->size();
    Real val = 0;
    for (uint i=0; i<dim; i++) {
      val += (*std_vec_)[i]*xval[i];
    }
    return val;
  }

  Real norm() const {
    Real val = 0;
    val = std::sqrt( dot(*this) );
    return val;
  }

  virtual Ptr<Vector<Real> > clone() const {
    return makePtr<StdVector>( makePtr<std::vector<Element>>(std_vec_->size(), static_cast<Element>(0)));
  }

  Ptr<const std::vector<Element> > getVector() const {
    return std_vec_;
  }

  Ptr<std::vector<Element> > getVector() {
    return std_vec_;
  }

  Ptr<Vector<Real> > basis( const int i ) const {

    ROL_TEST_FOR_EXCEPTION( i >= dimension() || i<0,
                                std::invalid_argument,
                                "Error: Basis index must be between 0 and vector dimension." );

    Ptr<Vector<Real>> e = clone();
    (*staticPtrCast<StdVector>(e)->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return static_cast<int>(std_vec_->size());
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
    uint dim  = std_vec_->size();
    for(uint i=0; i<dim; ++i) {
      (*std_vec_)[i] = f.apply((*std_vec_)[i]);
    }
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

    const StdVector & ex = dynamic_cast<const StdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    uint dim  = std_vec_->size();
    for (uint i=0; i<dim; i++) {
      (*std_vec_)[i] = f.apply((*std_vec_)[i],xval[i]);
    }

  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    uint dim  = std_vec_->size();
    for(uint i=0; i<dim; ++i) {
      r.reduce((*std_vec_)[i],result);
    }
    return result;
  }

  void setScalar( const Real C ) {
    uint dim = std_vec_->size();
    std_vec_->assign(dim,C);
  }

  virtual void print( std::ostream &outStream ) const {
    uint dim = std_vec_->size();
    for(uint i=0; i<dim; ++i) {
      outStream << (*std_vec_)[i] << " ";
    }
    outStream << std::endl;
  }

}; // class StdVector


} // namespace ROL

#endif
