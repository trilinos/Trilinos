// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_STDVECTOR_H
#define ROL_STDVECTOR_H

#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <initializer_list>
#include "ROL_Vector.hpp"

/** \class ROL::StdVector
    \brief Provides the std::vector implementation of the ROL::Vector interface.
*/


namespace ROL {

template <class Real, class Element=Real>
class StdVector : public Vector<Real> {

  using size_type = typename std::vector<Real>::size_type;

public:

  StdVector( const Ptr<std::vector<Element>> & std_vec ) : std_vec_(std_vec) {}

  StdVector( const int dim, const Element val=0.0 ) {
    std_vec_ = makePtr<std::vector<Element>>(dim,val);
  }

  StdVector( std::initializer_list<Element> ilist ) : 
    std_vec_( makePtr<std::vector<Element>>(ilist) ) {}

  Real& operator[] ( int i ) { return (*std_vec_)[i]; }
  const Real& operator[] ( int i ) const { return (*std_vec_)[i]; }

  void set( const Vector<Real> &x ) {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

    const StdVector &ex = static_cast<const StdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    std::copy(xval.begin(),xval.end(),std_vec_->begin());
  }

  void plus( const Vector<Real> &x ) {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

    const StdVector &ex = static_cast<const StdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    size_type dim  = std_vec_->size();
    for (size_type i=0; i<dim; i++) {
      (*std_vec_)[i] += xval[i];
    }
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

    const StdVector &ex = static_cast<const StdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    size_type dim  = std_vec_->size();
    for (size_type i=0; i<dim; i++) {
      (*std_vec_)[i] += alpha*xval[i];
    }
  }

  void scale( const Real alpha ) {
    for( auto& e : *std_vec_ ) e *= alpha;
//    size_type dim = std_vec_->size();
//    for (size_type i=0; i<dim; i++) {
//      (*std_vec_)[i] *= alpha;
//    }
  }

  virtual Real dot( const Vector<Real> &x ) const {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

    const StdVector& ex = static_cast<const StdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
//    size_type dim  = std_vec_->size();
//    Real val = 0;
//    for (size_type i=0; i<dim; i++) {
//      val += (*std_vec_)[i]*xval[i];
//    }
//    return val;
    return std::inner_product(std_vec_->begin(), std_vec_->end(), xval.begin(), Real(0));
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
    auto e = clone();
    (*staticPtrCast<StdVector>(e)->getVector())[i] = 1.0;
    return e;
  }

  int dimension() const {
    return static_cast<int>(std_vec_->size());
  }

  void applyUnary( const Elementwise::UnaryFunction<Real> &f ) {
//    size_type dim  = std_vec_->size();
//    for(size_type i=0; i<dim; ++i) {
//      (*std_vec_)[i] = f.apply((*std_vec_)[i]);
//    }
    for( auto& e : *std_vec_ ) e = f.apply(e);
  }

  void applyBinary( const Elementwise::BinaryFunction<Real> &f, const Vector<Real> &x ) {

    ROL_TEST_FOR_EXCEPTION( dimension() != x.dimension(),
                                std::invalid_argument,
                                "Error: Vectors must have the same dimension." );

    const StdVector & ex = static_cast<const StdVector&>(x);
    const std::vector<Element>& xval = *ex.getVector();
    size_type dim  = std_vec_->size();
    for (size_type i=0; i<dim; i++) {
      (*std_vec_)[i] = f.apply((*std_vec_)[i],xval[i]);
    }

  }

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const {
    Real result = r.initialValue();
    size_type dim  = std_vec_->size();
    for(size_type i=0; i<dim; ++i) {
      r.reduce((*std_vec_)[i],result);
    }
    return result;
  }

  void setScalar( const Real C ) {
    size_type dim = std_vec_->size();
    std_vec_->assign(dim,C);
  }

  void randomize( const Real l = 0.0, const Real u = 1.0 ) {
    Real a = (u-l);
    Real b = l;
//    Real x(0);
//    size_type dim = std_vec_->size();
//    for (size_type i=0; i<dim; ++i) {
//      x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
//      (*std_vec_)[i] = a*x + b;
//    }
    auto get_rand = [a,b]( Real& e ) { 
      auto x = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX); 
      e = a*x+b;
    };
    std::for_each( std_vec_->begin(), std_vec_->end(), get_rand );
  }

  virtual void print( std::ostream &outStream ) const {
//    size_type dim = std_vec_->size();
//    for(size_type i=0; i<dim; ++i) {
//      outStream << (*std_vec_)[i] << " ";
//    }
    for( auto e : *std_vec_ ) outStream << e << " ";
    outStream << std::endl;
  }

private:

  Ptr<std::vector<Element>>  std_vec_;

}; // class StdVector


} // namespace ROL

#endif
