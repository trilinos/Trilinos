
// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_ELEMENTWISE_VECTOR_H
#define ROL_ELEMENTWISE_VECTOR_H

#include "ROL_Vector.hpp"

/** @ingroup la_group
    \class ROL::ElementwiseVector
    \brief Intermediate abstract class which does not require users 
           implements plus, set, scale, axpy, norm, dot, or zero if 
           they implement the three elementwise functions: applyUnary,
           applyBinary, and reduce

           dot and norm are unweighted dot products and Euclidean norm 
           by default
*/

namespace ROL {

template< class Real>
class ElementwiseVector : public Vector<Real> {

public: 

  virtual ~ElementwiseVector() {}

  void plus( const Vector<Real> &x ) {
    this->applyBinary(Elementwise::Plus<Real>(),x);
  }

  void scale( const Real alpha ) {
    this->applyUnary(Elementwise::Scale<Real>(alpha));
  }

  virtual Real dot( const Vector<Real> &x ) const {
    ROL::Ptr<Vector<Real> > y = this->clone();
    y->set(*this);
    y->applyBinary(Elementwise::Multiply<Real>(),x);
    return y->reduce(Elementwise::ReductionSum<Real>());    
  }

  virtual Real norm() const {
      return std::sqrt(this->dot(*this));
  }

  void axpy( const Real alpha, const Vector<Real> &x ) {
    this->applyBinary(Elementwise::Axpy<Real>(alpha),x);
  }

  void zero() {
    this->applyUnary(Elementwise::Fill<Real>(Real(0)));
  }

  void set( const Vector<Real> &x ) {
    this->applyBinary(Elementwise::Set<Real>(),x);
  }

  // MUST overload these three functions
  virtual void applyUnary( const Elementwise::UnaryFunction<Real> &uf ) = 0;

  virtual void applyBinary( const Elementwise::BinaryFunction<Real> &bf,
                            const Vector<Real> &x ) = 0;

  virtual Real reduce( const Elementwise::ReductionOp<Real> &r ) const = 0;

}; // class ElementwiseVector


} // namespace ROL




#endif // ROL_ELEMENTWISE_VECTOR_H

