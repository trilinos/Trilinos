// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_DYADICOPERATOR_H
#define ROL_DYADICOPERATOR_H

#include "ROL_LinearOperator.hpp"


/** @ingroup func_group
    \class ROL::DyadicOperator
    \brief Interface to apply a dyadic operator to a vector
*/

namespace ROL {

// Implementation of a Dyadic operator x*y'
template<class Real> 
class DyadicOperator : public ROL::LinearOperator<Real> {
  
  typedef ROL::Vector<Real> V;

private:

  const ROL::Ptr<const V> x_;
  const ROL::Ptr<const V> y_;

public:
  
  DyadicOperator( const ROL::Ptr<const V> &x,
                  const ROL::Ptr<const V> &y ) : x_(x), y_(y) {}

  void apply( V &Hv, const V &v, Real &tol ) const {
    Hv.set(*x_);
    Hv.scale(v.dot(*y_));  
  }
    
  void applyInverse( V &Hv, const V &v, Real &tol ) const {

    ROL_TEST_FOR_EXCEPTION( true , std::logic_error, 
                                ">>> ERROR (ROL_DyadicOperator, applyInverse): "
                                "Not implemented."); 

  } 
}; // class DyadicOperator

} // namespace ROL




#endif // ROL_NULLOPERATOR_H

