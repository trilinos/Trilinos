// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_VECTOR_H
#include "ROL_Vector.hpp"
#else

#ifndef ROL_ELEMENTWISE_FUNCTION_H
#define ROL_ELEMENTWISE_FUNCTION_H

namespace ROL {

template<class Real>
class Vector;

namespace Elementwise {

// Interface class with function of a single argument
template<class Real>
class UnaryFunction {
public:
  virtual ~UnaryFunction() {}
  virtual Real apply( const Real &x ) const = 0;
};

// Wrap a generic function/functor of a single argument
template<class Real, class Func>
class UnaryFunctionWrapper : public UnaryFunction<Real> {
public:
  UnaryFunctionWrapper( const Func &f ) : f_(f) {}
  
  Real apply( const Real &x ) const {
    return f_(x);
  }
private:
  Func f_;
};

// Apply a Unary function to a ROL Vector
template<class Real>
void applyUnaryInPlace(ROL::Vector<Real> &x, const UnaryFunction<Real> &f) {
  x.applyUnary(f);
}

// Wrap a generic single argument function or functor and and apply it to a ROL Vector
template<class Real, class Func>
void applyUnaryInPlace(ROL::Vector<Real> &x, Func f) {
  UnaryFunctionWrapper<Real,Func> f_wrapped(f);
  x.applyUnary(f_wrapped);
}



// Interface class with function of two arguments
template<class Real>
class BinaryFunction {
public:
  virtual ~BinaryFunction() {}
  virtual Real apply( const Real &x, const Real &y ) const = 0;
};


// Wrap a generic function/functor of two arguments
template<class Real, class Func>
class BinaryFunctionWrapper : public BinaryFunction<Real> {
public:
  BinaryFunctionWrapper( const Func &f ) : f_(f) {}
  
  Real apply( const Real &x, const Real &y ) const {
    return f_(x,y); 
  }
private:
  Func f_;   
};

// Apply a Binary function to a ROL Vector
template<class Real>
void applyBinaryInPlace(ROL::Vector<Real> &x, const ROL::Vector<Real> &y, const BinaryFunction<Real> &f) {
  x.applyBinary(f,y);
}


// Wrap a generic two argument function or functor and and apply it to a ROL Vector
template<class Real, class Func>
void applyBinaryInPlace(ROL::Vector<Real> &x, const ROL::Vector<Real> &y, Func f) {
  BinaryFunctionWrapper<Real,Func> f_wrapped(f);
  x.applyBinary(f_wrapped,y);
}


} // namespace Elementwise 
} // namespace ROL

#include <ROL_Elementwise_Reduce.hpp>
#include <ROL_UnaryFunctions.hpp>
#include <ROL_BinaryFunctions.hpp>

#endif
#endif
