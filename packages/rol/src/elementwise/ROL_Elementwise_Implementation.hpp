
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

#ifndef ROL_ELEMENTWISE_IMPLEMENTATION_H
#define ROL_ELEMENTWISE_IMPLEMENTATION_H

#include <ROL_Vector.hpp>

namespace ROL {

template<class Real>
class Vector;

namespace Elementwise {

// Implementation of applyUnary nonmember function
template<class Real, class Func>
class ApplyUnary {
  typedef ROL::Vector<Real> V;
public:
  static void inPlace( V &x, Func f ) {
    UnaryFunctionWrapper<Real, Func> *f_wrapped =
      new UnaryFunctionWrapper<Real, Func> (f);

    try {
      x.applyUnary(f_wrapped);
    } catch (...) {
      delete f_wrapped;
      throw;
    }
    delete f_wrapped;
  }

};  

// Partial specialization for UnaryFunction<Real>*
template<class Real>
class ApplyUnary<Real, UnaryFunction<Real>* > {
  typedef ROL::Vector<Real> V;
public:
  static void inPlace ( V &x, const V &y, UnaryFunction<Real> *f ) {
    x.applyUnary(f);
  }
};


template<class Real, class Func>
void applyUnaryInPlace(ROL::Vector<Real> &x, Func f) {
  ApplyUnary<Real,Func>::inPlace(x,f);
}



// Implementation of applyBinary nonmember function
template<class Real, class Func>
class ApplyBinary {
  typedef ROL::Vector<Real> V;
public:
  static void inPlace( V &x, const V &y, Func f ) {
    BinaryFunctionWrapper<Real, Func> *f_wrapped =
      new BinaryFunctionWrapper<Real, Func> (f);

    try {
      x.applyBinary(f_wrapped, y);
    } catch (...) {
      delete f_wrapped;
      throw;
    }
    delete f_wrapped;
  }
};  

// Partial specialization for BinaryFunction<Real>*
template<class Real>
class ApplyBinary<Real, BinaryFunction<Real>* > {
  typedef ROL::Vector<Real> V;
public:
  static void inPlace ( V &x, const V &y, BinaryFunction<Real> *f ) {
    x.applyBinary(f,y);
  }
};

template<class Real, class Func>
void applyBinaryInPlace(ROL::Vector<Real> &x, const ROL::Vector<Real> &y, Func f) {
  ApplyBinary<Real,Func>::inPlace(x,y,f);
}



} // namespace Elementwise
} // namespace ROL

#endif
