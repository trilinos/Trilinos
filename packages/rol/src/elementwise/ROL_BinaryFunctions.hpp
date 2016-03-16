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

#ifndef ROL_BINARYFUNCTIONS_H
#define ROL_BINARYFUNCTIONS_H

#include "ROL_Elementwise_Function.hpp"

namespace ROL {
namespace Elementwise {

// Used to set every element in a vector to a specific value
template<class Real>
class Multiply : public BinaryFunction<Real> {
public:
  Multiply( ) { }

  Real apply( const Real &x, const Real &y ) const {
    return x*y;
  }  
}; // class Multiply

template<class Real>
class Plus : public BinaryFunction<Real> {
public:
  Real apply( const Real &x, const Real &y ) const {
    return x+y;
  }
};


// Used to set every element in a vector to a specific value
template<class Real>
class Divide : public BinaryFunction<Real> {
public:
  Divide( ) { }

  Real apply( const Real &x, const Real &y ) const {
    return x/y;
  }  
}; // class Divide


template<class Real> 
class Axpy : public BinaryFunction<Real> {
private:
  Real a_;
public:
  Axpy(Real a) : a_(a) {}
  Real apply( const Real &x, const Real &y ) const {
    return x+a_*y;
  }
};


template<class Real> 
class Aypx : public BinaryFunction<Real> {
private:
  Real a_;
public:
  Aypx(Real a) : a_(a) {}
  Real apply( const Real &x, const Real &y ) const {
    return a_*x+y;
  }
};


template<class Real> 
class Set : public BinaryFunction<Real> {
public:
  Real apply( const Real &x, const Real &y ) const {
    return y;
  }
};

template<class Real> 
class Lesser : public BinaryFunction<Real> {
public:
  Real apply( const Real &x, const Real &y ) const {
    return (x<y) ? x : y;
  }
};

template<class Real> 
class Greater : public BinaryFunction<Real> {
public:
  Real apply( const Real &x, const Real &y ) const {
    return (x>y) ? x : y;
  }
};


// Evaluate g(f(x,y))
template<class Real> 
class BinaryComposition : public BinaryFunction<Real> {

private:

  Teuchos::RCP<BinaryFunction<Real> > f_;
  Teuchos::RCP<UnaryFunction<Real> >  g_;

public:

  BinaryComposition( Teuchos::RCP<BinaryFunction<Real> > &f,
                     Teuchos::RCP<UnaryFunction<Real> > &g ) : f_(f), g_(g) {}
  Real apply( const Real &x, const Real &y ) const {
    return g_->apply(f_->apply(x,y));
  }

};


} // namespace Elementwise
} // namespace ROL




#endif // ROL_BINARYFUNCTIONS_H
