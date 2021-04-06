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
class DivideAndInvert : public BinaryFunction<Real> {
public:
  DivideAndInvert( ) { }

  Real apply( const Real &x, const Real &y ) const {
    return y/x;
  }  
}; // class DivideAndInvert



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
class Min : public BinaryFunction<Real> {
public:
  Min() {}
  Real apply(const Real &x, const Real &y) const {
    return std::min(x,y);
  }
};


template<class Real>
class Max : public BinaryFunction<Real> {
public:
  Max() {}
  Real apply(const Real &x, const Real &y) const {
    return std::max(x,y);
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

// Set x to one of two values based on whether y satisfies
// a comparative condition
template<class Real>
class ValueSet : public BinaryFunction<Real> {
private:
  const Real threshold_;
  const int option_;
  const Real c1_;
  const Real c2_;
public:
  static const int LESS_THAN    = 0;
  static const int EQUAL_TO     = 1;
  static const int GREATER_THAN = 2;
  ValueSet( const Real& threshold, const int option, const Real &c1=Real(1), const Real &c2=Real(0) ) :
    threshold_(threshold), option_(option), c1_(c1), c2_(c2) {}
 
  Real apply(const Real &x, const Real &y ) const {
    Real result(c2_);
    switch( option_ ) {
      case LESS_THAN:    { result = y <  threshold_ ? c1_ : c2_; break; }
      case EQUAL_TO:     { result = y == threshold_ ? c1_ : c2_; break; }
      case GREATER_THAN: { result = y >  threshold_ ? c1_ : c2_; break; }
    }
    return result;
  }
};


// Evaluate g(f(x,y))
template<class Real> 
class BinaryComposition : public BinaryFunction<Real> {

private:

  ROL::Ptr<BinaryFunction<Real> > f_;
  ROL::Ptr<UnaryFunction<Real> >  g_;

public:

  BinaryComposition( ROL::Ptr<BinaryFunction<Real> > &f,
                     ROL::Ptr<UnaryFunction<Real> > &g ) : f_(f), g_(g) {}
  Real apply( const Real &x, const Real &y ) const {
    return g_->apply(f_->apply(x,y));
  }

};


} // namespace Elementwise
} // namespace ROL




#endif // ROL_BINARYFUNCTIONS_H
