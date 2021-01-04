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

#pragma once
#ifndef ROL_BINARYFUNCTIONS_H
#define ROL_BINARYFUNCTIONS_H

namespace ROL {
namespace Elementwise {

template<class Real> 
class Axpy : public BinaryFunction<Real> {
public:
  Axpy(Real alpha) : alpha_(alpha) {}
  virtual ~Axpy() = default;

  Real apply( const Real &x, const Real &y ) const {
    return x+alpha_*y;
  }

  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

  Real get_alpha() const { return alpha_; } 

private:
  Real alpha_;
};

template<class Real> 
class Aypx : public BinaryFunction<Real> {
public:
  Aypx(Real alpha) : alpha_(alpha) {}
  virtual ~Aypx() = default;

  Real apply( const Real &x, const Real &y ) const {
    return alpha_*x+y;
  }

  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

  Real get_alpha() const { return alpha_; }

private:
  Real alpha_;
};



// Used to set every element in a vector to a specific value
template<class Real>
class Multiply : public BinaryFunction<Real> {
public:

  Multiply() = default;
  virtual ~Multiply() = default;

  Real apply( const Real &x, const Real &y ) const {
    return x*y;
  }

  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
  
}; // class Multiply

template<class Real>
class Plus : public BinaryFunction<Real> {
public:

  Plus() = default;
  virtual ~Plus() = default;

  Real apply( const Real &x, const Real &y ) const {
    return x+y;
  }

  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
};


// Used to set every element in a vector to a specific value
template<class Real>
class Divide : public BinaryFunction<Real> {
public:

  Divide() = default;
  virtual ~Divide() = default;

  Real apply( const Real &x, const Real &y ) const {
    return x/y;
  }
  
  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
}; // class Divide

template<class Real>
class DivideAndInvert : public BinaryFunction<Real> {
public:

  DivideAndInvert() = default;
  virtual ~DivideAndInvert() = default;

  Real apply( const Real &x, const Real &y ) const {
    return y/x;
  }  

  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
}; // class DivideAndInvert



template<class Real>
class Min : public BinaryFunction<Real> {
public:
  Min() = default;
  virtual ~Min() = default;
  Real apply(const Real &x, const Real &y) const {
    return std::min(x,y);
  }

  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
};


template<class Real>
class Max : public BinaryFunction<Real> {
public:
  Max() = default;
  virtual ~Max() = default;
  Real apply(const Real &x, const Real &y) const {
    return std::max(x,y);
  }

  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
};


template<class Real> 
class Set : public BinaryFunction<Real> {
public:
  Set() = default;
  virtual ~Set() = default;
  Real apply( const Real &x, const Real &y ) const {
    return y;
  }
  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
};

template<class Real> 
class Lesser : public BinaryFunction<Real> {
public:
  Lesser() = default;
  virtual ~Lesser() = default;
  Real apply( const Real &x, const Real &y ) const {
    return static_cast<Real>(x<y);
  }

  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
};

template<class Real> 
class Greater : public BinaryFunction<Real> {
public:
  Greater() = default;
  virtual ~Greater() = default;

  Real apply( const Real &x, const Real &y ) const {
    return static_cast<Real>(x>y);
  }

  void accept( typename BinaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
};

// Set x to one of two values based on whether y satisfies
// a comparative condition
template<class Real>
class ValueSet : public BinaryFunction<Real> {
public:

  static const int LESS_THAN    = 0;
  static const int EQUAL_TO     = 1;
  static const int GREATER_THAN = 2;

  ValueSet( Real threshold, int option, Real c1=Real(1), Real c2=Real(0) ) :
    threshold_(threshold), c1_(c1), c2_(c2), option_(option) {}
 
  Real apply( const Real &x, const Real &y ) const {
    Real result(c2_);
    switch( option_ ) {
      case LESS_THAN:    { result = y <  threshold_ ? c1_ : c2_; break; }
      case EQUAL_TO:     { result = y == threshold_ ? c1_ : c2_; break; }
      case GREATER_THAN: { result = y >  threshold_ ? c1_ : c2_; break; }
    }
    return result;
  }
  Real get_threshold() const { return threshold_; }
  Real get_c1() const { return c1_; }
  Real get_c2() const { return c2_; }
  int get_option() const { return option_; }

private:
  Real threshold_;
  Real c1_;
  Real c2_;
  int option_;
};

} // namespace Elementwise
} // namespace ROL




#endif // ROL_BINARYFUNCTIONS_H
