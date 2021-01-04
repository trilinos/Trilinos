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

#ifndef ROL_UNARYFUNCTIONS_H
#define ROL_UNARYFUNCTIONS_H

#include <ctime>
#include <random>
#include <chrono>

namespace ROL {
namespace Elementwise {

// Used to set every element in a vector to a specific value
template<typename Real>
class Fill : public UnaryFunction<Real> {
public:
  Fill( const Real& value ) : value_(value) {}
  Real apply( const Real& x ) const override{
    return value_;
  }  

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

  Real get_value() const { return value_; }

private:  
  Real value_;
}; // class Fill

// Used to shift every element in a vector by a specific value
template<typename Real>
class Shift : public UnaryFunction<Real> {
public:
  Shift( const Real& value ) : value_(value) {}
  Real apply( const Real& x ) const override {
    return x+value_;
  }

  Real get_value() const { return value_; }

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

private:
  Real value_;

}; // class Shift


// Get the elementwise reciprocal of a vector
template<typename Real> 
class Reciprocal : public UnaryFunction<Real> {
public:
  Real apply( const Real& x ) const override{
    return static_cast<Real>(1)/x;
  }  

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

}; // class Reciprocal

// Get the elementwise absolute value of a vector
template<typename Real>
class AbsoluteValue : public UnaryFunction<Real> {
public:
  Real apply( const Real& x ) const override {
    return std::abs(x); 
  }

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

};

template<typename Real>
class Sign : public Elementwise::UnaryFunction<Real> {
public:
  Sign() : zero_(0), one_(1) {}
  Real apply(const Real& x) const override{
    if(x==zero_) {
      return zero_;
    }
    else {
      return x>zero_ ? one_ : -one_;
    }
   }

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

private:
  Real zero_;
  Real one_;
};


// Compute the elementwise power of a vector
template<typename Real> 
class Power : public UnaryFunction<Real> {
public:
  Power( const Real& exponent ) : exponent_(exponent) {}

  Real apply( const Real& x ) const override {
    return std::pow(x,exponent_);
  } 

  Real get_exponent() const { return exponent_; }

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

private:
  Real exponent_;

}; // class Power


// Compute the elementwise square root of a vector
template<typename Real> 
class SquareRoot : public UnaryFunction<Real> {
public:
  SquareRoot( void ) {}

  Real apply( const Real& x ) const override {
    return std::sqrt(x);
  } 

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

}; // class Power


// Generate a normally distributed random number
// with mean mu and standard deviation sigma
template<typename Real> 
class NormalRandom : public UnaryFunction<Real> {
public:
  NormalRandom(const Real& mu = 0.0, const Real& sigma = 1.0,
               const unsigned &iseed = 0) {
    unsigned seed = iseed;
    if (seed == 0) {
      seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    gen_  = makePtr<std::mt19937_64>(seed);
    dist_ = makePtr<std::normal_distribution<Real>>(mu,sigma);
  }

  Real apply( const Real& x ) const override {
    return (*dist_)(*gen_);
  }

private:
  Ptr<std::mt19937_64>  gen_;
  Ptr<std::normal_distribution<Real>> dist_;

}; // class NormalRandom


// Generate a uniformly distributed random number
// between lower and upper
template<typename Real> 
class UniformlyRandom : public UnaryFunction<Real> {
public:
  UniformlyRandom( const Real& lower = 0.0, const Real& upper = 1.0) : 
    lower_(lower), upper_(upper) {
  }

  Real apply( const Real& x ) const override {
    return (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX)) * (upper_-lower_) + lower_;
  }

private:
  const Real lower_;
  const Real upper_;

}; // class UniformlyRandom

// Multiply element by a uniformly distributed random number
// between lower and upper
template<typename Real> 
class UniformlyRandomMultiply : public UnaryFunction<Real> {
private:
  const Real lower_;
  const Real upper_;

public:
  UniformlyRandomMultiply( const Real& lower = 0.0, const Real& upper = 1.0) : 
    lower_(lower), upper_(upper) {
  }

  Real apply( const Real& x ) const {
    return x*((static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX)) * (upper_-lower_) + lower_);
  }
}; // class UniformlyRandom



// Returns max(x,s) where s is the given scalar
template<typename Real>
class ThresholdUpper : public UnaryFunction<Real> {
public:
  ThresholdUpper( const Real threshold ) : 
    threshold_(threshold) {}

  Real apply( const Real& x ) const override {
    return std::max(threshold_,x);
  }

  Real get_threshold() const { return threshold_; }

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

private:
  const Real threshold_;
}; 

// Returns min(x,s) where s is the given scalar
template<typename Real>
class ThresholdLower : public UnaryFunction<Real> {
public:
  ThresholdLower( const Real threshold ) : 
    threshold_(threshold) {}

  Real apply( const Real& x ) const override {
    return std::min(threshold_,x);
  }

  Real get_threshold() const { return threshold_; }

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

private:
  const Real threshold_;

}; 

template<typename Real> 
class Scale : public UnaryFunction<Real> {
public:
  Scale( const Real value ) : value_(value) {}

  Real apply( const Real& x ) const override {
    return value_*x;
  }

  Real get_value() const { 
    return value_;
  }

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

private:
  Real value_;
};



template<typename Real>
class Logarithm : public UnaryFunction<Real> {
public:

  Real apply( const Real& x ) const {
    // To avoid circular dependency
    Real NINF = -0.1*std::abs(ROL::ScalarTraits<Real>::rmax()); 
    return (x>0) ? std::log(x) : NINF;
  }

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }

};


// Heaviside step function
template<typename Real>
class Heaviside : public UnaryFunction<Real> {
public:
 
  Real apply( const Real& x ) const {
    Real value = 0;
    if( x>0 ) {
      value = 1.0;
    } else if( x==0 ) {
      value = 0.5;
    } else {
      value = 0.0;
    }
    return value;
  }

  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
};

template<typename Real>
class Round : public UnaryFunction<Real> {
public:
  Round() {}
  Real apply(const Real& x) const {
    const Real half(0.5), fx = std::floor(x), cx = std::ceil(x);
    return (x-fx < half ? fx : cx); 
  }
  
  void accept( typename UnaryFunction<Real>::Visitor& visitor ) const override {
    visitor.visit( *this );
  }
};

} // namespace Elementwise
} // namespace ROL

#endif // ROL_UNARYFUNCTIONS_H
