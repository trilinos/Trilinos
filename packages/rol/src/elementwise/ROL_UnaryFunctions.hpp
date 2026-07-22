// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_UNARYFUNCTIONS_H
#define ROL_UNARYFUNCTIONS_H

#include "ROL_Types.hpp"
#include "ROL_Elementwise_Function.hpp"

#include <cstdlib>
#include <ctime>
#include <random>
#include <chrono>


namespace ROL {
namespace Elementwise {

// Used to set every element in a vector to a specific value
template<class Real>
class Fill : public UnaryFunction<Real> {
public:
  Fill( const Real &value ) : value_(value) {}
  Real apply( const Real &x ) const {
    return value_;
  }  
private:  
  Real value_;
}; // class Fill

// Used to shift every element in a vector by a specific value
template<class Real>
class Shift : public UnaryFunction<Real> {
private:
  Real value_;
public:
  Shift( const Real &value ) : value_(value) {}
  Real apply( const Real &x ) const {
    return x+value_;
  }
}; // class Shift


// Get the elementwise reciprocal of a vector
template<class Real> 
class Reciprocal : public UnaryFunction<Real> {
public:
  Real apply( const Real &x ) const {
    return static_cast<Real>(1)/x;
  }  
}; // class Reciprocal

// Get the elementwise absolute value of a vector
template<class Real>
class AbsoluteValue : public UnaryFunction<Real> {
public:
  Real apply( const Real &x ) const {
    return std::abs(x); 
  }

};

template<class Real>
class Sign : public Elementwise::UnaryFunction<Real> {
private:
  Real zero_;
  Real one_;
public:
  Sign() : zero_(0), one_(1) {}
  Real apply(const Real &x) const {
    if(x==zero_) {
      return zero_;
    }
    else {
      return x>zero_ ? one_ : -one_;
    }
   }
};


// Compute the elementwise power of a vector
template<class Real> 
class Power : public UnaryFunction<Real> {
private:
  Real exponent_;
public:
  Power( const Real &exponent ) : exponent_(exponent) {}

  Real apply( const Real &x ) const {
    return std::pow(x,exponent_);
  } 
}; // class Power


// Compute the elementwise square root of a vector
template<class Real> 
class SquareRoot : public UnaryFunction<Real> {
public:
  SquareRoot( void ) {}

  Real apply( const Real &x ) const {
    return std::sqrt(x);
  } 
}; // class Power


// Generate a normally distributed random number
// with mean mu and standard deviation sigma
template<class Real> 
class NormalRandom : public UnaryFunction<Real> {
private:
  Ptr<std::mt19937_64>  gen_;
  Ptr<std::normal_distribution<Real>> dist_;

public:
  NormalRandom(const Real &mu = 0.0, const Real &sigma = 1.0,
               const unsigned &iseed = 0) {
    unsigned seed = iseed;
    if (seed == 0) {
      seed = std::chrono::system_clock::now().time_since_epoch().count();
    }
    gen_  = makePtr<std::mt19937_64>(seed);
    dist_ = makePtr<std::normal_distribution<Real>>(mu,sigma);
  }

  Real apply( const Real &x ) const {
    return (*dist_)(*gen_);
  }
}; // class NormalRandom


// Generate a uniformly distributed random number
// between lower and upper
template<class Real> 
class UniformlyRandom : public UnaryFunction<Real> {
private:
  const Real lower_;
  const Real upper_;

public:
  UniformlyRandom( const Real &lower = 0.0, const Real &upper = 1.0) : 
    lower_(lower), upper_(upper) {
  }

  Real apply( const Real &x ) const {
    return (static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX)) * (upper_-lower_) + lower_;
  }
}; // class UniformlyRandom

// Multiply element by a uniformly distributed random number
// between lower and upper
template<class Real> 
class UniformlyRandomMultiply : public UnaryFunction<Real> {
private:
  const Real lower_;
  const Real upper_;

public:
  UniformlyRandomMultiply( const Real &lower = 0.0, const Real &upper = 1.0) : 
    lower_(lower), upper_(upper) {
  }

  Real apply( const Real &x ) const {
    return x*((static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX)) * (upper_-lower_) + lower_);
  }
}; // class UniformlyRandom



// Returns max(x,s) where s is the given scalar
template<class Real>
class ThresholdUpper : public UnaryFunction<Real> {

private:
  const Real threshold_;

public:
  ThresholdUpper( const Real threshold ) : 
    threshold_(threshold) {}

  Real apply( const Real &x ) const {
    return std::max(threshold_,x);
  }
}; 

// Returns min(x,s) where s is the given scalar
template<class Real>
class ThresholdLower : public UnaryFunction<Real> {

private:
  const Real threshold_;

public:
  ThresholdLower( const Real threshold ) : 
    threshold_(threshold) {}

  Real apply( const Real &x ) const {
    return std::min(threshold_,x);
  }
}; 


template<class Real> 
class Scale : public UnaryFunction<Real> {
private:
  Real value_;
public:
  Scale( const Real value ) : value_(value) {}
  Real apply( const Real &x ) const {
    return value_*x;
  }
};




template<class Real>
class Logarithm : public UnaryFunction<Real> {
public:

  Real apply( const Real &x ) const {
    // To avoid circular dependency
    Real NINF = -0.1*std::abs(ROL::ScalarTraits<Real>::rmax()); 
    return (x>0) ? std::log(x) : NINF;
  }

};


// Heaviside step function
template<class Real>
class Heaviside : public UnaryFunction<Real> {
public:
 
  Real apply( const Real &x ) const {
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

};

template<class Real>
class Round : public UnaryFunction<Real> {
public:
  Round(void) {}
  Real apply(const Real &x) const {
    const Real half(0.5), fx = std::floor(x), cx = std::ceil(x);
    return (x-fx < half ? fx : cx); 
  }
};

// Evaluate g(f(x))
template<class Real> 
class UnaryComposition : public UnaryFunction<Real> {

private:
  
  ROL::Ptr<UnaryFunction<Real> > f_;
  ROL::Ptr<UnaryFunction<Real> > g_; 
  
public:
  UnaryComposition( ROL::Ptr<UnaryFunction<Real> > &f,
                    ROL::Ptr<UnaryFunction<Real> > &g ) : f_(f), g_(g) {}
  Real apply( const Real &x ) const {
    return g_->apply(f_->apply(x));
  }

};



} // namespace Elementwise
} // namespace ROL

#endif // ROL_UNARYFUNCTIONS_H
