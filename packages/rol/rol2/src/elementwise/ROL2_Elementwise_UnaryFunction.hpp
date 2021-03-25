#pragma once
#ifndef ROL2_ELEMENTWISE_UNARYFUNCTION_HPP
#define ROL2_ELEMENTWISE_UNARYFUNCTION_HPP

#include <random>

namespace ROL2 {
namespace Elementwise {

// Forward Declaration of Derived classes 
template<typename> class AbsoluteValue;
template<typename> class Fill;
template<typename> class Heaviside;
template<typename> class Logarithm;
template<typename> class Power;
template<typename> class Reciprocal;
template<typename> class Round;
template<typename> class Scale;
template<typename> class Shift;
template<typename> class Sign;
template<typename> class SquareRoot;
template<typename> class ThresholdUpper;
template<typename> class ThresholdLower;

namespace Unary {

template<typename Real>
struct _Fill { 
  inline Real operator() ( Real x, Real v ) const noexcept { return v; } 
  using type = Fill<Real>;
};

template<typename Real>
struct _AbsoluteValue {
  inline Real operator() ( Real x, Real v ) const noexcept { return std::abs(x); };
  using type = AbsoluteValue<Real>;
};

template<typename Real>
struct _Heaviside {
  inline Real operator() ( Real x, Real v ) const noexcept { return x == 0 ? 0.5 : ( x>0 ? 1 : 0 ); };
  using type = Heaviside<Real>;
};

template<typename Real>
struct _Logarithm {
  inline Real operator() ( Real x, Real v ) const noexcept { return std::log(x); }
  using type = Logarithm<Real>;
};

template<typename Real> 
struct _Power {
  inline Real operator() ( Real x, Real v ) const noexcept { return std::pow(x,v); }
  using type = Power<Real>;
};


template<typename Real> 
struct _Reciprocal {
  inline Real operator() ( Real x, Real v ) const noexcept { return 1/v; }
  using type = Reciprocal<Real>;
};

template<typename Real> 
struct _Round {
  inline Real operator() ( Real x, Real v ) const noexcept { 
     auto fx = std::floor(x);
     auto cx = std::ceil(x);
     return 2*(x-fx) < 1 ? fx : cx;
   }
  using type = Round<Real>;
};

template<typename Real> 
struct _Scale {
  inline Real operator() ( Real x, Real v ) const noexcept { return v*x; }
  using type = Scale<Real>;
};

template<typename Real> 
struct _Shift {
  inline Real operator() ( Real x, Real v ) const noexcept { return v+x; }
  using type = Shift<Real>;
};

template<typename Real> 
struct _Sign {
  inline Real operator() ( Real x, Real v ) const noexcept { return x == 0 ? 0 : ( x>0 ? 1 : -1 ); };
  using type = Sign<Real>;
};

template<typename Real> 
struct _SquareRoot {
  inline Real operator() ( Real x, Real v ) const noexcept { return std::sqrt(x); }
  using type = SquareRoot<Real>;
};

template<typename Real> 
struct _ThresholdLower {
  inline Real operator() ( Real x, Real v ) const noexcept { return std::min(x,v); }
  using type = ThresholdLower<Real>;
};

template<typename Real> 
struct _ThresholdUpper {
  inline Real operator() ( Real x, Real v ) const noexcept { return std::max(x,v); }
  using type = ThresholdUpper<Real>;
};



} // namespace Unary

/** \class ROL2::Elementwise::UnaryFunction
 *  \brief Interface class with function of a single argument
 *
 */

template<typename Real>
using UnaryFunction = UFunction<Real,AbsoluteValue,
                                    Fill,
                                    Heaviside,
                                    Logarithm,
                                    Power,
                                    Reciprocal,
                                    Round,
                                    Scale,                                        
                                    Shift,
                                    Sign,
                                    SquareRoot,
                                    ThresholdUpper,
                                    ThresholdLower>;
namespace Unary {


///< TODO: Generalize Wrap to take, store, and call with arbitrary parameters
template<typename Real, template<typename> class F>
class Wrap : public UnaryFunction<Real> {
public:
  Wrap( Real value=0 ) : value_(value) {}

  void accept( typename UnaryFunction<Real>::visitor_type& visitor ) const override {
    visitor.visit( static_cast<const typename F<Real>::type&>(*this) );
  }

  Real apply( const Real& x ) const override { return f_(x,value_); };

  F<Real> get_function() const { return f_;     }
  Real    get_value()    const { return value_; }

private:
  F<Real> f_;
  Real value_;
};

} // namespace Unary


template<typename Real>
struct AbsoluteValue : public Unary::Wrap<Real,Unary::_AbsoluteValue> {};

template<typename Real>
struct Fill : public Unary::Wrap<Real,Unary::_Fill> {
  Fill( Real value ) : Unary::Wrap<Real,Unary::_Fill>::Wrap(value){}
};

template<typename Real> struct Heaviside : public Unary::Wrap<Real,Unary::_Heaviside> {};
template<typename Real> struct Logarithm : public Unary::Wrap<Real,Unary::_Logarithm> {};

template<typename Real>
struct Power : public Unary::Wrap<Real,Unary::_Power> {
  Power( Real value ) : Unary::Wrap<Real,Unary::_Power>::Wrap(value) {}
};

template<typename Real>
struct Reciprocal : public Unary::Wrap<Real,Unary::_Reciprocal> {};

template<typename Real> struct Round : public Unary::Wrap<Real,Unary::_Round> {};

template<typename Real>
struct Scale : public Unary::Wrap<Real,Unary::_Scale> {
  Scale( Real value ) : Unary::Wrap<Real,Unary::_Scale>::Wrap(value) {}
};

template<typename Real>
struct Shift : public Unary::Wrap<Real,Unary::_Shift> {
  Shift( Real value ) : Unary::Wrap<Real,Unary::_Shift>::Wrap(value) {}
};

template<typename Real> struct SquareRoot : public Unary::Wrap<Real,Unary::_SquareRoot> {};

template<typename Real>
struct Sign : public Unary::Wrap<Real,Unary::_Sign> {
  Sign( Real value ) : Unary::Wrap<Real,Unary::_Sign>::Wrap(value) {}
};

template<typename Real>
struct ThresholdLower : public Unary::Wrap<Real,Unary::_ThresholdLower> {
  ThresholdLower( Real value ) : Unary::Wrap<Real,Unary::_ThresholdLower>::Wrap(value) {}
};

template<typename Real>
struct ThresholdUpper : public Unary::Wrap<Real,Unary::_ThresholdUpper> {
  ThresholdUpper( Real value ) : Unary::Wrap<Real,Unary::_ThresholdUpper>::Wrap(value) {}
};

//---------------------------------------------------------------------------------------
// Visitor not yet supported:


template<typename Real>
class NormalRandom : public UnaryFunction<Real> {
public:
  NormalRandom( const Real&     mu = 0.0, 
                const Real&     sigma = 1.0,
                const unsigned& iseed = 0) {
    unsigned seed = iseed;
    if (seed == 0) seed = std::chrono::system_clock::now().time_since_epoch().count();
    gen_  = makePtr<std::mt19937_64>(seed);
    dist_ = makePtr<std::normal_distribution<Real>>(mu,sigma);
  }

  Real apply( const Real& x ) const override { return (*dist_)(*gen_); }

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
  Real lower_;
  Real upper_;
}; // class UniformlyRandom

// Multiply element by a uniformly distributed random number
// between lower and upper
template<typename Real>
class UniformlyRandomMultiply : public UnaryFunction<Real> {
public:
  UniformlyRandomMultiply( const Real& lower = 0.0, const Real& upper = 1.0) :
    lower_(lower), upper_(upper) {
  }

  Real apply( const Real& x ) const override {
    return x*((static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX)) * (upper_-lower_) + lower_);
  }

private:
  Real lower_;
  Real upper_;
}; // class UniformlyRandom

} // namespace Elementwise
} // namespace ROL2

#endif // ROL2_ELEMENTWISE_UNARYFUNCTION_HPP

