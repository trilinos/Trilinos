#pragma once
#ifndef ROL2_ELEMENTWISE_BINARYFUNCTION_HPP
#define ROL2_ELEMENTWISE_BINARYFUNCTION_HPP

namespace ROL2 {
namespace Elementwise {

template<typename Real> class Axpy;
template<typename Real> class Divide;
template<typename Real> class DivideAndInvert;
template<typename Real> class Greater;
template<typename Real> class Lesser;
template<typename Real> class Max;
template<typename Real> class Min;
template<typename Real> class Multiply;
template<typename Real> class Plus;
template<typename Real> class Set;
template<typename Real> class ValueSet;

namespace Binary {

template<typename Real>
struct _Divide { 
  inline Real operator() ( Real x, Real y ) const noexcept { return x/y; }
  using type = Divide<Real>;
};


template<typename Real>
struct _DivideAndInvert { 
  inline Real operator() ( Real x, Real y ) const noexcept { return y/x; }
  using type = DivideAndInvert<Real>;
};


template<typename Real>
struct _Greater { 
  inline Real operator() ( Real x, Real y ) const noexcept { return static_cast<Real>(x>y); }
  using type = Greater<Real>;
};


template<typename Real>
struct _Lesser { 
  inline Real operator() ( Real x, Real y ) const noexcept { return static_cast<Real>(x<y); }
  using type = Lesser<Real>;
};


template<typename Real>
struct _Max { 
  inline Real operator() ( Real x, Real y ) const noexcept { return x > y ? x : y; }
  using type = Max<Real>;
};


template<typename Real>
struct _Min { 
  inline Real operator() ( Real x, Real y ) const noexcept { return x < y ? x : y; }
  using type = Min<Real>;
};


template<typename Real>
struct _Multiply { 
  inline Real operator() ( Real x, Real y ) const noexcept { return x*y; }
  using type = Multiply<Real>;
};


template<typename Real>
struct _Plus { 
  inline Real operator() ( Real x, Real y ) const noexcept { return x+y; }
  using type = Plus<Real>;
};


template<typename Real>
struct _Set { 
  inline Real operator() ( Real x, Real y ) const noexcept { return y; }
  using type = Set<Real>;
};

} // namespace Binary

template<typename Real>
using BinaryFunction = Function<Real,Divide,
                                     DivideAndInvert,
                                     Greater,
                                     Lesser,
                                     Max,
                                     Min,
                                     Multiply,
                                     Plus,
                                     Set>;
namespace Binary {

template<typename Real, template<typename> class F>
class Wrap : public BinaryFunction<Real> {
public:
  void accept( typename BinaryFunction<Real>::visitor_type& visitor ) const override {
    visitor.visit( static_cast<const typename F<Real>::type&>(*this) );
  }

  Real apply( const Real& x, const Real& y ) const override { return f_(x,y); };

  F<Real> get_function() const { return f_; }

private:
  F<Real> f_;
}; // class Wrap

} // namespace Binary

template<typename Real> struct Divide          : public Binary::Wrap<Real,Binary::_Divide> {};
template<typename Real> struct DivideAndInvert : public Binary::Wrap<Real,Binary::_DivideAndInvert> {};
template<typename Real> struct Greater         : public Binary::Wrap<Real,Binary::_Greater> {};
template<typename Real> struct Lesser          : public Binary::Wrap<Real,Binary::_Lesser> {};
template<typename Real> struct Max             : public Binary::Wrap<Real,Binary::_Max> {};
template<typename Real> struct Min             : public Binary::Wrap<Real,Binary::_Min> {};
template<typename Real> struct Multiply        : public Binary::Wrap<Real,Binary::_Multiply> {};
template<typename Real> struct Plus            : public Binary::Wrap<Real,Binary::_Plus> {};
template<typename Real> struct Set             : public Binary::Wrap<Real,Binary::_Set> {};

//-------------------------------------------------------------------------------
// Visitor not yet supported:

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



} // namespace Elementwise
} // namespace ROL2

#endif //ROL2_ELEMENTWISE_BINARYFUNCTION_HPP

