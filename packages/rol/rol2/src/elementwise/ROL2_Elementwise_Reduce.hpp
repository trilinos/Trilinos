#pragma once
#ifndef ROL_ELEMENTWISE_REDUCE_HPP
#define ROL_ELEMENTWISE_REDUCE_HPP

namespace ROL2 {
namespace Elementwise {

template<typename> class EuclideanNormSquared;
template<typename> class ReductionAnd;
template<typename> class ReductionMax;
template<typename> class ReductionMin;
template<typename> class ReductionSum;

namespace Reduction {

template<typename Real, template<typename> class...DerivedTypes>
struct Op {
  using visitor_type = Visitor<DerivedTypes<Real>...>;
  virtual ~Op() = default;
  virtual void reduce( const Real& input, Real& output ) const = 0;
  virtual void reduce( const volatile Real& input, volatile Real& output ) const = 0;    
  virtual Real initialValue() const = 0;
  virtual void accept( visitor_type& ) const {}
};  

} // namespace Reduction 

template<typename Real>
using ReductionOp = Reduction::Op<Real,EuclideanNormSquared,
                                       ReductionAnd,
                                       ReductionMax,
                                       ReductionMin,
                                       ReductionSum>;

template<class Real>
class EuclideanNormSquared : public ReductionOp<Real> {
public:

  virtual ~EuclideanNormSquared() = default;

  void reduce( const Real& input, Real& output ) const override {
    output += input * input;
  }

  void reduce( const volatile Real& input, volatile Real& output ) const override {
    output += input * input;
  }

  Real initialValue() const override { return 0;  }

  void accept( typename ReductionOp<Real>::Visitor& visitor ) const override {
    visitor.visit(*this);
  }
}; 

template<typename Real>
class ReductionAnd : public ReductionOp<Real> {
public:
  virtual ~ReductionAnd() = default;

  void reduce( const Real& input, Real& output ) const override {
    output = (input*output)==0 ? 0.0 : 1.0;
  }

  void reduce( const volatile Real& input, volatile Real& output ) const override {
    output = (input*output)==0 ? 0.0 : 1.0;
  }

  Real initialValue() const override { return 1.0; }

  void accept( typename ReductionOp<Real>::Visitor& visitor ) const override {
    visitor.visit(*this);
  }
};


template<class Real>
class ReductionMax : public ReductionOp<Real> {
public:
  virtual ~ReductionMax() = default;

  ReductionMax() {
      ROL_TEST_FOR_EXCEPTION(!std::numeric_limits<Real>::is_specialized,
      std::logic_error,"ReductionMax() requires std::numeric_limits "
        "be specialized on supplied template parameter.\n" );
  }

  void reduce( const Real& input, Real& output ) const override {
    output = (input>output) ? input : output;
  }

  void reduce( const volatile Real& input, volatile Real& output ) const override {
    output = (input>output) ? input : output;
  }

  Real initialValue() const override { return std::numeric_limits<Real>::min(); }

  void accept( typename ReductionOp<Real>::Visitor& visitor ) const override {
    visitor.visit(*this);
  }
};

template<class Real>
class ReductionMin : public ReductionOp<Real> {
public:
  virtual ~ReductionMin() = default;

  ReductionMin() {
      ROL_TEST_FOR_EXCEPTION(!std::numeric_limits<Real>::is_specialized,
      std::logic_error,"ReductionMin() requires std::numeric_limits "
        "be specialized on supplied template parameter.\n" );
  }

  void reduce( const Real& input, Real& output ) const override {
    output = (input<output) ? input : output;
  }

  void reduce( const volatile Real& input, Real volatile& output ) const override {
    output = (input<output) ? input : output;
  }

  Real initialValue() const override { return std::numeric_limits<Real>::max(); } 

  void accept( typename ReductionOp<Real>::Visitor& visitor ) const override {
    visitor.visit(*this);
  }

};

template<class Real>
class ReductionSum : public ReductionOp<Real> {
public:
  virtual ~ReductionSum() = default;

  void reduce( const Real& input, Real& output ) const override {
    output = output + input;
  }

  void reduce( const volatile Real& input, volatile Real& output ) const override {
    output = output + input;
  }

  Real initialValue() const override { return 0; }

  void accept( typename ReductionOp<Real>::Visitor& visitor ) const override {
    visitor.visit(*this);
  }

};

} // namespace Elementwise
} // namespace ROL2

#endif //ROL_ELEMENTWISE_REDUCE_HPP

