// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MiniTensor_Solvers.h"

namespace minitensor {

//
// Define some nonlinear systems (NLS) to test nonlinear solution methods.
//

//
//
//
template<typename S, Index M = 1>
class SquareRoot : public Function_Base<SquareRoot<S, M>, S, M>
{
public:

  SquareRoot(S const c) : c_(c)
  {
  }

  static constexpr
  char const * const
  NAME{"Square Root"};

  using Base = Function_Base<SquareRoot<S, M>, S, M>;

  // Default value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = x(0) * x(0) - c_;

    return r;
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

private:
  S const
  c_{0.0};
};

//
//
//
template<typename S, Index M = 2>
class Quadratic : public Function_Base<Quadratic<S, M>, S, M>
{
public:
  Quadratic(S const a, S const b, S const c) :  a_(a), b_(b), c_(c)
  {
  }

  static constexpr
  char const * const
  NAME{"Quadratic"};

  using Base = Function_Base<Quadratic<S, M>, S, M>;

  // Default value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = 2.0 * c_ * (x(0) - a_);
    r(1) = 2.0 * c_ * (x(1) - b_);

    return r;
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

private:
  S const
  a_{0.0};

  S const
  b_{0.0};

  S const
  c_{0.0};
};

//
//
//
template<typename S, Index M = 2>
class Gaussian : public Function_Base<Gaussian<S, M>, S, M>
{
public:
  Gaussian(S const a, S const b, S const c) : a_(a), b_(b), c_(c)
  {
  }

  static constexpr
  char const * const
  NAME{"Inverted Gaussian"};

  using Base = Function_Base<Gaussian<S, M>, S, M>;

  // Default value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    T const
    xa = (x(0) - a_) * c_;

    T const
    xb = (x(1) - b_) * c_;

    T const
    e = std::exp(- xa * xa - xb * xb);

    r(0) = 2.0 * xa * e * c_ * c_;
    r(1) = 2.0 * xb * e * c_ * c_;

    return r;
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

private:
  S const
  a_{0.0};

  S const
  b_{0.0};

  S const
  c_{0.0};
};

//
//
//
template<typename S, Index M = 2>
class Banana : public Function_Base<Banana<S, M>, S, M>
{
public:

  Banana()
  {
  }

  static constexpr
  char const * const
  NAME{"Rosenbrock's Banana"};

  using Base = Function_Base<Banana<S, M>, S, M>;

  // Default value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = 2.0 * (x(0) - 1.0) + 400.0 * x(0) * (x(0) * x(0) - x(1));
    r(1) = 200.0 * (x(1) - x(0) * x(0));

    return r;
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
//
//
template<typename S, Index M = 2>
class Matyas : public Function_Base<Matyas<S, M>, S, M>
{
public:

  Matyas() {}

  static constexpr
  char const * const
  NAME{"Matyas"};

  using Base = Function_Base<Matyas<S, M>, S, M>;

  // Default value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = (13.0 * x(0) - 12.0 * x(1)) / 25.0;
    r(1) = (13.0 * x(1) - 12.0 * x(0)) / 25.0;

    return r;
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
//
//
template<typename S, Index M = 2>
class McCormick : public Function_Base<McCormick<S, M>, S, M>
{
public:

  McCormick() {}

  static constexpr
  char const * const
  NAME{"McCormick"};

  using Base = Function_Base<McCormick<S, M>, S, M>;

  // Default value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = std::cos(x(0) + x(1)) + 2.0 * x(0) - 2.0 * x(1) - 1.5;
    r(1) = std::cos(x(0) + x(1)) - 2.0 * x(0) + 2.0 * x(1) + 2.5;

    return r;
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
//
//
template<typename S, Index M = 2>
class StyblinskiTang : public Function_Base<StyblinskiTang<S, M>, S, M>
{
public:

  StyblinskiTang() {}

  static constexpr
  char const * const
  NAME{"Styblinski-Tang"};

  using Base = Function_Base<StyblinskiTang<S, M>, S, M>;

  // Default value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    return Base::value(*this, x);
  }

  // Explicit gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x) const
  {
    Index const
    dimension = x.get_dimension();

    assert(dimension == Base::DIMENSION);

    Vector<T, N>
    r(dimension);

    r(0) = 2.0 * x(0) * x(0) * x(0) - 16.0 * x(0) + 2.5;
    r(1) = 2.0 * x(1) * x(1) * x(1) - 16.0 * x(1) + 2.5;

    return r;
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Define some nonlinear functions (NLF) to test nonlinear optimization methods.
//

//
// Paraboloid of revolution
//
template<typename S, Index M = 2>
class Paraboloid : public Function_Base<Paraboloid<S, M>, S, M>
{
public:

  Paraboloid(S xc = 0.0, S yc = 0.0) : xc_(xc), yc_(yc)
  {
  }

  static constexpr
  char const * const
  NAME{"Paraboloid"};

  using Base = Function_Base<Paraboloid<S, M>, S, M>;

  // Explicit value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    assert(x.get_dimension() == Base::DIMENSION);

    T const
    a = x(0) - xc_;

    T const
    b = x(1) - yc_;

    T const
    f = a * a + b * b;

    return f;
  }

  // Default AD gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

private:
  S
  xc_{0.0};

  S
  yc_{0.0};
};

//
//
//
template<typename S, Index M = 2>
class Rosenbrock : public Function_Base<Rosenbrock<S, M>, S, M>
{
public:

  Rosenbrock(S a = 1.0, S b = 100.0) : a_(a), b_(b)
  {
  }

  static constexpr
  char const * const
  NAME{"Rosenbrock's Function 2D"};

  using Base = Function_Base<Rosenbrock<S, M>, S, M>;

  // Explicit value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & x)
  {
    T const
    a = (a_ - x(0));

    T const
    b = (x(1) - x(0) * x(0));

    return a * a + b_ * b * b;
  }

  // Default AD gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

private:
  S
  a_{1.0};

  S
  b_{100.0};
};

//
// Beale's function
//
template<typename S, Index M = 2>
class Beale : public Function_Base<Beale<S, M>, S, M>
{
public:

  Beale() {}

  static constexpr
  char const * const
  NAME{"Beale"};

  using Base = Function_Base<Beale<S, M>, S, M>;

  // Explicit value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
  {
    assert(X.get_dimension() == Base::DIMENSION);

    T const &
    x = X(0);

    T const &
    y = X(1);

    T const
    a = 1.5 - x + x * y;

    T const
    b = 2.25 - x + x * y * y;

    T const
    c = 2.625 - x + x * y * y * y;

    T const
    f = a * a + b * b + c * c;

    return f;
  }

  // Default AD gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Booth's function
//
template<typename S, Index M = 2>
class Booth : public Function_Base<Booth<S, M>, S, M>
{
public:

  Booth() {}

  static constexpr
  char const * const
  NAME{"Booth"};

  using Base = Function_Base<Booth<S, M>, S, M>;

  // Explicit value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
  {
    assert(X.get_dimension() == Base::DIMENSION);

    T const &
    x = X(0);

    T const &
    y = X(1);

    T const
    a = x + 2 * y - 7;

    T const
    b = 2 * x + y - 5;

    T const
    f = a * a + b * b;

    return f;
  }

  // Default AD gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Goldstein-Price function
//
template<typename S, Index M = 2>
class GoldsteinPrice : public Function_Base<GoldsteinPrice<S, M>, S, M>
{
public:

  GoldsteinPrice() {}

  static constexpr
  char const * const
  NAME{"Goldstein-Price"};

  using Base = Function_Base<GoldsteinPrice<S, M>, S, M>;

  // Explicit value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
  {
    assert(X.get_dimension() == Base::DIMENSION);

    T const &
    x = X(0);

    T const &
    y = X(1);

    T const
    a = x + y + 1;

    T const
    b = 19 - 14 * x + 3 * x * x - 14 * y + 6 * x * y + 3 * y * y;

    T const
    c = 2 * x - 3 * y;

    T const
    d = 18 - 32 * x + 12 * x * x + 48 * y - 36 * x * y + 27 * y * y;

    T const
    e = 1 + a * a * b;

    T const
    f = 30 + c * c * d;

    T const
    fn = e * f;

    return fn;
  }

  // Default AD gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Failure function to test failed mechanism
//
template<typename S, Index M = 1>
class Failure : public Function_Base<Failure<S, M>, S, M>
{
public:

  Failure() {}

  static constexpr
  char const * const
  NAME{"Failure"};

  using Base = Function_Base<Failure<S, M>, S, M>;

  // Explicit value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
  {
    // Set the flag to signal that an unrecoverable error happened.
    this->set_failed("Testing failure mechanism");

    T const
    fn = 0.0;

    return fn;
  }

  // Default AD gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Non-monotonic function to test monotonicity enforcement.
//
template<typename S, Index M = 1>
class Mesa : public Function_Base<Mesa<S, M>, S, M>
{
public:

  Mesa() {}

  static constexpr
  char const * const
  NAME{"Mesa"};

  using Base = Function_Base<Mesa<S, M>, S, M>;

  // Explicit value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
  {
    T const &
    x = X(0);

    T
    y = x * x;

    if (-1.0 <= x && x <= 1.0) {
      y = y + 100.0;
    }

    return y;
  }

  // Default AD gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Function to test boundedness or residual enforcement.
//
template<typename S, Index M = 1>
class Sigmoid : public Function_Base<Sigmoid<S, M>, S, M>
{
public:

  Sigmoid() {}

  static constexpr
  char const * const
  NAME{"Sigmoid"};

  using Base = Function_Base<Sigmoid<S, M>, S, M>;

  // Explicit value.
  template<typename T, Index N>
  T
  value(Vector<T, N> const & X)
  {
    T const &
    x = X(0);

    T const
    x2 = x * x;

    T const
    x4 = x2 * x2;

    T const
    x8 = x4 * x4;

    T const
    x16 = x8 * x8;

    T const
    x32 = x16 * x16;

    T
    y = x * x32;

    return y;
  }

  // Default AD gradient.
  template<typename T, Index N>
  Vector<T, N>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

  // Default AD hessian.
  template<typename T, Index N>
  Tensor<T, N>
  hessian(Vector<T, N> const & x)
  {
    return Base::hessian(*this, x);
  }

};

//
// Functions to test constraint interface.
//

//
// Identity
//
template<typename S, Index NC, Index NV>
class Identity : public Equality_Constraint<Identity<S, NC, NV>, S, NC, NV>
{
public:

  Identity() {}

  static constexpr
  char const * const
  NAME{"Identity Map"};

  using Base = Equality_Constraint<Identity<S, NC, NV>, S, NC, NV>;

  // Explicit value.
  template<typename T, Index N>
  Vector<T, NC>
  value(Vector<T, N> const & x)
  {
    assert(x.get_dimension() == NV);
    return x;
  }

  // Default AD gradient.
  template<typename T, Index N>
  Matrix<T, NC, NV>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }
};

//
// A nonlinear function
//
template<typename S, Index NC = 3, Index NV = 5>
class Nonlinear01 : public Equality_Constraint<Nonlinear01<S, NC, NV>, S, NC, NV>
{
public:

  Nonlinear01() {}

  static constexpr
  char const * const
  NAME{"Nonlinear 01"};

  using Base = Equality_Constraint<Nonlinear01<S, NC, NV>, S, NC, NV>;

  // Explicit value.
  template<typename T, Index N = 5>
  Vector<T, NC>
  value(Vector<T, N> const & x)
  {
    assert(x.get_dimension() == NV);

    Vector<T, NC>
    c(Filler::ZEROS);

    c(0) = dot(x, x) - 10.0;

    c(1) = x(1) * x(2) - 5.0 * x(3) * x(4);

    c(2) = x(0) * x(0) * x(0) + x(1) * x(1) * x(1) + 1.0;

    return c;
  }

  // Default AD gradient.
  template<typename T, Index N = 5>
  Matrix<T, NC, NV>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }
};

//
// Circumference feasible region
//
template<typename S, Index NC = 1, Index NV = 2>
class Circumference : public Equality_Constraint<Circumference<S, NC, NV>, S, NC, NV>
{
public:

  Circumference(S const r, S const xc = S(0.0), S const yc = S(0.0)) : r_(r)
  {
    c_(0) = xc;
    c_(1) = yc;
  }

  static constexpr
  char const * const
  NAME{"Circumference"};

  using Base = Equality_Constraint<Circumference<S, NC, NV>, S, NC, NV>;

  // Explicit value.
  template<typename T, Index N = 2>
  Vector<T, NC>
  value(Vector<T, N> const & x)
  {
    assert(x.get_dimension() == NV);

    Vector<T, NC>
    f(Filler::ZEROS);

    f(0) = r_ * r_ - norm_square(x - c_);

    return f;
  }

  // Default AD gradient.
  template<typename T, Index N = 2>
  Matrix<T, NC, NV>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

private:
  S
  r_{0.0};

  Vector<S, NV>
  c_;
};

//
// Circle feasible region
//
template<typename S, Index NC = 1, Index NV = 2>
class Circle : public Inequality_Constraint<Circle<S, NC, NV>, S, NC, NV>
{
public:

  Circle(S const r, S const xc = S(0.0), S const yc = S(0.0)) : r_(r)
  {
    c_(0) = xc;
    c_(1) = yc;
  }

  static constexpr
  char const * const
  NAME{"Circle constraint"};

  using Base = Inequality_Constraint<Circle<S, NC, NV>, S, NC, NV>;

  // Explicit value.
  template<typename T, Index N = 2>
  Vector<T, NC>
  value(Vector<T, N> const & x)
  {
    assert(x.get_dimension() == NV);

    Vector<T, NC>
    f(Filler::ZEROS);

    f(0) = r_ * r_ - norm_square(x - c_);

    return f;
  }

  // Default AD gradient.
  template<typename T, Index N = 2>
  Matrix<T, NC, NV>
  gradient(Vector<T, N> const & x)
  {
    return Base::gradient(*this, x);
  }

private:
  S
  r_{0.0};

  Vector<S, NV>
  c_;
};

} // namespace minitensor
